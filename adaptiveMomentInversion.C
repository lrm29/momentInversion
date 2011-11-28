/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application

    adaptiveMomentInversion

Description

    Test the implementation of the "Chebyshev algorithm". AKA "Wheeler's
    algorithm". Not sure what name would be appropriate.

    Written on Wednseday 17th November 2011
    By Laurence R. McGlashan
    l.mcglashan@sgi.com

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "cpuTime.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "scalarList.H"
#include "scalarMatrices.H"
#include "eigenSolvers.H"


using namespace Foam;


bool momentSetIsValid(const scalarList& momentList)
{
    label L = 0.5*(momentList.size() - 2);
    labelList pivotIndices(L + 1);
    scalarSquareMatrix hankelHadamardMatrix(L + 1, L + 1, 0.0);

    for (label m = 0; m < 2; ++m)
    {
        for (label row = 0; row < L + 1; ++row)
        {
            for (label col = 0; col < L + 1; ++col)
            {
                hankelHadamardMatrix[row][col] = momentList[row + col + m];
            }
        }

        LUDecompose(hankelHadamardMatrix, pivotIndices);

        // This checks whether the number of row evaluations is
        // odd or even. See numerical recipes.
        label d = 1;
        for (label size = 0; size < L + 1; ++size)
        {
            if (size != pivotIndices[size])
            {
                d *= -1.0;
            }
        }

        scalar sumDiagonal = d;
        for (label diag = 0; diag < L + 1; ++diag)
        {
            sumDiagonal *= hankelHadamardMatrix[diag][diag];
        }

        if (sumDiagonal < 0)
        {
            return false;
        }
    }

    return true;
}


label checkMoments
(
    const label N,
    const scalarList& moments,
    const scalar rMin
)
{
    if (moments[0] < 0)
    {
        FatalErrorIn("checkMoments")
            << "Number Density, m_0 < 0!"
            << abort(FatalError);
    }
    else if (moments[0] == 0)
    {
        WarningIn("checkMoments")
            << "Number Density, m_0 = 0. Returning w_i = p_i = 0."
            << endl;

        return 1;
    }
    else if (N == 1 || moments[0] < rMin)
    {
        WarningIn("checkMoments")
            << "Only one node specified or number density is too small." << nl
            << "Nodes = " << N << ", OR m_0(" << moments[0] << ") < " << rMin
            << endl;
        return 2;
    }

    return 0;
}


bool checkRealisability(const scalarList& b)
{
    if (min(b) < 0)
    {
        FatalErrorIn("checkRealisability")
            << "Non-realisable moments have been generated."
            << abort(FatalError);
    }
    return true;
}


void constructRecurrenceMatrix
(
    const label N,
    const scalarList& moments,

    scalarSquareMatrix& sig,
    scalarList& a,
    scalarList& b
)
{
    for (label i = 1; i < 2*N + 1; ++i)
    {
        sig[1][i] = moments[i-1];
    }

    a[0] = moments[1]/moments[0];
    b[0] = 0;

    for (label k = 2; k < N + 1; ++k)
    {
        for (label l = k; l < 2*N - k + 2; ++l)
        {
            sig[k][l] = sig[k-1][l+1]
                      - a[k-2]*sig[k-1][l]
                      - b[k-2]*sig[k-2][l];
        }
        a[k-1] = sig[k][k+1]/sig[k][k] - sig[k-1][k]/sig[k-1][k-1];
        b[k-1] = Foam::sqrt(sig[k][k]/sig[k-1][k-1]);
    }
}


scalarList chebyshev
(
    label& N,
    const scalarList& moments,
    const scalar cutoff,
    const scalar rMin
)
{
    const label oldN = N;

    label resultOfCheck = checkMoments(N, moments, rMin);

    if (resultOfCheck == 1)
    {
        return scalarList(2, 0.0);
    }
    else if (resultOfCheck == 2)
    {
        scalarList temp(2, moments[0]);
        temp[1] = moments[1]/moments[0];
        return temp;
    }

    scalarSquareMatrix sig(2*N + 1, 2*N + 1);
    scalarList a(N, 0.0);
    scalarList b(N, 0.0);

    constructRecurrenceMatrix(N, moments, sig, a, b);

    // Determine maximum number of quadrature nodes.
    for (label k = N; k >= 2; --k)
    {
        if (sig[k][k] <= cutoff)
        {
            N = k - 1;
            if (N == 1)
            {
                WarningIn("chebyshev")
                    << "Only one quadrature node required. Returning."
                    << endl;
                scalarList temp(2, moments[0]);
                temp[1] = moments[1]/moments[0];
                return temp;
            }
        }
    }

    // Calculate weights and positions with the new N.

    if (N != oldN)
    {
        constructRecurrenceMatrix(N, moments, sig, a, b);
    }

    checkRealisability(b);

    // Set up Jacobi matrix
    for (label n1 = N - 1; n1 >= 0; n1--)
    {
        if (n1 == 0)
        {
            WarningIn("chebyshev")
                << "Only one quadrature node required. Returning."
                << endl;
            N = 1;
            scalarList temp(2, moments[0]);
            temp[1] = moments[1]/moments[0];
            return temp;
        }

        scalarList WandA(2*(n1 + 1), 0.0);
        scalarList diagonal(n1 + 1, 0.0);
        scalarList subDiag(n1 + 1, 0.0);
        scalarSquareMatrix z(n1 + 1, n1 + 1, 0.0);

        for (label i = 0; i < n1 + 1; ++i)
        {
            diagonal[i] = a[i];
            subDiag[i] = b[i];
            z[i][i] = 1.0;
        }

        // Calculate eigenvalues and eigenvectors.
        QL(n1 + 1, diagonal, subDiag, z);

        scalarList w(n1 + 1, 0.0);
        for (label i = 0; i < n1 + 1; ++i)
        {
            w[i] = moments[0]*Foam::sqr(z[0][i]);
            WandA[i] = w[i];
            WandA[i + n1 + 1] = diagonal[i];
        }

        // @todo Put in some adaptive stuff here.
        scalarList dab(n1 + 1, 0.0);
        scalarList mab(n1 + 1, 0.0);
        for (label i = n1 + 1; i >= 2; i--)
        {
            //dab[i] = min(mag(
        }

        if (min(w)/max(w) > rMin)
        {
            N = n1 + 1;
            return WandA;
        }
        else
        {
            for (label wpI = 0; wpI < n1 + 1; ++wpI)
            {
                Info<< "Failed Inversion: " << scientific 
                    << "    w_" << wpI << ": " << WandA[wpI]
                    << ",   p_" << wpI << ": " << WandA[wpI + n1 + 1] << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append
    (
        "Number of Moments"
    );
    argList::validArgs.append
    (
        "Moments File"
    );
    argList args(argc, argv);

    cpuTime cpu;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Start: Adaptive Moment Inversion Algorithm." << nl << endl;

    // Number of moments
    const label M = args.argRead<label>(1);
    const fileName momentsFile = args.argRead<fileName>(2);

    // Number of nodes
    const label N = M/2;

    const scalar cutoff = 0.0;
    const scalar rMin = 1e-3;

    scalarList moments;

    IFstream file(momentsFile);
    while (file.good())
    {
        scalar temp = 0;
        if (file.read(temp)) moments.append(temp);
    }

    if (moments.size() != M)
    {
        FatalErrorIn("chebyshev")
            << "Size of moment list, " << moments.size()
            << ", does not equal number of moments, " << M
            << abort(FatalError);
    }

    Info<< "Number of moments: " << M << endl;
    Info<< "Initial number of quadrature nodes: " << N << endl;

    Info<< "Moments:" << endl;
    forAll(moments, mI)
    {
        Info<< "    m_" << mI << ": " << moments[mI] << endl;
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (!momentSetIsValid(moments))
    {
        WarningIn("chebyshev")
            << "Moments are not valid!"
            << endl;
    }

    label finalN = N;

    // Returns Weights and Abscissas.
    // Returns the number of quadrature points through finalN.
    scalarList WandA = chebyshev(finalN, moments, cutoff, rMin);

    if (finalN != N)
    {
        Info<< "Number of quadrature nodes has been changed from "
            << N << " to " << finalN << endl;
    }

    Info<< nl << "Weights(w) and Abscissas(p):" << endl;
    for (label wpI = 0; wpI < finalN; ++wpI)
    {
        Info<< scientific
            << "    w_" << wpI << ": " << WandA[wpI]
            << ",   p_" << wpI << ": " << WandA[wpI + finalN] << endl;
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << cpu.elapsedCpuTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;

}

// ************************************************************************* //
