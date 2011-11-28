/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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

Class
    eigenSolvers

Description
    Taken from the book Numerical Recipes in C++.

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef eigenSolvers_H
#define eigenSolvers_H

#include "scalar.H"
#include "scalarMatrices.H"

using namespace Foam;

//! Given a and b, calculates c using pythagoras.
Foam::scalar pythag
(
    const Foam::scalar& a,
    const Foam::scalar& b
)
{
    scalar absa=std::abs(a), absb=std::abs(b);

    return (absa > absb ? absa*std::sqrt(1.0+(absb/absa)*(absb/absa)) :
           (absb == 0.0 ? 0.0 : absb*std::sqrt(1.0+(absa/absb)*(absa/absb))));
}

//! Gives the sign of b to a.
Foam::scalar SIGN
(
    const Foam::scalar& a,
    const Foam::scalar& b
)
{
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

void QL
(
	const int& n,
    scalarList& d,
    scalarList& e,
    scalarSquareMatrix& z
)
{
    int m,l,iter,i,k;
    scalar s,r,p,g,f,dd,c,b;
    const scalar EPS=1e-30;

    for (i=1;i<n;i++) e[i-1]=e[i];
    e[n-1]=0.0;
    for (l=0;l<n;l++) {
        iter=0;
        do {
            for (m=l;m<n-1;m++) {
                dd=std::abs(d[m])+std::abs(d[m+1]);
                if (std::abs(e[m]) <= EPS*dd) break;
            }
            if (m != l) {
               if (iter++ == 30) throw("Too many iterations in tqli");
               g=(d[l+1]-d[l])/(2.0*e[l]);
               r=pythag(g,1.0);
               g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
               s=c=1.0;
               p=0.0;
               for (i=m-1;i>=l;i--) {
                   f=s*e[i];
                   b=c*e[i];
                   e[i+1]=(r=pythag(f,g));
                   if (r == 0.0) {
                       d[i+1] -= p;
                       e[m]=0.0;
                       break;
                   }
                   s=f/r;
                   c=g/r;
                   g=d[i+1]-p;
                   r=(d[i]-g)*s+2.0*c*b;
                   d[i+1]=g+(p=s*r);
                   g=c*r-b;
                   for (k=0;k<n;k++) {
                        f=z[k][i+1];
                        z[k][i+1]=s*z[k][i]+c*f;
                        z[k][i]=c*z[k][i]-s*f;
                   }
               }
               if (r == 0.0 && i >= l) continue;
               d[l] -= p;
               e[l]=g;
               e[m]=0.0;
            }
       } while (m != l);
   }
}

#endif
