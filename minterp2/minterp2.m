function [z] = minterp2(M, x, y, intmeth)
%function [z] = minterp2(M, x, y, intmeth)
%
% minterp2 :  very fast polynomial interpolation of gridded 
%             data (matrices). The interpolating polynomial
%             is of the general form:
%
%             P(x,y) =   (a0 + a1*x + a2*x^2 + a3*x^3 + ...)
%                      * (b0 + b1*y + b2*y^2 + b3*y^3 + ...)
%
% M :       a matrix to be interpolated (but see NOTE below)
% x,y :     vectors containing fractional column (x) and row (y)
%           indices at which M will be interpolated. Same convention
%           as in the MATLAB function 'interp2'.
% intmeth : (optional) a string selecting the interpolation method:
%           'linear', 'cubic', 'pchip' (or 'spline'), or 'lagrange'. 
%           Default is 'pchip'.
% z  :      vector with interpolation results; has the same shape
%           as x.
%
% NOTES:
% 1) The function has a mode that is useful when the same matrix must
% be interpolated many times. If the minterp2 function is called 
% with the matrix M as its sole argument, the matrix is stored 
% internally for subsequent interpolations:
%
%            minterp2(M);
%
% When minterp2 is afterwards called with an empty matrix argument,
% the previously stored matrix will be interpolated. For example,
%
%      z = minterp2([], x, y);
%
% interpolates a previously stored matrix at points (x(k),y(k)) using 
% pchip interpolation. 
%
% 2) The function returns NaN when the interpolation point is too
% close to the edge of the matrix of nodes for valid interpolation. 
%
% Supported interpolation methods are:
%
% linear  - bi-linear interpolation from map values of the 
%           four nearest neighboring nodes. Polynomial contains
%           the strictly linear terms and the 'x*y' term:
%           P(x,y) = a0 + a1*x + a2*y + a3*x*y
%
%                       o----o
%                       |  * |
%                       o----o
%
% cubic     - uses map values and all derivatives at 4 nearest neighboring 
%             points to calculate the 16 coefficients of a bi-cubic 
%             spline interpolation polynomial. Derivatives are computed 
%             from map data with central difference formulae.
%
%                  +----+----+----+
%                  |    |    |    |
%                  +----o----o----+
%                  |    |*   |    |
%                  +----o----o----+
%                  |    |    |    |
%                  +----+----+----+
%
% pchip     - piecewise cubic Hermite spline interpolation using 
%             the 16 nearest grid points. Slightly faster than cubic 
%             spline interpolation and also very good. DEFAULT.
%
% spline    - a synonym for 'pchip'
%
% lagrange  - cubic barycentric Lagrange interpolation using the 16 nearest
%             grid points (J.-P. Berrut and L. N. Trefethen,
%             "Barycentric Lagrange Interpolation", SIAM Review 46,
%             501-517, 2004)
%
% The cubic and pchip interpolations have a continuous first
% derivative.
%
% This function calls the MEX function 'minterp2', which must
% be compiled as follows:
%
% For Octave:  export CFLAGS="-O3 -march=native -fomit-frame-pointer"
%              mkoctfile --mex -s minterp2.c
%

% This software was developed at the National Institute of Standards and
% Technology by employees of the United States Federal Government in the
% course of their official duties. Pursuant to title 17 Section 105 of
% the United States Code this software is not subject to copyright
% protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its
% use by other parties, and makes no guarantees, expressed or implied,
% about its quality, reliability, or any other characteristic.
%
% Ulf Griesmann, NIST, October 2011 - September 2013

