function [mZ] = zern_eval(vSize,vCenter,radius,vPar,vPixSep,bAll)
% function [mZ] = zern_eval(vSize,vCenter,radius,vPar,vPixSep,bAll)
%
% Evaluation of a Zernike polynomial
%
% mZ      : array with values Zernike polynomial (defined according to 
%           ANSI Z80.28-2017)
%
% vSize   : [row,column] size of the image array to generate
% vCenter : [x,y] coordinates of the center of the Zernike normalization
%           disk relative to the image origin
% radius  : radius of the Zernike normalization disk (same units as
%           vCenter)
% vPar    : Zernike coefficients listed in ANSI standard order (groups of
%           increasing radial order, for each order from highest odd (sin) to
%           highest even (cos) term)
% vPixSep : (optional, default [1,1]) pixel separation in x and y
% bAll    : (optional, default 1) evaluate all radii, i.e. including radii
%           larger than one.
%
% Example: mZ = zern_eval([1001,1001],[501,501],500,[0;0;0;0;1],[1,1],0);
%
% The Zernike Polynomial is defined by:
%
%             ___k   ___n
%             \      \                                 sin(|m|*th)
%   Z(r,th) = /__    /__     a    * N    *  R   (r) *      or   
%              n=0    m=-n    n,m    n,m     n,m       cos(m*th)
%
% where cos is used for m > 0 and sin for m < 0. n equals the degree
% of the radial polynomial and m is the angular dependence (exponential)
% parameter. 
%
% The normalization constant N    is defined as:
%                             n,m  
%            2(n+1)   (1/2)
%   N    = ( ------- )      where d    = 1 if m=0 and 0 otherwise.
%    n,m     1+d                   m,0
%               m,0 
%
% The radial polynomial R   (r) equals 0 if n-|m| is odd, otherwise:
%                        n,m 
%
%                 (n-|m|)/2
%                 \          k             (n-k)!              n-2k
%   R   (r)  =    /__    (-1)  ------------------------------ r
%    n,m          k=0          k!((n+|m|)/2-k)!((n-|m|)/2-k)! 
%
%
% Examples
%
% Z#   n   m    Z
%                n,m
%
%  0   0   0    1                                 Piston
%  1   1  -1    2*r*sin(th)                       Tilt Y 
%  2   1   1    2*r*cos(th)                       Tilt X
%  3   2  -2    sqrt(6)*r^2*sin(2*th)             Astigmatism Y
%  4   2   0    sqrt(3)*(2*r^2-1)                 Power
%  5   2   2    sqrt(6)*r^2*cos(2*th)             Astigmatism X
%  6   3  -3    sqrt(8)*r^3*sin(3*th)             Trefoil Y
%  7   3  -1    sqrt(8)*(3*r^3-2*r)*sin(th)       Coma Y
%  8   3   1    sqrt(8)*(3*r^3-2*r)*cos(th)       Coma X
%  9   3   3    sqrt(8)*r^3*cos(3*th)             Trefoil X
%  10  4  -4    sqrt(10)*r^4*sin(4*th)            Tetrafoil Y
%  11  4  -2    sqrt(10)*(4*r^4-3*r^2)*sin(2*th)  Secondary astigmatism Y
%  12  4   0    sqrt(5)*(6*r^4-6*r^2+1)           Primary spherical
%  13  4   2    sqrt(10)*(4*r^4-3*r^2)*cos(2*th)  Secondary astigmatism X
%  14  4   4    sqrt(10)*r^4*cos(4*th)            Tetrafoil X
%
% This function implements the modified Kintner's method for the fast computation 
% of the Zernike terms [2].
%
% References:
%
% [1] ANSI Z8028-2017, (2017), "Ophthalmics - methods for reporting optical
%     aberrations of eyes," American National Standards Institute, New York, NY.
% [2] C. Singh, E. Walia, (2010), "Fast and numerically stable methods for the
%     computation of Zernike moments," Pattern Recognition 43(7):2497–2506.
% [3] C.W. Chong, P. Raveendran, R. Mukundan, (2003), "A comparative analysis
%     of algorithms for fast computation of Zernike moments," Pattern Recognition
%     36(3):731–742.
%
% See also zern_par_trans for the relation with other Zernike polynomial
% definitions

% This software was developed by employees of the National Institute
% of Standards and Technology (NIST), an agency of the Federal Government,
% and is being made available as a public service. Pursuant to title 17
% United States Code Section 105, works of NIST employees are not subject
% to copyright protection in the United States.  This software may be subject
% to foreign copyright.  Permission in the United States and in foreign
% countries, to the extent that NIST may hold copyright, to use, copy,
% modify, create derivative works, and distribute this software and its
% documentation without fee is hereby granted on a non-exclusive basis,
% provided that this notice and disclaimer of warranty appears in all copies.
%
% THIS SOFTWARE IS BEING PROVIDED FOR RESEARCH AND EDUCATIONAL PURPOSES ONLY.
% THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
% EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY
% THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM FROM
% INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO THE
% SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE.  IN NO EVENT
% SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT,
% INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM,
% OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY,
% CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS
% OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE
% OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.
%
% Author: Johannes Soons; NIST; Sep 2000 - Dec 2019 

% set defaults

  if nargin < 6, bAll = []; end
  if nargin < 5, vPixSep = []; end
  
  if isempty(bAll), bAll = 1; end  
  
  if isempty(vPixSep), vPixSep = 1; end
  if numel(vPixSep) == 1, vPixSep = [vPixSep,vPixSep]; end
  
% calculate radius and angle of each image point

  [mTh,mR] = zern_map2pol(vSize,vCenter,radius,vPixSep);
  
% calculate Zernike polynomial values

  if bAll == 1
    mZ = zern_eval_pol(mTh,mR,vPar);
  else
    vb = mR(:) <= 1;
    mZ = nan(vSize);
    mZ(vb) = zern_eval_pol(mTh(vb),mR(vb),vPar);
  end