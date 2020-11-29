function [vPar2] = zern_transform(vPar,vTrans,scale,rAngle)
% function [vPar2] = zern_transform(vPar,vTrans,scale,rAngle)
%
% Transform Zernike coefficients to describe a translation, rotation,
% and lateral scaling of the Zernike polynomial
%
% vPar2   : Transformed Zernike coefficients 
%
% vPar    : Zernike coefficients
% vTrans  : (optional, default [0,0]) [x,y] displacement center
%           (normalized by original Zernike normalization radius)
% scale   : (optional, default 1) scaling factor
% rAngle  : (optional, default 0) rotation angle (postive in counter
%           clockwise direction)
%
% The Zernike polynomial is transformed in the following manner:
%
%  - lateral scaling (stretching) of the polynomial relative
%    to the center of the Zernike normalization disk by a
%    factor scale
%  - rotating the polynomial by an angle rAngle around the center
%    of the Zernike normalization disk
%  - translating the polynomial by translation vTrans
%
% See also zern_subaperture for calculating the Zernike
% coefficients of a subaperture
%

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
% Author: Johannes Soons; NIST; Feb 2004 - Dec 2019 

% check input

  if nargin < 2, vTrans = []; end
  if nargin < 3, scale = []; end
  if nargin < 4, rAngle = []; end
  
  if isempty(vTrans), vTrans = [0,0]; end
  if isempty(scale), scale = 1; end
  if isempty(rAngle), rAngle = 0; end
  
% check input

  if numel(vTrans) ~= 2
    error(' Translation vector should contain two values');
  end
  
  if isempty(vPar), vPar2 = vPar; return; end
  
% calculate modified translation
%
% -T' = inv(S)*inv(R)*(T)

  vT(1) = cos(rAngle)*vTrans(1)+sin(rAngle)*vTrans(2);
  vT(2) = -sin(rAngle)*vTrans(1)+cos(rAngle)*vTrans(2);
  
  vT = vT/scale;
  
% calculate transformed Zernike parameters

  vPar2 = zern_subaperture(vPar,-vT,1/scale,-rAngle);