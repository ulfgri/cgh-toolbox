function [vPar,mapRes,mCov] = zern_estim(map,vCenter,radius,vjTerm,vPixSep,bNormal)
% function [vPar,mapRes,mCov] = zern_estim(map,vCenter,radius,vjTerm,vPixSep,bNormal)
%
% Estimate least-squares best-fit Zernike coefficients
% and residual for an image array
%
% See zern_eval for the definition of the Zernike terms
%
% vPar     : vector with estimated Zernike coefficients (in ANSI Z80.28-2017
%            order and augmented with zeros)
% mapRes   : residual values (map-Zernike)
% mCov     : estimated variance covariance matrix of the parameters
%            (NOT ordered nor padded with zeros. mCov(i,i) is the estimated
%            variance for the parameter vjTerm(i))
%
% map      : image array
% vCenter  : [x,y] coordinates of the center of the Zernike normalization
%            disk relative to the image origin
% radius   : radius of the Zernike normalization disk
% vjTerm   : vector with the indices of the Zernike terms to be estimated
%            (defined according to ANSI Z80.28-2017, with 0 indicating the
%            piston term)
% vPixSep  : (optional, default [1,1]) pixel separation in x and y
% bNormal  : (optional, default 0) perform fit under the assumption that
%            the Zernike polynomial terms are orthogonal over the 
%            domain of the points. This approximation may be useful
%            for cases where the number of Zernike coefficients is large.
%            The estimation will be performed in a slow but memory efficient
%            manner.


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

% set default values

  if nargin < 6, bNormal = []; end
  if nargin < 5, vPixSep = []; end
  
  if isempty(bNormal), bNormal = 0; end

  if isempty(vPixSep), vPixSep = 1; end
  if numel(vPixSep) == 1, vPixSep = [vPixSep,vPixSep]; end
  
% change to polar coordinates

  [mTh,mR] = zern_map2pol(size(map),vCenter,radius,vPixSep);
  
% estimate parameters

  if nargout > 2
    [vPar,mapRes,mCov] = zern_estim_pol(mTh(:),mR(:),map(:),vjTerm,bNormal);
  else
    [vPar,mapRes] = zern_estim_pol(mTh(:),mR(:),map(:),vjTerm,bNormal);
  end
  
  mapRes = reshape(mapRes,size(map));
