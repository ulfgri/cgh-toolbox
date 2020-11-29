function [vZ] = zern_eval_pos(mPos,vCenter,radius,vPar,bAll)
%function [vZ] = zern_eval_pos(mPos,vCenter,radius,vPar,bAll)
%
% Evaluates a Zernike Polynomial at points specified in  
% cartesian coordinates
%
% vZ      : values Zernike polynomial (defined according to ANSI
%           Z80.28-2017)
%
% mPos    : [x,y] coordinates of points, one per row
% vCenter : [x,y] coordinates of the center of the Zernike normalization disk
% radius  : radius of the Zernike normalization disk
% vPar    : Zernike coefficients listed in ANSI standard order (groups of
%           increasing radial order, for each order from highest odd (sin) to
%           highest even (cos) term)
% bAll    : (optional, default 1) evaluate all radii, i.e. including 
%           normalized radii larger than one.

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

% set defaults

  if nargin < 5, bAll = []; end
  if isempty(bAll), bAll = 1; end
  
% check input

  if isempty(mPos), vZ = []; return; end
  
% convert to polar coordinates  
  
  [vTh,vR] = zern_pos2pol(mPos,vCenter,radius);
  
% evaluate Zernike terms  
  
  if bAll == 1
    vZ = zern_eval(vTh,vR,vPar);      
  else
    vZ = nan(size(mPos,1),1);
    vb = vR <= 1;
    if any(vb)
      vZ(vb) = zern_eval(vTh(vb),vR(vb),vPar);
    end
  end