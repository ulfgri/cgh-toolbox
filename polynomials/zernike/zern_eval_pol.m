function [vZ] = zern_eval_pol(vTh,vR,vPar)
% function [vZ] = zern_eval_pol(vTh,vR,vPar)
%
% Evaluates a Zernike polynomial for points specified in 
% normalized polar coordinates
%
% vZ      : values Zernike polynomial (defined according to ANSI Z80.28-2017)
%
% vTh     : angular coordinates of points (counter clockwise starting at
%           the positive x-axis
% vR      : normalized radius coordinates of points
% vPar    : Zernike coefficients listed in ANSI standard order (groups of
%           increasing radial order, for each order from highest odd (sin) to
%           highest even (cos) term)
%
% see zern_eval for a definition of the Zernike terms
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
% Author: Johannes Soons; NIST; Feb 2004 - Dec 2019 

% check input

  if ~isequal(size(vTh),size(vR))
    error(' The angle array vTh and radius array vR do not have the same size');
  end
  
  if all(size(vPar) > 1)
    error(' The Zernike parameters should be specified as a vector');
  end
  
% restrict evaluation to valid terms

  vb = ~isnan(vTh(:)) & ~isnan(vR(:));
  
  if ~any(vb)
    vZ = nan(size(vTh));
    return;
  else
    vZ = zeros(size(vTh));      
    vZ(~vb) = NaN;
  end
  
% return 0 if vPar is empty for non-nan pixel location

  if all(vPar == 0) || isempty(vPar), return; end
  
% obtain radial and angular orders

  nm = zern_index_trans((0:(numel(vPar)-1))','#ansi','nm');
  
% set maximum angular order  
  
  mMax = max(abs(nm(:,2)));
  
% process angular orders one by one

  vRmm1 = [];
  vRmm2 = [];

  for m = 0:mMax
      
%   extract parameters relevant to this angular order

    vParM = vPar(abs(nm(:,2)) == m);
    
%   compute Zernike terms belonging to a single angular order

    if rem(m,2) == 0
      [vDum1,vRmm1] = zern_eval_pol_m(vTh(vb),vR(vb),vParM(:),m,vRmm1,0);
    else  
      [vDum1,vRmm2] = zern_eval_pol_m(vTh(vb),vR(vb),vParM(:),m,vRmm2,0);
    end
    
    vZ(vb) = vZ(vb)+vDum1;
  end