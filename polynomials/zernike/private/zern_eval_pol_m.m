function [vZ,vR22] = zern_eval_pol_m(vTh,vR,vPar,m,vRmm,bDes)
%function [vZ,vR22] = zern_eval_pol_m(vTh,vR,vPar,m,vRmm,bDes)
%
% zern_eval_pol_m: Internal function for the evaluation of Zernike terms
%   belonging to a single angular order
%
% vZ     : values Zernike polynomial (defined according to ANSI Z80.28-2017),
% vR22   : matrix or vector with R{m+2,m+2} radial polynomial values (see [2])
%
% vTh    : matrix or vector with angles for each point (defined relative to
%          the positive x -axis, positive in counter clockwise direction).
% vR     : matrix or vector with normalized radii for each point
% vPar   : Zernike coefficients for the specified angular order
%          listed in ANSI standard order (lowest radial order to highest
%          radial order; sin before cos for the same radial order)
% m      : absolute value angular order (frequency)
% vRmm   : (optional, default []) matrix or vector with R{m,m} radial
%          polynomial values (see R{q,q} in [2])
% bDes   : (optional, default 0) return a matrix where each column contains
%          the Zernike polynomial term for which vPar is nonzero.
%
% This function implements the modified Kintner's method for the fast computation 
% of the Zernike terms [2].
%
% References:
%
% [1] ANSI Z8028-2017, (2017), "Ophthalmics - methods for reporting optical
%     aberrations of eyes," American National Standards Institute, New York, NY.
% [2] C.W. Chong, P. Raveendran, R. Mukundan, (2003), "A comparative analysis
%     of algorithms for fast computation of Zernike moments," Pattern Recognition
%     36(3):731–742.
% [3] C. Singh, E. Walia, (2010), "Fast and numerically stable methods for the
%     computation of Zernike moments," Pattern Recognition 43(7):2497–2506.

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

  if nargin < 6, bDes = []; end
  if nargin < 5, vRmm = []; end
  
  if isempty(bDes), bDes = 0; end
  
% check input

  if ~isequal(size(vTh),size(vR))
    error('The angle array vTh and radius array vR do not have the same size');
  end
  
  if m < 0
    error(' m should be larger or equal to zero');
  end
  
% initialize

  if bDes == 0
    vZ = zeros(size(vR));
  else
    vZ = zeros(numel(vR),sum(vPar ~= 0));
  end
  
  vR22 = [];
  
  if all(vPar == 0), return; end
  
% set maximum radial order  

  if m == 0
    nMax = 2*numel(vPar)-2;
  else
      
%   enforce an even number of parameters

    if rem(numel(vPar),2) ~= 0
      vPar(end+1) = 0;
    end
    
    nMax = 2*numel(vPar)/2+m-2;
    
  end
 
% initialize

  mQ = zeros(numel(vR),2);
  
  ip = 1;
  id = 1;
  i1 = 1; 
  
  vSin = [];
  vCos = [];
  
% Calculate radial polynomial using modified Kintner's method [2]
  
  for p = m:2:nMax
      
%   mQ(:,i1) will contain R{p-2,m}
%   mQ(:,i2) will contain R{p-4,m}
%
%   increase i1 and i2 index values

    i1 = mod(i1,2)+1;
    i2 = mod(i1,2)+1;    
    
    switch p
      case m
          
%       calculate R{m,m) using equation 3.9 in [2]

        if isempty(vRmm)
          vRmm = vR(:).^m;
        end  
        
        mQ(:,i2) = vRmm;
        
      case m+2
          
%       calculate R{m+2,m} using equation 3.10 in [2]

        vR2 = vR(:).^2;  
        
        vR22 = vRmm.*vR2;
        
        mQ(:,i2) = p*vR22-(p-1)*mQ(:,i1);
        
      otherwise
          
%       calculate R{m+4,m} using equation 3.8 in [2]            
            
        k1 = (p+m)*(p-m)*(p-2)/2;
        k2 = 2*p*(p-1)*(p-2);
        k3 = -m^2*(p-1)-p*(p-1)*(p-2);
        k4 = -p*(p+m-2)*(p-m-2)/2;
        
        mQ(:,i2) = ((k2*vR2+k3).*mQ(:,i1)+k4*mQ(:,i2))./k1;
    end    
  
%   evaluate Zernike term (multiplication with normalization term
%   and sin or cos

    if m == 0
      if vPar(ip) ~= 0  
        if bDes == 1
            
%         add Zernike term to design matrix
  
          vZ(:,id) = mQ(:,i2)*sqrt(p+1);
          id = id+1;
        else
            
%         add Zernike term to previous sum

          vZ(:) = vZ(:)+vPar(ip)*mQ(:,i2)*sqrt(p+1);
        end
      end
      ip = ip+1;
    else
      if bDes == 1
        
%       add sin and cos terms to design matrix

        if vPar(ip) ~= 0
          if isempty(vSin), vSin = sin(m*vTh); end  
            
          vZ(:,id) = sqrt(2*(p+1))*mQ(:,i2).*vSin;
          id = id+1;
        end
        
        if vPar(ip+1) ~= 0
          if isempty(vCos), vCos = cos(m*vTh); end              
            
          vZ(:,id) = sqrt(2*(p+1))*mQ(:,i2).*vCos;
          id = id+1;
        end
      else
          
%       add sin and cos Zernike terms

        if vPar(ip) ~= 0
          if isempty(vSin), vSin = sin(m*vTh); end              
          
          vZ(:) = vZ(:)+vPar(ip)*sqrt(2*(p+1))*mQ(:,i2).*vSin;
        end
        
        if vPar(ip+1) ~= 0
          if isempty(vCos), vCos = cos(m*vTh); end              
          vZ(:) = vZ(:)+vPar(ip+1)*sqrt(2*(p+1))*mQ(:,i2).*vCos;
        end
      end
      ip = ip+2;
    end

%   determine if there is more work to do
    
    if all(vPar(ip:end) == 0)
      return;
    end
  end  
  
