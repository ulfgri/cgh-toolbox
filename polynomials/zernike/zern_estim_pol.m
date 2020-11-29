function [vPar,vRes,mCov] = zern_estim_pol(vTh,vR,vZ,vjTerm,bNormal)
% function [vPar,vRes,mCov] = zern_estim_pol(vTh,vR,vZ,vjTerm,bNormal)
%
% Estimate least-squares best-fit Zernike coefficients
% for points defined in normalized polar coordinates
%
% See zern_eval for the definition of the Zernike terms
%
% See zern_map2pol for the calculation of normalized polar coordinates
% for an image array
%
% vPar    : vector with estimated Zernike coefficients (in ANSI Z80.28-2017
%           order and augmented with zeros)
% vRes    : residual values (vZ-Zernike)
% mCov    : estimated variance covariance matrix of the parameters
%           (NOT ordered nor padded with zeros. mCov(i,i) is the estimated
%           variance for the parameter vjTerm(i))
%
% vTh     : angular coordinates of points (counter clockwise starting at
%           the positive x-axis
% vR      : normalized radius coordinates of points
% vZ      : surface values to fit
% vjTerm  : vector with the indices of the Zernike terms (defined according
%           to ANSI Z80.28-2017, with 0 indicating the piston term)
% bNormal : (optional, default 0) perform fit under the assumption that
%           the Zernike polynomial terms are orthogonal over the 
%           domain of the points. This approximation may be useful
%           for cases where the number of Zernike coefficients is large.
%           The estimation will be performed in a slow but memory efficient
%           manner.

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
% Version: 2.0
% Author: Johannes Soons; NIST; Feb 2004 - Sept 2019 
% Review: Johannes Soons; NIST; Sept 2019 
% Status: TESTING
% ---------------------------------------------------------

% set default values

  if nargin < 5, bNormal = []; end
  if isempty(bNormal), bNormal = 0; end
  
% check input

  if ~isequal(size(vTh),size(vR))
    error(' The angle array vTh and radius array vR do not have the same size');
  end

  if ~isequal(size(vTh),size(vZ))
    error(' The angle array vTh and data array vZ do not have the same size');
  end
  
  if isempty(vjTerm)
    vPar = [];
    mCov = [];
    vRes = vZ;
    return;
  end
  
  if numel(unique(vjTerm(:))) ~= numel(vjTerm)
    error(' vjTerm contains duplicate entries');
  end
  
  vb = (~isnan(vTh(:)) & ~isnan(vR(:))) & ~isnan(vZ(:));
  
  if ~any(vb)
    vRes = nan(size(vZ));
    vPar = nan(max(vjTerm(:))+1,1);
    
    if nargout > 2, mCov = nan(numel(vjTerm)); end
    return;
  end
  
% initialize  

  vRes = vZ;
  vRes(~vb) = NaN;
  
  if nargout > 2, mCov = zeros(numel(vjTerm)); end
 
% construct expanded parameter vector with all relevant 
% radial and angular orders

  vPar = zeros(max(vjTerm(:))+1,1);
  
% set parameter vector with 1 for parameters to be estimated  
  
  vParDes = vPar;
  vParDes(vjTerm(:)+1) = 1;
  
% set reverse index ( i in vParDes to j in vjTerm )

  viParDes = zeros(numel(vParDes),1);
  viParDes(vjTerm(:)+1) = (1:numel(vjTerm))';
  
% set properties of terms to estimate  

  nm = zern_index_trans((1:numel(vPar))'-1,'#ansi','nm');  
  mMax = max(nm(:,2));
  
  if bNormal == 1
      
%   memory efficient term by term estimation asssuming orthogonality

%   process angular orders one by one

    for m = 0:mMax
      
%     identify parameters relevant to this angular order

      vi = find(abs(nm(:,2)) == m);
      
      if ~isempty(vi) && any(vParDes(vi) ~= 0)
          
%       set relevant radial orders

        vn = unique(nm(vi(vParDes(vi) ~= 0),1));
      
        for i = 1:numel(vn)
          
%         get terms belonging to this radial order          
          
          vj = find(nm(vi,1) == vn(i));
        
%         set parameter vector for Zernike terms with angular order m,
%         where 1 indicates that the term will be estimated
        
          vParM = zeros(numel(vi),1);
          vParM(vj) = vParDes(vi(vj));
        
%         compute relevant columns design matrix

          mX = zern_eval_pol_m(vTh(vb),vR(vb),vParM(:),m,[],1);
          
          k = 0;
        
          for j = 1:numel(vj)
              
            if vParM(vj(j)) ~= 0

%             estimate coefficients                
                 
              k = k+1;
            
              vPar(vi(vj(j))) = mX(:,k)\vZ(vb);
              
%             update residuals
          
              vRes(vb) = vRes(vb)-mX(:,k)*vPar(vi(vj(j)));
              
%             set entry for covariance matrix

              p = viParDes(vi(vj(j)));

              mCov(p,p) = 1./sum(mX(:,k).^2);
            end
          end  
        end
      end
    end
    
  else
      
%   full least squares solution

    mX = zeros(sum(vb),numel(vjTerm));
    
%   process angular orders one by one

    vRmm1 = [];
    vRmm2 = [];

    for m = 0:mMax
      
%     extract parameters relevant to this angular order

      vi = find(abs(nm(:,2)) == m);
      
%     set parameter vector for Zernike terms with angular order m,
%     where 1 indicates that the term will be estimated
      
      vParM = vParDes(vi);
      
%     compute design matrix colomns for this angular order

      vbp = vParM ~= 0;
      
      if rem(m,2) == 0
        [mX(:,viParDes(vi(vbp))),vRmm1] = zern_eval_pol_m(vTh(vb),vR(vb),vParM(:),m,vRmm1,1);
      else  
        [mX(:,viParDes(vi(vbp))),vRmm2] = zern_eval_pol_m(vTh(vb),vR(vb),vParM(:),m,vRmm2,1);
      end
    end
    
%   solve system

    vPar(vjTerm(:)+1) = mX\vZ(vb);
    
%   set covariance matrix

    if nargout > 2
      mCov = inv(mX'*mX);
    end
    
%   calculate residuals

    vRes(vb) = vRes(vb)-mX*vPar(vjTerm(:)+1);
  end    
  
% muliply coveriance matrix with estimated variance residuals

  if nargout > 2
    s2 = sum(vRes(vb).^2)./(sum(vb)-numel(vjTerm));
    mCov = mCov*s2;
  end