function [vpar,rmse,spar,tpar,ppar,mcov] = fdh_estim(map,mr,ma,nm,weights,salgor,snoise)
% function [vpar,rmse,spar,tpar,ppar,mcov] = fdh_estim(map,mr,ma,nm,weights,salgor,snoise)
%
% Estimates coefficients of a decomposition in Free Disk Harmonic
% (FDH) polynomials
%
% map     : array with heights (column vector, matrix, or three dimensional array). If the number
%           of dimensions of map exceeds that of mr or ma, map is interpreted as a stack of height
%           arrays. Each column k of vpar will contain the estimated FDH parameters of map(:,:,k)
%           or map(:,k).
% mr      : matrix or column vector with normalized radii at each point of the image
% ma      : matrix or column vector with angles at each point of the image
%           (defined relative to y, x corresponds to pi/2)
% nm      : Radial and azimuthal order for each disk harmonic term (one term per row)
% vpar    : matrix or column vector with estimated parameters (padded with zeros) 
% weights : optional, array having the same size as mr or image used to weigh
%           each observation (e.g., the intensity image). The weights are assumed to 
%           equal the inverse of the noise variance at each point. 
% salgor  : optional, algorithm used to solve the least squares problem
%           'qr'    : Solution by QR decomposition of the regression matrix using
%                     Householder transformations. This is a numerically accurate method but very 
%                     memory intensive. Use for many FDH terms may yield memory problems.
%           'normal': DEFAULT. This algorithm explicitly calculates the normal equations (X'X) b = X'y
%                     before solving them. The algorithm does not assume orthogonality but
%                     is numerically inferior to the QR procedure. This option is fast and
%                     not memory intensive, allowing for many FDH terms to be estimated.
%           'orthogonal' : This algorithm assumes orthogonality of the FDH terms. In practice
%                     this is only an approximation due to the discretization of the image.
%                     The algorithm explicitly calculates the relevant elements of the normal
%                     equations. This option is fast and not memory intensive, allowing for many
%                     FDH terms to be estimated. 
% snoise  : optional, estimate of the standard deviation of uncorrelated measurement noise.
%           If present, this standard deviation is used to estimate the statistical properties of 
%           each estimated parameter. If snoise is not defined, the root mean squared residual
%           error is used instead, which is not an accurate measure if the model is incomplete.
% rmse    : root mean squared residual error
% spar    : estimated standard deviation for each estimated parameter
% tpar    : t-statistic for each estimated parameter (ratio of parameter value to its standard deviation)
%           As the model is close to being orthogonal, this parameter can be used to evaluate the 
%           significance of individual parameters and subsets of parameters.
% ppar    : probability that an estimated parameter is significant (has a value unequal to zero), 
%           assuming Gaussian uncorrelated noise.
% mcov    : Covariance matrix estimated paramers
%  
% The weights are implemented to yield the following problem:
%
%   Minimize the sum of  weights .* (map - fitted_map).^2
%
% References:
%
% N. M. Milton and M. Lloyd-Hart, "Disk Harmonic Functions for Adaptive Optics Simulations,"
% in Adaptive Optics: Analysis and Methods/Computational Optical Sensing and Imaging/Information
% Photonics/Signal Recovery and Synthesis Topical Meetings on CD-ROM, Technical Digest
% (Optical Society of America, 2005), paper AWA3 .
%
% S. C. Verrall and R. Kakarala, "Disk-harmonic coefficients for invariant
% pattern recognition," J. Opt. Soc. Am. A 15, 389-401 (1998)

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Version 1.0
% Author: Johannes Soons; NIST; 2006 
% Review: Johannes Soons; NIST; 2006 
% Status: Need to modify statistical output for weights
% ---------------------------------------------------------

    % set default values
    if (nargin < 7), snoise = []; end
    if (nargin < 6), salgor = []; end
    if (nargin < 5), weights = []; end

    if isempty(salgor), salgor = 'normal'; end

    % solve according to desired routine
    switch lower(salgor)
      case 'qr'
        [vpar,rmse,spar,tpar,ppar,mcov] = fdh_estim_qr(map,mr,ma,nm,weights,snoise);
      case 'normal'
        [vpar,rmse,spar,tpar,ppar,mcov] = fdh_estim_normal(map,mr,ma,nm,weights,snoise);
      case 'orthogonal'
        [vpar,rmse,spar,tpar,ppar,mcov] = fdh_estim_ortho(map,mr,ma,nm,weights,snoise);
      otherwise
        error([salgor,' is not a valid algorithm option']);
    end

end


%------------------------------------------------------------------------------------

function [vpar,rmse,spar,tpar,ppar,mcov] = fdh_estim_qr(map,mr,ma,mnm,weights,snoise)
%
% estimates disk harmonic decomposition with QR decomposition
%

    % determine relevant size parameters and nature of the problem
    [ni,nj] = size(mr);
    nm = numel(map)/(ni*nj);

    if isempty(weights) 
        weights = ones(ni*nj,1); 
    end
    nw = numel(weights)/(ni*nj);
    
    % limit evaluation to non NaN values
    vi = find((mr <= 1) & ~isnan(mr));

    % calculate regression matrix
    np = length(mnm(:,1));
    mX = zeros(length(vi),np);
  
    for i = 1:np
        mX(:,i) = fdh_eval(mr(vi),ma(vi),mnm(i,:),1);
    end

    % solve system sequentially for each case to avoid memory problems
    vpar = zeros(np,nm);
    no   = zeros(nm,1);
    mcov = zeros(np,np,nm);
  
    for i = 1:nm

        % consider only non NaN values
        vii = find(~isnan(map((i-1)*ni*nj+vi)));
        vim = (i-1)*ni*nj+vi(vii);
        viw = (i-1)*ni*nj*(nw > 1)+vi(vii);
        no(i) = length(vim);
        
        % estimate parameters and key statistics (note: covariance matrix has to
        % be recalculated due to possible variations in weights and NaN observations
        mapw = map(vim).*(weights(viw).^(1/2));
        mXw = gtimes(mX(vii,:),weights(viw).^(1/2));
        
        vpar(:,i) = mXw\mapw;
        mcov(:,:,i) = inv(mXw'*mXw);
        sse(i,1) = mapw'*mapw-vpar(:,i)'*mXw'*mapw;        
        mse(i,1) = sse(i,1)/(no(i)-np);
    end

    % calculate statistical properties
    rmse = mse.^(1/2);
    
    if isempty(snoise)
        snoise = rmse;
    else
        if length(snoise) < nm
            snoise = repmat(snoise,nm,1);
        end
    end
  
    spar = zeros(np,nm);
    tpar = zeros(np,nm);
    ppar = zeros(np,nm);
    
    for i = 1:nm
        rmse(i) = sqrt(mse(i)); 
        spar(:,i) = snoise(i)*(diag(mcov(:,:,i))).^(1/2);
        vi = spar(:,i) ~= 0;
        tpar(vi,i) = abs(vpar(vi,i)./spar(vi,i));
        tpar(~vi,i) = inf;    
        
        for j = 1:np
            ppar(j,i) = tcdf(tpar(j,i),no(i)-np);
        end
    end
end


%------------------------------------------------------------------------------------

function [vpar,rmse,spar,tpar,ppar,mcov] = fdh_estim_ortho(map,mr,ma,mnm,weights,snoise)
%
% fdh_estim_ortho: estimates FDH coefficients assuming orthogonality of the
%                  FDH terms.
%
% note: don't worry about scaling and centering the regression matrix
%       columns as the FDH terms are already scaled

    % determine relevant size parameters and nature of the problem
    [ni,nj] = size(mr);
    nm = numel(map)/(ni*nj);

    if isempty(weights) 
        weights = ones(ni*nj,1); 
    end
    nw = numel(weights)/(ni*nj);

    np = length(mnm(:,1));

    % limit evaluation to non NaN values
    vi = find((mr <= 1) & ~isnan(mr));

    % initialize
    vpar = zeros(np,nm);
    mcov = zeros(np,np,nm);
    sse = zeros(nm,1);
    mse = zeros(nm,1);
    no  = zeros(nm,1);

    for ip = 1:np

        % get FDH term
        vdes = fdh_eval(mr(vi),ma(vi),mnm(ip,:),1);

        for i = 1:nm

            % consider only non NaN values
            vii = find(~isnan(map((i-1)*ni*nj+vi)));
            vim = (i-1)*ni*nj+vi(vii);
            viw = (i-1)*ni*nj*(nw > 1)+vi(vii);

            % estimate parameters and key statistics (note: covariance matrix has to
            % be recalculated due to possible variations in weights and NaN observations
            X = vdes(vii)'*(vdes(vii).*weights(viw));
            Y = vdes(vii)'*(map(vim).*weights(viw));  
            
            vpar(ip,i) = Y/X;
            mcov(ip,ip,i) = X.^(-1);          
          
            if ip == 1
                sse(i) = sse(i)+map(vim)'*(map(vim).*weights(viw));
                no(i) = length(vim);
            end
          
            sse(i) = sse(i) - 2*vpar(ip,i)*Y+vpar(ip,i)^2*X;
        end
         
        mse = sse./(no-np);
    end

    % calculate statistical properties
    rmse = mse.^(1/2);
  
    if isempty(snoise)
        snoise = rmse;
    else
        if length(snoise) < nm
            snoise = repmat(snoise,nm,1);
        end
    end
  
    spar = zeros(np,nm);
    tpar = zeros(np,nm);
    ppar = zeros(np,nm);
  
    for i = 1:nm
        spar(:,i) = snoise(i)*(diag(mcov(:,:,i))).^(1/2);
        vi = spar(:,i) ~= 0;
        tpar(vi,i) = abs(vpar(vi,i)./spar(vi,i));
        tpar(~vi,i) = inf;    

        for j = 1:np
            ppar(j,i) = tcdf(tpar(j,i),no(i)-np);
        end
    end
end


%------------------------------------------------------------------------------------

function [vpar,rmse,spar,tpar,ppar,mcov] = fdh_estim_normal(map,mr,ma,mnm,weights,snoise)
%
% Estimates FDH coefficients by explicit calculation
% of the normal equations before solving them. This procedure
% is numerically inferior to the QR decomposition of the 
% regression matrix.
%
% note: don't worry about scaling and centering the regression matrix
%       columns as the FDH terms are already scaled

    % set maximum number of elements in the regression matrix to avoid memory problems
    nmax = 1E7;

    % initialize
    [ni,nj] = size(mr);
    nm = numel(map)/(ni*nj);

    if isempty(weights) 
        weights = ones(ni*nj,1); 
    end
    nw = numel(weights)/(ni*nj);

    % limit evaluation to non NaN values
    vi = find((mr <= 1) & ~isnan(mr));
    no = length(vi);

    % calculate regression matrix
    np = length(mnm(:,1));
  
    mX = zeros(np,np,nm);
    vY = zeros(np,nm);

    % estimate number of sub problems
    ns = ceil((no*np)/nmax);
    
    for is = 1:ns

        % set observations in this patch
        vis = (is:ns:no)';

        % calculate regression matrix
        mdes = zeros(length(vis),np);
        for ip = 1:np
            mdes(:,ip) = fdh_eval(mr(vi(vis)),ma(vi(vis)),mnm(ip,:),1);
        end

        for i = 1:nm

            % consider only non NaN values
            vii = find(~isnan(map((i-1)*ni*nj+vi(vis))));
            vim = (i-1)*ni*nj+vi(vis(vii));
            viw = (i-1)*ni*nj*(nw > 1)+vi(vis(vii));

            % Calculate contribution to system (note: covariance matrix has to be recalculated
            % due to possible variations in weights and NaN observations
            mX(:,:,i) = mX(:,:,i) + ...
                        bsxfun(@times,mdes(vii,:),weights(viw).^(1/2))' * ...
                        bsxfun(@times,mdes(vii,:),weights(viw).^(1/2));
            
            vY(:,i) = vY(:,i) + ...
                      bsxfun(@times,mdes(vii,:),weights(viw).^(1/2))' * ...
                      bsxfun(@times,map(vim),weights(viw).^(1/2));
        end
    end

    vpar = zeros(np,nm);
    no = zeros(nm,1);
    sse = zeros(nm,1);
    mse = zeros(nm,1);

    for i = 1:nm

        % solve normal equations
        vpar(:,i) = mX(:,:,i)\vY(:,i);

        % return statistics
        vii = find(~isnan(map((i-1)*ni*nj+vi)));
        vim = (i-1)*ni*nj+vi(vii);
        viw = (i-1)*ni*nj*(nw > 1)+vi(vii);
        no(i) = length(vim);
    
        mcov(:,:,i) = inv(mX(:,:,i));
        sse(i) = map(vim)'*(map(vim).*weights(viw).^(1/2))-...
                 2*vpar(:,i)'*vY(:,i)+vpar(:,i)'*mX(:,:,i)*vpar(:,i);
        mse(i) = sse(i)/(length(vim)-np);
    end

    % calculate statistical properties
    rmse = mse.^(1/2);

    if isempty(snoise)
        snoise = rmse;
    else
        if length(snoise) < nm
            snoise = repmat(snoise,nm,1);
        end
    end
  
    spar = zeros(np,nm);
    tpar = zeros(np,nm);
    ppar = zeros(np,nm);

    for i = 1:nm
        spar(:,i) = snoise(i)*(diag(mcov(:,:,i))).^(1/2);
        vi = spar(:,i) ~= 0;
        tpar(vi,i) = abs(vpar(vi,i)./spar(vi,i));
        tpar(~vi,i) = inf;    

        for j = 1:np
            ppar(j,i) = tcdf(tpar(j,i),no(i)-np);
        end
    end

end
