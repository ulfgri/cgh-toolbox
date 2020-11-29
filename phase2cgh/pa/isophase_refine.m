function [ciso,piso] = isophase_refine(fphase,par,aiso,apha,tol,brac,verbose=false)
%function [ciso,piso] = isophase_refine(fphase,par,aiso,aph,tol,brac,verbose=false)
%
% Adds additional points to the isophases between the vertex points
% found in the initial pilot approximation and then refines all points
% using a root finding algorithm to obtain a polygon approximation of
% the isophases that satisfy the tolerance requirements.
%
% INPUT
% fphase :  function handle of phase function fphase(x,y,par)
% par :     structure with phase function parameters
% aiso :    cell array with approximate phase function polynomials
% apha :    phase values of the approximate polynomials
% tol :     tolerances for isophase polynomials
% brac :    size of bracket for the refinement of approximate isophase points
% verbose:  (optional) a logical; display progress information. Default is 'false';
%
% OUTPUT
% ciso :    cell array with isophase contours
% piso :    vector with corresponding phase values
    
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Short version: This software is in the PUBLIC DOMAIN.
%
% Initial version: Ulf Griesmann, NIST, July 2019

    % check arguments
    if nargin < 5
        error('function must have at least five arguments.');
    end

    if verbose
        fprintf('   Refining isophase curves .. ');
        fflush(stdout);
    end
    
    % init output
    ciso = cell(numel(aiso),1); 
    piso = zeros(numel(aiso),1);
    
   % set the parameters for vfzero
    opts = optimset('TolX',tol.vertex, 'MaxIter',tol.maxit);    
        
    % loop over isophase curves
    for k = 1:numel(aiso)

        % approximate vertices of isophase
        xy = aiso{k};
        
        % STEP 1: determine the precise locations of the initial vertices
        % ---------------------------------------------------------------
        D = xy(2:end,:) - xy(1:end-1,:);     % vectors between initial vertices
        N = [D(:,2),-D(:,1)]; 
        N = N./vecnorm(N,2,2);               % unit normal vectors at xy(n,:)
        N(end+1,:) = N(end,:);               % use same for last vertex
        
        X1 = xy - 0.5*brac*N;                % isophase search brackets
        X2 = xy + 0.5*brac*N;
        
        xyi = vfmatch(fphase,par,X1,X2,repmat(apha(k),size(X1,1),1),opts);

        % STEP 2: determine where additional vertices are needed
        % ------------------------------------------------------
        D = xyi(2:end,:) - xyi(1:end-1,:);   % vectors between refined vertices
        d = vecnorm(D,2,2);                  % distances between vertices
        mxy = xyi(1:end-1,:) + 0.5 * D;      % mid-points between vertices
        
        % Estimate the isophase sag at the mid-points between
        % vertices and the additional number of vertices needed
        [~,~,~,~,~,isocur] = vphderiv(fphase,par,mxy,tol.hderiv);
        R = abs(1./(isocur+eps));            % isophase curvature @ midpoints
        sag = R - 0.5 * sqrt(4*R.^2 - d.^2); % isophase curve sags @ midpoints
        S = 2 * sqrt(tol.segdev.*(2*R - tol.segdev)); % allowable secant length
        n = find(S >= d);                    % no additional vertices needed here,
        S(n) = d(n);                         % so keep existing distances
        nadd = ceil(d./S) - 1;               % number of additional vertices

        % STEP 3: add additional vertices through linear interpolation, calculate
        % bracketing points for root finding, and find the vertices
        % -------------------------------------------------------------------
        if any(nadd)

            N = [D(:,2),-D(:,1)];
            N = N./vecnorm(N,2,2);           % unit normal vectors at xy(n,:)

            ia = find(nadd);                 % additional vertices needed here
            ia = repelem(ia,nadd(ia));       % create additional vertices ...
            ii = cell2mat(arrayfun(@(x)[1:x]',nadd,'UniformOutput',false));
            xadd = xyi(ia,:) + ii.* D(ia,:)./(nadd(ia)+1);
            
            X1 = xadd - 2*sag(ia).*N(ia,:);  % brackets, +/- sag should
            X2 = xadd + 2*sag(ia).*N(ia,:);  % now suffice ... use +/-  2*sag
        
            % find the additional isophase vertices within tol.vertex
            xya = vfmatch(fphase,par,X1,X2,repmat(apha(k),size(X1,1),1),opts);
            
            % splice the additional vertices inbetween the initial ones
            xy = zeros(size(xyi,1)+size(xya,1),2);
            ixi = 1 + [0:length(nadd)]' + [0;cumsum(nadd)];
            xy(ixi,:) = xyi;
            ixa = true(size(xy,1),1); 
            ixa(ixi) = false;
            xy(ixa,:) = xya;
            
        end
           
        ciso{k} = xy;
        piso(k) = apha(k);  % just copy them

        % make noise
        if verbose && ~mod(k,25)
            fprintf('%d .. ', k); 
            fflush(stdout);
        end
    
    end % for
              
    if verbose
        fprintf('\n   %d isophase curves refined.\n', numel(ciso));
    end

end
