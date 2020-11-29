function [eixy,pval,kidx] = isophase_line_xing(fphase,par,xa,xb,vphase,tol,grid)
%function [eixy,pval,kidx] = isophase_line_xing(fphase,par,xa,xb,vphase,tol,grid)
%
% Finds locations on a straight line, usually a tile edge, where the phase
% function 'fphase' assumes the values 2*pi*k + vphase(m), k = ...-2,-1,0,1,2,...
%
% INPUT
% fphase :   function handle of phase function, fphase(x,y,par)
% par :      phase function parameters
% xa :       coordinate pair xa = [x,y] at start of edge
% xb :       last point of the edge
% vphase :   phase offset values
% tol :      tolerances (see: phase2cgh_tile.m)
% grid :     number of samples on edge for initial evaluation
%
% OUTPUT
% eixy :     a nx2 array with coordinates of intersections, one per row
% pval :     a vector with phase values for each of the intersections in 'pxy'.
% kidx :     a vector with phase indices k for each of the intersections

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
% Initial version: Ulf Griesmann, NIST, September 2018
    
    % initialize output arrays
    eixy = pval = kidx = [];

    % set the parameters for fzero
    opts = optimset('TolX',tol.vertex, 'MaxIter',tol.maxit);
    
    % sample the phase function on the line
    X = linspace(xa, xb, grid)';          % one coordinate pair per row
    P = fphase(X(:,1),X(:,2),par);        % phase values on edge
    Pmax = max(P); Pmin = min(P);         % phase extrema

    % locatate any phase discontinuities on the line
    didx = find_discon_grid(P,pi);

    % find k(m) such that 2*pi*k + vphase(m) > Pmin & < Pmax on the edge
    kmin = ceil (0.5 * (Pmin - vphase)/pi);
    kmax = floor(0.5 * (Pmax - vphase)/pi);

    % loop over phase offset values and find all isophase curve intersections
    for m = 1:numel(vphase)
        
        % check if the line has any intersections with isophase curves
        if (2*pi*kmax(m) + vphase(m) > Pmax) & ...
           (2*pi*kmin(m) + vphase(m) < Pmin)
            continue
        end
        
        % unique phase values at the isophase-line intersections,
        K = [kmin(m):kmax(m)];          % k at intersections
        phtarget = 2*pi*K + vphase(m);  % phase values at intersections
        
        % find isophase-line intersection brackets for all isophases
        % i.e. zero-crossings of P - (2*pi*k + vphase)
        PK = repmat(P,1,numel(K)) - phtarget;
        [zci,~] = find( diff(sign(PK),1,1) ); % zci is index into vectors X,P
        x1 = X(zci,:);                        % x1,x2 bracket isophase locations
        x2 = X(zci+1,:);
        
        % several intersections can have the same phase value. 
        % Determine the target phase for all intersections
        if size(x1,1) > length(phtarget)
            P1 = fphase(x1(:,1),x1(:,2),par);    % phase brackets
            P2 = fphase(x2(:,1),x2(:,2),par);
            Pa = min([P1,P2],[],2);
            Pb = max([P1,P2],[],2);
            PT = repmat(phtarget,size(x1,1),1);  % phase value table
            KT = repmat(K,size(x1,1),1);         % corresponding k values
            tidx = phtarget>=Pa & phtarget<=Pb;
            phtarget = PT(tidx)';
            K = KT(tidx)';
        end

        % remove any brackets that include a phase discontinuity
        if ~isempty(didx)
            nod = true(size(x1,1),1);
            for k=1:length(didx)
                nod &= any(abs(x1-X(didx(k),:)) > 5*eps, 2);
            end
            x1 = x1(nod,:);
            x2 = x2(nod,:);
            phtarget = phtarget(nod);
            K = K(nod);
        end

        % calculate precise isophase-edge intersections
        xy = vfmatch(fphase,par,x1,x2,phtarget',opts);
        
        % return results
        eixy = vertcat(eixy,xy);         % locations of intersections on the line
        pval = vertcat(pval,phtarget');  % phase values at intersections
        kidx = vertcat(kidx,K');         % phase index values (k) at intersections
        
    end
    
end
