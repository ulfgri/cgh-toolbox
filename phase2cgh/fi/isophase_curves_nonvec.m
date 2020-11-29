function [ciso,piso,kiso,eiso] = isophase_curves_nonvec(fphase,par,tile,eixy,pval,kval,pedg,tol,verbose)
%function [ciso,piso,kiso,eiso] = isophase_curves_nonvec(fphase,par,tile,eixy,pval,kval,pedg,tol,verbose)
%
% finds approximations to isophase curves by polygons that have vertices on the 
% isophase and connect two points with equal phase values on the edges of a tile.
%
% INPUT
% fphase :   function handle of phase function, fphase(x,y,par)
% par :      phase function parameters
% tile :     polygon circumscribing a rectangular tile area
% eixy :     a nx2 matrix with intersection coordinates, one per row
% pval :     vector of phase values at the intersections
% kval :     vector with phase indices k at the intersections
% pedg :     vector with edge number for each of the intersections
% tol:       structure containing permissible error bounds. All tolerances are
%            defined as a FRACTION OF THE LENGTH UNIT used for x and y coordinates (usually
%            um). See 'phase2cgh_tile' for documentation.
% verbose :  (optional) a logical. Display progress information when set to 'true'.
%            Primarily useful for debugging. Default is 'false'.
%
% OUTPUT
% ciso :     cell array with polygon approximations to isophase lines
% piso :     phase values of the isophase lines. piso(k) is phase of cai{k}.
% kiso :     corresponding phase indices for the phase values in 'piso'.
% eiso :     nx2 matrix with pairs of edge numbers at the start and at the
%            end of polygons in 'cai'. One row per polygon.
    
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
% ---------------------------------------------------------
    
    % estimated number of polygon vertices
    numvert = 1000;  % 
    
    % some input checking
    if nargin < 8
        error('function ''isophase'' must have at least eight input arguments.');
    end
    if nargin < 9, verbose = false; end

    % set the parameters for fzero
    opts = optimset('TolX',tol.vertex, 'MaxIter',tol.maxit);
 
    % initialize output
    ciso = {};
    piso = kiso = eiso = [];
    cidx = 0;  % number of elements in arrays [cpke]iso

    % nothing to do in the case of a tile without edge intersections
    if isempty(eixy)
        return
    end
    
    % indices of edge intersections that need to be traced
    aei = true(size(kval));
    aidx = 1;  % number of current available edge intersection
        
    % trace isophases accross the tile
    while true
    
        % pick the next edge intersection that needs connecting
        while ~aei(aidx)
            aidx = aidx + 1;
        end

        % initialize next isophase polygon
        nxy = zeros(numvert,2);   % preallocate new polygon
        curr_xy = eixy(aidx,:);   % point of departure
        nidx = 1;                 % number of vertices in new polygon
        nxy(nidx,:) = curr_xy;    % first point of polygon
            
        % assemble an isophase polygon
        while true
            
            % gradient and isophase curvature at current location
            [Px,Py,~,~,~,isocur] = phderiv(fphase,par,curr_xy,tol.hderiv);

            % local coordinate system @ curr_xy: unit tangent T, unit normal N
            [T,N] = normtang(Px,Py);
                
            % at the initial edge of the tile decide on the correct sign of
            % the tangent T for following the isophase INTO the tile area
            if nidx == 1 
                direct = isodirect(curr_xy,T,N,tile,tol.vertex);
            end
            
            % predict the location of the next vertex on the isophase
            % using the isophase osculating circle at the current location
            R = abs(1 / (isocur + eps));             % radius of osculating circle    
            S2 = tol.segment*(8*R - 4*tol.segment);  % square of secant length
            nu = 0.5*S2/R;                           % coordinates (tau,nu) of next
            tau = sqrt(S2 - nu^2);                   % point in moving frame T,N
            pred_xy = curr_xy + direct*tau*T + sign(isocur)*nu*N;
                                
            % check if the predicted point is still within the tile area.
            % if not, proceed to finding the closest edge intersection
            if ~isinpolygon(tile,pred_xy)
                break
            end

            % calculate the exact location of the vertex
            [x1,x2] = find_bracket(fphase,par,pred_xy,pval(aidx),nu*N,tol.maxbs);
            if isempty(x1)
                error('failed to find a valid bracketing in ''isophase_tile''.');
            end
            next_xy = fmatch(fphase,par,x1,x2,pval(aidx),opts);
            
            % add new vertex to polygon
            nidx = nidx + 1;
            nxy(nidx,:) = next_xy;
                
            % go look for the next vertex
            curr_xy = next_xy;
                
        end
            
        % find out where the isophase leaves the tile
        [edge_xy,edge_num,iidx] = edge_intersect(curr_xy,pred_xy,tile,eixy);
            
        % edge intersection is terminal point in new polygon
        nidx = nidx + 1;
        nxy(nidx,:) = edge_xy;
        nxy = nxy(1:nidx,:);  % truncate to actual number of vertices
        
        % add new polygon to output
        cidx = cidx + 1;
        ciso{cidx} = nxy;
        piso(cidx) = pval(aidx);
        kiso(cidx) = kval(aidx);
        eiso(cidx,:) = [pedg(aidx),edge_num]; % initial edge,terminal edge

        % signal some progress
        if verbose
            if ~mod(cidx,10) 
                fprintf('.'); 
                if is_octave(), fflush(stdout()); end 
            end
            if ~mod(cidx,50)
                fprintf(' %d ', cidx);
                if is_octave(), fflush(stdout()); end 
            end
        end

        % strike start & end points from the edge intersection list
        aei([aidx,iidx]) = false;
        
        if ~any(aei) % all edge intersections are connected by polygons
            break
        end
                
    end
    
    if verbose
        fprintf('\n');
    end
        
end
