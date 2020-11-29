function [clciso,eiso] = isophase_clip(fphase,par,ciso,piso,tile,tol,brac,verbose=false)
%function [clciso,eiso] = isophase_clip(fphase,par,ciso,piso,tile,tol,brac,verbose=false)
%
% Isophase curves that extend beyond the borders of a tile are clipped 
% to the tile boundary.
%
% INPUT
% fphase :   function handle of phase function fphase(x,y,par)
% par :      structure with phase function parameters
% ciso :     cell array with Mx2 matrices of isophase curve vertices
% piso :     vector with corresponding phase values
% tile :     a closed polygon circumscribing the tile area
% tol :      structure with tolerances. See 'cghparset' for details
% brac :     size (range) of bracket for the refinement of isophase points
% verbose:   (optional) a logical; display progress information. Default is 'false';
%
% OUTPUT
% clciso :   cell array with isophase curves that begin and end on tile edges
% eiso :     nx2 matrix with pairs of edge numbers at the start and at the
%            end of polygons in 'cliso'. One row per polygon.
    
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Short version: This software is in the Public Domain.
%
% Initial version: Ulf Griesmann, NIST, July 2019

    if verbose
        fprintf('   clip approximate isophase curves to tile edges.\n');
    end
    
    % set the parameters for vfzero
    opts = optimset('TolX',tol.vertex, 'MaxIter',tol.maxit);
    
    % clip the isophase curves outside the tile boundary and 
    % determine on which side the isophases cross the tile
    [clciso,pidx,ned] = clip_to_tile(tile, ciso);
    piso = piso(pidx);

    % collect tile edge intersections for refinement
    tei = zeros(2*numel(clciso),2);  % isophase - tile edge intersections
    for k=1:numel(clciso)
        tei(2*k-1,:) = clciso{k}(1,:);   % first vertex of polygon
        tei(2*k  ,:) = clciso{k}(end,:); % last vertex of polygon
    end
    
    % isophase search direction vectors along tile boundaries
    S = [1,0; 0,1; -1,0; 0,-1];
    
    % build tables of bracketing points and target phase values
    Pt = reshape(repmat(piso,1,2)',2*length(piso),1); % [start,end,start,end,...]
    ned = reshape(ned', numel(ned),1);                % same shape ...
    X1 = X2 = zeros(size(tei));
    for k=1:size(tei,1)
        X1(k,:) = tei(k,:) - 0.5*brac*S(ned(k),:);
        X2(k,:) = tei(k,:) + 0.5*brac*S(ned(k),:);
    end

    % refine the intersections on the tile boundaries
    rtei = vfmatch(fphase,par,X1,X2,Pt,opts);
    
    % Vertices near a corner may be outside the tile after refinement
    rned = onedges(tile, rtei, tol.vertex);
    if any(~rned) % found at least one
        
        idx = find(rned==0);
        tc = tile(1,:) + 0.5*[tile(2,1)-tile(1,1), tile(3,2)-tile(2,2)]; % tile center
        X1 = X2 = Px = [];
        
        for k = idx'
            
            % find the correct edges
            switch ned(k)
              case 1 
                if rtei(k,1) < tc(1) % find tile quadrant
                    edg = 4;
                else
                    edg = 2;
                end
              case 2
                if rtei(k,2) < tc(2)
                    edg = 1;
                else
                    edg = 3;
                end
              case 3
                if rtei(k,1) > tc(1)
                    edg = 2;
                else
                    edg = 4;
                end
              case 4
                if rtei(k,2) > tc(2)
                    edg = 3;
                else
                    edg = 4;
                end
 
            end
               
            % new brackets for wayward corner vertices
            X1(end+1,:) = tei(k,:) - 0.5*brac*S(edg,:);
            X2(end+1,:) = tei(k,:) + 0.5*brac*S(edg,:);
            Px(end+1) = Pt(k);            
            rned(k) = edg; % correct the edge number
        end

        % re-optimize vertices on the correct edges
        rteix = vfmatch(fphase,par,X1,X2,Px,opts);

        % update any re-optimized vertices
        for k = 1:length(Px)
            rtei(idx(k),:) = rteix(k,:);
        end
    end
    
    % update the first and last vertex of each polygon
    for k=1:numel(clciso)
        clciso{k}(1,:)   = rtei(2*k-1,:); % update first vertex
        clciso{k}(end,:) = rtei(2*k,  :); % update last vertex
    end
    eiso = [rned(1:2:end-1),rned(2:2:end)];

    % orient polygons so that they begin at the intersection that is
    % closest to the tile origin; collect polygon start points and edge 
    % numbers
    for k=1:numel(clciso)
        if edge_point_order(eiso(k,1),clciso{k}(1,:),eiso(k,2),clciso{k}(end,:)) < 0
            clciso{k} = polyrev(clciso{k});     % reverse polygon orientation
            eiso(k,:) = eiso(k,[2,1]);          % swap edge numbers
        end
    end
    
    if verbose
        fprintf('   %d isophase curves clipped\n', numel(clciso));
    end

end


%------------------------------------------------------------------

function [clpo,pidx,nedg] = clip_to_tile(tile,po)
%
% Trims polygons that cross a tile to the edges of the tile
%
% INPUT
% tile : a 5x2 matrix with a polygon circumscribing a square tile
%        area
% po :   a cell array of polygons, i.e. nx2 matrices with vertices,
%        one per row. The polygons are sorted such that they run 
%        from the edge with the smaller edge number to the edge
%        with the larger one.
%
% OUTPUT
% clpo : a cell array of polygons that end on the tile edge
% pidx : a vector that contains for each clpo{k} the index of the 
%        original polygon in cell array po. clpo can contain fewer 
%        polygons than po.
% nedg : a nx2 matrix of tile edges at which the intersections occur
%
% NOTE: this function does not implement a general purpose polyline
% clipping algorithm because it is not required here. We know that
% all polygons begin and end outside the tile area, and that they
% have two intersections with the tile, except when the
% intersections coincide with tile corner or the polygon is wholly
% outside the tile.

    clpo = {};
    pidx = nedg = [];
    idx = 1;
    
    % clip polygons
    for k = 1:numel(po)
        
        % calculate tile-polygon intersections
        [poi,nseg] = polyinters(tile,po{k});
        
        % polygons with fewer than two intersections are outside
        % the tile
        if isempty(poi) || size(poi,1) < 2
            continue
        end

        % nseg(:,2) contains the intersections of po{k} with the tile edges
        [~,ixs] = sort(nseg(:,2));  % sort them
        poi = poi(ixs,:);
        nseg = nseg(ixs,:);
        
        % create clipped polygon(s)
        for m = 1:2:size(nseg,1)
            clpo{idx} = vertcat(poi(m,:),po{k}(nseg(m,2)+1:nseg(m+1,2),:),poi(m+1,:));
            nedg(idx,:) = [nseg(m,1),nseg(m+1,1)]; % tile edge intersections
            pidx(idx) = k;                         % store index of original polygon
            idx += 1;
        end % for m

    end % for k
    
end
