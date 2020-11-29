function [ciso,eiso] = isophase_curves(fphase,par,tile,eixy,pval,kval,nedg,tol,verbose=false)
%function [ciso,eiso] = isophase_curves(fphase,par,tile,eixy,pval,kval,nedg,tol,verbose=false)
%
% Finds approximations to isophase curves by polygons that have vertices on the 
% isophase and either connect two points with equal phase values on the edges
% of a tile, or trace an isophase curve to a phase discontinuity.
%
% INPUT
% fphase :   function handle of phase function, fphase(x,y,par)
% par :      phase function parameters
% tile :     polygon circumscribing a rectangular tile area
% eixy :     a nx2 matrix with edge intersection coordinates, one per row
% pval :     vector of phase values at the tile edge intersections
% kval :     vector with corresponding phase indices k at the intersections
% nedg :     vector with edge number (1..4) for each of the intersections
% tol:       structure containing permissible error bounds. See 'cghparset'.
% verbose :  (optional) a logical. Display progress information when set to 'true'.
%            Primarily useful for debugging. Default is 'false'.
%
% OUTPUT
% ciso :     cell array with polygon approximations to isophase lines
% eiso :     nx2 matrix with pairs of edge numbers at the start and at the
%            end of polygons in 'ciso'. One row per polygon.
    
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
%
% Initial version: Ulf Griesmann, NIST, July 2019
% ---------------------------------------------------------

    if verbose
        fprintf('\n');
    end
    
    % number of polygon vertices for pre-allocation
    numvert = 1000;
    
    % some input checking
    if nargin < 9
        error('function ''isophase_curves'' must have at least 9 input arguments.');
    end

    % set the parameters for fzero
    opts = optimset('TolX',tol.vertex, 'MaxIter',tol.maxit);
 
    % initialize output
    ciso = {};
    eiso = [];
    cidx = 0;  % number of elements in arrays [cpke]iso

    % nothing to do in the case of a tile without edge intersections
    if isempty(eixy)
        return
    end

    % indices of edge intersections that have been connected by
    % isophase lines
    cei = false(size(pval));          % Connected Edge Intersections
    
    % determine the initial vertices of the polygons to trace
    aidx = isophase_vertex_init(pval,kval);
    upval = pval(aidx);    
    if verbose
        fprintf('   Calculating %d isophase curves...\n', length(upval));
    end
    
    % trace isophases accross the tile
    npol = length(upval);             % number of polygons
    nxy = cell(npol,1);
    nxy(1:npol) = zeros(numvert,2);   % preallocate isophase polygons
    isac = true(npol,1);              % logical array to flag active polygons
    
    % copy the points of departure from the list of all edge intersections
    curr_xy = eixy(aidx,:);

    % mark the starting points in 'cei', the master list of edge intersections
    cei(aidx) = true;
        
    % the initial edge intesections are the polygon starting vertices
    for k=1:npol
        nxy{k}(1,:) = curr_xy(k,:);
    end
    nvtx = ones(npol,1);              % number of vertices in new polygons
    
    % pre-allocate variables for intermediate vertices
    next_xy = pred_xy = zeros(length(upval),2);
        
    % assemble the set of isophase polygonz
    nstp = 0; % step count
    
    while true % steps along the isophase curves

        nstp += 1;

        % gradient and isophase curvatures at current location
        [Px,Py,~,~,~,isocur] = vphderiv(fphase,par,curr_xy(isac,:),tol.hderiv);

        % local coordinate systems @ curr_xy: unit tangent T, unit normal N
        [T,N] = normtang(Px,Py);

        % at the initial edges of the tile decide on the correct sign of
        % the tangent T for following the isophase INTO the tile area.
        if nstp == 1 
            direct = isodirect(curr_xy,T,N,tile,tol.vertex);
        end
        
        % predict the locations of the next set of vertices on the isophases
        % using the isophase osculating circles at the current locations
        [pred_xy(isac,:),B] = ...
            predict_xy(curr_xy(isac,:),isocur,T,N,direct(isac),tol);
        
        % calculate the exact locations of the newly predicted vertices
        [x1,x2] = find_bracket(fphase,par,pred_xy(isac,:),upval(isac),B,tol);
        next_xy(isac,:) = vfmatch(fphase,par,x1,x2,upval(isac),opts);
        
        % which new vertices are within the tile area ?
        iidx = isinpolygon(tile,next_xy(isac,:));  % inside index
        ia = find(isac); % for index lookup
        
        % store any new predicted vertices that are outside the tile
        % and turn off the tracing of those isophases
        if any(~iidx)
            for k=1:length(iidx)
                if ~iidx(k) && isac(ia(k))
                    nvtx(ia(k)) += 1;
                    nxy{ia(k)}(nvtx(ia(k)),:) = next_xy(ia(k),:);
                    nxy{ia(k)} = nxy{ia(k)}(1:nvtx(ia(k)),:); % truncate polygon
                    isac(ia(k)) = false;                      % turn off isophase tracing
                end
            end
        end

        % done isophase tracing when all polygons terminate outside the tile area.
        if all(~isac)
            break
        end

        % add the new vertices to the polygons that are still within the tile
        for k=1:length(iidx)
            if iidx(k)
                nvtx(ia(k)) += 1;
                nxy{ia(k)}(nvtx(ia(k)),:) = next_xy(ia(k),:);
            end
        end
        
        % go look for the next set of vertices on the isophases
        curr_xy = next_xy;
        
        % make some noise
        if verbose % print number of steps / active isophases
            if ~mod(nstp,20)
                fprintf('%d/%d..',nstp,sum(isac)); 
                fflush(stdout());
            end
        end
        
    end % while ...
    
    % At this point all isophase lines have been traced accross the tile and
    % the polygons terminate outside the tile. Look up the edge intersections 
    ins_xy = out_xy = zeros(npol,2);  % last points in polygons
    for k=1:npol
        ins_xy(k,:) = nxy{k}(end-1,:);
        out_xy(k,:) = nxy{k}(end,:);
    end
    [edge_xy,edge_num,eidx] = edge_intersect(ins_xy,out_xy,tile,eixy,nedg,tol);
    
    % the edge intersections replace the terminal polygon vertices
    for k = 1:npol
        nxy{k}(nvtx(k),:) = edge_xy(k,:);
    end
    
    % also mark them in the master list 'cei'
    cei(eidx) = true;
    
    % return traced polygons
    ciso = nxy;
    eiso = [nedg(aidx),edge_num]; % initial edge,terminal edge

    if verbose
        fprintf('\n     %d isophase contours completed.\n', numel(ciso));
        fflush(stdout());
    end

    % Now check if any unconnected intersections remain. This
    % should never happen.
    if any(~cei)
        uei = find(~cei);
        fprintf('     %d unconnected edge intersections\n', length(uei));
        error('this should never happen - likely a bug');
    end

end
