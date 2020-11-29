function [ciso,eiso] = boundary_curves(fphase,par,tile,eixy,pval,kval,nedg,tol,sing,verbose=false)
%function [ciso,eiso] = boundary_curves(fphase,par,tile,eixy,pval,kval,nedg,tol,sing,verbose=false)
%
% Finds approximations to phase boundary curves by polygons for phase functions
% with phase discontinuities due to a Hilbert phase term. The isophase polygons have 
% vertices on the isophase, and they connect points with equal phase values on 
% the edges +- 2*pi the number discontinuities traversed.
%
% INPUT
% fphase :   function handle of phase function, fphase(x,y,par)
% par :      phase function parameters
% tile :     polygon circumscribing a rectangular tile area
% eixy :     a nx2 matrix with edge intersection coordinates, one per row
% pval :     vector of phase values at the edge intersections
% kval :     vector with corresponding phase indices k at the intersections
% nedg :     vector with edge number (1..4) for each of the intersections
% tol:       structure containing permissible error bounds. See 'cghparset'.
% sing:      a structure with information about a phase singularity on the tile edge.
% verbose :  (optional) a logical. Display progress information when set to 'true'.
%            Default is 'false'.
%
% OUTPUT
% ciso :     cell array with polygon approximations to phase boundary lines
% eiso :     nx2 matrix with pairs of edge numbers at the start and at the
%            end of polygons in 'ciso'. One row per polygon.
%
% NOTE: the term "isophase" is still used for phase boundaries even though "isophase" 
% curves no longer have a constant phase value because the phase changes whenever
% a curve traverses a phase discontinuity.
    
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
% Initial version: Ulf Griesmann, NIST, September 2019
% ---------------------------------------------------------

    if verbose
        fprintf('\n');
    end
    
    % number of polygon vertices for pre-allocation
    numvert = 1000;
    
    % some input checking
    if nargin < 9
        error('function ''boundary_curves'' must have at least 9 input arguments.');
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

    % indices of isophase-edge intersections that have been connected by
    % isophase lines
    cei = false(size(pval));          % Connected Edge Intersections

    % Boundary curve following starts again at edges that have a subset of unique
    % phase values. Unlike in the case of phase functions without discontinuities
    % this will in general be a smaller subset of phase boundary curves because
    % they do not have a constant phase values. This results in a larger number
    % of duplicate phase values at the tile edge intersections.
    
    % Start tracing boundary curves from a sub-set of edge intersections
    % 'aidx' are indices of actively traced polygons in 'cei'
    aidx = boundary_vertex_init(pval,kval,eixy,nedg,tile,sing,tol);
    upval = pval(aidx);
    if verbose
        fprintf('   Calculating %d phase boundary curves...\n', length(upval));
    end
        
    % trace boundary curves accross the tile
    npol = length(upval);             % number of polygons
    nxy = cell(npol,1);
    nxy(1:npol) = zeros(numvert,2);   % preallocate polygons
    isac = true(npol,1);              % logical array to flag active polygons

    % vector to keep track of boundary curves that terminate at the singularity
    tats = zeros(npol,1);

    % copy the points of departure from the list of all edge intersections
    curr_xy = eixy(aidx,:);

    % mark the starting points in 'cei', the master list of edge intersections
    cei(aidx) = true;
        
    % the initial edge intesections are the polygon starting vertices
    for k=1:npol
        nxy{k}(1,:) = curr_xy(k,:);
    end
    nvtx = ones(npol,1);              % number of vertices in new polygons
    
    % pre-allocate variables for intermediate vertices and phase values
    next_xy = pred_xy = zeros(length(upval),2);
    next_ph = zeros(length(upval),1);
    
    % assemble the set of phase boundary polygonz
    nstp = 0; % step count
    
    while true % steps along the curves

        nstp += 1;

        % gradient and isoline curvatures at current locations
        [Px,Py,~,~,~,isocur] = vwphderiv(fphase,par,curr_xy(isac,:),tol.hderiv);

        % local coordinate systems @ curr_xy: unit tangent T, unit normal N
        [T,N] = normtang(Px,Py);
        
        % at the initial edges of the tile decide on the correct sign of
        % the tangent T for following the isophase INTO the tile area.
        if nstp == 1 
            direct = isodirect(curr_xy,T,N,tile,tol.vertex);
        end
        
        % predict the locations of the next set of vertices on the isophases
        % ------------------------------------------------------------------
        [pred_xy(isac,:),B] = ...
            predict_xy(curr_xy(isac,:),isocur,T,N,direct(isac),tol);

        % check if we are (too) close to the singularity
        % ---------------------------------------------------------------
        isd = [];
        
        if ~isempty(sing)

            % check if predicted vertex is close to singularity
            odis = 1e10*ones(size(isac));
            ldis = zeros(size(isac));
            [odis(isac),ldis(isac)] = point_line_dist(curr_xy(isac,:),pred_xy(isac,:),sing.pos);
            idx_near = odis < sing.era;
            if any(idx_near) % stay clear of singularity to have room for bracketing
                isd = find(idx_near & (ldis(idx_near) < tol.seglen));
            end
        
            % check if any polygons are progressing into a singularity
            idx_near = vecnorm(pred_xy(isac,:)-sing.pos,2,2) < sing.era;
            if any(idx_near)
                isd = find(idx_near);
            end
        
            % terminate any boundary lines that have reached the singularity
            if ~isempty(isd)
                for k=1:length(isd)
                    nvtx(isd(k)) += 1;
                    nxy{isd(k)}(nvtx(isd(k)),:) = sing.pos;      % last vertex is singularity
                    nxy{isd(k)} = nxy{isd(k)}(1:nvtx(isd(k)),:); % truncate polygon
                    tats(isd(k)) = onedges(tile,sing.pos,tol.vertex);
                    isac(isd(k)) = false;                        % turn off isophase tracing
                end
            
                % done boundary tracing when all polygons have terminated
                if all(~isac)
                    break
                end
            
                % predicted vertices for revised set of active boundary curves
                [Px,Py,~,~,~,isocur] = vwphderiv(fphase,par,curr_xy(isac,:),tol.hderiv);
                [T,N] = normtang(Px,Py);
                [pred_xy(isac,:),B] = ...
                    predict_xy(curr_xy(isac,:),isocur,T,N,direct(isac),tol);
            end
        end
        
        % check if any phase boundaries are crossing a discontinuity
        % ----------------------------------------------------------
        curr_ph = pred_ph = zeros(size(isac));
        curr_ph(isac) = fphase(curr_xy(isac,:)(:,1),curr_xy(isac,:)(:,2),par);
        pred_ph(isac) = fphase(pred_xy(isac,:)(:,1),pred_xy(isac,:)(:,2),par);
        indi = abs(curr_ph - pred_ph) > pi;  % check if phase value changes
        if any(indi) % adjust boundary phases by multiples of 2*pi
            upval(indi) += 2*pi * round( 0.5*(pred_ph(indi) - curr_ph(indi))/pi );
        end
        
        % calculate the exact locations of the predicted vertices
        % -----------------------------------------------------------------        
        [x1,x2] = find_bracket(fphase,par,pred_xy(isac,:),upval(isac),B,tol);
        next_xy(isac,:) = vfmatch(fphase,par,x1,x2,upval(isac),opts);
        
        % which new vertices are within the tile area ?
        iidx = isinpolygon(tile,next_xy(isac,:));  % inside index
        ia = find(isac);
        
        % store any new vertices that are outside the tile and turn off the
        % tracing of those boundaries
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

        % done boundary tracing when all polygons have terminated
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
    
    % At this point all boundary lines have been traced accross the tile.
    % Look up the edge intersections of the phase boundary curves.
    ins_xy = out_xy = zeros(sum(~tats),2);  % last two points in polygons
    m = 1;
    for k=1:npol
        if ~tats(k)
            ins_xy(m,:) = nxy{k}(end-1,:);
            out_xy(m,:) = nxy{k}(end,:);
            m += 1;
        end
    end
    [edge_xy,toe_num,eidx] = edge_intersect(ins_xy,out_xy,tile,eixy,nedg,tol);

    % the edge intersections replace the terminal polygon vertices
    m = 1;
    for k = 1:npol
        if ~tats(k)
            nxy{k}(nvtx(k),:) = edge_xy(m,:);
            m += 1;
        end
    end
    
    % also mark them in the master list 'cei'
    cei(eidx) = true;

    % return traced polygons and edge number data
    ciso = nxy;
    edge_num = tats;              % edge numbers at singularity
    edge_num(~tats) = toe_num;    % terminate on edge number
    eiso = [nedg(aidx),edge_num]; % initial edge,terminal edge

    % this batch is done
    if verbose
        fprintf('\n     %d phase boundary contours completed.\n',numel(ciso));
        fflush(stdout());
    end

    % Now check if there are any remaining unconnected edge intersections.
    % 'boundary_curves' is called recursively on the remaining intersections
    if any(~cei)
    
        % indices of remaining intersections
        ridx = find(~cei);
        [rciso,reiso] = ...
            boundary_curves(fphase,par,tile, ...
                            eixy(ridx,:),pval(ridx),kval(ridx),nedg(ridx), ...
                            tol,sing,verbose);

        % add return values to the output
        ciso = vertcat(ciso,rciso);
        eiso = vertcat(eiso,reiso);

    end

    % remove any degenerate polygons that have a length of two and
    % both vertices are identical. These can occur when the start
    % vertex of a polygon is a tile corner.
    % not sure this is needed any longer (U.G., 24 Nov 2019)
    good = true(numel(ciso),1);
    for k = 1:numel(ciso)
        if size(ciso{k},1) == 2 && ...
           all( abs(ciso{k}(1,:)-ciso{k}(2,:)<5*eps) )
           good(k) = false;
        end
    end
    if ~all(good)
        ciso = ciso(good);
        eiso = eiso(good,:);
    end
    
end

%-- la fin --------------------------------------------------------------------
