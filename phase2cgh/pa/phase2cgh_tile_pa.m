function cap = phase2cgh_tile_pa(fphase,par,llc,urc,vphase,tol,grid,verbose=false)
% function cap = phase2cgh_tile_pa(fphase,par,llc,urc,vphase,tol,grid,verbose=false)
%  
% Calculates the layout of a computer generated hologram (CGH) for a
% given phase function in a rectangular area using a pilot approximation
% approach in which an initial estimate of the isophase lines is made
% with a marching squares algorithm followed by refinement of the
% isophase curves. The CGH layout is described by polygons that
% approximate the true phase contours. The polygons are boundaries
% of areas with phase values within intervals 2*pi * k + vphase,
% where k is an arbitrary but fixed integer for each polygon and vphase
% is a pair of boundary phase values.
%
% NOTE: This function is meant to be called by the 'phase2cgh' function.
%
% INPUT
% fphase:    function handle to user defined function describing the 2D phase
%            distribution (in radians)
% par:       structure with parameters describing the phase distribution
%            (see user defined function describing phase distribution for details)
%            The phase function will be called as fphase(x,y,par), where x,y are
%            vectors with coordinate pairs.
% llc:       [x,y] row vector, coordinates of the lower left corner of the tile area
% urc:       [x,y] row vector, coordinates of the upper right corner of the tile area
% vphase:    (optional) boundaries for the phase values within the polygon area;
%            vphase(1) < vphase(2). Default is [0,pi].
% grid :     a structure with sampling resolution on tile edge and area
% tol:       structure containing permissible error bounds. All tolerances are defined
%            as a FRACTION OF THE LENGTH UNIT used for x and y coordinates (usually
%            um). See 'cghparset' for default values.
% verbose:   (optional) a logical; display progress information. Default is 'false';
%
% OUTPUT
% cap:       cell array with closed polygons

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
% Ulf Griesmann, NIST, July 2019
%
% This function is based on an earlier implementation by 
% Johannes A. Soons, NIST
% ---------------------------------------------------------
    
    % --------------------------------------------------------------------
    % STEP 1: find an initial (pilot) approximation to the isophase lines
    % using the 'contourc' function on a relatively coarse grid. The
    % initial search is performed for an area that is slightly larger
    % than the specified tile area to ensure that the isophase lines
    % cross the tile boundaries.
    % --------------------------------------------------------------------

    % enlarge the tile area by one grid spacing to ensure that the
    % approximate isophases cross the tile boundaries even after
    % vertex refinement. The number of grid points remains unchanged.
    dx = (urc(1)-llc(1))/(grid.area(1)-1); % grid spacing
    dy = (urc(2)-llc(2))/(grid.area(2)-1);
    ellc = llc - 0.5*[dx,dy];         % extend tile area by one
    eurc = urc + 0.5*[dx,dy];         % grid spacing
    
    % find approximate isophase curves
    [aiso,apha,kiso] = isophase_approx(fphase,par,ellc,eurc,vphase,grid.area,verbose);
    if isempty(aiso) % tile has no isophase curves
        
        % check phase in one of the tile corners
        tph = fphase(llc(1),llc(2),par);
        if (tph > vphase(1)) && (tph < vphase(2)) % tile is filled
            cap = tile;
        else
            cap = [];
        end
        
        pok = [];
        return
    end

    % --------------------------------------------------------------------
    % STEP 2: We want to avoid closed isophase lines within a tile
    % area. An isophase curve is assumed to be closed when all its
    % vertices are within the area of the tile. While this is not a
    % conclusive test that the true isophase is closed, sub-dividing
    % the tile will cause no harm even if the isophase is not closed.
    % --------------------------------------------------------------------

    tile = [llc(1),llc(2);... % edge 1  % lower left corner
            urc(1),llc(2);... % edge 2
            urc(1),urc(2);... % edge 3
            llc(1),urc(2);... % edge 4
            llc(1),llc(2)];   % must be closed polygon
    
    isoctr = closed_isophase_check(tile,aiso,tol.vertex);
    
    if ~isempty(isoctr) % partition the tile at centroid
            
        if verbose
            fprintf('Calling ''phase2cgh_tile_pa'' recursively on partitioned tile.\n');
        end

        [tile_a,tile_b] = partile(tile, isoctr);
        cap_a = phase2cgh_tile_pa(fphase,par,tile_a(1,:),tile_a(3,:),...
                                  vphase,tol,grid,verbose);
        cap_b = phase2cgh_tile_pa(fphase,par,tile_b(1,:),tile_b(3,:),...
                                  vphase,tol,grid,verbose);
        cap = vertcat(cap_a,cap_b);
        
        return  % our work is done
        
    end

    % --------------------------------------------------------------------
    % STEP 3: At this point no closed isophase lines are found. A
    % sufficiently close polynomial approximation of an isophase
    % contour requires additional points between the vertex points of
    % the initial approximation. The number of additional vertices is
    % calculated from the radius of the isophase osculating circle at
    % the mid-points between two vertices and the allowable deviation
    % of the polygon from the true isophase curve.
    % --------------------------------------------------------------------

    brac = 0.5*max([dx,dy]); % initial size of search bracket: ~1/2 cell size
    [ciso,piso] = isophase_refine(fphase,par,aiso,apha,tol,brac,verbose);

    % --------------------------------------------------------------------
    % STEP 4: The end-vertices of the refined isophase curves are
    % likely to be located near, but not on the edges of the tile
    % area, and we want to replace them with vertices that are located
    % precisely on the edges. Since all isophase lines cross the tile
    % edges, this can be done by clipping the isophase lines to the
    % original tile edges followed by refinement of the edge crossings.
    % --------------------------------------------------------------------
    
    [ciso,eiso] = isophase_clip(fphase,par,ciso,piso,tile,tol,brac,verbose);
    if isempty(ciso)     % tile has no isophase curves
        cap = pok = [];
        return
    end

    %---------------------------------------------------------------------------
    % STEP 5: assemble isophase curves with adjacent phase indices into closed
    % polygons that circumscribe areas with phase values in the phase intervals
    % 2*pi*k + vphase, k = ...-2,-1,0,1,2,...
    %--------------------------------------------------------------------------

    cap = isophase_fill(fphase,par,ciso,eiso,vphase,tol,tile,verbose);
    
end % finis terra
