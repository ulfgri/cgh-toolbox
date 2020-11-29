function [cap] = phase2cgh_tile_fi(fphase,par,llc,urc,vphase,tol,grid,verbose=false)
% function [cap] = phase2cgh_tile_fi(fphase,par,llc,urc,vphase,tol,grid,verbose=false)
%
% Calculates the layout of a computer generated hologram (CGH) for a
% given phase function in a rectangular area using an isophase contour
% following algorithm. The CGH layout is described by polygons that
% approximate the true phase contours. The polygons are boundaries of
% areas with phase values within intervals 2*pi * k + vphase, where k
% is an arbitrary but fixed integer for each polygon and vphase is a
% pair of boundary phase values.
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
% vphase:    boundaries for the phase values within the polygon area
% tol:       structure containing permissible error bounds. All tolerances are
%            defined as a FRACTION OF THE LENGTH UNIT used for x
%            and y coordinates (usually um). See 'cghparset' for default values
% grid:      structure with number of samples for phase
%            function evaluation. See 'cghparset'
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
% Ulf Griesmann, NIST, September 2018
% ---------------------------------------------------------

    pkg load optim

    % Construct polygon circumscribing the tile. Tile edges are 
    % numbered as follows:
    %
    %          y ^
    %            |   3
    %            +<------+
    %            |       ^
    %          4 |       | 2
    %            v       |
    %            o------>+ ---> x
    %                1
    %
    % Tile vertices are numbered in counter-clockwise direction
    % beginning with the vertex at the lower left corner of the
    % tile (o). Tiles without edge intersections are handled after
    % checking for a phase function extremum.
    %--------------------------------------------------------------------------
    
    tile = [llc(1),llc(2);... % edge 1  % lower left corner
            urc(1),llc(2);... % edge 2
            urc(1),urc(2);... % edge 3
            llc(1),urc(2);... % edge 4
            llc(1),llc(2)];   % must be closed polygon
    
    %--------------------------------------------------------------------------
    % STEP 1: We first want to find out if the tile are contains any
    % closed isophase curves. If yes, the tile is partitioned into two
    % sub-tiles at the location of the extremum and 'phase2cgh_tile_fi'
    % is called recursively on both sub-tiles.  For the detection of
    % closed isophase curves we sample the function on a coarse grid
    % and use the 'contourc' to detect closed curves.
    %--------------------------------------------------------------------------
        
    % calculate isophase curve approximations (is fast)
    aiso = isophase_approx(fphase,par,llc,urc,vphase,grid.area,verbose);
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

    
    % check for presence of closed isophase curves
    isoctr = closed_isophase_check(tile,aiso,tol.vertex);

    if ~isempty(isoctr) % partition the tile at centroid

        if verbose
            fprintf('Calling ''phase2cgh_tile_fi'' recursively on partitioned tile.\n');
        end

        [tile_a,tile_b] = partile(tile, isoctr);
        cap_a = phase2cgh_tile_fi(fphase,par,tile_a(1,:),tile_a(3,:),...
                                  vphase,tol,grid,verbose);
        cap_b = phase2cgh_tile_fi(fphase,par,tile_b(1,:),tile_b(3,:),...
                                  vphase,tol,grid,verbose);
        cap = vertcat(cap_a,cap_b);

        return  % our work is done

    end
    
    %---------------------------------------------------------------------------
    % STEP 2: find all locations on the tile edges where isophase curves cross
    % an edge. 
    %---------------------------------------------------------------------------

    [eixy,pval,kval,nedg] = isophase_edges(fphase,par,tile,vphase, ...
                                           tol,grid.edge,verbose);

    %---------------------------------------------------------------------------
    % STEP 3: no closed isophase curves were found in the tile. The isophases
    % are now traced accross the tile area starting from the edges, i.e. pairs
    % of points with equal phase values, or index k, on the edges are connected
    % with an isophase line that is approximated by a polygon.
    %--------------------------------------------------------------------------

    [ciso,eiso] = isophase_curves(fphase,par,tile,eixy,pval,kval,nedg, ...
                                  tol,verbose); 

    %---------------------------------------------------------------------------
    % STEP 4: assemble isophase curves with adjacent phase indices into closed
    % polygons that circumscribe areas with phase values in the phase intervals
    % 2*pi*k + vphase, k = ...-2,-1,0,1,2,...
    %--------------------------------------------------------------------------

    cap = isophase_fill(fphase,par,ciso,eiso,vphase,tol,tile,verbose);
    
end % finis terra
