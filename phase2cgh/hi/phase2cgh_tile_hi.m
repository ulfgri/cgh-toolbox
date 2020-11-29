function [cap] = phase2cgh_tile_hi(fphase,par,llc,urc,vphase,tol,grid,sing,verbose=false)
% function [cap] = phase2cgh_tile_hi(fphase,par,llc,urc,vphase,tol,grid,sing,verbose=false)
%
% Calculates the layout of a computer generated hologram (CGH) for a
% phase function WITH HILBERT TERM in a rectangular area using an
% phase boundary following algorithm. The CGH layout is described by
% polygons that approximate the true phase contours. The polygons are
% boundaries of areas with phase values within intervals 2*pi * k +
% vphase, where k is an arbitrary but fixed integer for each polygon
% and vphase is a pair of boundary phase values.
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
%            vphase(1) < vphase(2).
% tol:       structure containing permissible error bounds. 
%            See 'cghparset' for default values
% grid:      structure with number of samples for phase function evaluation.
%               grid.area :  sampling grid for tile area.
%               grid.edge :  sampling grid for tile edge.
% sing:      a structure array with information about phase singularities.
%            See 'cghparset' for details
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
    % STEP 1: Check if any of the singularities is on a tile edge, and which
    % one. 
    %--------------------------------------------------------------------------
        
    sxy = reshape([sing(:).pos],2,numel(sing))';  % one per row
    isoe = find( onedges(tile,sxy,2*eps) );       % singularities on a tile edge
    if isempty(isoe)
        sing = [];
    elseif length(isoe) == 1
        sing = sing(isoe);
    else
        error('more than one singularity on tile edges - redesign tiling.');
    end    
        
    %---------------------------------------------------------------------------
    % STEP 2: find all locations on the tile edges where isophase segments cross
    % an edge. 
    %---------------------------------------------------------------------------

    [eixy,pval,kval,nedg] = boundary_edges(fphase,par,tile,vphase, ...
                                           tol,grid.edge,sing,verbose);

    %---------------------------------------------------------------------------
    % STEP 3: The phase boundaries are now traced accross the tile area starting 
    % from the edges.
    %--------------------------------------------------------------------------

    [ciso,eiso] = boundary_curves(fphase,par,tile,eixy,pval,kval,nedg, ...
                                  tol,sing,verbose);

    %---------------------------------------------------------------------------
    % STEP 4: assemble isophase curves with adjacent phase indices into closed
    % polygons that circumscribe areas with phase values in the phase intervals
    % 2*pi*k + vphase, k = ...-2,-1,0,1,2,...
    %--------------------------------------------------------------------------

    cap = boundary_fill(fphase,par,ciso,eiso,vphase,tol,tile,sing,verbose);

end % finis terra
