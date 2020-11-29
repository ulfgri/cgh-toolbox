function cic = closed_isophase_check(tile,aiso,tol)
%function cic = closed_isophase_check(tile,aiso,tol)
%
% Checks if a closed isophase line exists within a tile area. A
% closed polygon is one that has no vertex on any edge of the tile
% and all vertices are strictly inside the tile area.
%
% INPUT
% tile :  a 5x2 matrix defining a polygon that circumscribes a
%         rectangular area. All edge intersections are assumed 
%         to be orthogonal (no checking !)
% aiso :  a cell array with polygons
% tol :   tolerance for the placement of a vertex to be considered
%         on the edge of the tile
%
% OUTPUT
% cic :   a 1x2 vector with the coordinates of the centroid of 
%         the closed isophase curve. Empty, when no closed isophase
%         curve exists.

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

    cic = [];
    
    for k=1:numel(aiso)

        % first eliminate polygons with vertices outside the tile
        % area ('pa' generates these)
        if any( ~isinpolygon(tile, aiso{k}) )
            continue
            
        % for other algorithms, isophase curves have 2 intersections
        elseif sum( onedges(tile,aiso{k},tol) > 0 ) < 2 % numbers of vertices on edges
            cic = polycentr(aiso(k));                   % polygon centroid
            break                                       % found one
        end
        
    end

end
