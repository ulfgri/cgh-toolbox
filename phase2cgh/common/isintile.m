function lti = isintile(tile,xy,tol=0)
%function lti = isintile(tile,xy,tol=0)
%
% Checks if one or more points are inside a rectangular tile area
% within a specified tolerance.
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
% tile (o).
%
% INPUT
% tile :  a 5x2 matrix defining a polygon that circumscribes a
%         rectangular area. All edge intersections are assumed 
%         to be orthogonal (no checking !)
% xy :    a nx2 matrix with plane coordinates, one per row
% tol :   (optional) a positive scalar. points that are on the edge 
%         within the specified tolerance are considered to be inside 
%         the tile. Default is 0.
%
% OUTPUT
% lti :   a nx1 logical vector. 'true' if the corresponding xy(n,:)
%         is inside the tile area. Points on tile edges are
%         considered to be outside the tile by default.

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

    tol = abs(tol);
    
    llc = tile(1,:);
    urc = tile(3,:);
    
    lti = (xy(:,1) - llc(1) >= -tol) & (urc(1) - xy(:,1) >= -tol) & ...
          (xy(:,2) - llc(2) >= -tol) & (urc(2) - xy(:,2) >= -tol);

end
