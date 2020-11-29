function nedge = onedges(tile, xy, tol)
%function nedge = onedges(tile, xy, tol)
%
% Determines if a point within a tile is located on any of the 
% tile edges within a specified tolerance. 
%          y ^
%            |   3
%            +<------+
%            |       ^
%          4 |       | 2
%            v       |
%            o------>+ ---> x
%                1
%
% Tile vertices are numbered in counter-clockwise direction beginning with
% the vertex at the lower left corner of the tile (o). Tiles without edge
% intersections are handled after checking for a phase function extremum.
%
% INPUT
% tile :  a 5x2 matrix of vertices defining a rectangular tile
% xy :    a Mx2 matrix with points, one per row
% tol :   tolerance for the point-on-edge condition
%
% OUTPUT
% nedge : nedge(k) is the number of the edge on which point xy(k)
%         is located, or 0 if the point is not on any edge

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
% ---------------------------------------------------------

    % translate tile make (0,0) the tile origin
    V = xy   - tile(1,:);
    T = tile - tile(1,:);
    
    % orthogonal distances of xy from the tile edges
    D(:,1) = abs(V(:,2));          % distance from edge 1 
    D(:,2) = abs(T(2,1) - V(:,1)); % distance from edge 2
    D(:,3) = abs(T(3,2) - V(:,2)); % distance from edge 3   
    D(:,4) = abs(V(:,1));          % distance from edge 4

    % create a table with edge numbers
    E = repmat([1,2,3,4],size(D,1),1);
    
    % filter out points that are not on any edge
    E(D > tol) = 0;
    
    % if a point is in a corner it is on the edge with the larger index
    nedge = max(E,[],2);

    % exclude any points that are outside the tile
    % T(2,1): tile width; T(3,2): tile height
    nedge( max([D(:,2),D(:,4)],[],2) > T(2,1) | ...
           max([D(:,1),D(:,3)],[],2) > T(3,2) ) = 0;
    
end

