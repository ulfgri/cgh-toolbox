function [dis,ori,cit] = edge_dist(X,E,x,e,tile,tol=[])
%function [dis,ori,cit] = edge_dist(edg1,xy1,edg2,xy2,tile,tol=[])
% 
% Calculates the distances along the edge of a tile between points 
% on the edges of a tile. The distance is the SHORTER of the two paths
% that connect the first point to the second. The largest possible 
% distance is half of the tile circumference. The orientation of the
% shortest path is also returned.
% 
%          y ^
%            |   3
%            +<----*-+
%            |       ^
%          4 |       | 2
%            v       |
%            o--*--->+ ---> x
%                1
%
% INPUT
% X,E :      nx2 matrix with edge points and vector with corresponding
%            edge numbers E
% x,e :      coordinates of a point on the tile edge e
% tile :     polygon circumscribing a rectangular tile
% tol:       (optional) if present, points must be on the tile
%            edges within the specified tolerance
%
% OUTPUT
% dis :      vector with distances along the tile edge between point (x,e)
%            and all points in (X,E)
% ori :      orientations of the shortest path from the first to the
%            second point. Equal to 1 if path is CCW (positive), equal
%            to -1 if the path runs CW (negative).
% cit :      circumference of the tile

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
% Ulf Griesmann, NIST, August-September 2019

    % check input
    if nargin < 5
        error('function requires 5 arguments.');
    end

    % check on-edgeness
    if ~isempty(tol) && any( ~onedges(tile,vertcat(X,x),tol) )
        error('one or more points are not on the tile edge.')
    end

    % tile edge lengths
    Dx = vecnorm(tile(2,:) - tile(1,:));
    Dy = vecnorm(tile(3,:) - tile(2,:));
    
    % distances around tile
    Dt = [Dx,Dy,Dx,Dy];
    cit = sum(Dt);          % tile circumference must be returned
        
    % edge distances (not vectorized but comprehensible ...)
    dis = ori = zeros(size(X,1),1);
    for m = 1:size(X,1)
        [dis(m),ori(m)] = edge_dist_single(X(m,:),E(m),x,e,tile,Dt);
    end
    
end % function


%--------------------------------------------------------------------    
    
function [dis,ori] = edge_dist_single(xy1,edg1,xy2,edg2,tile,Dt)
% 
% Distance of two edge points
%
% INPUT
% xy1,xy2 :  coordinates of two points on two tile edges
% e1,e2 :    edge numbers of the two points
% tile :     nx2 matrix with vertices of a polygon circumscribing
%            the area of a tile. One vertex per row.
% Dt :       lengths of the tile edges
%
% OUTPUT
% dis :      distance along the tile edge between two points
% ori :      orientation of the shortest path from the first to the
%            second point. Equal to 1 if path is CCW (positive), equal
%            to -1 if the path runs CW (negative).
    
    if edg1 == edg2 % both points are on the same edge
      
        dis = vecnorm(xy2 - xy1);
        if any(xy2 - xy1) < 0
            ori = -1;
        else
            ori = 1;
        end
        
    else

        % total distance of involved edges
        eidx = enxt = edg1;
        while true            
            enxt = mod(enxt,4) + 1;
            eidx = [eidx,enxt];
            if enxt == edg2
                break
            end
        end

        D = sum(Dt(eidx)) - exdist(edg1,xy1,'prev') - exdist(edg2,xy2,'next');
        dis = min(D,sum(Dt)-D); % pick the smaller of two ways
        
        if dis == D  % we went round the tile in CCW direction
            ori = 1;
        else
            ori = -1;
        end
        
    end
    
    function d = exdist(edg,xy,direct)
    % calculates the distance from a point on a tile edge
    % either to the next tile corner ('next') or to the 
    % previous tile corner ('prev')
 
        % distance from the previous tile corner
        dp = vecnorm(tile(edg,:) - xy);
        
        switch direct
          case 'prev'
            d = dp;
          case 'next'
            d = abs(Dt(edg) - dp);
        end
    end
    
end % function
