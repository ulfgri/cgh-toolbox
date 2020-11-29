function ridx = edge_sort(xy, edg)
%function ridx = edge_sort(xy, edg)
%
% Sorts points on the edges of a tile in counterclockwise
% order. The points in xy are assumed to be on the edges of
% a rectangular tile, and the tile edges are numbered as follows:
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
% The origin of the tile is its lower left corner
%
% INPUT
% xy :   coordinates of points on edges, one point per row
% edg :  edges on which the points are located
%
% OUTPUT
% ridx :  indices of sorted coordinates
    
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Ulf Griesmann, NIST, September 2018

    ridx = [];
    idx = [1:length(edg)]';

    % group points according to edges
    ie1 = edg==1;
    ie2 = edg==2;    
    ie3 = edg==3;
    ie4 = edg==4;

    % edge 1: sort x-values ascending
    if any(ie1)
        X = xy(ie1,:);
        [~,sidx] = sort(X(:,1),'ascend');
        ridx = vertcat(ridx,idx(ie1)(sidx));
    end
    
    % edge 2: sort y-values ascending
    if any(ie2)
        X = xy(ie2,:);
        [~,sidx] = sort(X(:,2),'ascend');
        ridx = vertcat(ridx,idx(ie2)(sidx));
    end
    
    % edge 3: sort x-values descending
    if any(ie3)
        X = xy(ie3,:);
        [~,sidx] = sort(X(:,1),'descend');
        ridx = vertcat(ridx,idx(ie3)(sidx));
    end
    
    % edge 4: sort y-values descending
    if any(ie4)
        X = xy(ie4,:);
        [~,sidx] = sort(X(:,2),'descend');
        ridx = vertcat(ridx,idx(ie4)(sidx));
    end
    
end
