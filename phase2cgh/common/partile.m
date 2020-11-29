function [sta,stb] = partile(tile, pit, direct=[])
%function [sta,stb] = partile(tile, pit, direct=[])
%
% Divides a rectangular tile through a vertical or horizontal cut at a
% specified point insdide the tile such that the resulting partition
% has the largest possible squareness.
%
% INPUT
% tile :   a polygon describing the tile, one vertex per row,
%          tile(1,:) is the lower left corner.
% pit :    coordinate pair [x,y] of a point inside the tile
% direct:  (optional) direction of the cut; either 'vert' or 'horz'.
%          When the direction is not specified, the tile will be 
%          divided into tiles of minimal squareness.
%
% OUTPUT
% sta,stb : polygons describing two subtiles of the input tile. The
%           polygons are normalized such that the first vertex is
%           the lower left corner.

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
    
    % first check if point is inside the tile
    if ~isinpolygon(tile,pit)
        error('point ''pit'' is not inside the tile.');
    end
    
    % tiles resulting from horizontal cut
    hsta = [tile(1,:); ...
            tile(2,:); ...
            [tile(2,1),pit(2)]; ...
            [tile(1,1),pit(2)]; ...
            tile(1,:)];
    hstb = [[tile(2,1),pit(2)]; ...
            [tile(1,1),pit(2)]; ...
            tile(3,:); ...
            tile(4,:); ...
            [tile(2,1),pit(2)]];
    
    % tiles resulting from vertical cut
    vsta = [tile(1,:); ...
            [pit(1),tile(1,2)]; ...
            [pit(1),tile(3,2)]; ...
            tile(4,:); ...
            tile(1,:)];
    vstb = [[pit(1),tile(1,2)]; ...
            tile(2,:); ...
            tile(3,:); ...
            [pit(1),tile(3,2)]; ...
            [pit(1),tile(1,2)]];

    % return partitioned tiles
    if isempty(direct)
        
        % minimal squareness of horizontal/vertical cuts
        msqh = min(squareness(hsta),squareness(hstb));
        msqv = min(squareness(vsta),squareness(vstb));
        
        % return the partition with larger squareness
        if msqh > msqv
            sta = hsta; 
            stb = hstb;
        else
            sta = vsta;
            stb = vstb;
        end
        
    else
        
        switch direct
          case 'horz'
            sta = hsta; 
            stb = hstb;            
          case 'vert'
            sta = vsta;
            stb = vstb;
          otherwise
            error('unknown tile partition direction.');
        end
        
    end
end

%---------------------------------------------------------------------------------


function S = squareness(tile)
%
% calculates the squareness of a rectangular tile which is defined
% here as the length of the short side divided by the length of the
% long side. The squareness is 1 for a square tile and 0 for a line.
%
% INPUT
% tile :  a polygon describing the tile, one vertex per row,
%         first vertex is the lower left corner.
%
% OUTPUT
% S :     tile squareness
    
    % lengths of two consecutive tile edges
    La = vecnorm(tile(2,:) - tile(1,:));
    Lb = vecnorm(tile(3,:) - tile(2,:));
    
    % squareness index
    S = min([La,Lb]) / max([La,Lb]);
    
end

%---------------------------------------------------------------------------------
