function xyb = point_between_vertices(xy1,e1,xy2,e2,tile,dist)
%
% Finds a point xyb between two phase boundary - edge intersections
% 1) xyb is located on the line connecting two points xy1 and xy2
% 2) it is a specified distance from xy1
%
% INPUT
% xy1,e1 :  nx2 matrix with coordinates and vector with
%           corresponding edge numbers of first edge intersection
% xy2,e2 :  same for second edge intersection
%           it is assumed that e2 >= e1
% tile :    tile polygon
% dist:     distance from polygon at which xyb are located
%
% OUTPUT
% xyb :     coordinates of points between edge intersections

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
% Ulf Griesmann, NIST, September 2018 - September 2019

    % some input checking
    numv = length(e1);
    if any([size(xy1,1),size(xy2,1),length(e2)] ~= numv)
        error('incompatible input arguments');
    end
    
    % pre-alloc output
    xyb = zeros(numv,2);

    for k = 1:numv

        if e1(k) == e2(k)
        
            % use point close to start vertex on the line xy2 - xy1
            xyd = xy2(k,:) - xy1(k,:);
    
        else
        
            % use a point close to the polygon in the direction of a tile corner
            con = end_connector(e1(k),e2(k),tile); % not vectorizable
            xyd = con(1,:) - xy1(k,:);
            
        end
        
        xyb(k,:) = xy1(k,:) + dist * xyd/vecnorm(xyd);
        
    end
end
