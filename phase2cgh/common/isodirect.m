function direct = isodirect(xy,T,N,tile,h)
%function direct = isodirect(xy,T,N,tile,h)
%
% Determines the direction in which an isophase line must be followed
% from the edge of a tile to the interior.
%
% INPUT
% xy :     nx2 matrix with points on the edges of the tile
% T,N :    nx2 matrices with unit tangent and normal vectors, one
%          for each point in xy
% tile :   a polygon circumscribing the tile area
% h :      vector with steps along the direction of the isophase
%          tangent, or a scalar with a stepsize that applies to all points.
%
% OUTPUT
% direct:  direction of the isophase tangent that leads INTO the 
%          the tile area

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
% Initial version: Ulf Griesmann, NIST, September 2018
    
    % take a step in positive direction
    xyp = xy + h .* T;
    
    % return direction factors
    direct = ones(size(xy,1),1);
    direct(~isinpolygon(tile,xyp)) = -1;
    
end
