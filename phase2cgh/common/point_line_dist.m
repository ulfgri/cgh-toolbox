function [odi,ldi] = point_line_dist(x1,x2,xp)
%function [odi,ldi] = point_line_dist(x1,x2,xp)
%
% Calculates the shortest distance of points from straight lines
%
% INPUT
% x1,x2 :  nx2 matrices with coordinates of point pairs defining a line
% xp :     either a point or n points in the plane
%
% OUTPUT
% odi :    orthogonal distances from point(s) xp to lines
% ldi :    distances x1 to point(s) xp

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
% Ulf Griesmann, NIST, September 2019

    % replicate xp if it's the same for all lines
    if size(x1,1) > 1 && size(xp,1) == 1
        xp = repmat(xp,size(x1,1),1);
    end
    
    % angle cosines between lines x1-x2 and x1-xp
    cosa = dot(x2-x1,xp-x1,2) ./ (vecnorm(x2-x1,2,2).*vecnorm(xp-x1,2,2));
    
    % distances
    ldi = vecnorm(xp-x1,2,2);
    odi = sqrt(1-cosa.^2) .* ldi; % is always > 0
    
end % function
