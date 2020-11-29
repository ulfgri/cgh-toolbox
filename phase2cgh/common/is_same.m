function same = is_same(xy1,xy2)
%function same = is_same(xy1,xy2)
% 
%  checks if two points are the same
%
% INPUT
% xy1,xy2 :  1x2 vectors with coordinates of a point
%
% OUTPUT
% same :     logical, true if both points are the same

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Initial version: Ulf Griesmann, NIST, September 2018

    same = all(abs(xy1-xy2) < 2*eps);
    
end
