function [didx] = find_discon_grid(P,trsh)
%function [didx] = find_discon_grid(P,trsh)
%
% Finds a discontinuity in a vector (grid) of phase values, typically
% sampled along a line. 
%
% INPUT
% P :       a vector with phase values
% trsh :    minimum difference between consecutive phase values
%           that is identified as a discontinuity
%
% OUTPUT
% didx :    vector of indices at which discontinuities occor in P
%           The indices point to the elements in P BEFORE the 
%           discontinuity.

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Initial version: Ulf Griesmann, NIST, July 2019

    if iscolumn(P), P = P'; end
    
    % calculate modulus of derivative
    dP = abs(diff(P));
    
    % find locations of phase jumps
    didx = find(dP > trsh);
    
end
