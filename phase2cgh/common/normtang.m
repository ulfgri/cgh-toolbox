function [T,N] = normtang(Px,Py)
%function [T,N] = normtang(Px,Py)
% 
% Calculates unit normal and unit tangent vectors from the phase
% gradients at points on a set of isophase curves
%
% INPUT
% Px,Py :  vectors with components of the phase gradient vector
%
% OUTPUT
% T :     nx2 matrix with the unit tangent vectors, one per row
% N :     nx2 matrix with the unit normal vectors, one per row
%         T,N are oriented such that the exterior product T /\ N > 0.

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

    N = [Px,Py] ./ sqrt(Px.*Px + Py.*Py); % same directions as gradients
    T = [N(:,2),-N(:,1)];                 % tangents are orthogonal to normals
    
    % we want T,N with positve wedge products: T /\ N > 0
    wp_neg = T(:,1) .* N(:,2) - T(:,2) .* N(:,1) < 0;
    T(wp_neg) = -T(wp_neg);
    
end
