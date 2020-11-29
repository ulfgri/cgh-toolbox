function [vph] = phase_even_asphere(vx,vy,spar)
% function [vph] = phase_even_asphere(vx,vy,spar)
%
% Phase function for an even (rotationally invariant) asphere that can
% be described with the following equation
%
%                     phi(rho) = rho^2 / (R + sqrt(R^2-(1+k)rho^2)) + 
%                                a1 + a2 * rho^2 + a3 * rho^4 + .....
%
% NOTE: All dimensions must have the same units
%
% INPUT
% vx:      a vector with desired x coordinates
% vy:      a vector with desired y coordinates
% spar:    structure with parameters describing the phase function
%            spar.R      : vertex radius in the same units as vx, vy
%            spar.k      : conic constant
%            spar.a      : a vector with polynomial coefficients
%                          describing the phase function.
%            spar.vcen   : center coordinate of the zone plate
%
% OUTPUT
% vph:     phase in radians

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.

    % make column vectors
    if isrow(vx), vx = vx'; end
    if isrow(vy), vy = vy'; end

    % function is defined relative to center    
    if ~isfield(spar,'vcen')
        spar.vcen = [0,0];
    end
    vx = vx - spar.vcen(1);
    vy = vy - spar.vcen(2);
  
    % calculate phase
    rho2 = (vx.^2 + vy.^2);  % rho^2 is column vector
    vph = rho2 ./ (spar.R + sqrt(spar.R^2 - (1+spar.k).*rho2)) + ...
          fast_poly_eval(spar.a, rho2);

end
