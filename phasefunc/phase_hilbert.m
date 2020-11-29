function [ph] = phase_hilbert(x,y,par)
% function [ph] = phase_hilbert(x,y,par)
%
% Hilbert phase function with one or more terms. Can be added to a 
% Fresnel (or other) phase function to create zone lenses with 
% topological charge (helical wave/phase fronts).
%
% INPUT
% x:      vector with x coordinates
% y:      vector with y coordinates
% par:    structure array with parameters for the Hilber phase terms
%            par(k).p      :  topological charge
%            par(k).hcen   :  center coordinate (phase singularity)
%            par(k).aoff   :  angular offset in radians
%
% OUTPUT
% ph:     column vector with (discontinuous !) phase in radians

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Author: Ulf Griesmann, NIST, January 2013; updated, August 2019

    % want column vectors
    if isrow(x), x = x'; end
    if isrow(y), y = y'; end

    ph = zeros(size(x));
    
    for k = 1:length(par)
    
        % function is defined relative to center    
        xev = x - par(k).hcen(1);
        yev = y - par(k).hcen(2);

        % calculate phase
        ph += arg( exp(i * par(k).p * (par(k).aoff + atan2(yev,xev)) ) );
  
    end
end
