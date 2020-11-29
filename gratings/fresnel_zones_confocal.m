function [rad] = fresnel_zones_confocal(rmax, fleno, fleni, lambda);
%function [rad] = fresnel_zones_confocal(rmax, fleno, fleni, lambda);
%
% Calculates the transition radii for a confocal Fresnel zone
% plate. The substrate of the zone plate is assumed to be so thin and
% any spherical aberration due to the substrate is neglected.
%
% INPUT (all dimensions and lengths must use the same unit)
% rmax   : maximum radius of the zone plate
% fleno  : flocal length (distance from zone plate center to focus)
%          on the object side.
% fleni  : flocal length (distance from zone plate center to focus)
%          on the image side.
% lambda : wavelength of light.
%
% OUTPUT
% rad    : a column vector with radii of transitions between opaque and
%          transparent zones.

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Short version: This software is in the PUBLIC DOMAIN.
% Ulf Griesmann, NIST, February 2012

    % check arguments
    if nargin < 3
        error('missing arguments');
    end

    % number of boundaries within rmax
    N = ceil( 4 * (sqrt(flen^2 + rmax^2) - flen) / lambda );
    
    % generate list of radii 
    m = [1:N];
    rad = sqrt( 0.5 * m * lambda .* (flen + m * lambda/8) )';

end


