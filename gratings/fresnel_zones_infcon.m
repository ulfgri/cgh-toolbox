function [rad] = fresnel_zones_infcon(rmax, flen, lambda);
%function [rad] = fresnel_zones_infcon(rmax, flen, lambda);
%
% Calculates the transition radii for a Fresnel zone plate with one
% infinite conjugate (the "textbook" Fresnel zone plate). 
%
% INPUT
% rmax   : maximum radius of the zone plate
% flen   : flocal length (distance from zone plate center to focus)
%          in the same units as the radii.
% lambda : wavelength of light in the same units as the focal length
%          and radii
%
% OUTPUT
% rad    : a column vector with radii of transitions between opaque and
%          transparent zones

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Short version: This software is in the PUBLIC DOMAIN.
% Ulf Griesmann, NIST, February 2008

    % check arguments
    if nargin < 3
        error('missing arguments');
    end

    % number of boundaries within rmax
    N = ceil( 2 * (sqrt(flen^2 + rmax^2) - flen) / lambda );

    % generate list of radii 
    m = [1:N];
    rad = sqrt( m * lambda .* (flen + m * lambda/4) )';
    
end
