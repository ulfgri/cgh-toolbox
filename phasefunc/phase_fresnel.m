function [vph] = phase_fresnel(vx,vy,spar)
% function [vph] = phase_fresnel(vx,vy,spar)
%
% Calculates the 2D spatial phase distribution (in radians) of a
% Fresnel zone plate with collimated illumination.
%
% INPUT
% vx:      desired x coordinates
% vy:      desired y coordinates
% spar:    structure with parameters of the phase function
%            spar.rfocus :      primary focal length of the zone plate
%            spar.fcen :        center coordinate of the zone plate
%            spar.wavelength :  design wavelength
%
% OUTPUT
% vph:     phase in radians

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Version: 1.0
% Author: Johannes Soons; NIST; Feb 2008 - Sep 2011 
% Review: Johannes Soons; NIST; Sep 2011 
% Status: OK
% ---------------------------------------------------------

    % function is defined relative to center  
    vx = vx - spar.fcen(1);
    vy = vy - spar.fcen(2);
    
    if isrow(vx), vx = vx'; end
    if isrow(vy), vy = vy'; end
    
    vph = 2*pi*(sqrt(spar.rfocus^2 + vx.^2 + vy.^2) - spar.rfocus)/spar.wavelength;
  
end
