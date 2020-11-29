function [vph] = phase_map(vx,vy,spar)
% function [vph] = phase_map(vx,vy,spar)
%
% Calculate the 2D spatial phase distribution (in radians) of an
% arbitrary phase function that was sampled on a grid.
%
% INPUT
% vx:      desired x coordinates
% vy:      desired y coordinates
% spar:    structure with parameters describing the phase
%          distribution function phi(x,y)
%            spar.map         : map to be interpolated. Can be
%                               empty when map was stored by
%                               minterp2 first. 
%            spar.mapsize     : size(map,1) of the phase map array
%            spar.vcenm       : pixel coordinates of the center (0,0)
%            spar.vpixsep     : pixel seperation (dx,dy) of the phase map in x and y
%            spar.intmeth     : interpolation method
%                                  'linear' - bilinear interpolation
%                                  'cubic'  - bicubic interplation
%                                  'pchip'  - piecewise cubic
%                                             Hermite interpolation (best)
%
% OUTPUT
% vph:     phase in radians
%
% NOTE: the map to be interpolated must first be loaded into
%       minterp2. See documentation for 'minterp2'

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Initial version, Ulf Griesmann, NIST, October 2011

    % translate position into pixel coordinates
    vx = vx./spar.vpixsep(1)+spar.vcenm(1);
    vy = vy./spar.vpixsep(2)+spar.vcenm(2);    
    
    % interpolate
    vph = minterp2(spar.map, vx, spar.mapsize+1-vy, spar.intmeth);
    
end
