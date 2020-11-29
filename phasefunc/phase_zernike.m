function [vph] = phase_zernike(vx,vy,spar)
% function [vph] = phase_zernike(vx,vy,spar)
%
% Phase function for an optic that is modeled by a Zernike polynomial
%
% NOTE: All dimensions must have the same units
%
% INPUT
% vx:      a vector with desired x coordinates
% vy:      a vector with desired y coordinates
% spar:    structure with parameters describing the phase function
%            spar.zc         : a vector with Zernike polynomial
%                              coefficients describing a Zernike
%                              phase map.
%                              Coefficients are in OST order and
%                              for the OST coordinate definition.
%            spar.vcen        : center of the zone plate
%            spar.rrad        : Normalization radius for the unit
%                               circle.
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

    % Zernike coefficients must be a column vector
    if isrow(spar.zc), spar.zc = spar.zc'; end

    % function is defined relative to center    
    vx = vx - spar.vcen(1);
    vy = vy - spar.vcen(2);
    
    % calculate phase
    vph = zern_eval_pos([vx,vy], spar.vcen, spar.rrad, spar.zc, 1, 1);

end
