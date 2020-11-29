function [vph] = phase_zernike_radial(vx,vy,spar)
% function [vph] = phase_zernike_radial(vx,vy,spar)
%
% Phase function for an optic that is described by a Zernike
% polynomial with radial terms only.
%
% NOTE: All dimensions must have the same units
%
% INPUT
% vx:      a vector with desired x coordinates
% vy:      a vector with desired y coordinates
% vph:     phase in radians
% spar:    structure with parameters describing the phase function
%            spar.zc    : a vector with the RADIAL Zernike polynomial
%                         coefficients. Coefficients are in OST(*) order
%                         and for the OST coordinate definition. A
%                         maximum of 20 Zernike terms are allowed.
%            spar.vcen  : center of the zone plate
%            spar.rrad  : Normalization radius for the unit circle
%
% OUTPUT
% vph:     phase in radians
%
% (*) OST stands for D. Malacara, "Optical Shop Testing", Wiley &
%     Sons. Same definition of Zernike polynomials as in Born and Wolf.

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
    nrt = length(spar.zc);      % number of radial terms
    if nrt > 20
        error('no more than 20 radial Zernike terms are supported.');
    end
    
    % function is defined relative to center    
    vx = vx - spar.vcen(1);
    vy = vy - spar.vcen(2);

    % actual indices of radial terms in the OST Zernike polynomials
    idxtab = [1,5,13,25,41,61,85,113,145,181,221,265,313,365,421,481,545,613,685,761];
    zc = zeros(idxtab(nrt),1);    % all Zernike terms
    zc(idxtab(1:nrt)) = spar.zc;  % copy radial terms into correct position

    % calculate phase
    vph = zern_eval_pos([vx,vy], spar.vcen, spar.rrad, zc, 1, 1);

end
