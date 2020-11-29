function ph = phase_freshilb(x,y,par)
%function ph = phase_freshilb(x,y,par)
%
% Phase function for spiral zone plates, a product of the Fresnel
% and Hilbert complex phase functions
%
% INPUT
% x:       vector with x coordinates
% y:       vector with y coordinates
% par:     structure with parameters of the phase function(s)
%            par.fresnel : parameters of the Fresnel phase term (see: phase_fresnel)
%            par.hilbert : parameters of the Hilbert phase term (see: phase_hilbert)
%
% OUTPUT
% ph:     phase in radians

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Initial version, Ulf GRIESMANN, NIST, July 2019
    
    ph = phase_hilbert(x,y,par.hilbert) + phase_fresnel(x,y,par.fresnel);
    
end
