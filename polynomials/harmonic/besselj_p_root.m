function vz = besselj_p_root(nu,zmin,zmax)
%function vz = besselj_p_root(nu,zmin,zmax)
%
% find roots of the derivative of the Bessel function of
% the first kind
%
% INPUT
% nu   : order of the Bessel function
% zmin : smallest possible root value
% zmax : largest possible root value
    
% OUTPUT
% vz   : roots
%

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Version: 1.0
% Author: Johannes Soons; NIST; Sep 2006 
% Review: Johannes Soons; NIST; Sep 2006 
% Status: OK
% ---------------------------------------------------------

    opts = optimset('TolX',eps, 'MaxIter',250);

    % find zero crossings
    if zmin > zmax 
        dum = zmin; 
        zmin = zmax; 
        zmax = dum; 
    end
    dz = 0.05;
    vze = (zmin:dz:zmax)';
    vi  = find(abs(diff(sign(besselj_p(nu,vze)) == 1)) == 1);
    
    % refine zeros
    vz = vfzero(@(z)besselj_p(nu,z),vze(vi)+[-2,3]*dz,opts);
  
end
