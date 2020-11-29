function fdh_init(nn=100, nm=100)
%function fdh_init(nn=100, nm=100)
%
% Generate spatial frequencies and normalization constants
% for Free Disk Harmonic functions
%
% INPUT
% nn :  maximum order for radial terms.
% nm :  maximum order for azimuthal terms.
%
% OUTPUT
% Data will be written to a file 'fdh.mat'

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Loosely based on original code by Dick Brunson, 12/1/97
%
% Version: 1.0
% Author: Johannes Soons; NIST; Sep 2006 
% Review: Johannes Soons; NIST; Sep 2006 
% Status: OK
% removed unfinished +/-1 normlization. Ulf Griesmann, Sep 2019
% ---------------------------------------------------------

    pkg load optim  % for vfzero
    
    % initialize arrays    
    fdh.mf = zeros(nn,nm+1);
    fdh.ma = zeros(nn,nm+1);
    vn = (1:nn);

    fprintf('\nPre-computing Disk Harmonic parameters (%d radial, %d azimuthal terms)\n',nn,nm);
    fprintf('processing order: '); fflush(stdout());
    
    for m = 0:nm

        if m && ~mod(m,10)
            fprintf('%d .. ',m);
            fflush(stdout());
        end

        % obtain roots of derivative of Bessel function of the first kind
        vz = besselj_p_root(m,m,nn*pi*2);
      
        % check validity of roots (roots should be spaced pi apart at infinity)
        if abs(diff(vz)-pi) > 0.5 
            error('Irregular root sequence for order %g',m)
        end
        if length(vz) <= nn 
            error('Missing roots for order %g',n);
        end
        
        vz = vz(1:nn);
        
        fdh.ma(:,m+1) = sqrt(1./((1-(m./vz).^2).*besselj(m,vz).^2));
        fdh.mf(:,m+1) = vz;
    
  end % for m ...

  fprintf('\n\n');
  
  save -binary 'fdh.mat' fdh
    
end
