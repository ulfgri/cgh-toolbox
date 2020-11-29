function [map] = fdh_eval(mr,ma,nm,vp,intmeth)
% function [map] = fdh_eval(mr,ma,nm,vp,intmeth)
%
% Evaluates a set of Free Disk-Harmonic (FDH) terms
%
% mr      : matrix or column vector with normalized radii for each point
% ma      : matrix or column vector with angles for each point
%           (defined relative to y, x corresponds to pi/2)
% nm      : Radial and azimuthal order for each disk harmonic term
%           (one term per row)
% vp      : coefficient of an FDH term
% map     : sum of all specified FDH terms
% intmeth : (optional, default 'spline') interpolation method
%           'none'    - no interpolation of Bessel functions (slow)
%           'nearest' - nearest neighbor interpolation
%           'linear'  - linear interpolation
%           'spline'  - spline interpolation
%
% References:
%
% N. M. Milton and M. Lloyd-Hart, "Disk Harmonic Functions for Adaptive Optics Simulations,"
% in Adaptive Optics: Analysis and Methods/Computational Optical Sensing and Imaging/Information
% Photonics/Signal Recovery and Synthesis Topical Meetings on CD-ROM, Technical Digest
% (Optical Society of America, 2005), paper AWA3 .
%
% S. C. Verrall and R. Kakarala, "Disk-harmonic coefficients for invariant
% pattern recognition," J. Opt. Soc. Am. A 15, 389-401 (1998)

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

    persistent fdh;

    if nargin < 5, intmeth = ''; end
    if isempty(intmeth), intmeth = 'spline'; end

    % initialize
    map = zeros(size(mr))+NaN;
  
    if isempty(vp) || isempty(nm) 
        return 
    end    
    if length(nm(:,1)) ~= length(vp) 
        error('Size of coefficient vector ''vp'' is incompatible with matrix order ''nm''.'); 
    end

    % set points to evaluate and initialize 
    vi = find((mr(:) <= 1) & ~isnan(mr(:)));
    map(vi) = 0;

    % obtain pre-calculated fdh parameters
    if isempty(fdh)
        if isempty( locate('fdh.mat') )
            fdh_init;
        end
        load fdh.mat fdh;
    end
  
    [nmaxa,mmax] = size(fdh.mf);
    mmax = mmax-1;
  
    % ensure fdh parameters are present        
    if any((nm(:,1) > nmax) | nm(:,1) < 0) 
        error(' Specified radial order cannot be processed, check fdh_init.m'); 
    end
    if any(abs(nm(:,2) > mmax)) 
        error(' Specified azimuthal order cannot be processed, check fdh_init.m'); 
    end  

    if strcmp(lower(intmeth),'none')

        for i = 1:length(nm(:,1))

            if nm(i,1) == 0
                if nm(i,2) ~= 0
                    error(' Specified disk harmonic term does not exist n = %g, m = %g',nm(i,1),nm(i,2));
                else
                    map(vi) = map(vi)+vp(i);
                end
            else

                if nm(i,2) < 0
                    map(vi) = map(vi)+vp(i)*sqrt(2)*fdh.ma(nm(i,1),-nm(i,2)+1)*...
                              besselj(-nm(i,2),fdh.mf(nm(i,1),nm(i,2)+1)*mr(vi)).*...
                              sin(-nm(i,2)*ma(vi));
                else
                    if nm(i,2) == 0
                        map(vi) = map(vi)+vp(i)*fdh.ma(nm(i,1),1)*...
                                  besselj(0,fdh.mf(nm(i,1),1)*mr(vi));
                    else
                        map(vi) = map(vi)+vp(i)*sqrt(2)*fdh.ma(nm(i,1),nm(i,2)+1)*...
                                  besselj(nm(i,2),fdh.mf(nm(i,1),nm(i,2)+1)*mr(vi)).*...
                                  cos(nm(i,2)*ma(vi));
                    end
                end
            end
        end
        
    else
        
        vr = (0:0.002:1.006)';

        for i = 1:length(nm(:,1))

            if nm(i,1) == 0
                if nm(i,2) ~= 0
                    error(' Specified disk harmonic term does not exist n = %g, m = %g',nm(i,1),nm(i,2));
                else
                    map(vi) = map(vi)+vp(i);
                end
            else

                % speed up evaluation bessel function using interpolation
                vrr = vr*fdh.mf(nm(i,1),abs(nm(i,2))+1);
                vb = besselj(abs(nm(i,2)),vrr);     
                if nm(i,2) < 0
                    map(vi) = map(vi)+vp(i)*sqrt(2)*fdh.ma(nm(i,1),-nm(i,2)+1)*...
                              interp1(vrr,vb,fdh.mf(nm(i,1),-nm(i,2)+1)*mr(vi),lower(intmeth)).*...
                              sin(-nm(i,2)*ma(vi));
                else
                    if nm(i,2) == 0,
                        map(vi) = map(vi)+vp(i)*fdh.ma(nm(i,1),1)*...
                                  interp1(vrr,vb,fdh.mf(nm(i,1),1)*mr(vi),lower(intmeth));
                    else
                        map(vi) = map(vi)+vp(i)*sqrt(2)*fdh.ma(nm(i,1),nm(i,2)+1)*...
                                  interp1(vrr,vb,fdh.mf(nm(i,1),nm(i,2)+1)*mr(vi),lower(intmeth)).*...
                                  cos(nm(i,2)*ma(vi));
                    end
                end
            end
        end
    end
  
end
