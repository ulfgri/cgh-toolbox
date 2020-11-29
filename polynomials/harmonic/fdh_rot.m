function [nm2,vpar2] = fdh_rot(nm,vpar,rangle)
% function [nm2,vpar2] = fdh_rot(nm,vpar,rangle)
%
% Calculates the effect of a rotation over angle rangle on the 
% coeffients vpar of a Disk Harmonic decomposition
%
% INPUT
% nm      : Radial and azimuthal order for each disk harmonic term
%           (one term per row)
% vpar    : Parameters of the Disk Harmonic decomposition.  The set of 
%           parameters will be expanded to ensure symmetry of cos and sin terms
% rangle  : rotation angle in radian (clockwise rotation is positive!)
%
% OUTPUT
% nm2 :     corresponding coefficients for the rotated decomposition
% vpar2

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Version 1.0
% Author: Johannes Soons; NIST; 2006 
% Review: Johannes Soons; NIST; 2006 
% Status: OK
% ---------------------------------------------------------

    vpar2 = zeros(size(vpar));
    nm2 = nm;

    vi = nm(:,2) == 0;
    vpar2(vi) = vpar(vi);
    nm2(vi,:) = nm(vi,:);
  
    vi = find(nm(:,2) ~= 0);
  
    for i = 1:length(vi)

        i1 = find(nm2(:,1) == nm(vi(i),1) & nm2(:,2) == nm(vi(i),2));
        i2 = find(nm2(:,1) == nm(vi(i),1) & nm2(:,2) == -nm(vi(i),2));
    
        % pad parameters if required
        if isempty(i2)
            vpar2 = [vpar2;0];
            nm2 = [nm2;[nm(vi(i),1),-nm(vi(i),2)]];
            i2 = length(nm2(:,1));
        end

        % rotate
        %
        % sin(l*(phi-a)) = sin(l*phi) cos(l*a) - cos(l*phi) sin(l*a)
        % cos(l*(phi-a)) = cos(l*phi) cos(l*a) + sin(l*phi) sin(l*a)        
        if nm(vi(i),2) < 0
            vpar2(i1) = vpar2(i1)+vpar(vi(i))*cos(abs(nm(vi(i),2))*rangle);
            vpar2(i2) = vpar2(i2)-vpar(vi(i))*sin(abs(nm(vi(i),2))*rangle);
        else
            vpar2(i1) = vpar2(i1)+vpar(vi(i))*cos(abs(nm(vi(i),2))*rangle);
            vpar2(i2) = vpar2(i2)+vpar(vi(i))*sin(abs(nm(vi(i),2))*rangle);
        end
    end

    % sort parameters
    [nm2,vi] = unique(nm2,'rows');
    vpar2 = vpar2(vi);
  
end
