function [xyd] = find_discon_line(fphase,par,xya,xyb,pha,dtol)
%function [xyd] = find_discon_line(fphase,par,xya,xyb,pha,dtol)
%
% Finds the locations in the x-y plane at which straight lines meet
% phase discontinuities using a bisection algorithm.
%
% INPUT
% fphase :  phase function handle fphase(x,y,par)
% par :     phase function parameters
% xya :     coordinates of points on the line BEFORE discontinuity 
% xyb :     coordinates of points on the line AFTER discontinuity 
% pha :     vector with phase values at points xya
% dtol :    spatial tolerance of the discontinuity intersection location
%
% OUTPUT
% xyd :     nx2 matrix with locations at which the lines meet the
%           phase discontinuities on the same side of the
%           discontinuity as the points xya

% NOTE: The current interval is bisected into parts that are not
% precisely equal. This is done to avoid evaluation of the phase
% function at exactly the phase discontinuity, where it is
% undefined, in some situations. 
    
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Initial version: Ulf Griesmann, NIST, August 2019

% input sanity check
    if ~all([size(xya,1),size(xyb,1)] == length(pha))
        error('incompatible matrix dimensions on input.')
    end

    % unit vector from xya to xyb
    D = (xyb - xya);
    d = vecnorm(D,2,2);               % distance
    D = D./d;                         % unit direction

    % Find the discontinuity along the line xya + delta*D
    delta = zeros(size(d));
    xyc = xya;                        % current location

    while true

        % take a bisection step in D direction
        d = 0.4999*d;                 % bisect distance d (see NOTE)
        delta += d;                   % and add to current value of delta
        xyd = xya + delta .* D;       % in x-y plane

        % check if we have reached the point of marginal return
        if all(vecnorm(xyd-xyc,2,2)<dtol)
            break
        end

        % evaluate phase at next point, check for phase jump
        phd = fphase(xyd(:,1),xyd(:,2),par);
        jmp = abs(phd - pha) > pi;    % should be 2*pi
    
        % in case of phase jump unwind last bisection step
        if any(jmp)            
            delta(jmp) -= d(jmp);
            xyd(jmp,:) = xyc(jmp,:);  % restore previous value
        end

        % remember location after bisection step
        xyc = xyd;

    end    

    xyd = xyc;  % makes sure we stay on THIS side of disontinuity

end
