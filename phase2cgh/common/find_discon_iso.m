function [xyd] = find_discon_iso(fphase,par,xya,ica,T,N,direct,ph,tol)
%function [xyd] = find_discon_iso(fphase,par,xya,ica,T,N,direct,ph,tol)
%
% Finds the locations in the x-y plane at which isophase curves meet
% phase discontinuities using a bisection algorithm that seeks the 
% discontinuity intersection along a 2nd-order approximation of the 
% isophase curves.
%
% INPUT
% fphase :  phase function handle fphase(x,y,par)
% par :     phase function parameters
% xya :     coordinates of points on an isophase line BEFORE discontinuity 
% ica :     isophase curvatures at vertices xya
% T,N :     moving Frenet coordinate system vectors at xya
% direct :  direction in which the isophase curve is being followed
% ph :      vector with phase values of the isophases
% tol :     structure with tolerances
%
% OUTPUT
% xyd :     nx2 matrix with locations at which isophase curves meet
%           the phase discontinuities.

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
    if ~all([size(xya,1),length(ica),size(T,1),size(N,1),length(direct)] == ...
            length(ph))
        error('incompatible matrix dimensions on input.')
    end
    
    % osculating circle properties at points xya
    R = abs(1./(ica + eps));                  % radii of osculating circles
    S2 = tol.segdev*(8*R - 4*tol.segdev);     % square of secant lengths
    Tau = sqrt(S2 .*(1 - 0.25*S2./(R.*R)) );  % value of tau at next vertex 

    % Find the discontinuity along the arc connecting vertices xya to the
    % predicted vertices on the other side of the discontinuity. The
    % discontinuity is located by a bisection algorithm along coordinate
    % tau, with tau in the interval [0,Tau].
    tau = zeros(size(Tau));               % actual value(s) of tau
    xyc = xya;                            % current location

    while true

        % take a bisection step in T direction
        Tau = 0.5*Tau;                    % bisect Tau
        tau += Tau;                       % and add to current value of tau
        nu = R - sqrt(R.*R - tau.*tau);   % (tau,nu) is point on oscul. circle
        xyd = xya + direct.*tau.*T + sign(ica).*nu.*N; % in x-y plane

        % check if we have reached the point of marginal return
        if all(vecnorm(xyd-xyc,2,2)<tol.discon)
            break
        end

        % evaluate phase at next point, check for phase jump
        phd = fphase(xyd(:,1),xyd(:,2),par);
        jmp = abs(phd - ph) > pi;         % point is not exactly on isophase ... 
    
        % in case of phase jump unwind last bisection step
        if any(jmp)            
            tau(jmp) -= Tau(jmp);
            xyd(jmp,:) = xyc(jmp,:);      % restore previous value
        end

        % remember location after bisection step
        xyc = xyd;

    end    

    xyd = xyc;  % makes sure we stay on THIS side of disontinuity

end
