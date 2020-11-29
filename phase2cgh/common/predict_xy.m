function [pxy,B] = predict_xy(cxy,isocur,T,N,direct,tol)
%
% Calculates the approximate locations of the next vertices based on
% the osculating circles (2nd order approximations) of the isophases
% at the current vertex locations.
%
% INPUT
% cxy :   coordinates of vertices on isophases, one per row
% isocur: isophase curvatures at the vertices
% T,N :   isophase tangent and normal vectors at cxy
% direct: direction along isophase (along or against T)
% tol :   structure with permissible tolerances. See 'cghparset'.
%
% OUTPUT
% pxy :   coordinates of the predicted vertices
% B :     bracketing vector; normal vectors at pxy with lengths that are
%         equal to the estimated distance, along the normal, between pxy 
%         and the isophase.

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
%
% Initial version: Ulf Griesmann, NIST, July-October 2019

    R = abs(1./(isocur + eps));                   % radii of osculating circles
    S2 = tol.segdev*(8*R - 4*tol.segdev);         % square of secant lengths
    if ~isempty(tol.seglen)                       % set maximum allowable
        S2(S2>tol.seglen^2) = tol.seglen^2;       % segment length 
    end
    nu = 0.5 * S2 ./ R;                           % coordinates (tau,nu) of next
    tau = sqrt(S2 - nu.^2);                       % point in moving frames T,N
    pxy = cxy + direct.*tau.*T + sign(isocur).*nu.*N;
    
    % use Frenet-Serret Eqn to estimate N at pxy
    Np = N - 2 * asin(sqrt(S2) ./ (2*R)) .* T;
    
    % bracketing vector
    B = nu .* Np;
    
end
