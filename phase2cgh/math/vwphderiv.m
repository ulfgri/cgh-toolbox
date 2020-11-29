function [Px,Py,Pxx,Pyy,Pxy,isocur] = vwphderiv(pfun,par,xy,h)
%function [Px,Py,Pxx,Pyy,Pxy,isocur] = vwphderiv(pfun,par,xy,h)
%
% Vectorized form of 'wphderiv'. Calculates first and second
% derivatives of a wrapped, non-continuous phase function using a
% symmetric secant approximation, and the isophase curvatures.
%
% INPUT
% pfun :  function handle of the phase function
% par :   a structure with parameters passed to the phase
%         function: pfun(x,y,par)
% xy :    nx2 matrix with coordinates at which to calculate the
%         derivatives, one per row.
% h :     coordinate increment for derivative calculations
%
% OUTPUT
% Px,Py :       vectors with components of the first partial
%               derivatives at locations in xy
% Pxx,Pyy,Pxy : vectors with second partial derivatives and
%               cross-derivative at locations in xy
% isocur :      vector with isophase curvatures at locations in xy

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
%
% Initial version: Ulf Griesmann, NIST, July 2019

    if nargin < 4
        error('function ''vwphderiv'' requires 4 input arguments.');
    end

    x = xy(:,1);
    y = xy(:,2);
    
    if nargout <= 2 % only first derivatives are requested
        
        % evaluate phase function
        X = [x+h; x-h;   x;   x; x];   
        Y = [  y;   y; y+h; y-h; y];  
        
        Ph = reshape(exp(i*pfun(X,Y,par)), length(x), 5);

        % first partial derivatives
        Z  = Ph(:,5); % Z = exp(i * phase)
        Zx = 0.5 * (Ph(:,1) - Ph(:,2)) / h;
        Zy = 0.5 * (Ph(:,3) - Ph(:,4)) / h;
    
        Px = real(-i * Zx ./ Z);
        Py = real(-i * Zy ./ Z);
        
    else
        
        %      1    2    3    4    5    6    7    8    9
        X = [x+h; x-h;   x;   x;   x; x+h; x+h; x-h; x-h];
        Y = [  y;   y; y+h; y-h;   y; y+h; y-h; y+h; y-h];
        
        Ph = reshape(exp(i*pfun(X,Y,par)),length(x),9);

        % first partial derivatives
        Z  = Ph(:,5); % Z = exp(i * phase)
        Zx = 0.5 * (Ph(:,1) - Ph(:,2)) / h;
        Zy = 0.5 * (Ph(:,3) - Ph(:,4)) / h;    
        Px = real(-i * Zx ./ Z);
        Py = real(-i * Zy ./ Z);

        % second partial derivatives
        h2 = h*h;
        Zxx = (Ph(:,1) + Ph(:,2) - 2*Ph(:,5)) / h2;
        Zyy = (Ph(:,3) + Ph(:,4) - 2*Ph(:,5)) / h2;
        Zxy = 0.25*(Ph(:,6) - Ph(:,7) - Ph(:,8) + Ph(:,9)) / h2;
        Pxx = real(-i * (Zxx./Z + Px.*Px));
        Pyy = real(-i * (Zyy./Z + Py.*Py));
        Pxy = real(-i * (Zxy./Z + Px.*Py));
        
        % isophase curvatures @ (x,y)
        isocur = (2*Pxy.*Px.*Py - Pxx.*Py.*Py - Pyy.*Px.*Px) ./ (Px.*Px + Py.*Py).^1.5;

    end

end
