function [Px,Py,Pxx,Pyy,Pxy,isocur] = wphderiv(pfun,par,xy,h)
%function [Px,Py,Pxx,Pyy,Pxy,isocur] = wphderiv(pfun,par,xy,h)
%
% Calculates the first derivatives of a wrapped, non-continuous phase
% function using a symmetric secant approximation.
%
% INPUT
% pfun :  function handle of the phase function
% par :   a structure with parameters passed to the phase
%         function: pfun(x,y,par)
% xy :    coordinate where phase function is evaluated
% h :     coordinate increment derivative calculation
%
% OUTPUT
% Px,Py :      components of the phase function gradient
% Pxx,Pyy,Pxy: second partial derivatives at (x,y)
% isocur :     the curvature of the isophase curve at (x,y),
%              i.e., the inverse of the radius of the osculating
%              circle. NOT THE PHASE SURFACE CURVATURE.

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
        error('function ''wphderiv'' requires 4 input arguments.');
    end

    % evaluate phase function
    x = xy(1); 
    y = xy(2);
    
    if nargout <= 2 % first derivatives only
        
        X = [x+h; x-h;   x;   x; x];   
        Y = [  y;   y; y+h; y-h; y];   
        
        Ph = exp(i*pfun(X,Y,par));

        % calculate partial derivatives
        Z  = Ph(5); % Z = exp(i * phase)
        Zx = 0.5 * (Ph(1) - Ph(2)) / h;
        Zy = 0.5 * (Ph(3) - Ph(4)) / h;
    
        Px = real(-i * Zx / Z);
        Py = real(-i * Zy / Z);
        
    else
        
        %      1    2    3    4    5    6    7    8    9
        X = [x+h; x-h;   x;   x;   x; x+h; x+h; x-h; x-h];
        Y = [  y;   y; y+h; y-h;   y; y+h; y-h; y+h; y-h];
        
        Ph = exp(i*pfun(X,Y,par));

        % first partial derivatives
        Z  = Ph(5); % Z = exp(i * phase)
        Zx = 0.5 * (Ph(1) - Ph(2)) / h;
        Zy = 0.5 * (Ph(3) - Ph(4)) / h;
        Px = real(-i * Zx / Z);
        Py = real(-i * Zy / Z);

        % second partial derivatives
        h2 = h*h;
        Zxx = (Ph(1) + Ph(2) - 2*Ph(5)) / h2;
        Zyy = (Ph(3) + Ph(4) - 2*Ph(5)) / h2;
        Zxy = 0.25 * (Ph(6) - Ph(7) - Ph(8) + Ph(9) ) / h2;
        Pxx = real(-i * (Zxx/Z + Px*Px));
        Pyy = real(-i * (Zyy/Z + Py*Py));
        Pxy = real(-i * (Zxy/Z + Px*Py));
        
        % isophase curvatures @ (x,y)
        isocur = (2*Pxy*Px*Py - Pxx*Py*Py - Pyy*Px*Px) / (Px*Px + Py*Py)^1.5;
        
    end

end
