function [Px,Py,Pxx,Pyy,Pxy,isocur] = vphderiv(pfun,par,xy,h)
%function [Px,Py,Pxx,Pyy,Pxy,isocur] = vphderiv(pfun,par,xy,h)
%
% Vectorized form of 'phderiv' to estimate assorted derivatives of a
% phase function pfun(x,y,par) at multiple locations (x,y) using a
% symmetric secant approximation, and estimates the curvatures of the
% isophase curves at the (x,y).
%
% INPUT
% pfun :    function handle of a scalar phase function. 
% par :     a structure with parameters passed to the phase
%           function pfun(x,y,par).
% xy :      a Mx2 matrix with coordinate pairs (x,y), one per row, 
%           at which to calculate the derivative(s) of the phase function.
% h :       a scalar, small coordinate increment for calculating 
%           the derivatives.
%
% OUTPUT
% Px, Py :     components of the gradient vectors at (x,y)
% Pxx,Pyy,Pxy: second partial derivativess at (x,y)
% isocur :     the curvatures of the isophase curves at (x,y), 
%              i.e., the inverse of the radii of the osculating
%              circles. NOT THE PHASE SURFACE CURVATURES.

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
% Initial version: Ulf Griesmann, NIST, September 2018

    if nargin < 4
        error('function ''vphderiv'' requires 4 input arguments.');
    end
    
    % Evaluate the phase function. All function calls are combined into a
    % single evaluation of the phase function for multiple locations
    % to take advantage of vectorization in the phase function
    % implementation.

    x = xy(:,1);
    y = xy(:,2);
    N = size(xy,1);

    if nargout <= 2 % only gradient is requested
        
        X = [x+h; x-h;   x;   x];   
        Y = [  y;   y; y+h; y-h];   
        
        Ph = reshape(pfun(X,Y,par),N,4);
    
        % first partial derivatives at (x,y)
        Px = 0.5 * (Ph(:,1) - Ph(:,2)) / h;
        Py = 0.5 * (Ph(:,3) - Ph(:,4)) / h;
        
    else

        %      1    2    3    4    5    6    7    8    9
        X = [x+h; x-h;   x;   x;   x; x+h; x+h; x-h; x-h];   
        Y = [  y;   y; y+h; y-h;   y; y+h; y-h; y+h; y-h];   
        
        Ph = reshape(pfun(X,Y,par),N,9);

        % first and second partial derivatives at (x,y)
        h2 = h^2;
        Px = 0.5 * (Ph(:,1) - Ph(:,2)) / h;
        Py = 0.5 * (Ph(:,3) - Ph(:,4)) / h;
        Pxx = (Ph(:,1) + Ph(:,2) - 2*Ph(:,5)) / h2;
        Pyy = (Ph(:,3) + Ph(:,4) - 2*Ph(:,5)) / h2;
        Pxy = 0.25 * (Ph(:,6) - Ph(:,7) - Ph(:,8) + Ph(:,9) ) / h2;
                
        % curvatures of isophases @ (x,y) (quite a beautiful relationship !)
        isocur = (2*Pxy.*Px.*Py - Pxx.*Py.*Py - Pyy.*Px.*Px) ./ (Px.*Px + Py.*Py).^1.5;
        
    end
    
end
