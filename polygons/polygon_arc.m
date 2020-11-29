function [ar] = polygon_arc(xa, xb, C, nv)
%function [ar] = polygon_arc(xa, xb, C, nv)
%
% Generates a polygon arc between two Cartesian points in 
% positive direction
%
% INPUT
% xa :  start point (Cartesian coordinates)
% xb :  end point (Cartesian coordinates)
% C  :  center of arc
% nv :  number of vertices
%
% OUTPUT
% ar :  nv x 2 matrix of vertex points on the arc
%

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
% Ulf Griesmann, NIST, September 2018
% ---------------------------------------------------------

     % shift to (0,0), transform to polar coordinates
     xsa = xa - C;
     [aa,r] = cart2pol(xsa(1), xsa(2));
     if aa < 0, aa += 2*pi; end
     
     xsb = xb - C;
     [ab,r] = cart2pol(xsb(1), xsb(2));
     if ab < 0, ab += 2*pi; end
    
     % generate vertices in polar coord and transform back
     if abs(ab - aa) < 5*eps
         ab = aa + 2*pi; % full circle
     end
     A = linspace(aa,ab,nv)';
     R = repmat(r,nv,1);
     [x,y] = pol2cart(A,R);
     
     % shift back
     ar = [x,y] + C;
    
end
