function order = edge_point_order(e1,x1,e2,x2)
%function order = edge_point_order(e1,x1,e2,x2)
%
% For two points on a tile edge the function determines
% which point is ahead of the other. The edges can include
% interior edges at discontinuities (5) and singularities (6).
%
% INPUT
% e1,e2 :  edge numbers of two points on the edges
% x1,x2 :  coordinate of two points on the edges
%
% OUTPUT
% order :  1 if origin distance of x2 > origin distance of x1
%          0 if the points coincide (they shouldn't)
%         -1 if origin distance of x2 < origin distance of x1

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Short version: This software is in the PUBLIC DOMAIN.
% Ulf Griesmann, NIST, September 2018

    if e1 == e2  % same edge, compare coordinates
        switch e1
            
          case 1  % on edge 1 compare x
            if x2(1) > x1(1)
                order = 1;
            elseif x2(1) < x1(1)
                order = -1;
            else
                order = 0;
            end
            
          case 2  % on edge 2 compare y
            if x2(2) > x1(2)
                order = 1;
            elseif x2(2) < x1(2)
                order = -1;
            else
                order = 0;
            end
            
          case 3  % on edge 3 compare -x
            if x2(1) < x1(1)
                order = 1;
            elseif x2(1) > x1(1)
                order = -1;
            else
                order = 0;
            end
            
          case 4  % on edge 4 compare -y
            if x2(2) < x1(2)
                order = 1;
            elseif x2(2) > x1(2)
                order = -1;
            else
                order = 0;
            end
            
          otherwise
            error('edge number out of range 1-4.');
        end
    else                 % different edge, compare edge #'s
        if e2 > e1
            order = 1;
        else
            order = -1;
        end
    end

end
