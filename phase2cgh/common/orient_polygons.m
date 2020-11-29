function [corp,eorp] = orient_polygons(ciso,eiso)
%function [corp,eorp] = orient_polygons(ciso,eiso)
%
% Orients polygons such that their start vertices
% are closer to the tile origin than their end vertices
%
% INPUT
% ciso :     cell array with polygon approximations to isophase lines
% eiso :     nx2 matrix with pairs of edge numbers at the start and at the
%            end of polygons in 'ciso'. One row per polygon.
%
% OUTPUT
% corp :     cell array with oriented polygons
% eorp :     edge numbers of oriented polygons

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
% Ulf Griesmann, NIST, October 2019
% ---------------------------------------------------------

    corp = ciso;
    eorp = eiso;
    
    for k=1:numel(ciso)
        xy = ciso{k};
        if edge_point_order(eiso(k,1),xy(1,:),eiso(k,2),xy(end,:)) < 0
            corp{k} = polyrev(xy);
            eorp(k,:) = eiso(k,[2,1]);
        end
    end

end
