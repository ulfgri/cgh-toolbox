function [xyd] = find_discon_edges(fphase,par,tile,trsh,nse,dtol)
%function [xyd] = find_discon_edges(fphase,par,tile,trsh,nse,dtol)
%
% Finds the locations on the edges of a tile at which a line
% discontinuity intersects the edges (always a pair of intersections)
%
% INPUT
% fphase:  function handle to user defined function describing the 2D phase
%          distribution (in radians)
% par:     structure with parameters describing the phase distribution
%          (see user defined function describing phase distribution for details)
% tile :   a polygon circumscribing a rectangular tile
% trsh :   minimum difference between consecutive phase values that is
%          identified as a discontinuity
% nse :    number of samples on the edges
% dtol :   discontinuity location tolerance
%
% OUTPUT
% xyd :    2x2 matrix with tile-discontinuity intersections, one per row

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Initial version: Ulf Griesmann, NIST, September 2019

    xyd = xa = xb = [];
    
    for k=1:4
        xy = linspace(tile(k,:),tile(k+1,:),nse)';
        P = fphase(xy(:,1),xy(:,2),par);
        didx = find_discon_grid(P,trsh);
        
        if ~isempty(didx) % find discontinuity from both directions and average
            xa = vertcat(xy(didx,:),xy(didx+1,:)); % xa,xb bracket discontinuity
            xb = vertcat(xy(didx+1,:),xy(didx,:));
            xyd(end+1,:) = mean( find_discon_line(fphase,par,xa,xb,[P(didx);P(didx+1)],dtol) );
        end
        
    end
end
