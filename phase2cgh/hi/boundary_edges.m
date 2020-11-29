function [eixy,pval,kval,nedg] = boundary_edges(fphase,par,tile,vphase,tol,grid,sing,verbose)
%function [eixy,pval,kval,nedg] = boundary_edges(fphase,par,tile,vphase,tol,grid,sing,verbose)
%
% NOT INTENDED TO BE CALLED BY USERS DIRECTLY
% Finds points on the edges of a rectangular area where a function
% fphase has the values 2*pi*k+vphase. 
%
% INPUT
% fphase:    function handle to user defined function describing the 2D phase
%            distribution (in radians)
% par:       structure with parameters describing the phase distribution
%            (see user defined function describing phase distribution for details)
%            The phase function will be called as fphase(x,y,par), where x,y are
%            vectors with coordinate pairs.
% tile :     a closed polygon circumscribing a rectangular tile area. The
%            polygon is a mx2 matrix with one vertex per row.
% vphase:    (optional) boundaries for the phase values within the polygon area;
%            vphase(1) < vphase(2).
% tol:       structure containing permissible error bounds. See 'cghparset' for details.
% grid:      number of samples on the tile edges used for phase function evaluation.
% sing:      a structure with information about a phase singularity on a tile edge.
%            Can be empty. See 'cghparset' for details
% verbose:   (optional) a logical; display progress information. Default is 'false';
%
% OUTPUT
% eixy :     a nx2 matrix with intersection coordinates, one per row
% pval :     vector of phase values at the intersections
% kval :     vector with phase indices k at the intersections
% nedg :     vector with edge number for each of the intersections

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
    
    % init output arguments
    eixy = pval = kval = nedg = [];
        
    % find edge intersections
    for ke = 1:4

        [xy,ph,kv] = isophase_line_xing(fphase,par,tile(ke,:),tile(ke+1,:), ...
                                        vphase,tol,grid);

        if ~isempty(xy)
            
            % sort points along the edge
            ne = repmat(ke,size(xy,1),1);
            sidx = edge_sort(xy, ne);
            xy = xy(sidx,:);
            ph = ph(sidx);
            kv = kv(sidx);        
            
            % last intersection can't be the tile corner
            if all(abs(xy(end,:)-tile(ke+1,:)) < 5*eps)
                nin = size(xy,1)-1;
                xy = xy(1:nin,:);
                ph = ph(1:nin);
                kv = kv(1:nin);
                ne = ne(1:nin);
            end
        
            % append to arrays
            eixy = vertcat(eixy, xy);
            pval = vertcat(pval, ph); 
            kval = vertcat(kval, kv); 
            nedg = vertcat(nedg, ne);
    
        end
        
        if verbose
            fprintf('     %d phase boundary - edge %d crossings\n', size(xy,1),ke);
        end

    end

    % return if we don't have any intersections
    if isempty(eixy)
        return
    end

    if verbose
        fprintf('     %d crossings in total\n', size(eixy,1));
    end

    % the algorithm that detects boundary-edge intersections can also return 
    % exceptional points, which must be removed from the list of isophase-edge
    % intersections, because they are not vertices of an isophase.
    [px,py] = vwphderiv(fphase,par,eixy,tol.hderiv);
    idz = (px == 0) & (py == 0);
    if any(idz)
        eixy = eixy(~idz,:);
        pval = pval(~idz);
        kval = kval(~idz);
        nedg = nedg(~idz);
    end

    % make sure that none of the intersections coincides with a singular
    % phase point. Phase boundary lines can terminate at a singular point,
    % but they cannot start at the singular point because the start direction
    % would be undefined.
    if ~isempty(sing)
        ids = all( abs(eixy-sing.pos) < tol.vertex, 2);
        if any(ids)
            eixy = eixy(~ids,:);
            pval = pval(~ids);
            kval = kval(~ids);
            nedg = nedg(~ids);
        end
    end
    
end
