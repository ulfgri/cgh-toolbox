function [sidx] = boundary_vertex_init(pval,kval,eixy,nedg,tile,sing,tol)
%function [sidx] = boundary_vertex_init(pval,kval,eixy,nedg,tile,sing,tol)
%
% Determines a set of initial vertices of the phase phase curves
% on the tile edge such that duplicate tracing of isophase curves
% is avoided.
%
% INPUT
% pval :  vector of phase values at the tile edge intersections
% kval :  corresponding phase indices
% eixy :  nx2 matrix with SORTED edge intersection coordinates
% nedg :  corresponding edge number of intersections
% tile :  polygon circumscribing tile area
% sing :  structure with singularity information
% tol :   structure with tolerances
%
% OUTPUT
% sidx :  indices of phase boundary curve start intersections

% NOTE: the algorithm assumes that phase boundary lines do not cross and
%       that each boundary curve has two vertices on the tile edge
    
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Initial version: Ulf Griesmann, NIST, October 2019
    
    % it can happen that there is only one boundary line left that
    % terminates on the singularity
    if length(pval) == 1
        sidx = 1;
        return
    end
    
    % find out where along the edge the singularity is located and
    % which intersection FOLLOWS the singularity
    if ~isempty(sing)                             % if present on this tile
        sedg = onedges(tile,sing.pos,tol.vertex); % edge with singularity
        eixy = vertcat(eixy,sing.pos);            % add singularity to intersections
        nedg = vertcat(nedg,sedg);                % and edges
        eixy = eixy(edge_sort(eixy,nedg),:);      % sort expanded intersections
    
        % intersection that follows the singularity in ORIGINAL eixy
        ixfs = find(vecnorm(eixy-sing.pos,2,2) < 2*eps);

        % rotate pval,kval CCW such that ixfs --> 1
        pidx = shift([1:length(pval)], -ixfs+1);
        pval = pval(pidx);
        kval = kval(pidx);
    end
    
    % divide pval into vphase(1) and vphase(2) groups; 1=low, 2=high
    poff = pval - 2*pi*kval;      % phase offsets
    pmin = min(poff);
    pmax = max(poff);
    [~,plev] = histc(poff,[pmin,0.5*(pmin+pmax),pmax]);
    plev(plev==3) = 2; % see 'histc' documentation
    
    % true if polygon start vertex
    lix = false(length(pval),1);
    
    % find adjacent edge intersections with identical phase levels
    duplicate_phase = find( ~diff(plev',1) ); % need a row vector
    if isempty(duplicate_phase)
        
        % this can happen when 'boundary_curves' was called
        % recursively and only a few boundary curves that 
        % terminate on the singularity remain to be traced.
        lix(:) = true;

    else
    
        % mark start vertices
        prv = 1;        % previous index
        isi = true;     % flag - is start index
        
        for cur = unique([duplicate_phase,length(lix)]) % avoid duplication of last index
            if isi
                lix(prv:cur) = true;
            end
            isi = ~isi;
            prv = cur+1;
        end
        
    end

    % rotate 1 --> ixfs and return indices of start vertices
    if ~isempty(sing)
        lix = shift(lix, ixfs-1);
    end

    % return start indices
    sidx = find(lix);
    
end
