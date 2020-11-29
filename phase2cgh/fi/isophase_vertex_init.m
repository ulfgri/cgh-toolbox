function [sidx] = isophase_vertex_init(pval,kval)
%function [sidx] = isophase_vertex_init(pval,kval)
%
% Determines a set of initial vertices of the isophase curves on
% the tile edge such that duplicate tracing of isophase curves is 
% avoided.
%
% INPUT
% pval :  vector of phase values at the tile edge intersections
% kval :  corresponding phase indices
%
% OUTPUT
% sidx :  indices of isophase curve start intersections

% NOTE: the algorithm assumes that isophase lines do not cross and
%       that each isophase curve has two vertices on the tile edge
    
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
        error('edge intersections cannot all be start vertices');
    end
    
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

    % consistency check (each isophase has two edge intersections !)
    if 2*sum(lix) ~= length(pval)
        error('2 * #start vertices ~= #edge intersections');
    end
    
    % return indices of start vertices
    sidx = find(lix);
    
end
