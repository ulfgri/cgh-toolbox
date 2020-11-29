function [vi] = tile_neighbor(stile,it)
%function [vi] = tile_neighbor(vsize,it)
%
% Returns the indices of tiles in a collection of tiles that
% surround a particular tile (including in diagonal direction)
%
% stile : structure array with tile information
%            stile.llc : [x,y] coordinates of the lower left tile corner
%            stile.urc : [x,y] coordinates of the upper right tile corner
% it :    index number of tile whose neighbors are to be found
% vi:     vector with indices of the tile neighbors

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic. 
%
% Version 2.0
% Author: Johannes Soons; NIST; Oct 2012 
% Review: Johannes Soons; NIST; Oct 2012 
% Status: OK
% ---------------------------------------------------------

    llc = [stile.llc];
    llc = reshape(llc,2,numel(stile))';
    urc = [stile.urc];
    urc = reshape(urc,2,numel(stile))';
  
    % check left edge, interval check: s2 < e1 and s1 < e2
    vi = find(abs(stile(it).llc(1)-urc(:,1)) < eps &...
              (llc(:,2) <= stile(it).urc(2)+eps &...
               stile(it).llc(2)-eps <= urc(:,2)));
          
    % check right edge, interval check: s2 < e1 and s1 < e2
    vi = [vi;find(abs(stile(it).urc(1)-llc(:,1)) < eps &...
                  (llc(:,2) <= stile(it).urc(2)+eps &...
                   stile(it).llc(2)-eps <= urc(:,2)))];

    % check bottom edge, interval check: s2 < e1 and s1 < e2
    vi = [vi;find(abs(stile(it).llc(2)-urc(:,2)) < eps &...
                  (llc(:,1) <= stile(it).urc(1)+eps &...
                   stile(it).llc(1)-eps <= urc(:,1)))];
                
    % check top edge, interval check: s2 < e1 and s1 < e2
    vi = [vi;find(abs(stile(it).urc(2)-llc(:,2)) < eps &...
                  (llc(:,1) <= stile(it).urc(1)+eps &...
                   stile(it).llc(1)-eps <= urc(:,1)))];
            
end

  