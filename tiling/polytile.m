function [stile] = polytile(bnd_out, bnd_inn, tilpar, plotpar=[]);
%function [stile] = polytile(bnd_out, bnd_inn, tilpar, plotpar=[]);
%
% Calculates a tiling made of rectangular tiles for an area
% circumscribed by a polygon such that the tiling covers the
% polygon. The polygon can have holes that are described by
% one or more inner polygons.
%
% INPUT:
% bnd_out :   outer boundary of the area as a n x 2 matrix of
%             vertices.
% bnd_inn :   (Optional) inner boundary for the corresponding outer 
%             boundary as a n x 2 matrix of vertices or a cell
%             array of inner boundaries. 
% tilpar :    tile parameters. Either a vector with tile size [width,height]
%             or a structure with parameters that define the tiling:
%               tilpar.tsize  :    tile size [width,height]
%               tilpar.offset :    tile offset from center. Default is [0,0]. 
%             Square tiles can be specified using a scalar tile size.
% plotpar:    optional structure with plot parameters
%                plotpar.lwidth : line width          (default: 1)
%                plotpar.fntsiz : font size in points (default: 12)
%                plotpar.withnum: plot tile numbers   (default: false)
%             No tiling is plotted if this parameter is not supplied, or empty.
%
% OUTPUT:
% stile : structure array with tile information
%   stile(k).llc :  [x,y] coordinates of the lower left corner of the tile area
%   stile(k).urc :  [x,y] coordinates of the upper right corner of the tile area

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Ulf Griesmann, NIST, June 2015

    % check arguments
    if nargin < 3
        error('polytile :  too few arguments.');
    end

    if ~isempty(bnd_inn) && ~iscell(bnd_inn) 
        bnd_inn = {bnd_inn}; 
    end

    % assemble tile parameters
    if ~isstruct(tilpar)  % for backward compatibility
        tilpar = struct('tsize',tilpar);
    end
    if ~isfield(tilpar,'offset'), tilpar.offset = [0,0]; end
    if isscalar(tilpar.tsize)
        tilpar.tsize = [tilpar.tsize,tilpar.tsize];
    end
    if isscalar(tilpar.offset)
        tilpar.offset = [tilpar.offset,tilpar.offset];
    end
        
    % center of outer polygon
    cop = polycentr({bnd_out});
    
    % find bounding box for the outer boundary that is an integer multiple of tsize
    % bounding box: lower left corner [x0,y0], upper right corner [x1,y1]
    sb = sort(bnd_out);
    x0 = sb(1,1);
    y0 = sb(1,2);
    x1 = sb(end,1);
    y1 = sb(end,2);
    
    % lower left corner of tiled area
    llxy = cop + tilpar.offset - tilpar.tsize .* ceil((cop - [x0,y0])./tilpar.tsize);
    
    % calculate tiling of rectangle (bottom left tile corners)
    xx = [llxy(1):tilpar.tsize(1):x1]';
    yy = [llxy(2):tilpar.tsize(2):y1]';
    llcc = [kron(ones(length(yy),1),xx), ...    % generate bottom left corner
            kron(yy,ones(length(xx),1))];       % coordinates of all tiles in rectangle
    
    % check if the tiles are wholly or partially within the boundary
    keep = false(1,size(llcc,1));               % remember which tiles to keep
    T = [0,0; tilpar.tsize(1),0; tilpar.tsize; 0,tilpar.tsize(2); 0,0];  % reference tile
    
    for k=1:size(llcc,1)
                      
        tile = llcc(k,:) + T;                   % tile k polygon
        
        % check if tile is fully inside any of the inner polygons or
        % fully outside the outer polygon
        if all(~isinpolygon(bnd_out,tile))
            continue  % don't keep this tile
        end
        if ~isempty(bnd_inn) && any(cellfun(@(p)all( isinpolygon(p,tile) ),bnd_inn))
            continue
        end

        % tile is part of the tiling
        keep(k) = true;
        
    end
        
    % discard tiles outside boundaries
    llcc = llcc(keep,:);
    urcc = llcc + tilpar.tsize;
    
    % pre-allocate output structure array and fill it
    stile = repmat(struct('llc',[0,0], 'urc',[0,0]),size(llcc,1),1);
    for k=1:size(llcc,1)
        stile(k).llc = llcc(k,:);
        stile(k).urc = urcc(k,:);
    end

    % plot if desired
    if ~isempty(plotpar)
        
        if ~isfield(plotpar,'lwidth'),  plotpar.lwidth = 1;  end;
        if ~isfield(plotpar,'fntsiz'),  plotpar.fntsiz = 12; end;
        if ~isfield(plotpar,'withnum'), plotpar.withnum = false; end;
        
        fprintf('\n  Number of tiles: %d\n\n', size(llcc,1));
        
        figure
        plot(bnd_out(:,1),bnd_out(:,2),'r','linewidth',plotpar.lwidth);
        hold on
        if ~isempty(bnd_inn)
            for b=1:numel(bnd_inn)
                bnd = bnd_inn{b};
                plot(bnd(:,1),bnd(:,2),'r-','linewidth',plotpar.lwidth);
            end
        end

        % plot tiles
        for k=1:size(llcc,1)
            tile = llcc(k,:) + T;
            plot(tile(:,1),tile(:,2),'b','linewidth',plotpar.lwidth);
            if plotpar.withnum
                txtpos = llcc(k,:) + [0.1,0.25].*tilpar.tsize; 
                text(txtpos(1),txtpos(2),sprintf('%d',k),'FontSize',plotpar.fntsiz);
            end
        end
        axis('equal');
        grid on
        title('Polygon Tiling');
        hold off
    end
    
end
