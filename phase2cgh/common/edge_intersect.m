function [xy,edgnum,eidx] = edge_intersect(xy_ins,xy_out,tile,eixy,nedg,tol)
%function [xy,edgnum,eidx] = edge_intersect(xy_ins,xy_out,tile,eixy,nedg,tol)
%
% Looks up the locations at which the isophase polygons intersect the
% tile edges after the isophase following algorithm has reached the
% first point outside the tile. Since the locations of all edge
% intersections have already been calculated, this is now just a
% matter of looking up the correct intersection.
%
% NOTE: this function uses the edge distance to find the nearest
% intersection without recourse to the phase of the isophase curves
%
% INPUT
% xy_ins :  nx2 matrix with last isophase vertices INSIDE the tile
%           before crossing an edge. One row per isophase
% xy_out :  nx2 matrix with first isophase vertices OUTSIDE the
%           tile after crossing an edge
% tile :    polygon circumscribing a rectangular tile area
% eixy :    coordinates of all edge intersections, one per row
% nedg :    edge numbers for the intersections in 'eixy'
% tol:      structure containing permissible error bounds. See 'cghparset'.
%
% OUTPUT
% xy :      nx2 matrix with coordinates of the exact edge intersections
% edgnum :  vector with the numbers of the intersected edges
% eidx :    index of the edge intersections in the list 'eixy' of intersections

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Initial version: Ulf Griesmann, NIST, September 2018

    % some sanity checking
    nv = size(xy_ins,1);   % number of vertices
    if nv ~= size(xy_out,1)
        error('numbers of inside and outside vertices must be equal.'); 
    end
    
    % allocate output, etc.
    edgnum = eidx = zeros(nv,1);
    dist = zeros(length(nedg),1);
        
    % find the intersections of the isophase secants and the tile. The
    % distance from the secant intersection to the true intersection
    % should be no larger than about tol.segdev if the polygon
    % vertices are spaced correctly.
    for k = 1:size(xy_ins,1)

        % find the secant intersection
        [xysi,nseg] = polyinters(tile,[xy_ins(k,:);xy_out(k,:)]);

        % if no intersection was found, something is seriously wrong
        if isempty(xysi)
            fprintf('\nno edge intersection was found for polygon %d\n', k);
            fprintf('input arguments are likely invalid or inconsistent.\n\n');
            error('aborted');
        end

        % close to a tile corner it can happen that the corner is traversed
        % by a single polygon segment and xy_ins is a point on the tile edge. 
        % This situation can result in several complications:
        % 1) 'polyinters' returns two intersections. One of them is xy_ins, 
        %    and it must be discarded
        % 2) in rare cases, when the terminal intersection is very close to
        %    the edge on which xy_ins is located, the intersection is found on
        %    the same edge as xy_ins, and close to xy_ins (this is a consequence
        %    of the quadratic interpolation used for vertex prediction). 
        %    Since xy_ins and xy_out MUST be on different tile edges this case 
        %    is handled separately. FIXME: use circle-line intersection
        % 3) xy_ins is a tile corner and 'polyinters' returns two identical 
        %    intersections that are the tile corner
        
        % check if xy_ins is on a tile edge
        on_edge = onedges(tile,xy_ins(k,:),5*eps);
        if on_edge
        
            % first handle case of two edge intersections; pick the correct one
            if size(xysi,1) > 1
                ixo = all(abs(xysi - xy_ins(k,:)) < tol.vertex,2);
                if all(ixo) % case 3
                    xysi = xysi(1,:); % pick one
                    nseg = nseg(1,:);
                elseif all(~ixo)
                    error('failed to find edge intersection');
                else        % case 1
                    xysi = xysi(~ixo,:);
                    nseg = nseg(~ixo,:);
                end
            
            % check if the edge intersection was erroneously found on the initial edge
            elseif on_edge == nseg(1,1)
                [xysi,nseg(1,1)] = find_nearest_inters(xysi,nseg(1,1),eixy,nedg,tile);
                 
            else
                fprintf('\nunrecognized type of corner intersection for polygon %d\n\n', k);
                error('aborted');
            end
        end
        
        % store the number of the intersected tile edge
        edgnum(k) = nseg(1,1);
              
        % find the closest edge intersection among all intersections
        % on the intersected edge
        dist(:) = realmax(); % largest possible number
        same = nedg==nseg(1,1);
        dist(same) = vecnorm(eixy(same,:) - xysi,2,2);
        eidx(k) = find(dist == min(dist));   

    end % for

    % return exact isophase-tile intersections
    xy = eixy(eidx,:);
    
end


%------------------------------------------------------------------

function [xyo,eo] = find_nearest_inters(xyi,ei,eixy,nedg,tile)
%
% finds the nearest intersection on the tile edge that is NOT 
% on the same edge as a specified intersection
% 
% INPUT
% xyi,ei :  a phase boundary intersection on the tile edge
% eixy :    all known phase intersections, NOT necessarily ordered
% nedg :    the corresponding edge numbers
% tile :    polygon circumscribing a tile
%
% OUTPUT
% xyo,eo :  the nearest intersection on an edge other than ei
    
    % can't assume that the edge intersections were sorted
    ix = edge_sort(eixy,nedg);
    sixy = eixy(ix,:);
    sedg = nedg(ix);
    
    % find edges with 'legitimate' intersections, em < ei < ep
    em = ei - 1;
    if em == 0
        em = 4;
    end
    ep = ei + 1;
    if ep == 5
        ep = 1;
    end
    
    % intersection candidates
    ix = sedg==em;
    mixy = sixy(ix,:)(end,:); % last one on the preceeding edge
    ix = sedg==ep;
    pixy = sixy(ix,:)(1,:);   % first one on the following edge
    
    % pick the one that is nearer
    dist = edge_dist([mixy;pixy],[em;ep],xyi,ei,tile);
    if dist(1) < dist(2)
        xyo = mixy;
        eo = em;
    else
        xyo = pixy;
        eo = ep;
    end
end
