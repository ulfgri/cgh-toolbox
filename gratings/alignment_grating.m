function [carp] = alignment_grating(alpha, lambda, D, outline, aperture=[])
%function [carp] = alignment_grating(alpha, lambda, D, outline, aperture=[])
%
% Returns a Littrow type alignment grating for aligning a flat
% substrate in a collimated beam. The grating has 50% duty cycle.
%
% INPUT
% alpha   :  angle of incidence in DEGREES.
% lambda  :  wavelength of light.
% D       :  unit normal vector orthogonal to the grating lines 
%            describing the orientation of the grating in the x-y plane.
% outline :  a polygon or cell array of SIMPLE polygons describing the 
%            outline of the grating in the same units as argument 
%            'lambda'.
% aperture:  (optional) a structure array containing aperture polygons and a
%            corresponding Boolean operation that will be applied to the 
%            output polygons:
%            aperture(k).poly : clipping polygon(s); see 'polybool.m'
%            aperture(k).oper : clipping operation, 'and', 'or', etc.
%            This argument is useful e.g. to create grating outlines with 
%            holes
%
% OUTPUT
% carp :     a cell array of polygons describing the grating lines

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Short version: This software is in the PUBLIC DOMAIN.
% Ulf Griesmann, NIST, December 2012
% made much more efficient, U.G, August 2014

    % check arguments
    if nargin < 4
        error('function requires 4 arguments.');
    end

    if ~iscell(outline)
        outline = {outline};
    end

    if ~isempty(aperture) && ~isstruct(aperture)
        error('argument ''aperture'' must be a structure');
    end

    % make sure D points to the upper half plane
    if atan2(D(2),D(1)) < 0
        D = -D;
    end
    
    % enclosing rectangle
    A = acos(D(1));                          % angle between D and x-axis
    bbox = oriented_bbox(outline, pi/2 - A); % make D point North
    
    % calculate grating in rectangle
    gp = 0.5*lambda/sind(alpha);        % grating period
    hp = gp/2;                          % half pitch
    nlin = rounde( norm(bbox(4,:)-bbox(1,:)) / hp);
    incr = bsxfun(@times, hp*repmat([0:nlin-1]',1,2), D);
    L = incr + bbox(1,:); % points on the left
    R = incr + bbox(2,:); % points to the right
    
    % sort points into polygons
    G = zeros(2*size(L,1),2);
    G([1:4:size(G,1)],:) = L([1:2:size(L,1)],:);
    G([4:4:size(G,1)],:) = L([2:2:size(L,1)],:);
    G([2:4:size(G,1)],:) = R([1:2:size(L,1)],:);
    G([3:4:size(G,1)],:) = R([2:2:size(L,1)],:);
    
    % convert to cell array of polygons, 5 points per polygon
    carp = mat2cell(G, repmat(4,1,nlin/2), 2);

    % clip to outline
    carp = polybool(carp,outline,'and');
   
    % clip to aperture(s)
    if ~isempty(aperture)
        for k = 1:numel(aperture)
            carp = cellfun(@(P)polybool(P,aperture(k).poly,aperture(k).oper), ...
                           carp,'UniformOutput',false);
            carp = flatten_cells(carp); % eliminate any nested cell arrays
        end
    end
    
end


%-------------------------------------------------------------

function [B] = oriented_bbox(P, ang)
%
% calculates a bounding box in a user-defined orientation
%
% P :    polygon or cell array of polygons (nx2 matrices)
% ang :  rotation angle of bounding box in radians  
% B :    4x2 matrix of corner coordinates, one per row

    if ~iscell(P), P = {P}; end

    % rotate polygons by ang around [0,0]
    R = cell(size(P));
    for k = 1:length(P)
        R{k} = ctra_trotz(P{k}, ang, [0,0]); 
    end
    
    % calculate aligned bounding box
    abb = aligned_bbox(R);
    
    % rotate back
    B = ctra_trotz(abb, -ang, [0,0]);
   
end


%-------------------------------------------------------------

function [B] = aligned_bbox(P)
%
% calculates axis-aligned bounding box polygons
%
% P :  polygon or cell array of polygons (nx2 matrices)
% B :  4x2 matrix of corner coordinates, one per row

    if ~iscell(P), P = {P}; end
   
    % extreme coordinates
    X = zeros(length(P),2);
    Y = zeros(length(P),2);
    
    % find minimal and maximal coordinate values
    for k = 1:length(P)
        X(k,1) = min(P{k}(:,1));
        X(k,2) = max(P{k}(:,1));
        Y(k,1) = min(P{k}(:,2));
        Y(k,2) = max(P{k}(:,2));
    end
    
    % for all polygons
    xmin = min(X(:));
    xmax = max(X(:));
    ymin = min(Y(:));
    ymax = max(Y(:));
    
    % bounding box
    B = [xmin,ymin; xmax,ymin; xmax,ymax; xmin,ymax]; 
    
end


%-------------------------------------------------------------

function [mpos] = ctra_trotz(mpos,a,vcen)
%
% Rotate points around a specified center
%
% mpos   : position of the fiducials (x,y) in pixel coordinates.
% a      : angle of rotation
% vcen   : center of rotation
%

    % rotate
    mpos = bsxfun(@plus, rotz(bsxfun(@minus,mpos,vcen),a), vcen);

end


%-------------------------------------------------------------

function [rmp] = rotz(mp,az)
%function [rmp] = rotz(mp,az)
%
% Transformation of points by rotation around the z-axis
%
% INPUT
% mp   : n2x matrix with points to be transformed (one per row)
% az   : angle(s) around z in radians
%
% OUTPUT
% rmp : transformed points
%
% Simplified version of 'ctra_rotz' by Johannes Soons, NIST, Dec 2001

    if isempty(mp) | isempty(az) 
        error('missing argument(s).');
    end

    rmp = zeros(size(mp));
    rmp(:,1) = cos(az).*mp(:,1)-sin(az).*mp(:,2);
    rmp(:,2) = sin(az).*mp(:,1)+cos(az).*mp(:,2);

end
  

%-------------------------------------------------------------

function [enum] = rounde(num)
%function [enum] = rounde(num)
%
% Rounds the argument(s) to the nearest even number(s)
%
% num :   input argument matrix
% enum :  input rounded to the nearest even number(s)
%
% Ulf Griesmann, August 2014

    enum = 2 * floor(num/2 + 0.5);

end
