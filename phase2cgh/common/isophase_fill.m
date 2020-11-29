function [cap] = isophase_fill(fphase,par,ciso,eiso,vphase,tol,tile,verbose=false)
%function [cap] = isophase_fill(fphase,par,ciso,eiso,vphase,tol,tile,verbose=false)
%
% Fills (paints) a set of isophase approximation polygons by assembling them
% into simple polygons that enclose phase values in the semi-open intervals 
% [2*pi*k + vphase(1), 2*pi*k + vphase(1)) where k = ...-2,-1,0,1,2,...   
%
% INPUT
% fphase:    function handle to user defined function describing the 2D phase
%            distribution (in radians)
% par:       structure with parameters describing the phase distribution
%            (see user defined function describing phase distribution for details)
% ciso :     cell array with polygon approximations to isophase lines
% eiso :     nx2 matrix with pairs of edge numbers at the start and at the
%            end of polygons in 'ciso'. One row per polygon.
% vphase:    boundaries for the phase values within the polygon area;
%            vphase(1) < vphase(2).
% tol:       structure containing permissible error bounds. See 'cghparset'
% tile :     a closed polygon circumscribing a rectangular tile area. The
%            polygon is a 5x2 matrix with one vertex per row.
% verbose :  (optional) a logical. Display progress information when set to 'true'.
%            Default is 'false'.
%
% OUTPUT
% cap:       cell array with closed polygons
%
% NOTE: This function ONLY relies on the spatial (topological) relationship between
%       boundary curves. The fundamental premises of the boundary filling algorithm
%       are that phase boundary curves DO NOT CROSS and that the boundary-edge
%       intersections can be ordered on the tile edge.
    
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
%
% Ulf Griesmann, NIST, September 2018 - September 2019
% ---------------------------------------------------------

    % input checking
    if nargin < 7
        error('function requires 7 input arguments.');
    end

    if verbose
        fprintf('   Assembling %d isophase curves into closed polygons ...',numel(ciso));
        fflush(stdout());
    end

    % initialize output
    cap = {};
    pok = [];

    % handle the case of a tile without any edge-boundary intersections
    if isempty(ciso)

        % calculate phase at tile center
        tc = polycentr({tile});
        ptc = fphase(tc(1),tc(2),par);

        % check if tile is inside or outside phase interval
        if isinrange(ptc, vphase)
            cap{1} = tile;
        end

        return
    end

    % orient all polygons such that their initial vertex is closer to the tile
    % origin than the terminal vertex
    [ciso,eiso] = orient_polygons(ciso,eiso);

    % sort the isophase curves according to the distance of their
    % start vertices from the tile origin. Makes debugging easier.
    % Polygon end vertices are sorted along with the start vertices
    NP = numel(ciso);
    VS = zeros(NP,2);
    for k=1:numel(ciso)
        VS(k,:) = ciso{k}(1,:);    % start (initial) vertices
    end
    ivs = edge_sort(VS,eiso(:,1));
    VS = VS(ivs,:);                % needed below
    ciso = ciso(ivs);
    eiso = eiso(ivs,:);
    
    % sort ALL edge intersections for lookup of neighboring intersections
    VT = zeros(NP,2);
    for k=1:numel(ciso)
        VT(k,:) = ciso{k}(end,:);  % corresponding end (terminal) vertices
    end
    VA = vertcat(VS,VT);
    EA = vertcat(eiso(:,1),eiso(:,2));
    ias = edge_sort(VA,EA);
    VA = VA(ias,:);                % ALL intersections, sorted 
    
    % logical array to keep track of available (unconnected) polygons
    connected = false(numel(ciso),1);
    
    % index of current polygon in array of polygons
    icp = 0;

    % main loop
    while true
        
        if all(connected)
            break
        end

        % fetch next available (unconnected) polygon
        [pol,E,ip] = next_poly();

        % may be the last remaining polygon
        if sum(~connected) == 1

            pol = connect_last(pol,E);
            
            add_to_completed(pol);
            connected(ip) = true;
            break
        
        end

        % connect the polygon to the appropriate neighbor
        pol = connect_to_neighbor(pol,E,1);
        
        % and add the result to the output polygons
        add_to_completed(pol);
        connected(ip) = true;
        
    end % while

    % make sure all polygons have CCW orientation
    cap = polycw(cap,'ccw');
    
    if verbose
        fprintf(' Done\n');
        fprintf('   %d closed polygons\n\n',numel(cap));
    end

    
    % -----------------------------------------------
    % nested auxiliary functions for polygon handling
    % -----------------------------------------------

    % connects polygon with the neighbor such that phase is in [vphase(1),vphase(2))
    % ------------------------------------------------------------------------------
    function [polo] = connect_to_neighbor(pol,epol,crec)
        % pol :  polygon to be connected to next polygon
        % epol : edge numbers of polygon
        % polo : connected polygon
        % crec : counts the number of recursive calls to self
            
        if crec > tol.maxpol  % too many attempts to connect polygon
            fprintf('\ntoo many calls to isophase_fill>connect_to_neighbor\n');
            fprintf('possible error in isophase lines\n');
            error('isophase filling aborted');
        end
        
        % look up the neighbor polygons at the terminal end
        [xym,em,ivm,lsm, xyp,ep,ivp,lsp]= find_neighbors(pol(end,:));
        
        % check which neighbor should be connected to the current polygon
        xybm = point_between_vertices(xym,em,pol(end,:),epol(2),tile,tol.hderiv);
        xybp = point_between_vertices(pol(end,:),epol(2),xyp,ep,tile,tol.hderiv);
        
        % connect current polygon with the preceeding polygon ?
        if isinrange(fphase(xybm(1),xybm(2),par),vphase)
            if ivm == icp  % self
                polo = polyrev(pol);
                ecp = E([2,1]);  % edges of connected polygons
            else
                if lsm     % next vertex is start vertex
                    polo = vertcat(polyrev(ciso{ivm}), ...
                                   end_connector(eiso(ivm,1),epol(2),tile),polyrev(pol));
                    ecp = [eiso(ivm,2),epol(1)];
                else       % next vertex is terminal vertex
                    polo = vertcat(ciso{ivm},end_connector(eiso(ivm,2),epol(2),tile),...
                                  polyrev(pol));
                    ecp = [eiso(ivm,1),epol(1)];
                end
                connected(ivm) = true;
            end
            
        % connect current polygon with the next polygon ?
        elseif isinrange(fphase(xybp(1),xybp(2),par),vphase)        
            if ivp == icp % self
                polo = pol;
                ecp = E;
            else
                if lsp    % next vertex is start vertex
                    polo = vertcat(pol,end_connector(epol(2),eiso(ivp,1),tile), ...
                                   ciso{ivp});
                    ecp = [epol(1),eiso(ivp,2)];
                else      % next vertex is terminal vertex
                    polo = vertcat(pol,end_connector(epol(2),eiso(ivp,2),tile), ...
                                  polyrev(ciso{ivp}));
                    ecp = [epol(1),eiso(ivp,1)];
                end
                connected(ivp) = true;
            end
            
        else
            fprintf('\n');
            error('no polygon with correct phase range found');
        end
        
        % reverse combined polygon so that it starts on the lower edge
        polo = polyrev(polo);
        ecp = ecp([2,1]);
        
        % find the m-neighbor of the connected polygon's terminal end
        xym = find_neighbors(polo(end,:));
        
        % close polygon if start end is a neighbor of the terminal end
        % if not, find another polygon to connect
        if is_same(xym,polo(1,:))
            polo = vertcat(end_connector(ecp(1),ecp(2),tile),polyrev(polo));
            return
        else
            polo = connect_to_neighbor(polo,ecp,crec+1);
        end
        
    end
    
    % connect the last polygon
    % ------------------------
    function [polo] = connect_last(pol,epol)
        xyb = point_between_vertices(pol(end,:),epol(2),pol(1,:),epol(1),tile,tol.hderiv);
        if isinrange(fphase(xyb(1),xyb(2),par),vphase)
            polo = vertcat(pol,end_connector(epol(2),epol(1),tile));
        else
            polo = vertcat(polyrev(pol),end_connector(epol(1),epol(2),tile));
        end
    end

    % advance to next unconnected polygon and return it
    % -------------------------------------------------
    function [xy,E,ipo] = next_poly()
        if all(connected)
            xy = E = [];
            return
        end
        while true
            icp += 1;
            if icp > NP
                icp = 1;
            end
            if ~connected(icp)
                break
            end
        end
        xy = ciso{icp};
        E = eiso(icp,:);
        ipo = icp;
    end
    
    % adds polygon to output cell array
    % ---------------------------------
    function add_to_completed(pol)
        pol(end+1,:) = pol(1,:);
        cap{end+1} = pol;
    end

    % look up neighbors of a polygon vertex on the tile edge
    % ------------------------------------------------------
    function [xym,em,ivm,lsm, xyp,ep,ivp,lsp] = find_neighbors(xy)
        % xy :   a polygon vertex
        % xym :  preceeding polygon vertex
        % em :   edge number of xym
        % ivm :  index of preceeding polygon 
        % lsm :  logical, true if start vertex
        % xyp :  following polygon vertex
        % ep :   edge number of xyp
        % ivp :  index of following polygon
        % lsp :  logical, true if start vertex
        idx = find_vertex(VA,xy);          % look up vertex
        ix = min(idx) - 1;                 % find preceeding edge vertex
        if ~ix, ix = 2*NP; end
        xym = VA(ix,:);

        if nargout > 1
            ix = max(idx) + 1;             % find following edge vertex
            if ix > 2*NP, ix = 1; end
            xyp = VA(ix,:);
        
            [ivm,lsm,em] = find_tail_index(xym);        
            [ivp,lsp,ep] = find_tail_index(xyp);
        end
    end
    
    % find polygon index of an initial or terminal vertex point (tail)
    % ----------------------------------------------------------------
    function [idx,lsv,edg] = find_tail_index(xy)
        ix = find_vertex(VS,xy);      % try start vertices
        if ~isempty(ix)
            idx = ix;
        else                          % must be terminal
            idx = find_vertex(VT,xy);
        end
        if is_same(xy,ciso{idx}(1,:)) % start vertex
            lsv = true;
            edg = eiso(idx,1);
        else
            lsv = false;
            edg = eiso(idx,2);
        end
    end
    
end


%---------------------------------------------------------------------------------

function [idx] = find_vertex(V,xy)
%
% Finds the location of a specific vertex in a list of vertices
%
% INPUT
% V :   nx2 matrix with polygon vertices, one per row
% xy :  1x2 matrix with vertex coordinates
%
% OUTPUT
% idx : index for which P(idx,:) == xy

    if size(xy,1) > 1
        error('isophase curves have shared vertices');
    end
    
    idx = find(vecnorm(V-xy,2,2) < 2*eps);

end
        
%- Ende ------------------------------------------------------------------------
