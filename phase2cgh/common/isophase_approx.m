function [aiso,apha,kiso] = isophase_approx(fphase,par,llc,urc,vphase,grid,verbose=false)
%function [aiso,apha,kiso] = isophase_approx(fphase,par,llc,urc,vphase,grid,verbose=false)
%
% Finds initial approximations of isophase lines in a tile area using the
% 'contourc' function of Octave.
%
% INPUT
% fphase :  function handle of phase function fphase(x,y,par)
% par :     structure with phase function parameters
% llc:      [x,y] row vector, coordinates of the lower left corner of the tile area
% urc:      [x,y] row vector, coordinates of the upper right corner of the tile area
% vphase:   (optional) boundaries for the phase values within the polygon area;
%           vphase(1) < vphase(2).
% grid:     (optional) vector [nx,ny] with initial number of phase function evaluations 
%           in x- and y-direction, or a scalar nx = ny = grid.
% verbose:  (optional) a logical; display progress information. Default is 'false'.
%
% OUTPUT
% aiso :    cell array with approximate isophase contours
% apha :    vector with corresponding phase values
% kiso :    a cell array with phase indices for the phase values in 'apha', kiso{1}
%           for vphase(1) and kiso{2} for vphase(2).
    
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
% Initial version: Ulf Griesmann, NIST, July 2019

    if verbose
        fprintf('   Calculating approximate isophases for tile\n');
        fprintf('     LLC: %.3f, %.3f\n', llc(1),llc(2));   
        fprintf('     URC: %.3f, %.3f\n', urc(1),urc(2));
    end
        
    % evaluate the phase function on the grid
    xx = linspace(llc(1),urc(1),grid(1));
    yy = linspace(llc(2),urc(2),grid(2));
    [X,Y] = meshgrid(xx,yy);
    P = fphase(X(:),Y(:),par);
    P = reshape(P,size(X));

    % determine the range phase indices and phase values in the tile area
    Pmin = min(P(:)); Pmax = max(P(:));
    
    % find k(m) such that 2*pi*k + vphase(m) > Pmin & < Pmax
    kmin = ceil (0.5 * (Pmin - vphase)/pi);
    kmax = floor(0.5 * (Pmax - vphase)/pi);
    K1 = [kmin(1):kmax(1)]';
    K2 = [kmin(2):kmax(2)]';
    plev = [2*pi*K1 + vphase(1); 2*pi*K2 + vphase(2)];  % all phase levels in range
    kiso = [K1;K2];                                     % corresponding k-values
    [plev,idx] = sort(plev);                            % sort phase levels
    kiso = kiso(idx);                                   % and phase indices
   
    % calculate contour lines at these phase levels
    [pcon,cph] = contourc(X,Y,P,plev);
    
    % parse out the approximate isophase contours in useable form
    % aiso{k} : Mx2 matrices with polygon vertices
    % apha(k) : phase value of aiso{k}
    k = m = 1; aiso = {}; apha = [];
    while m < size(pcon,2)
        aiso{k} = pcon(:,m+1:m+pcon(2,m))';
        apha(k)  = pcon(1,m);
        k = k + 1;
        m = m + 1 + pcon(2,m);
    end    

    % the marching squares algorithm used by the 'contourc' function produces
    % vertices that are spaced quite unevenly. Eliminate vertices that are very 
    % closely spaced. A fine grid may also produce too many vertices.
    % FIXME: the following does not remove all superfluous vertices
    gs = ceil(min(abs(xx(2)-xx(1)),abs(yy(2)-yy(1)))); % grid spacing
    for p = 1:numel(aiso)
        xy = aiso{p};
        dist = vecnorm(xy(2:end,:)-xy(1:end-1,:),2,2);             % distances
        keep = [dist > gs | logical(mod(1:length(dist),2)');true]; % always keep one
        aiso{p} = xy(keep,:);
    end
    
    if verbose
        fprintf('     %d approximate isophase curves found.\n', length(apha));
    end
end
