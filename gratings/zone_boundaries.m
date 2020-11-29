function [rzb,dpdr] = zone_boundaries(phi,par,rl,ru,vph=[0,pi],tol=[],grd=10000)
%function [rzb,dpdr] = zone_boundaries(phi,par,rl,ru,vph=[0,pi],tol=[],grd=10000)
%
% Finds the radial zone boundaries phi(r) = 2*pi*k + pb, k = 0,1,2,3,
% ... for a given RADIAL phase function in the interval
% [rl,ru]. Optionally, the derivative dphi/dr will be calculated at
% the zone boundaries.
%
% INPUT
% phi :  function handle for a radial phase function. phi is of the
%        general form phi(x,y,par) and must have radial symmetry
% rl :   radius lower bound
% ru :   radius upper bound
% vph :  a 1x2 vector with phase boundaries. Default is [0,pi]
% tol :  tolerances for zone boundary finding
%           tol.vertex :  tolerance for zone boundary
%           tol.maxit  :  max number of iterations
%           tol.hderiv :  increment for derivative calculation
% grd :   radial sampling grid. Default is 10000
%
% OUTPUT
% rzb :    a column vector with radial zone boundaries
% dpdr :  (Optional) derivative dphi/dr|_{zb} at zone boundaries

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Short version: This software is in the PUBLIC DOMAIN.
% Ulf Griesmann, NIST, October 2012, updated October 2019

    pkg load optim
    
    % check arguments
    if nargin < 3
        error('expecting 3 arguments.');
    end
   
    if any([rl,ru] < 0)
        error('interval boundaries rl,ru must be positive.');
    end
    
    if isempty(tol)
        tol = struct('vertex',1e-5, 'maxit',100, 'hderiv',15e-3);
    end
    
    % set the parameters for fzero
    opts = optimset('TolX',tol.vertex, 'MaxIter',tol.maxit);
    
    % calculate zone boundaries along x-axis (could also use y-axis)
    bxy = isophase_line_xing(phi,par,[rl,0],[ru,0],vph,tol,grd);
    [rzb,sidx] = sort(bxy(:,1)); % y-values are all 0
    
    % calculate derivatives at boundaries
    if nargout > 1       
        dpdr = vphderiv(phi,par,bxy(sidx,:),tol.hderiv);
    end
    
end
