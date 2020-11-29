function [cmp] = circle2polygon(mcen,vrad,rtol,vpos1,vpos2,bverbose)
% function [cmp] = circle2polygon(mcen,vrad,rtol,vpos1,vpos2,bverbose)
%
% Generates a set of polygons that approximate a set of circular
% apertures (holes). If desired, the polygons describe the
% intersection with a rectangular domain.
%
% mcen:      [x,y] coordinates of each circular aperture
% vrad:      radius of each aperture
% rtol:      (optional, default 1E-3) maximum allowed distance between polygon
%            vertex and circular aperture
% vpos1:     (optional, default []) [x,y] coordinates of the lower left corner
%            of the rectangular domain
% vpos2:     (optional, default []) [x,y] coordinates of the upper right corner
%            of the rectangular domain
% cmp:       cell array with polygons
% bverbose:  (optional, default 0) display progress
%
% NOTE: It is assumed that there is no overlap between the circles

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Author: Johannes Soons; NIST; Oct 2012 
% Review: Johannes Soons; NIST; Oct 2012 
% Status: TESTING
% ---------------------------------------------------------

  if nargin < 6, bverbose = []; end
  if nargin < 5, vpos2 = []; end
  if nargin < 4, vpos1 = []; end  
  if nargin < 3, rtol = []; end
  if isempty(rtol), rtol = 1E-3; end
  if isempty(bverbose), bverbose = 0; end
  
  if size(mcen,1) ~= length(vrad)
    error('Number of centers does not match number of radii');
  end
  
  if isempty(vpos1) ~= isempty(vpos2)
    if isempty(vpos1)
      vpos1 = [-1E99,-1E99];
    else
      vpos2 = [1E99,1E99];
    end
  end

% switch algorithm depending on whether clipping is required

  if isempty(vpos1)

    cmp = cell(size(mcen,1),1);          
      
    for ic = 1:size(mcen,1)
      cmp{ic} = gen_polygon(mcen(ic,:),vrad(ic),rtol);

      if bverbose == 1 && mod(ic,10000) == 0
        disp(sprintf(' Number of polygons: %d',ic));
      end
    end
    
  else
      
%   generate domain

    mdom = [vpos1(1),vpos1(2);...
            vpos2(1),vpos1(2);...
            vpos2(1),vpos2(2);...
            vpos1(1),vpos2(2)];
        
    cmp = {};

    for ic = 1:size(mcen,1)
        
%     check if we need to generate the circle

      if (mcen(ic,1)+vrad(ic) > vpos1(1) &&...
          mcen(ic,1)-vrad(ic) < vpos2(1)) &&...
         (mcen(ic,2)+vrad(ic) > vpos1(2) &&...
          mcen(ic,2)-vrad(ic) < vpos2(2))

        mpol = gen_polygon(mcen(ic,:),vrad(ic),rtol);
      
%       determine whether polygon needs to be clipped

        if (mcen(ic,1)+vrad(ic) > vpos2(1) ||...
            mcen(ic,1)-vrad(ic) < vpos1(1)) ||...
           (mcen(ic,2)+vrad(ic) > vpos2(2) ||...
            mcen(ic,2)-vrad(ic) < vpos1(2))
      
          cmpi = polybool(mpol,mdom,'and');
          
          if ~isempty(cmpi)
            cmp(end+1) = cmpi;
          end
        else
          cmp{end+1} = mpol;
        end
    
        if bverbose == 1 && mod(ic,10000) == 0
          disp(sprintf(' Number of polygons: %d',ic));
        end
      end      
    end
  end
end
  
function mpol = gen_polygon(vcen,rrad,rtol)
%
% generate polygon
% 
% vcen : center
% rrad : radius
% rtol : tolerance value

% generate angular interval

  da = sqrt(6*(1-rrad^2/(rtol+rrad)^2));
  da = 2*pi/ceil(2*pi/da);
    
% adjust radius to ensure equal area

  rrad = sqrt(rrad^2*da/sin(da));

% generate polygon

  va = (0:da:(2*pi-da-eps))'+rand(1)*da;
  mpol = [vcen(1)+rrad*cos(va),vcen(2)+rrad*sin(va)];
    
end
