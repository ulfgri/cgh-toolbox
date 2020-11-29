function [vph] = phase_poly2(vx,vy,spar)
% function [vph] = phase_poly2(vx,vy,spar)
%
% Returns the spatial phase distribution (in radians) modeled as a
% bivariate polynomial
%
% INPUT
% vx:      desired x coordinates
% vy:      desired y coordinates
% spar:    structure with parameters describing the phase distribution function
%            spar.mpol : polynomial parameters (index 1 is x power, starting
%                        from zero, index 2 is y power, starting from zero)
%
%             E.g.:      spar.mpol = [1,2,3;4,5,6;7,8,9] corresponds to
%                        vph = 1 +     2*y +     3*y^2 + ...
%                              4*x +   5*x*y +   6*x*y^2 + ...
%                              7*x^2 + 8*x^2*y + 9*x^2*y^2
%
% OUTPUT
% vph:     phase field in radians

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Version: 1.0
% Author: Johannes Soons; NIST; Feb 2008 - Sep 2011
% Review: Johannes Soons; NIST; Sep 2011 
% Status: OK
% ---------------------------------------------------------

  vph = polyval2d(vx,vy,spar.mpol);

end


%----------------------------------------------------------------------------

function [vp] = polyval2d(vx,vy,pc)
% function [vp] = polyval2d(vx,vy,pc)
%
% Evaluates a bivariate polynomial
%
% INPUT
% vx:      desired x coordinates
% vy:      desired y coordinates
% pc:      polynomial coefficients (index 1 is x power, starting
%          from zero, index 2 is y power, starting from zero)
% OUTPUT
% vp:      vector with polynomial values
%
% EXAMPLE
%   pc = [1,2,3;4,5,6;7,8,9] corresponds to
%     vp = 1 +     2*y +     3*y^2 + ...
%          4*x +   5*x*y +   6*x*y^2 + ...
%          7*x^2 + 8*x^2*y + 9*x^2*y^2
%  

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Version: 1.0
% Author: Johannes Soons; NIST; Feb 2008 - Sep 2011
% Review: Johannes Soons; NIST; Sep 2011 
% Status: OK
% ---------------------------------------------------------

  if isrow(vx), vx = vx'; end
  if isrow(vy), vy = vy'; end

  vpx = ones(length(vx),1);   

  nx = size(pc,1);

  vp = zeros(max(length(vx),length(vy)),1);  
  
  for i = 1:1:nx

    ny = max([1,max(find(pc(i,:) ~= 0))]);
     
    if ny <= 1
      vpy = repmat(pc(i,ny),length(vy),1);
    else
      vpy = pc(i,ny)*vy+pc(i,ny-1);
      for j = (ny-2):-1:1,
        vpy = vpy.*vy+pc(i,j);
      end
    end  
  
    vp = vp+vpy.*vpx;
    
    if i < nx, vpx = vpx.*vx; end
  end

end
