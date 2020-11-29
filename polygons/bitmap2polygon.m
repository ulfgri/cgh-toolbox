function [cmp] = bitmap2polygon(map,rtol,bedge,bopt)
% function [cmp] = bitmap2polygon(map,rtol,bedge,bopt)
%
% Convert outline(s) of a binary image into a polygon
%
% map     : binary image
% rtol    : (optional, default 1) maximum allowed distance of image edge to polygon
% bedge   : (optional, default 1) enforce map edge to be represented without error
% bopt    : (optional, default 0) find the polygon with the fewest edges given a
%           starting point (time intensive)%                   
% cmp     : cell array with polygons

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Author: Johannes Soons; NIST; Feb 2008 
% Review: Johannes Soons; NIST; Feb 2008 
% Status: OK
% ---------------------------------------------------------

  if nargin < 4, bopt = []; end
  if nargin < 3, bedge = []; end
  if nargin < 2, rtol = []; end
  
  if isempty(bopt), bopt = 0; end
  if isempty(bedge), bedge = 1; end
  if isempty(rtol), rtol = 1; end

% identify pixels on edge

  mk = 1/8*[-1,-1,-1;-1,8,-1;-1,-1,-1];
  mape = conv2(map,mk,'same');
  mape = (mape > 0);

  ip = 0;
  
  [nr,nc] = size(mape);

  while any(mape(:))
    
%   walk along the edge of a new polygon

    vi = find(mape(:) == 1);
    [vpe(1),vpe(2)] = ind2sub(size(mape),vi(1));
    mpe = [];
    while ~isnan(vpe(1,1)),
      mpe(end+1,:) = vpe;
      vpe = nextpixel(mape,mpe(end,:));
      mape(mpe(end,1),mpe(end,2)) = 0;
    end
    
    mpe(:,3) = rtol;
    
    if bedge == 1
      mpe((mpe(:,1) == 1 | mpe(:,1) == nr) |...
          (mpe(:,2) == 1 | mpe(:,2) == nc),3) = eps;
    end
    
    if size(mpe,1) > 2

%     start new polygon
      
      ip = ip+1;
      
      if bopt == 0
      
%       identify edges through sequential local optimization 

        ic = 1;
        ie = 0;

        while ic < (size(mpe,1)-1)

%         start new edge

          ie = ie+1;
          ics = ic;
      
          cmp{ip}(ie,:) = mpe(ic,1:2);
      
%         next pixel has to be valid

          ic = ic+1;
      
          rmax = 0;
          rtolval = mpe(ic,3);

          while rmax <= rtolval & ic < size(mpe,1)

%           try to update solution

            ic = ic+1;

            vdir = vec_uni(bsxfun(@plus,mpe(ic,1:2),-cmp{ip}(ie,:)));
        
            rmax = max(abs(-vdir(2)*(mpe((ics+1):(ic-1),1)-cmp{ip}(ie,1))+...
                            vdir(1)*(mpe((ics+1):(ic-1),2)-cmp{ip}(ie,2))));
            
            if mpe(ic,3) < rtolval, rtolval = mpe(ic,3); end            
          end
          ic = ic-1; 
        end
      
      else

%       identify edges through global optimization 

%       1) obtain range for every vertex

        np = size(mpe,1);
        ns = 1;
        mpe(:,4) = ns;
        rmax = 0;
        vi = (1:np)';

        while ~isempty(vi),
          ns = ns+1;
          mdir(vi,:) = vec_uni(bsxfun(@plus,mpe(mod(vi+ns-1,np)+1,1:2),-mpe(vi,1:2)));

          it = 1;  
          while it < ns-1 & ~isempty(vi),
            vi = vi(abs(-mdir(vi,2).*(mpe(mod(vi+it-1,np)+1,1)-mpe(vi,1))+...
                         mdir(vi,1).*(mpe(mod(vi+it-1,np)+1,2)-mpe(vi,2))) <= mpe(vi,3));
            it = it+1;
          end
      
          if ~isempty(vi),
            mpe(vi,4) = ns;
          end
        end

%       2) find shortest path

        vd = zeros(np+1,1)+inf;
        vp = zeros(np+1,1)+NaN;

        vd(1) = 0;

        for i = 1:np,
          for j = 1:mpe(i,4),
            if (i+j <= np+1) & vd(i+j) > vd(i)+1,
              vd(i+j) = vd(i)+1;
              vp(i+j) = i;
            end
          end
        end

%       extract optimal sequence

        vi = np+1;

        while vi(end) > 1,
          vi(end+1) = vp(vi(end));
        end

        cmp{ip} = mpe(vi(end:-1:2),1:2);
      end    

%     convert to pixel coordinates

      cmp{ip} = trans_sub2pix(size(mape),cmp{ip});

%     close contour

      cmp{ip}(end+1,:) = cmp{ip}(1,:);
    end
  end  
end


%--------------------------------------------------------------------

function [vnij] = nextpixel(map,vij)
% function [vnij] = nextpixel(map,vij)
%
% nextpixel: Find a valid neighbouring pixel in coordinate directions
%
  [nr,nc] = size(map);

  if vij(2) < nc & map(vij(1),vij(2)+1) == 1, vnij = vij+[0,1]; return; end
  if vij(1) > 1  & map(vij(1)-1,vij(2)) == 1, vnij = vij+[-1,0]; return; end
  if vij(2) > 1  & map(vij(1),vij(2)-1) == 1, vnij = vij+[0,-1]; return; end
  if vij(1) < nr & map(vij(1)+1,vij(2)) == 1, vnij = vij+[1,0]; return; end

  vnij = [NaN,NaN];

end


%--------------------------------------------------------------------

function [mv] = vec_uni(mv)
%function [mv] = vec_uni(mv)
%
% vec_uni: Returns the unit direction vector of each row

    if isempty(mv),
        return
    end

    vl = vecnorm(mv,2,2);

    if ~isempty(find(vl < eps)) 
        error('Vector with zero length'); 
    end

    mv = bsxfun(@times,mv,1./vl);

end


%--------------------------------------------------------------------

function [mpix] = trans_sub2pix(vsize,mij,mjj)
%function [mpix] = trans_sub2pix(vsize,mij,mjj)
%
% trans_sub2pix: convert matrix row and column numbers of matrix element A(i,j)
%                to the pixel coordinates (x,y)
%
% mpix  : [x,y] pixel coordinates
% vsize : size of the matrix (or the number of rows)
% mij   : the indices (i,j) of the corresponding matrix element A(i,j)
% mjj   : if specified, mjj contains the column indices and mij the row indices
%
% A note on coordinate systems:
% 
% Image data are displayed in the same orientation as used by e.g.
% NASA's 'fv' software to display FITS image files. 
%
% For an NxM matrix the matrix elements (i,j) are displayed at
% coordinates (x,y) where :
%
%    x = j  (!) 
%    y = N - i + 1
%
% Likewise, a pixel at coordinates (x,y) represents a matrix
% element A(i,j) with indices
%
%    j = x
%    i = N - y + 1
%
% Note that array element (N,1) has pixel coordinates (1,1) and will
% be displayed at the bottom left of the screen. Array element (N,M)
% will be displayed at the bottom right of the screen

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic. 
%
% Author: Johannes Soons & Ulf Griesmann; NIST; Mar 2003 
% Review: Johannes Soons; NIST; Jun 2006 
% Status: OK
% ---------------------------------------------------------

    if nargin < 2,
        error('missing argument(s)');
    end
  
    if isempty(mij)
        mpix = [];
        return
    end
  
    if nargin == 2
        mpix = [mij(:,2),vsize(1) + 1 - mij(:,1)];
    else
        mpix = [mjj,vsize(1) + 1 - mij];
    end
    
end
