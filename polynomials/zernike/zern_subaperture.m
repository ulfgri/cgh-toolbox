function [vPar2] = zern_subaperture(vPar,vTrans,scale,rAngle)
% function [vPar2] = zern_subaperture(vPar,vTrans,scale,rAngle)
%
% Transforms Zernike coefficients to those of a sub-aperture of the
% unit disk
%
% The transformed Zernike coefficients describe the original Zernike
% polynomial in a new coordinate system obtained by:
%
%  - translating the origin of the coordinate system (center of 
%    the Zernike normalization disk) by a translation vTrans
%    (normalized by the original Zernike normalization radius)
%  - rotating the translated coordinate frame around its origin 
%    by an angle rAngle
%  - scaling the x and y coordinates of the new coordinate system
%    by a factor scale
%
% In other words, the transformed Zernike coefficients describe
% the original Zernike polynomial in a subaperture with origin vTrans,
% radius scale*radius_original, and angle rAngle
%
% This corresponds to the following transformation of the
% Zernike polynomial relative to a fixed coordinate system
%
%  - translating the wavefront by translation -vTrans
%  - rotating the wavefront by an angle -rAngle around 
%    the center of the fixed coordinate system
%  - lateral scaling (stretching) of the polynomial relative to
%    the center of the fixed coordinate system by a
%    factor 1/scale
%
% see zern_transform for translating, rotating, and scaling the Zernike polynomial
%
% vPar2   : Transformed Zernike coefficients 
%
% vPar    : Zernike coefficients
% vTrans  : [x,y] displacement center (normalized by original pupil diameter)
% scale   : (optional, default 1) scaling factor (new pupil diameter divided by old one)
% rAngle  : (optional, default 0) rotation angle (postive in counter clockwise direction)
%
% The algorithm and code used is described in:
%
% Linda Lundström and Peter Unsbo, "Transformation of Zernike coefficients: scaled,
% translated, and rotated wavefronts with circular and elliptical pupils,"
% J. Opt. Soc. Am. A, Vol. 24, No. 3, 569-577, 2007.

% This software was developed by employees of the National Institute
% of Standards and Technology (NIST), an agency of the Federal Government,
% and is being made available as a public service. Pursuant to title 17
% United States Code Section 105, works of NIST employees are not subject
% to copyright protection in the United States.  This software may be subject
% to foreign copyright.  Permission in the United States and in foreign
% countries, to the extent that NIST may hold copyright, to use, copy,
% modify, create derivative works, and distribute this software and its
% documentation without fee is hereby granted on a non-exclusive basis,
% provided that this notice and disclaimer of warranty appears in all copies.
%
% THIS SOFTWARE IS BEING PROVIDED FOR RESEARCH AND EDUCATIONAL PURPOSES ONLY.
% THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
% EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY
% THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM FROM
% INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO THE
% SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE.  IN NO EVENT
% SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT,
% INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM,
% OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY,
% CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS
% OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE
% OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.

% check input

  if nargin < 2, vTrans = []; end
  if nargin < 3, scale = []; end
  if nargin < 4, rAngle = []; end
  
  if isempty(vTrans), vTrans = [0,0]; end
  if isempty(scale), scale = 1; end
  if isempty(rAngle), rAngle = 0; end
  
% check input

  if numel(vTrans) ~= 2
    error(' Translation vector should contain two values');
  end
  
  if isempty(vPar), vPar2 = vPar; return; end
  
  nParOrig = numel(vPar);
  
% express translation in polar coordinates

  [transAngle,transRadius] = cart2pol(vTrans(1),vTrans(2));
  
% expand parameters to cover all terms up to order n  
  
  jMax = length(vPar)-1;
  nMax = ceil((-3+sqrt(9+8*jMax))/2);
  jMax = nMax*(nMax+3)/2;
  
  vPar((end+1):(jMax+1)) = 0;

% Matrix P transforms from standard to Campbell order  
  
  P = zeros(jMax+1);
  
% Matrix N contains the normalization coefficients
 
  N = zeros(jMax+1);
  
% Matrix R contains the coefficients of the radial polynomials  
  
  R = zeros(jMax+1); 

% vParc is a complex representation of vPar

  vParc = zeros(jMax+1,1);
  
  counter = 1;
  
  for m = -nMax:nMax
      
%   Meridional indexes  

    for n = abs(m):2:nMax
        
%     radial indexes

%     set ANSI index of current (n,m) term

      jnm = (m+n*(n+2))/2;
      
      P(counter,jnm+1) = 1;
      N(counter,counter) = sqrt(n+1);

      for s=0:(n-abs(m))/2
          
        R(counter-s,counter)=(-1)^s*factorial(n-s) / (factorial(s)*factorial((n+m)/2-s)*...
          factorial((n-m)/2-s));
      end
      
      if m < 0
        vParc(jnm+1)= (vPar((-m+n*(n+2))/2+1)+1i*vPar(jnm+1)) /sqrt(2);
      elseif m == 0
        vParc(jnm+1) = vPar(jnm+1);
      else
        vParc(jnm+1) = (vPar(jnm+1)-1i*vPar((-m+n*(n+2))/2+1))/sqrt(2);
      end
      
      counter=counter+1;
    end
  end

% Coordinate-transfer matrix  
  
  ETA = [];
  
  for m = -nMax:nMax
    for n =abs(m):2:nMax
      ETA = [ETA, P*(transform(n,m,jMax,scale,transRadius,transAngle,rAngle))];
    end
  end

  vParc2 = (R*N*P)\(ETA*R*N*P*vParc);
  
% convert complex Zernike terms   
  
  vPar2 = zeros(jMax+1,1);
  
  for m = -nMax:nMax
    for n =abs(m):2:nMax
      jnm = (m+n*(n+2))/2;
      if m < 0
        vPar2(jnm+1) = imag(vParc2(jnm+1)-vParc2((-m+n*(n+2))/2+1))/sqrt(2);
      elseif m ==0
        vPar2(jnm+1) = real(vParc2(jnm+1));
      else    
        vPar2(jnm+1) = real(vParc2(jnm+1)+vParc2((-m+n*(n+2))/2+1))/sqrt(2);
      end
    end
  end
  
% try to reduce number of parameters

  vi = find(vPar2(:) ~= 0,1,'last');
  vPar2 = vPar2(1:max(max(vi),nParOrig));
  
function Eta = transform(n,m,jMax,scale,transRadius,transAngle,rAngle)

% Returns coefficients for transforming a ro^n*exp(i*m*theta)-term into '-terms

  Eta = zeros(jMax+1,1);

  for p = 0:((n+m)/2)
    for q= 0:((n-m)/2)
      nnew = n-p-q;
      mnew = m-p+q;
      jnm = (mnew+nnew*(nnew+2))/2;
      
      Eta(floor(jnm+1)) = Eta(floor(jnm+1)) + nchoosek((n+m)/2,p)*...
        nchoosek((n-m)/2,q) * scale^(n-p-q)*transRadius^(p+q)*exp(1i*((p-q)*(transAngle-rAngle)+m*rAngle));
    end
  end
