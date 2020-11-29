function [mjTerm,vScale] = zern_index_trans(mjTerm,source,target)
% function [mjTerm,vScale] = zern_index_trans(mjTerm,source,target)
%
% Translates a Zernike index from one definition to another
%
% mjTerm : index (each row contains the description of one term)
% vScale : scale factor to apply to target Zernike coefficient
%          to obtain the same Zernike polynomial value (normalization and sign)
%
% source : string indicating type of source index
% target : string indicating type of target index
%
% source and target can have the following values
%
% #ansi  : number of the term according to the ANSI Z80.28-2017 standard
% #ost   : number of the term according to Opical Shop Testing, 2nd edition
% #uoa   : number of the term according to the University of Arizona convention (Fringe)
% #noll  : number of the term according to the Noll convention (also used by Zemax)
% #oslo  : number of the term according to the OSLO convention (similar to
%          University of Arizona, but with clockwise angle convention
%          starting at Y+)
%  
% The relation between the various conventions is as follows
%
%             #ANSI   n    m    #ost   #uoa  #noll  #oslo
%
% Normalized    Y     Y    Y      N      N     Y      N
%
% Index         0     0    0      1      0     1      0
%               1     1   -1      3      2     3      1
%               2     1    1      2      1     2      2
%               3     2   -2      4      5     5      5
%               4     2    0      5      3     4      3
%               5     2    2   (-)6      4     6   (-)4
%               6     3   -3  (-)10     10     9   (-)9
%               7     3   -1      9      7     7      6
%               8     3    1      8      6     8      7
%               9     3    3   (-)7      9    10  (-)10
%              10     4   -4  (-)11     17    15  (-)17
%              11     4   -2     12     12    13     12
%              12     4    0     13      8    11      8
%              13     4    2  (-)14     11    12  (-)11 
%              14     4    4     15     16    14     16
%
% Where n equals the degree of the radial polynomial and m the angular
% frequency (exponential) parameter.  The normalization factor ensures
% that the RMS over the unit circle equals 1.
%
% References: 
%
% [1] ANSI Z8028-2017, (2017), "Ophthalmics - methods for reporting optical
%     aberrations of eyes," American National Standards Institute, New York, NY.
% [2] L.N. Thibos, R.A. Applegate, J.T. Schwiegerling, R. Webb R, 
%     Members OVST, (2002), "Standards for reporting the optical aberrations of
%     eyes," Journal of Refractive Surgery, 18:S652–S660.
% [3] D. Malacara (Ed), (1992), "Optical Shop Testing," John Wiley & Sons,
%     2nd Edition.
% [4] D. Malacara (Ed), (2007), "Optical Shop Testing," John Wiley & Sons,
%     3rd Edition.
% [5] R.J. Noll, (1976), "Zernike polynomials and atmospheric turbulence,"
%     Journal of the Optical Society of America, 66(3):207–211.
% [6] J. Loomis, (1978), "A computer program for analysis of interferometric
%     data," Optical Interferograms - Reduction and Interpretation,
%     ASTM STP 666 (American Society for Testing and Materials (ASTM)), pp 71–86.
% [7] J.C. Wyant, K. Creath, (1992), "Basic wavefront aberration theory for
%     optical metrology," Applied Optics and Optical Engineering, Vol. XI, pp 2–53.

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
%
% Any mention of commercial products is for information only; it does not imply
% recommendation or endorsement by NIST.
%
% Author: Johannes Soons; NIST; Feb 2004 - Dec 2019 

% check input and initialize scale

  if isempty(mjTerm)
    vScale = [];
    return;
  end

  if strcmpi(source,'nm')

    if size(mjTerm,2) ~= 2
      error('Each row of mjTerm must have 2 elements');
    end 
      
    vScale = ones(size(mjTerm,1),1);
  else
      
    if size(mjTerm,2) ~= 1
      error('Each row of mjTerm must have 1 element');
    end 
      
    vScale = ones(size(mjTerm,1),1);
  end

% check if work needs to be done

  if strcmpi(source,target), return; end
  
% translate source into nm

  switch(lower(source))
    case 'nm'
        
%     check for errors        

      if (any(isnan(mjTerm(:,1))) || any(mjTerm(:,1) < 0)) ||...
           any(fix(mjTerm(:,1)) ~= mjTerm(:,1))
        error(' n should be an integer, larger or equal to zero');
      end
      
      if (any(isnan(mjTerm(:,2))) || any(abs(mjTerm(:,2)) > mjTerm(:,1))) ||...
           any(fix(mjTerm(:,2)) ~= mjTerm(:,2))
        error(' m should be an integer, smaller or equal to n');
      end
      
      if any(rem(sum(abs(mjTerm),2),2) ~= 0)
        error(' The sum of the radial order n and angular frequency |m| should be an even number');
      end
      
%     set [n,m] values      
        
      nm = mjTerm;  

    case '#ansi'
        
%     check for errors        
        
      if (any(isnan(mjTerm)) || any(mjTerm < 0)) ||...
           any(fix(mjTerm) ~= mjTerm)
        error(' Zernike index should be an integer, larger or equal to zero');
      end
      
%     set [n,m] values 
%
%     Equation (5) of [2]
      
      nm(:,1) = ceil((-3+sqrt(9+8*mjTerm))/2);
      
%     Equation (6) of [2]      

      nm(:,2) = 2*mjTerm-nm(:,1).*(nm(:,1)+2);      

    case '#noll'
        
%     check for errors        
        
      if (any(isnan(mjTerm)) || any(mjTerm < 1)) ||...
           any(fix(mjTerm) ~= mjTerm)
        error(' Zernike index should be an integer, larger than zero');
      end
      
%     set [n,m] values

%     n can be obtained using Equation (5) of [2], modified for
%     the offset of the index value
        
      nm(:,1) = ceil((-3+sqrt(1+8*mjTerm))/2);
      
%     obtain the offset within the respective n group      
      
      nm(:,2) = (mjTerm-1)-(nm(:,1)+1).*nm(:,1)./2;
      
%     calculate |m| from the offset
%      
%     1) if n is even and the offset is odd: |m| = offset+1
%     2) if n is odd and the offset is even: |m| = offset+1      

      nm(:,2) = nm(:,2)+(rem(nm(:,1),2) ~= rem(nm(:,2),2)); 
      
%     set sign for m
%
%     1) if j is even: positive m (cosine term)
%     2) if j is odd: negative m (sinus term)
      
      nm(rem(mjTerm,2) ~= 0,2) = -nm(rem(mjTerm,2) ~= 0,2);
      
    case '#ost'
        
%     check for errors        
        
      if (any(isnan(mjTerm)) || any(mjTerm < 1)) ||...
           any(fix(mjTerm) ~= mjTerm)
        error(' Zernike index should be an integer, larger than zero');
      end
      
%     set [n,m] values

%     n can be obtained using Equation (5) of [2], modified for
%     the offset of the index value
        
      nm(:,1) = ceil((-3+sqrt(1+8*mjTerm))/2);
      
%     for #ost m > 0 corresponds to sin. However, for equal n
%     the sin terms (m > 0) are numbered first. Therefore
%     the ANSI numbering scheme holds
%
%     Equation (6) of [2], modified for offset index value      

      nm(:,2) = 2*(mjTerm-1)-nm(:,1).*(nm(:,1)+2);      
      
%     apply modified angle definition

      [nm,vScale] = zern_index_trans_angle(nm,vScale);
      
%     remove normalization from parameter value 
      
      vScale = vScale./zern_norm(nm(:,1),nm(:,2));
      
    case '#uoa'
 
%     check for errors        
        
      if (any(isnan(mjTerm)) || any(mjTerm < 0)) ||...
        any(fix(mjTerm) ~= mjTerm)
          error(' Zernike index should be an integer, larger or equal to zero');
      end
      
%     set group number n+|m|

      nm(:,1) = 2*ceil(-1+sqrt(1+mjTerm));
      
%     obtain offset relative to highest index in group

      nm(:,2) = 1/4*nm(:,1).^2+nm(:,1)-mjTerm;
      
%     set angular order (odd offset yields sin term)

      nm(rem(nm(:,2),2) ~= 0,2) = -nm(rem(nm(:,2),2) ~= 0,2);
      nm(:,2) = sign(nm(:,2)).*ceil(abs(nm(:,2))./2);
      
%     set radial order

      nm(:,1) = nm(:,1)-abs(nm(:,2));
      
%     remove normalization from parameter value 
      
      vScale = vScale./zern_norm(nm(:,1),nm(:,2));
    
    case '#oslo'
        
%     This convention is similar to that of UOA, except that the
%     angle is defined relative to the Y-axis in clockwise direction
 
%     check for errors        
        
      if (any(isnan(mjTerm)) || any(mjTerm < 0)) ||...
        any(fix(mjTerm) ~= mjTerm)
          error(' Zernike index should be an integer, larger or equal to zero');
      end
      
%     set group number n+|m|

      nm(:,1) = 2*ceil(-1+sqrt(1+mjTerm));
      
%     obtain offset relative to highest index in group

      nm(:,2) = 1/4*nm(:,1).^2+nm(:,1)-mjTerm;
      
%     set angular order (odd offset yields sin term)

      nm(rem(nm(:,2),2) ~= 0,2) = -nm(rem(nm(:,2),2) ~= 0,2);
      nm(:,2) = sign(nm(:,2)).*ceil(abs(nm(:,2))./2);
      
%     set radial order

      nm(:,1) = nm(:,1)-abs(nm(:,2));
      
%     apply modified angle definition

      [nm,vScale] = zern_index_trans_angle(nm,vScale);
      
%     remove normalization from parameter value 
      
      vScale = vScale./zern_norm(nm(:,1),nm(:,2));
    
    otherwise
      
      error('Unknown source type');
  end
  
% translate nm into target

  switch(lower(target))
    case 'nm'
      
      mjTerm = nm;

    case '#ansi'
        
%     Equation (4) from reference [2]        

      mjTerm = (nm(:,1).*(nm(:,1)+2)+nm(:,2))./2;
      
    case '#noll'
        
%     set index to number of terms used by previous orders
%     plus 1 for the starting index

      mjTerm = (nm(:,1).*(nm(:,1)+1))./2+1;
      
%     add maximum offset for angular order

      mjTerm = mjTerm+abs(nm(:,2));
      
%     add corrections
%
%     if m == 0                    : no correctiom
%     if mjTerm is even and m > 0  : no correction
%     if mjTerm is even and m < 0  : subtract 1 (cos -> sin)
%     if mjTerm is odd and m > 0   : subtract 1 (sin -> cos)
%     if mjTerm is odd and m < 0   : no correction

      mjTerm = mjTerm-(rem(mjTerm,2) == 0 & nm(:,2) < 0)-...
                      (rem(mjTerm,2) ~= 0 & nm(:,2) > 0);
     
    case '#ost'

%     apply normalization
      
      vScale = vScale.*zern_norm(nm(:,1),nm(:,2));
        
%     apply modified angle definition

      [nm,vScale] = zern_index_trans_angle(nm,vScale);
      
%     for #ost m > 0 corresponds to sin. However, for equal n
%     the sin terms (m > 0) are numbered first. Therefore
%     the ANSI numbering scheme holds

%     calculate index with offset of 1. Equation (4) from reference [2]        

      mjTerm = (nm(:,1).*(nm(:,1)+2)+nm(:,2))./2+1;

    case '#uoa'

%     apply normalization 
      
      vScale = vScale.*zern_norm(nm(:,1),nm(:,2));
      
%     set highest index for n+|m| group

      mjTerm = 1/4*(nm(:,1)+abs(nm(:,2))).^2+...
       (nm(:,1)+abs(nm(:,2)));
   
%     subtract maximum offset for angular order

      mjTerm = mjTerm-2*abs(nm(:,2));
      
%     add corrections
%
%     if m >= 0 : no correction
%     if m < 0  : add 1 (cos -> sin)
      
      mjTerm(nm(:,2) < 0) = mjTerm(nm(:,2) < 0)+1;

    case '#oslo'
        
%     This convention is similar to that of UOA, except that the
%     angle is defined relative to the Y-axis in clockwise direction

%     apply normalization 
      
      vScale = vScale.*zern_norm(nm(:,1),nm(:,2));
      
%     apply modified angle definition

      [nm,vScale] = zern_index_trans_angle(nm,vScale);
      
%     set highest index for n+|m| group

      mjTerm = 1/4*(nm(:,1)+abs(nm(:,2))).^2+...
       (nm(:,1)+abs(nm(:,2)));
   
%     subtract maximum offset for angular order

      mjTerm = mjTerm-2*abs(nm(:,2));
      
%     add corrections
%
%     if m >= 0 : no correction
%     if m < 0  : add 1 (cos -> sin)
      
      mjTerm(nm(:,2) < 0) = mjTerm(nm(:,2) < 0)+1;

    otherwise
      
      error('Unknown target type');
  end

function [nm2,vScale] = zern_index_trans_angle(nm,vScale)
% function [nm2,vScale] = zern_index_trans_angle(nm,vScale)
%
% zern_index_trans_angle: Translates nm index values from one
%   angle definition to another
%
% nm     : n and m numbers (one pair per row)
% vScale : scale factor to be applied to the target parameter (sign)
%
% angle' = pi/2-angle
%
  if isempty(vScale)
    vScale = ones(size(nm,1),1);
  end 
  
  nm2 = nm;
  
% sinus terms (m < 0)
%  
% sin(|m|*angle') = sin(|m|*(pi/2-angle))
%   = sin(|m|*pi/2-|m|*angle)
%   = sin(|m|*pi/2)*cos(|m|*angle)-cos(|m|*pi/2)*sin(|m|*angle)
%   = -sin(|m|*angle) for mod(|m|,4) == 0
%      cos(|m|*angle) for mod(|m|,4) == 1
%      sin(|m|*angle) for mod(|m|,4) == 2
%     -cos(|m|*angle) for mod(|m|,4) == 3

  vb = nm(:,2) < 0;
  if any(vb)
     
%   sin -> cos      
      
    vConvert = [1,-1,1,-1]';
    nm2(vb,2) = nm(vb,2).*vConvert(mod(abs(nm(vb,2)),4)+1);
    
%   sign conversion

    vConvert = [-1,1,1,-1]';
    vScale(vb) = vScale(vb).*vConvert(mod(abs(nm(vb,2)),4)+1);
  end
  
% cosine terms (m > 0)
%  
% cos(m*angle') = cos(m*(pi/2-angle))
%   = cos(m*pi/2-m*angle)
%   = cos(m*pi/2)*cos(m*angle)+sin(m*pi/2)*sin(m*angle)
%   =  cos(m*angle) for mod(m,4) == 0
%      sin(m*angle) for mod(m,4) == 1
%     -cos(m*angle) for mod(m,4) == 2
%     -sin(m*angle) for mod(m,4) == 3

  vb = nm(:,2) > 0;
  if any(vb)
     
%   cos -> sin      
      
    vConvert = [1,-1,1,-1]';
    nm2(vb,2) = nm(vb,2).*vConvert(mod(nm(vb,2),4)+1);
    
%   sign conversion

    vConvert = [1,1,-1,-1]';
    vScale(vb) = vScale(vb).*vConvert(mod(nm(vb,2),4)+1);
  end
