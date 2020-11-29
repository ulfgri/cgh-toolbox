function [vPar2] = zern_par_trans(vPar,source,target)
% function [vPar2] = zern_par_trans(vPar,source,target)
%
% Translates a Zernike parameter vector from one one definition to
% another. The returned parameter vector is padded with zeros.
%
% vPar2  : vector with target parameters
%
% vPar   : vector with source parameters
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

  vPar = vPar(:);
  
  nParOrig = numel(vPar);
  
% set index vector  

  switch lower(source)      
    case {'#ost','#noll'}
      vjTerm = (1:numel(vPar))';
    case {'#uoa','#ansi','#oslo'}
      vjTerm = (0:(numel(vPar)-1))';
    otherwise
      error(' Unrecognized source format')
  end

% translate index and obtain scale factor

  [vjTerm,vScale] = zern_index_trans(vjTerm,source,target);
  
% resort parameter vector and apply scaling

  switch lower(target)
    case {'#uoa','#ansi','#oslo'}
      vPar2 = zeros(max(vjTerm)+1,1);
      vPar2(vjTerm+1) = vPar.*vScale;
    case {'#ost','#noll'}      
      vPar2 = zeros(max(vjTerm),1);
      vPar2(vjTerm) = vPar.*vScale;
    otherwise
      error(' Unrecognized target format')
  end
  
% try to reduce size

  vi = find(vPar2(:) ~= 0,1,'last');
  vPar2 = vPar2(1:max(max(vi),nParOrig));  
