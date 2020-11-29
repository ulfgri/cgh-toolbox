function [vjTerm] = zern_terms_symmetry(vjTerm,property)
% function [vjTerm] = zern_terms_symmetry(vjTerm,property)
%
% Selects from a vector with the ANSI Z8028-2017 indices for Zernike
% terms those with specified symmetry properties.
%
% See zern_eval for a definition of the Zernike terms.
%
% vjTerm   : vector with ANSI indices of Zernike terms
% property : desired property
%
%   'asym_x'     : asymmetric in X
%   'sym_x'      : symmetric in X
%   'asym_y'     : asymmetric in Y
%   'sym_y'      : symmetric in Y
%   'asym_rot'   : rotation asymmetric
%   'sym_rot'    : rotation symmetric
%   'asym_point' : point asymmetric
%   'sym_point'  : point symmetric

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
% Author: Ulf Griesmann & Johannes Soons; NIST; Feb 2004 - Dec 2019 

  if isempty(vjTerm), return; end

% get nm terms  
  
  nm = zern_index_trans(vjTerm,'#ansi','nm');
  
% extract indices

  switch lower(property)
    case 'asym_x'
        
%     asymmetric terms in X

%     sin(k*theta) = symmetric around Y if k = odd
%     cos(k*theta) = symmetric around Y if k = even

      nm = zern_index_trans(vjTerm,'#ansi','nm');
      vb = (nm(:,2) < 0 & rem(nm(:,2),2) == 0) |...
           (nm(:,2) > 0 & rem(nm(:,2),2) ~= 0);

    case 'sym_x'
        
%     symmetric terms in X

%     sin(k*theta) = symmetric around Y if k = odd
%     cos(k*theta) = symmetric around Y if k = even

      vb = (nm(:,2) < 0 & rem(nm(:,2),2) ~= 0) |...
           (nm(:,2) >= 0 & rem(nm(:,2),2) == 0);
       
    case 'asym_y'
        
%     asymmetric terms in Y

%     sin(k*theta) = asymmetric around X
%     cos(k*theta) = symmetric around X

      nm = zern_index_trans(vjTerm,'#ansi','nm');
      vb = nm(:,2) < 0;

    case 'sym_y'
        
%     symmetric terms in X

%     sin(k*theta) = asymmetric around X
%     cos(k*theta) = symmetric around X

      vb = nm(:,2) >= 0;
       
    case 'asym_rot'
        
%     rotation asymmetric terms

      vb = nm(:,2) ~= 0;
      
    case 'sym_rot'
        
%     rotation symmetric terms

      vb = nm(:,2) == 0;
      
    case 'asym_point'
        
%     asymmetric point terms

%     sin(k*theta) = symmetric if k = even
%     cos(k*theta) = symmetric if k = even

      nm = zern_index_trans(vjTerm,'#ansi','nm');
      vb = (nm(:,2) < 0 & rem(nm(:,2),2) ~= 0) |...
           (nm(:,2) > 0 & rem(nm(:,2),2) ~= 0);

    case 'sym_point'
        
%     symmetric terms in X

%     sin(k*theta) = symmetric if k = even
%     cos(k*theta) = symmetric if k = even

      vb = (nm(:,2) < 0 & rem(nm(:,2),2) == 0) |...
           (nm(:,2) >= 0 & rem(nm(:,2),2) == 0);
      
    otherwise
      error('Invalid property');
  end

  vjTerm = vjTerm(vb);