function [vjTerm] = zern_terms(varargin)
% function [vjTerm] = zern_terms(varargin)
%
% Returns a vector with the ANSI indices of the Zernike terms
% corresponding to the list of argument strings.
%
% The following argument strings are recognized:
%
% 'constant', 'const' or  'con'   :  constant (or piston) term
% 'piston'   or  'pis'            :  constant (or piston) term
% 'tilt'     or  'tlt'            :  the tilt terms
% 'power'    or  'pwr'            :  the focus shift (power) term
% 'astigmatism, 'astig' or  'ast' :  2nd order astigmatism terms
% 'coma'     or  'com'            :  3rd order coma terms
% 'spherical, 'spher' or 'sph'    :  the primary spherical aberration term
% 'trefoil'  or  'tre'            :  the trefoil terms
%
% See zern_eval for a definition of the Zernike terms.
%
% vjTerm: vector with ANSI Z80.28-2017 indices of the Zernike terms (sorted)
%
% EXAMPLE:
%
%      znum = zern_terms('pwr','con','tlt','sph');
%
%      returns a vector with indices of the Zernike terms for
%      constant, tilt, focus shift (power), and spherical aberration

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

% initialize

  vjTerm = [];

% translate arguments

  for k=1:nargin
    switch lower(varargin{k})
     
      case {'con', 'const', 'piston', 'pis', 'constant'}
         vjTerm(end+1) = 0;
	
      case {'tlt', 'tilt'}
         vjTerm(end+1) = 1;
         vjTerm(end+1) = 2;
	
      case {'pwr', 'power'}
         vjTerm(end+1) = 4;
       
      case {'ast', 'astig', 'astigmatism'}
         vjTerm(end+1) = 3;
         vjTerm(end+1) = 5;
       
      case {'com', 'coma'}
         vjTerm(end+1) = 7;
         vjTerm(end+1) = 8;
       
      case {'sph', 'spher', 'spherical'}
         vjTerm(end+1) = 12;

      case {'tre', 'trefoil'}
         vjTerm(end+1) = 6;
         vjTerm(end+1) = 9;         
       
      otherwise
         error([varargin{k},' is not a valid description for a Zernike term.'])
    end
  end
     
% return sorted result
  
  vjTerm = sort(vjTerm);