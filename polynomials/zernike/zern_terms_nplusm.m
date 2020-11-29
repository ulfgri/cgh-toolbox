function [vjTerm] = zern_terms_nplusm(nplusm)
% function [vjTerm] = zern_terms_nplusm(nplusm)
%
% Returns a vector with the ANSI Z80.28-2017 indices of the Zernike
% terms that have a specified sum of the radial order n and angular
% frequency |m|.
%
% vjTerm: vector with ANSI Z8028-2017 indices of the Zernike terms (not sorted)
%
% nplusm: value(s) for the n+|m| Zernike terms of interest
%
% see zern_eval for a definition of the Zernike terms
%
% EXAMPLE:
%
%   vjTerm = zern_terms_nplusm(2) : paraxial terms (tilt and defocus)
%   vjTerm = zern_terms_nplusm(4) : third order aberrations
%   vjTerm = zern_terms_nplusm(6) : fifth order aberrations
%

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
% Author: Johannes Soons; NIST; Dec 2019 

% initialize

  vjTerm = [];
  
  for i = 1:numel(nplusm)
      
%   check input      
      
    if fix(nplusm(i)) ~= nplusm(i) || ...
      (nplusm(i) < 0 || rem(nplusm(i),2) ~= 0)
  
      error(' The n plus |m| parameter should be an even number greater or equal to 0');
    end
    
%   obtain fringe/UoA index range

    if nplusm(i) == 0
      i1 = 0;
      i2 = 0;
    else
      i1 = 1/4*(nplusm(i)-2)^2+(nplusm(i)-2)+1;
      i2 = 1/4*nplusm(i)^2+nplusm(i);
    end
    
    vjTerm = [vjTerm;(i1:i2)'];
  end
  
% translate to ANSI index  
  
  vjTerm = zern_index_trans(vjTerm,'#uoa','#ansi');
  
% sort results

  vjTerm = sort(vjTerm);