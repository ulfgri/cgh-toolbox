function [zTerm] = zern_symbolic(vjTerm,bCart,bDisplay)
% function [zTerm] = zern_symbolic(vjTerm,bCart,bDisplay)
%
% Generates the symbolic definition of an ANSI Z8028-2017
% Zernike polynomial term in either polar or cartesian coordinates.
%
% zTerm    : Symbolic variable with Zernike term 
%
% vjTerm   : ANSI Z8028-2017 index number(s) of the Zernike term
% bCart    : (optional, default 0) display the results in cartesian
%            coordinates normalized by the Zernike normalization radius) 
% bDisplay : (optional, default 1) whether to display the result
%
% See zern_eval for a definition of the Zernike terms
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
% Author: Johannes Soons; NIST; Nov 2001 - Dec 2019

% set defaults

  if nargin < 2, bCart = []; end
  if nargin < 3, bDisplay = []; end
  
  if isempty(bCart), bCart = 0; end
  if isempty(bDisplay), bDisplay = 1; end  
  
% initialize  

  syms term r theta x y;  
  
% extract n and m terms  
  
  nm = zern_index_trans(vjTerm,'#ansi','nm');  

  for i = 1:length(vjTerm)
%
    n = nm(i,1);
    m = nm(i,2);
%    
%   construct term
%
    term = 0;
    
%   radial polynomial    

    for k = 0:((n-abs(m))/2)
      p = (-1)^k * prod(1:(n-k))/(prod(1:k)*prod(1:((n+abs(m))/2-k))*prod(1:((n-abs(m))/2-k)));
      term = term + p * r^(n-2*k);
    end

%   normalization    
    
    term = term*sqrt((2*(n+1))/(1+(m==0)));      
    
%   sin and cos terms    
    
    if m < 0
      term = term * sin(abs(m)*theta);
    else
      term = term * cos(m*theta);
    end
    
%   convert to cartesian coordinates    

    if bCart
      term = expand(term);
      term = subs(term,cos(theta),x/r);
      term = subs(term,sin(theta),y/r);      
      term = subs(term,r,sqrt(x^2+y^2));      
    end

    term = simplify(term);
    
%   display    
     
    if bDisplay
      disp(sprintf(' Zernike term %g, n = %g, m = %g:',vjTerm(i),n,m));
      pretty(term);
      disp(' ');
    end
    
    if nargout == 1 
        zTerm(i) = term; 
    end
    
  end % for
