function [vPar2] = zern_mirror_line(vPar,rAngle)
% function [vPar2] = zern_mirror_line(vPar,rAngle)
%
% Calculates the effect on Zernike coefficients of a mirror operation
% at a line through the center of the Zernike unit disk.
%
% vPar2   : modified Zernike coefficients. The set of parameters is
%           typically expanded 
%
% vPar    : Zernike coefficients (see zern_eval for definition)
%           
% rAngle  : angle mirror line

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
% Author: Johannes Soons; NIST; Feb 2004 - Dec 2019 

% determine if there is anyting to do

  if isempty(vPar) || all(vPar == 0)
    vPar2 = vPar;
    return;
  end
      
% pad vPar to ensure symmetry of cos and sin terms

  nParOrig = numel(vPar);

  nm = zern_index_trans((0:(numel(vPar)-1))','#ansi','nm');
  
  mAdd = (nm(end,1)-nm(end,2))/2;

  if mAdd > 0
    vPar(end+(1:mAdd)) = 0;
    nm = [nm;[repmat(nm(end,1),mAdd,1),nm(end,2)+2*(1:mAdd)']];
  end  
  
% initialize

  vPar2 = zeros(size(vPar));
  
% mirror around angle a:
%
% th' = 2*a-th;
%  
% sin(m*(2*a-th)) = sin(2*m*a) cos(m*th) - cos(2*m*a) sin(m*th)
%                 = - cos(2*m*a) sin(m*th) + sin(2*m*a) cos(m*th)
%
% cos(m*(2*a-th)) = cos(2*m*a) cos(m*th) + sin(2*m*a) sin(m*th)
%
% process sin and cos terms

  vi = find(nm(:,2) ~= 0);
  
  for i = 1:length(vi)
      
%   modify original [sin,cos] term      
    
    vPar2(vi(i)) = vPar2(vi(i))+vPar(vi(i))*sign(nm(vi(i),2))*...
      cos(2*nm(vi(i),2)*rAngle);
  
%   add reverse [cos,sin] term  
    
    j = zern_index_trans([nm(vi(i),1),-nm(vi(i),2)],'nm','#ansi');
  
    vPar2(j+1) = vPar2(j+1)+vPar(vi(i))*sin(2*abs(nm(vi(i),2))*rAngle);
  end
    
% process radial terms  
  
  vi = find(nm(:,2) == 0);
  if ~isempty(vi)
    vPar2(vi) = vPar(vi);
  end
  
% try to reduce size

  vi = find(vPar2(:) ~= 0,1,'last');
  vPar2 = vPar2(1:max(max(vi),nParOrig));
