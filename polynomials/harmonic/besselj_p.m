function db = besselj_p(nu,z)
%function db = besselj_p(nu,z)
%
% First derivative of a Bessel function of the first kind
%
% nu: order of the Bessel function
% z : function argument
%
% see help besselj for details

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Version: 1.0
% Author: Johannes Soons; NIST; Sep 2006 
% Review: Johannes Soons; NIST; Sep 2006 
% Status: OK
% ---------------------------------------------------------

  db = 1/2*(besselj(nu-1,z)-besselj(nu+1,z));
  
end
