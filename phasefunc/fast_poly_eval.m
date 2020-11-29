function [y] = fast_poly_eval(a, x)
%function [y] = fast_poly_eval(a, x)
%
% Fast evaluation of a polynomial
%    
%         y = a1 + a2*x + a3*x^2 + ....
%
% INPUT
% a :  a vector with polynomial coefficients
% x :  a vector with points at which the polynomial is
%      evaluated
%
% OUTPUT
% y :  a vector with polynomial values

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Author: Ulf Griesmann, June 2013
