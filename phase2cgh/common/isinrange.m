function [inrng,k] = isinrange(pval, vphase)
%function [inrng,k] = isinrange(pval, vphase)
%
% Determines if 'pval' is in any of the semi-open intervals
% [2*pi*k + vphase(1),2*pi*k + vphase(1))
%
% INPUT
% pval :    a vector with phase value(s)
% vphase :  a 1x2 matrix with limits of a phase range [lowest,highest]
%
% OUTPUT
% inrng :   true if 2*pi*k + vphase(1) <= val <= 2*pi*k + vphase(2)
% k :       the values of k corresponding to 'pval'

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
%
% Ulf Griesmann, NIST, September 2018 - September 2019

    % calculate k
    k = floor(0.5*(pval - vphase(1)) / pi);
    
    % check range
    inrng = (pval > 2*pi*k + vphase(1)) & (pval < 2*pi*k + vphase(2));
    
end
