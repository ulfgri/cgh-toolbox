function disptime(tins)
%function disptime(tins)
%
% Displays a time interval as a well-formatted string
%
% INPUT
% tins :  time in seconds, e.g. the output of the 'toc' or
%         'cputime' funcions.
    
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Ulf Griesmann, NIST, June 2015

    tind = tins / 86400;  % convert to time in days
    fprintf('\nElapsed time: %s\n\n', datestr(tind, 'HH:MM:SS.FFF'));
    
end
