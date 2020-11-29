function polywrite(fname, cap);
% function polywrite(fname, cap);
% 
% Writes polygon vertex data to a binary file
%
% INPUT
% fname :  file name
% cap :    cell array of polygons

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Version: 1.0
% Author: Ulf Griesmann; NIST; Mar 2008
% Review: Ulf Griesmann; NIST; Mar 2008
% Status: OK

    % check arguments
    if nargin < 2
        error('missing argument(s)');
    end

    % open file
    pf = fopen(fname, 'w');
    if pf < 0
        error('polywrite :  could not open file.');
    end

    % write polygon data
    for k = 1:numel(cap)
        fwrite(pf, size(cap{k},1), 'int');
        fwrite(pf, cap{k}, 'double');
    end

    % close file
    fclose(pf);
    
end
