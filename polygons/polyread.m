function [cap] = polyread(fname);
% function [cap] = polyread(fname);
%
% Read a list of polygons that were stored with the polywrite
% function.
%
% INPUT
% fname :  name of file with polygon data
%
% OUTPUT   
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

    % open file
    pf = fopen(fname, 'r');
    if pf < 0
        error('polyread :  could not open file');
    end
    
    % read all polygons
    cap = {};
    while ~feof(pf)
        num = fread(pf, 1, 'int');
        cap{end+1} = fread(pf, [num,2], 'double');
    end
    
    % close file
    fclose(pf);

end
