function [nm] = fdh_upto_nm(n, m)
% function [nm] = fdh_upto_nm(n, m)
%
% Returns a table with harmonic function terms
% up to radial order n and angular order m
%
% INPUT
% n :   highest radial order
% m :   highest angular order
%
% OUTPUT
% nm :  table with indices for all terms up to radial order n 
%       and angular order m (ready for feeding into fdh_* functions)

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
% 
% Version: 1.0
% Author: Ulf Griesmann; NIST; Sep 2006 
% Review: Ulf Griesmann; NIST; Sep 2006 
% Status: OK
% ---------------------------------------------------------

    if nargin < 2 
        error('missing argument(s)'); 
    end

    nm = zeros(1+n*(2*m+1),2);

    ic = 1;
    for i = 1:n
        ic = ic+1;
        nm(ic,:) = [i,0];
        for j = 1:m
            ic = ic+1;        
            nm(ic,:) = [i,-j];  
            ic = ic+1;        
            nm(ic,:) = [i,j];
        end
    end

end
