function openblas_threads(n)
%function openblas_threads(n)
%
% Checks the value of the environment variable OPENBLAS_NUM_THREADS.   
% The function reads the environment variable and compares its value
% to the function argument. A mismatch triggers an error.
% 
% n :  (OPTIONAL) required number of threads. A value of 0 means that
%      OPENBLAS_NUM_THREADS is expected to be unset. If the function
%      is called without argument, the value of OPENBLAS_NUM_THREADS
%      is displayed.
% 
% NOTE: for parallel processing applications using the 'parcellfun'
%       function, the number of threads should be set to 1 to 
%       avoid thread contention.

% Ulf Griesmann, NIST, September 2014

    % do nothing if we are running on Windows
    if ispc
        return
    end
    
    % get environment
    ont = getenv('OPENBLAS_NUM_THREADS');
    if ~isempty(ont)
        ont = str2num(ont);
    end

    % check argument 
    if nargin < 1
        
        if isempty(ont)
            fprintf('OPENBLAS_NUM_THREADS not found\n');
        else
            fprintf('OPENBLAS_NUM_THREADS = %d\n', ont);
        end

    else
       
        if n == 0
            if ~isempty(ont)
                error('OPENBLAS_NUM_THREADS should not be set\n');
            end
        else
            if n ~= ont
                error(sprintf('OPENBLAS_NUM_THREADS = %d; required %d threads\n', ont, n));
            end
        end
        
    end

end
