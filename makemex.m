%
% Octave script to make .mex files on Windows
%

% don't run on Linux
if isunix
    fprintf('\n>>> Use the makemex-octave shell script to compile mex functions\n\n');
    return
end

% check if we are running on MATLAB
if ~(exist('OCTAVE_VERSION')==5)

    fprintf('\n>>> This toolbox requires Octave ...\n\n');

else % good news: we are on Octave with gcc

    fprintf('Compiling mex functions for Octave on Windows ...\n');

    setenv('CFLAGS', '-O2 -fomit-frame-pointer -march=native -mtune=native');
    setenv('CXXFLAGS', '-O2 -fomit-frame-pointer -march=native -mtune=native');

    cd phasefunc
    mex fast_poly_eval.c
    
    cd ../minterp2
    mex minterp2.c

    cd ..
    fprintf('Done.\n');

end
