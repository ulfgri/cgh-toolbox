function cgh_version()
%function cgh_version()
%
% Displays the version of the CGH toolbox and the version
% of the interpreter running the software

% Ulf Griesmann, NIST, October 2019

    % toolbox version and version date
    tb_version = '177';
    tb_date = '2020 November 27';
    
    ltb_ver = [tb_version, '  (', tb_date, ')'];
    if exist('OCTAVE_VERSION')
        interpreter = 'Octave';
    else
        interpreter = 'Matlab';
    end
    lsw_ver = [interpreter,' ',version()];
    
    fprintf('\n-----------------------------------------------\n');
    fprintf('Interpreter version : %s\n', lsw_ver);
    fprintf('CGH Toolbox version : %s\n', ltb_ver);
    fprintf('-----------------------------------------------\n\n');

end
