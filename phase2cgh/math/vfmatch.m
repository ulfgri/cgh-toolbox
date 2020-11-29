function [x,fval] = vfmatch(func,par,x1,x2,ftarget,opts)
%function [x,fval] = vfmatch(func,par,x1,x2,ftarget,opts)
%
% A vectorized form of the 'fmatch' function. Finds the locations
% along lines in a 2D domain where a user-defined function matches
% specified target values. The lines are defined by two points on each
% line and variables lambda which parametrize the positions along the
% lines:
%
%      x(lambda) = x1 + lambda .* (x2 - x1)
%
% This function seeks values of lambda such that
%
%      func(x(lambda),y(lambda),...) = ftarget
%
% INPUT
% func :     function handle of a user-defined function func(x,y,par)
% par :      structure with function parameters passed to func(x,y,par)
% x1,x2 :    Mx2 matrices with 2d coordinates [x,y] that define a line in the plane
%            and bracket the target value along the line connecting
%            x1 and x2. One coordinate pair per row.
% ftarget :  a vector with target values for the user defined function.
%            func(x(lambda),...) - ftarget lambda must have a zero crossing 
%            along the line connecting x1 and x2.
% opts :     options for 'vfzero'. These are set by the calling
%            function. See 'vfzero' help text for details.
%
% OUTPUT
% x :        a Mx2 matrix with locations along the lines from x1 to x2 where 
%            func(x,...) == ftarget
% fval :     a vector with function values at the final values of lambda
%
% NOTE: this function makes use of the 'vfzero' function that is
% provided by the 'optim' package in Octave.
    
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
%
% Initial version, Ulf Griesmann, NIST, July 2019

    % perform some input checking, but not much
    if nargin < 6
        error('function requires six input arguments.');
    end

    % call the 'vfzero' solver
    [lambda_final,fval,info] = ...
        vfzero(@(lambda)(func_line(func,par,x1,x2,lambda)-ftarget), ...
               repmat([0,1],length(ftarget),1), opts);

    if info <= 0
        fprintf('\nCall to ''vfzero'' returned with warning in function vfmatch.\n');
        fprintf('>>> info = %d\n', info);
    end

    % return the solution
    x = x1 + lambda_final .* (x2 - x1);
        
end

function fval = func_line(func,par,x1,x2,lambda)
% We need a wrapper for the phase function because anonymous 
% functions can only have one executable statement

    X = x1 + lambda .* (x2 - x1);
    fval = func(X(:,1),X(:,2),par);
    
end
