function [x,fval] = fmatch(func,par,x1,x2,ftarget,opts)
%function [x,fval] = fmatch(func,par,x1,x2,ftarget,opts)
%
% Finds the location along a line in a 2D domain where a user-defined
% function matches a specified target value. The line is defined by
% two points on the line and a variable, lambda, which parametrizes
% the position along the line:
%
%      x(lambda) = x1 + lambda * (x2 - x1)
%
% This function seeks a value of lambda such that
%
%      func(x(lambda),y(lambda),...) = ftarget
%
% INPUT
% func :     function handle of a user-defined function func(x,y,par)
% par :      structure with function parameters passed to func(x,y,par)
% x1,x2 :    2D coordinates [x,y] that define a line in the plane
%            and bracket the target value along the line connecting x1 and x2
% ftarget :  a target value for the user defined function.
%            func(x(lambda),...) - ftarget lambda must have a zero crossing 
%            along the line connecting x1 and x2.
% opts :     options for 'fzero'. These are set by the calling
%            function. See 'fzero' help text for details.
%
% OUTPUT
% x :        location along the line from x1 to x2 where 
%            func(x,...) == ftarget
% fval :     function value at the final value of lambda
    
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
% Initial version, Ulf Griesmann, NIST, September 2018

    % perform some input checking, but not much
    if nargin < 6
        error('function requires six input arguments.');
    end

    % call the 'fzero' solver
    [lambda_final,fval,info] = ...
        fzero(@(lambda)(func_line(func,par,x1,x2,lambda)-ftarget),[0,1],opts);

    if info <= 0
        fprintf('\nCall to ''fzero'' returned with warning in function fmatch.\n');
        fprintf('>>> info = %d\n', info);
        fprintf('>>> ftarget = %f\n', ftarget);
        fprintf('>>> fval = %f\n', fval);
        fprintf('>>> x1 = [%f,%f]; x2 = [%f,%f]\n\n', x1(1),x1(2),x2(1),x2(2));
    end

    % return the solution
    x = x1 + lambda_final * (x2 - x1);
        
end

function fval = func_line(func,par,x1,x2,lambda)
% We need a wrapper for the phase function because anonymous 
% functions can only have one executable statement in Octave/Matlab.

    x = x1 + lambda * (x2 - x1); % here x1,x2 are coordinate pairs
    fval = func(x(1),x(2),par);
    
end
