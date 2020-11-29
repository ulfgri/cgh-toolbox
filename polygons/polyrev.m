function [xyo] = polyrev(xyi);
% function [xyo] = polyrev(xyi);
%
% Reverses the order of the vertices in a polygon
%
% INPUT
% xyi :  input polygon; a nx2 matrix with polygon vertices, one per row
%
% OUTPUT   
% xyo :  output polygon; same shape as input polygon

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Ulf GRIESMANN, NIST, August 2019
    
    xyo = xyi(end:-1:1,:);

end
