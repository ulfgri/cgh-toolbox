function [con] = end_connector(Ea, Eb, tile)
%function [con] = end_connector(Ea, Eb, tile)
%
% Connects points on two edges of a tile via one or more tile
% corners in counter-clockwise (positive) direction
%
% INPUT
% Ea, Eb :  edge numbers of two polygon ends to be joined
% tile :    nx2 matrix with polygon circumscribing the tile area
%
% OUTPUT
% con :     nx2 matrix with connecting points (can be empty when Ea = Eb)
 
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
%
% Ulf Griesmann, NIST, August 2019

    % init output
    con = [];
    
    % if both edges are the same we are done
    if Ea == Eb
        return
    end

    % connect points through corners
    ecurr = Ea;
    while true
        
        % next edge
        enext = mod(ecurr,4) + 1;

        % add another corner between E1 and E2
        con = vertcat(con, tile(enext,:));

        % advance
        if enext == Eb
            break
        else
            ecurr = enext;
        end
        
    end    
end
