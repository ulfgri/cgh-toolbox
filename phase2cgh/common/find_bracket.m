function [x1,x2] = find_bracket(fphase,par,x,piso,B,tol)
%function [x1,x2] = find_bracket(fphase,par,x,piso,B,tol)
%
% Finds bracketings of isophase curves along the direction defined by
% vectors D. Brackets are [x-k*D,x+k*D], k = [1.5,2:M].^1.7
%
% INPUT
% fphase :   function handle of phase function, fphase(x,y,par)
% par :      phase function parameters
% x :        nx2 matrix with estimated coordinates of isophase crossings, one per row
% piso :     vector with corresponding phase values of the isophase curves
% B :        nx2 matrix with vectors defining the direction of the bracketing
%            and the size of the initial bracketing
% tol :      structure with tolerances. See function 'cghparset'
%
% OUTPUT
% x1,x2 :    nx2 matrices with coordinates of two points on either side of an 
%            isophase curve, one per row
        
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%  
% Short version: This software is in the PUBLIC DOMAIN.
% Initial version: Ulf Griesmann, NIST, September 2018

    % some sanity checking
    sx = size(x,1);
    if ~all([size(piso,1),size(B,1)] == sx)
        error('arguments ''x,B,piso'' must have equal length.');
    end

    % preallocate matrices for bracket values
    ph1 = ph2 = zeros(sx,1);
    x1 = x2 = zeros(sx,2);

    % 'true' when a bracket was found
    bdone = bdisc = false(sx,1);
    inb = find(~bdone);
    
    % evaluate phase for bracket points (make 1st one larger)
    for k = [sqrt(2),2:tol.maxbs].^2  % non-lin reduces number of iterations ?

        % define brackets
        x1(inb,:) = x(inb,:) - k*B(inb,:);
        x2(inb,:) = x(inb,:) + k*B(inb,:);

        % evaluate phase at locations without valid brackets
        ph1(inb) = fphase(x1(inb,1),x1(inb,2),par);
        ph2(inb) = fphase(x2(inb,1),x2(inb,2),par);

        % check if (more) valid brackets were found, update 
        bdone = bdone | sign(ph1-piso) .* sign(ph2-piso) < 0;
        
        % check if a phase discontinuity was encountered
        ind = find(abs(ph2 - ph1) > pi); % would be 2*pi
        if ~isempty(ind)
            bdone(ind) = true;   % don't keep looking here
            bdisc(ind) = true;
        end

        % remaining un-bracketed points
        inb = find(~bdone);
        
        if all(bdone)
            ind = find(bdisc);
            break
        end        
    end
        
    if ~isempty(inb)   % bracket search has failed outright
        ph1 = fphase(x1(inb,:)(:,1),x1(inb,:)(:,2),par);
        ph2 = fphase(x2(inb,:)(:,1),x2(inb,:)(:,2),par);
        phv = piso(inb);
        for k=1:length(inb)
            fprintf('\n   Failed bracket search at index: %d\n', inb(k) );
            fprintf('        x1 = (%.3f, %.3f); ph1 = %f\n', x1(inb(k),:), ph1(k) );  
            fprintf('        x2 = (%.3f, %.3f); ph2 = %f\n', x2(inb(k),:), ph2(k) );  
            fprintf('   pred xy = (%.3f, %.3f); pht = %f\n\n', x(inb(k),:),phv(k) );  
        end
        fprintf('   Possible remedies: REDUCE tol.seglen and/or INCREASE tol.maxbs\n\n');
        error('failed to find a valid bracketing in ''find_bracket''.');
    end
    
    if ~isempty(ind)   % too close to discontinuity, seek a tight bracket
        ph1 = fphase(x1(ind,:)(:,1),x1(ind,:)(:,2),par); % phase at x1,x2
        ph2 = fphase(x2(ind,:)(:,1),x2(ind,:)(:,2),par);
        phv = piso(ind);                                 % target phase values
        for k=1:length(ind)
            if abs(ph1(k)-phv(k)) > pi                   % x1 is on the wrong side
                xya = x2(ind(k),:);
                xyb = x1(ind(k),:);
                x1(ind(k),:) = find_discon_line(fphase,par,xya,xyb,phv(k),tol.discon);
            elseif abs(ph2(k)-phv(k)) > pi               % x2 is on the wrong side
                xya = x1(ind(k),:);
                xyb = x2(ind(k),:);
                x2(ind(k),:) = find_discon_line(fphase,par,xya,xyb,phv(k),tol.discon);
            else
                error('line x1-x2 does not cross discontinuity in ''find_bracket''.');
            end
        end
    end

end
