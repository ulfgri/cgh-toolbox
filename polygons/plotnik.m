function plotnik(xy,tile=[])
% plotnik(xy,tile=[])
%
% Plots a cell array of polygons for debugging(*) purposes. Start vertices
% are indicated by a green dot. Optionally, a tile polygon can be plotted.
%
% INPUT
% xy :  cell array of polygons
% tile : polygon circumscribing a rectangular tile
%
% (*) Плотник: carpenter

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Ulf Griesmann, NIST, July 2019
    
    figure
    hold on
    
    if ~iscell(xy), xy = {xy}; end
    
    % plot tile
    if ~isempty(tile)
        plot(tile(:,1),tile(:,2),'b-', 'linewidth',2);
    end
    
    % plot polygons
    for k=1:numel(xy)
        plot(xy{k}(:,1),xy{k}(:,2),'b-');
        plot(xy{k}(:,1),xy{k}(:,2),'r.', 'markersize',12);
    end

    % green dots on initial polygon vertices
    pvi = zeros(numel(xy),2);
    for k=1:numel(xy)
        pvi(k,:) = xy{k}(1,:);
    end
    plot(pvi(:,1),pvi(:,2), 'g.', 'markersize',12);
    
    % grid on
    axis equal
end
