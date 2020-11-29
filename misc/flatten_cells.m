function carr = flatten_cells(carr)
%
% Turns a nested cell array (cell array of cell arrays of ...) into a flat 
% cell array. 
%

% A clever little function, found on the Web - author unknown, June 2013

    try
        carr = cellfun(@flatten_cells,carr,'UniformOutput',0);
        if any(cellfun(@iscell,carr))
            carr = [carr{:}];
        end
    catch
        % reached non-cell, pass through unchanged
    end
    
end
