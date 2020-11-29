function [cap] = phase2cgh(stile,fphase,phapar,cghpar,ncpu=[],verbose=0)
%function [cap] = phase2cgh(stile,fphase,phapar,cghpar,ncpu=[],verbose=0)
% 
% Generates the layout of a computer generated hologram, defined by
% polygons, from a phase function.
%
% IMPORTANT: - all dimensions must have the same units.
%            - use 'cghparset' to create the 'cghpar' structure
%
% INPUT:
% stile    :   structure array with CGH area tiling
%                   stile(k).llc : [x,y], lower left corner of tile k
%                   stile(k).urc : [x,y], upper right corner of tile k
% fphase :     function handle of user defined function describing the 2D phase
%              distribution fphase(x,y,phapar). MUST RETURN PHASE IN RADIANS. 
% phapar :     structure with the parameters of the phase function
%              that is passed to the phase function fphase(x,y,phapar).
%              see user defined phase functions for details.
% cghpar :     structure with parameters controlling the generation
%              of CGH polygons. See 'cghparset' for default values.
% ncpu :       (Optional) number of processors to use for the calculation. This 
%              argument is (currently) ignored when the function is used with Octave
%              on Windows. On Linux, the number of available processors, as returned 
%              by the 'nproc()' function, is used by default.
% verbose   :  (Optional) verbosity level: 0 is quiet, 1 prints a message for each tile
%              being processed, 2 prints progess messages at all stages of polygon 
%              generation. Default is 0.
%
% OUTPUT:
% cap :        cell array of polygons describing the CGH layout

% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Ulf Griesmann, NIST, October 2018
%---------------------------------------------------------------------------

    % load Octave optimization package for vfzero
    pkg load optim;
    
    % check and validate input arguments
    if nargin < 4
        error('phase2cgh: function requires at least four arguments.');
    end

    % supply default parameters
    if isempty(ncpu)
      ncpu = nproc('all');
    end
    
    % set up a cell array with parameters, one structure per tile
    % first the parameters that are shared by all tiles
    tpar.phfunc = fphase;                    % phase function handle
    tpar.par_phase = phapar;                 % phase function parameters
    tpar.par_cgh = cghpar;                   % CGH computation parameters
    tpar.verbose = verbose;                  % verbosity

    % duplicate a couple to make them easier to access in 'process_tile'
    tpar.vphase = cghpar.vphase;           % phase interval
    tpar.aperture = cghpar.aperture;       % aperture polygon(s)

    % then add parameters that are tile specific
    parlist = cell(1,numel(stile));
    for k = 1:numel(parlist)
        tpar.num = k;
        tpar.llc = stile(k).llc;
        tpar.urc = stile(k).urc;
        parlist{k} = tpar;
    end

    % calculate CGH polygons in each tile
    if ispc || ncpu == 1
        cap = cellfun(@process_tile, parlist, 'UniformOutput',false);
    else
        cap = parcellfun(ncpu, @process_tile, parlist, 'UniformOutput',false);
    end

    cap = flatten_cells(cap); % eliminate any nested cell arrays
    
end


%---------------------------------------------------------------------------

function [cap] = process_tile(tpar)
%
% Auxiliary function used to distribute data to multiple processes
% that can be executed on differenct processors.
%
% INPUT
% tpar :  parameters of the tile
%           tpar.num :       tile number
%           tpar.phfunc :    function handle of phase function
%           tpar.par_phase : structure with phase function parameters
%           tpar.llc :       lower left corner of tile
%           tpar.urc :       upper right corner of tile
%           tpar.par_cgh :   tolerances and parameters for CGH generation
%           tpar.algorithm:  algorithm to use for CGH computation
%           tpar.verbose :   verbosity level
%           tpar.debug :     trap errors or not
%
% OUTPUT
% cap :   cell array of polygons

    tverbose = false;  
    if tpar.verbose > 0  
        fprintf('   Processing tile %d (%s)\n', tpar.num, tpar.par_cgh.algorithm);
        if tpar.verbose == 2
            tverbose = true;
        end
    end
    
    if tpar.par_cgh.debug
        cap = dispatch_tile_processor(tpar.par_cgh.algorithm, tpar.phfunc, tpar.par_phase, ...
                                      tpar.llc, tpar.urc, ...
                                      tpar.par_cgh.vphase, tpar.par_cgh.tol, ...
                                      tpar.par_cgh.grid, tverbose, tpar.par_cgh.sing);
    else
        try
            cap = dispatch_tile_processor(tpar.par_cgh.algorithm, tpar.phfunc, tpar.par_phase, ...
                                          tpar.llc, tpar.urc, ...
                                          tpar.par_cgh.vphase, tpar.par_cgh.tol, ...
                                          tpar.par_cgh.grid, tverbose, tpar.par_cgh.sing);
        catch
            cap = {};
            fprintf('\n   >>> ERROR PROCESSING TILE %d\n\n',tpar.num);
        end
    end

    % clip the resulting polygons to the aperture(s)
    % NOTE: each polygon is clipped individually to avoid problems due to the
    % clipper library merging adjacent polygons and creating hole polygons.
    % Reason: the GDSII model of polygons is not the same as the polygon model 
    % of the Clipper library.
    if ~isempty(tpar.aperture) && ~isempty(cap)
        for k = 1:numel(tpar.aperture)
            cap = cellfun(@(P)polybool(P,tpar.aperture(k).poly,tpar.aperture(k).oper), ...
                          cap,'UniformOutput',false);
            cap = flatten_cells(cap); % eliminate any nested cell arrays
        end
    end

end

    
%---------------------------------------------------------------------------

function [cap] = dispatch_tile_processor(alg,fphase,ppar,llc,urc,vphase,tol,grid,verbose,sing)
%
% call algorithm-specific tile processing function

    switch alg
        
      case 'pa'
        cap = phase2cgh_tile_pa(fphase, ppar, llc, urc, vphase, tol, grid, verbose);
        
      case 'fi'
        cap = phase2cgh_tile_fi(fphase, ppar, llc, urc, vphase, tol, grid, verbose);
        
      case 'hi'
        cap = phase2cgh_tile_hi(fphase, ppar, llc, urc, vphase, tol, grid, sing, verbose);
    
    end
end
    
