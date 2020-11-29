function cghpar = cghparset(varargin)
%function cghpar = cghparset('par1',val1, 'par2',val2, ...)
%function cghpar = cghparset(oldpar, 'par1',val1, 'par2',val2, ...)
%
% Creates or modifies a structure with parameters and tolerances for
% computer generated hologram (CGH) layout calculations and provides
% parameter default values.
%
% INPUT
% oldpar :   (optional) an existing structure with parameters to be modified
% par/val :  (optional) list of name/value pairs consisting of a string
%            specifying the parameter followed by its value. If the argument
%            list is empty, default parameters are returned. The following 
%            parameters can be set:
%                 
%            algorithm: selects default algorithm for CGH computation
%                          'pa' :  pilot approximation
%                          'fi' :  isophase following algorithm (default)
%                          'hi' :  phase functions with Hilbert phase terms
%            vphase :   1x2 vector with phase boundaries of polygon areas that enclose
%                       phase values between vphase(1) and vphase(2). Default is [0,pi).
%            tol :      a structure with tolerances for the isophase finding
%                           tol.vertex : tolerance for polygon vertex positions
%                                        Default is 1e-4 length units
%                           tol.segdev:  maximum allowable deviation of segment from contour 
%                                        Default is 5e-3 length units.
%                           tol.seglen:  If present, limits the length of polgon segments to
%                                        the specfied maximum. Needed for isophase lines
%                                        curvature sign changes. Default is [].
%                           tol.maxit :  maximum number of iterations for root finding .
%                                        Default is 150.
%                           tol.maxbs :  maximum number of search steps for finding
%                                        isophase bracketing points. Default is 32.
%                           tol.hderiv : increment for calculation of derivatives 
%                                        Default is 15e-3 length units.
%                           tol.discon : tolerance for discontinuity location
%                                        Default is 1e-5 length units.
%                           tol.maxpol : maximum number of polygon segments to connect.
%                                        Default is 10.
%            grid :     a structure with sampling resolution on tile edge and area
%                            grid.edge : sampling resolution on tile edges. 
%                                        Default is 10000. 
%                            grid.area : sampling resolution in tile area. 
%                                        Default is [100,100]. For the 'pa' algorithm the
%                                        grid can be specfied as a vector [nx,ny].
%            sing:      a structure array with information about a phase singularity that is
%                       needed for phase functions with a Hilbert phase term. Only needed
%                       for the 'hi' algorithm. A tile can contain ONLY ONE singularity,
%                       the singularity must be located on the tile edge. Default is [].
%                            sing(k).pos :  the location of a singularity in the x-y plane
%                            sing(k).era :  exclusion radius around singularity. Isophase 
%                                           following terminates when the distance to the 
%                                           singularity is smaller than sing.era.
%            aperture:  a structure array containing aperture polygons and a
%                       corresponding Boolean operation that will be applied to the 
%                       output polygons:
%                             aperture(k).poly : clipping polygon(s);
%                                                see 'polybool.m'
%                             aperture(k).oper : clipping operation, 
%                                                'and', 'or', etc.
%                       Clipping operations are applied successively, 
%                       THE ORDER CAN MATTER. Default is [].
%            debug:     If set to 'true', errors in the calculations on a CGH tile are not 
%                       trapped and the Octave thread displays the error cause and terminates. 
%                       Default is 'false'.
%
% OUTPUT
% cghpar :     structure with parameters and tolerances
    
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 Section 105 of the United States Code this software is not
% subject to copyright protection and is in the Public Domain. This software is an
% experimental system. NIST assumes no responsibility whatsoever for its use by other
% parties, and makes no guarantees, expressed or implied, about its quality, reliability,
% or any other characteristic.
%
% Ulf Griesmann, NIST, October 2018

    % some default values
    EDGE = 10000;      % number of samples for isoline-edge intersections
    AREA = [100,100];  % number of samples for coarse isophase calculation

    % set output defaults
    if numel(varargin) > 0 && isstruct(varargin{1})
        cghpar = varargin{1};
        sidx = 2;
    else
        cghpar = struct('vphase',[0,pi], ...
                        'tol',struct('vertex',1e-4, ...
                                     'segdev',5e-3, ...
                                     'seglen',[], ...
                                     'hderiv',15e-3, ...
                                     'maxbs',32, ...
                                     'maxit',150, ...
                                     'discon',1e-5, ...
                                     'maxpol',10), ...
                        'grid',struct('edge',EDGE, 'area',AREA), ... 
                        'sing',[], ...
                        'aperture',[], ...
                        'algorithm','fi', ...
                        'debug',false);
        sidx = 1;
    end
    
    % replace defaults from name/value argument pairs
    for k = sidx:2:numel(varargin)
        
        if ~ischar(varargin{k})
            error('error in keyword / value pair');
        end
        
        switch varargin{k}            
          case {'tol','grid'} % structures with variable numbers of fields
                S = varargin{k+1};
                F = fieldnames(S);
                for m = 1:numel(F)
                    cghpar.(varargin{k}).(F{m}) = S.(F{m});
                end
          otherwise
              cghpar.(varargin{k}) = varargin{k+1};
        end
        
    end

    % a few sanity and spelling checks
    allnam = {'vphase','tol','grid','sing','aperture','algorithm','debug'};
    actnam = fieldnames(cghpar);
    correct = ismember(actnam,allnam);
    if any(~correct)
        fprintf('   unrecognized or mis-spelled CGH parameter(s):\n   >>> ');
        fprintf('%s   ', actnam{~correct}); fprintf('\n');
        error('parameter name error.');
    end
    tolnam = {'vertex','segdev','seglen','hderiv','maxbs','maxit','discon','maxpol'};
    actnam = fieldnames(cghpar.tol);
    correct = ismember(actnam,tolnam);
    if any(~correct)
        fprintf('   unrecognized or mis-spelled CGH tolerance(s):\n   >>> ');
        fprintf('%s   ', actnam{~correct}); fprintf('\n');
        error('tolerance name error.');
    end
    
    % check if grid is fully specified
    if ~isfield(cghpar.grid,'area'), cghpar.grid.area = AREA; end
    if ~isfield(cghpar.grid,'edge'), cghpar.grid.edge = EDGE; end
    
    if numel(cghpar.vphase) < 2
        error('argument ''vphase'' must be a 1x2 vector');
    end
    
    if ~ismember(cghpar.algorithm, {'pa','fi','hi'})
        error('an unknown or unimplemented algorithm was specified.');
    end
    
    if strcmp(cghpar.algorithm, 'hi') && isempty(cghpar.sing)
        error('a singularity structure is required for the ''hi'' algorithm.');
    end

end
