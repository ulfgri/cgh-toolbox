% Zernike functions
% Version 2.1 11-Dec-2019
%
% NOTE: This toolbox is based on Zernike polynomials defined in the
%       ANSI Z80.28-2017 standard (see zern_eval for definition).
%       Functions are available to translate parameters from one
%       definition to another (see e.g., zern_par_trans)
%
% Definition and conversion
% -------------------------
%   zern_par_trans       - Translates a Zernike parameter vector from one definition
%                          to another.
%   zern_index_trans     - Translates a Zernike index from one definition to another
%   zern_terms           - Returns a vector with the ANSI indices of common Zernike
%                          terms (e.g., focus)
%   zern_terms_nplusm    - Returns a vector with the ANSI indices of the Zernike
%                          terms corresponding to a specified sum of the radial
%                          order n and angular frequency |m|.
%   zern_terms_symmetry  - Extract indices of Zernike terms with specified symmetry
%                          properties
%   zern_terms_upto_nm   - Returns a vector with the ANSI indices of the Zernike
%                          terms up to a radial order n and up to angular
%                          order |m|
%   zern_symbolic        - Generate the symbolic definition of an ANSI Zernike
%                          polynomial term.
%   zern_norm            - Calculate the ANSI normalization component of a specified
%                          Zernike term
%
% Evaluation
% ----------
%   zern_eval            - Evaluation of a Zernike polynomial for an image array
%   zern_eval_pol        - Evaluation of a Zernike polynomial at points specified
%                          in normalized polar coordinates
%   zern_eval_pos        - Evaluation of a Zernike Polynomial at points specified in  
%                          cartesian coordinates
%   zern_map2pol         - Calculate angles and normalized radii for each element of
%                          an image array
%   zern_pos2pol         - Calculate angles and normalized radii for points
%                          specified in cartesian coordinates
%
% Estimation
% ----------
%   zern_estim           - Estimate least-squares best-fit Zernike coefficients
%                          and residual image for an image array
%   zern_estim_pol       - Estimate least-squares best-fit Zernike coefficients
%                          and residuals for data specified in normalized polar
%                          coordinates
%
% Transformations
% ---------------
%   zern_transform       - Calculate Zernike coefficients after applying
%                          a translation, scaling, and rotation
%   zern_subaperture     - Calculate Zernike coefficients of a sub aperture
%   zern_mirror_line     - Calculate effect of a mirror operation over a line
%                          on the Zernike coefficients
%   zern_mirror_point    - Calculate effect of a mirror operation over the
%                          center of the Zernike normalization disk on the
%                          Zernike coefficients
%   zern_rot             - Calculate effect of a rotation on the Zernike coeffients
%
% Auxiliary
% ---------
%   private/sub2pos      - user defined function that converts an image array
%                          index [i,j] into a position [x,y]
%

% This software was developed by employees of the National Institute
% of Standards and Technology (NIST), an agency of the Federal Government,
% and is being made available as a public service. Pursuant to title 17
% United States Code Section 105, works of NIST employees are not subject
% to copyright protection in the United States.  This software may be subject
% to foreign copyright.  Permission in the United States and in foreign
% countries, to the extent that NIST may hold copyright, to use, copy,
% modify, create derivative works, and distribute this software and its
% documentation without fee is hereby granted on a non-exclusive basis,
% provided that this notice and disclaimer of warranty appears in all copies.
%
% THIS SOFTWARE IS BEING PROVIDED FOR RESEARCH AND EDUCATIONAL PURPOSES ONLY.
% THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
% EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY
% THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM FROM
% INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO THE
% SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE.  IN NO EVENT
% SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT,
% INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM,
% OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY,
% CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS
% OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE
% OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.
%
% Author: Johannes Soons & Ulf Griesmann; NIST; Feb 2000 - Dec 2019 
