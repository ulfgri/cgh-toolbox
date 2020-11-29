Octave Toolbox for Computer-Generated Holograms (CGH)
=====================================================

Ulf GRIESMANN, 2020

A toolbox for [Octave](https://octave.org) to compute the layouts of
computer-generated holograms from scalar phase functions describing
the optical action of holograms in optical imaging systems. The
optical function of a hologram is generally modeled with optical
ray-tracing software and it can be encapsulated by a scalar optical
phase function $\phi$(x,y). The toolbox converts phase functions into
equivalent binary or multi-level holograms. The algorithm used in this
toolbox takes advantage of the relationship between the local
derivatives of phase functions and the local geometry (curvature) of
isophase lines. The isophase-following algorithm is easily extended
to phase functions with singularities and discontinuities that occur
in some optical applications. [A full description of the toolbox
algorithms has been published.](https://nvlpubs.nist.gov/nistpubs/jres/125/jres.125.024.pdf).
The [GDSII Toolbox](https://github.com/ulfgri/gdsii-toolbox) can be
used to create hologram layouts suitable for fabrication with
direct-write lithography systems. The CGH toolbox also
includes a family of functions for the robust and efficient estimation
and evaluation of Zernike polynomials, which are widely used in optical
applications.


Copyright
---------
The functions in this toolbox are in the **Public Domain**. For details
see the file LICENSE in the root directory of the toolbox.


Documentation
-------------
Additional documentation is available on:

https://sites.google.com/site/ulfgri/numerical/cgh-toolbox

in a user manual that describes the installation of the toolbox
and its use through a set of detailed examples.


Software Dependencies
---------------------
The CGH toolbox makes use of the [Octave Parallel Package](https://octave.sourceforge.io/parallel)
to speed up the computation of holograms on computers with multiple
CPU cores. The Parallel Package can now be used on all common
operating systems. In addition, the CGH toolbox requires the
[Octave Optimization Package](https://octave.sourceforge.io/optim). Symbolic
computations of Zernike polynomials require the
[Octave Symbolic Package](https://octave.sourceforge.io/symbolic) and an
installation of [Python](https://www.python.org) with
[SymPy](https://www.sympy.org/en/index.html).


Compiling
---------
The CGH toolbox contains several MEX functions that must be compiled
with a C compiler before the toolbox can be used.

For Octave on Linux, the mex functions are compiled by executing

$ ./makemex-octave

at the shell prompt. In Octave on Windows the mex functions are
compiled by changing to the ./cgh-toolbox directory and running

>> makemex

at the Octave command prompt.


Help
----
If you find a bug in the software, please send a message to
ulf.griesmann@nist.gov or ulfgri@gmail.com and I will try to correct it.
