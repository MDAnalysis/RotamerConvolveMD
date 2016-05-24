# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# Convolve MTSS rotamers with MD trajectory.
# Copyright (c) 2011-2013 Philip Fowler, Oliver Beckstein
# Published under the GNU Public Licence, version 2 (or higher)
#
# Includes a rotamer library for MTSS at 298 K by Gunnar Jeschke,
# which is published under the same licence by permission.
"""\
:mod:`rotcon` --- Convolve Rotamers
===================================

:Author:    Philip W Fowler, Oliver Beckstein
:Year:      2013
:Licence:   GNU Public Licence, version 2 (or higher)
:Copyright: © 2013 Philip W Fowler, Oliver Beckstein
:Citation:  LS Stelzl, PW Fowler, MSP Sansom, O Beckstein. J Mol Biol
            426 (2014), 735-751. doi: 10.1016/j.jmb.2013.10.024

This package contains the *MTSL rotamer library R1A_298K* provided by
`Gunnar Jeschke`_, which is also published under the GPL with his
permission.

Summary
=======

This package analyses molecular dynamics trajectories or
conformational ensembles in terms of spin-label distances as probed in
double electron-electron resonance (DEER) experiments. The spin labels
are fitted on trajectories and the spin label mobility is taken into
account using a rotamer library.


Background
==========

Double electron electron spin resonance (DEER) is an EPR technique for
measuring distances between two spin labels that have been covalently
attached to a protein. Two cysteine residues are introduced into the
protein and subsequently labelled. The positions are chosen to report
on the expected conformational change. A commonly used spin label is
(1-oxyl-2,2,5,5-tetramethylpyrroline-3-methyl)-methanethiosulfonate
(MTSL). MTSL has a linker with five rotatable bonds and is therefore
very flexible. The distance distributions between the two spin labels
measured by experiments are typically broad and often multi-modal. The
distributions are therefore a convolution of the flexibility of the
MTSL spin label and the conformational spread of the proteins in the
sample. To ensure that we compared like with like we developed a
method that

1. maps rotamer libraries of the MTSL spin label onto each position,

2. discards those rotamers that sterically clash with the protein
   (typically distances <2 Å) and

3. calculates all (weighted) distance pairs between the remaining
   rotamers and

4. thereby estimates a distance distribution for that structure.

The code was written in Python using the MDAnalysis_ library
[Michaud-Agrawal2011]_ and a published rotamer library for MTSL
[Polyhach2011]_. It is available for download from the MDAnalysis
website, https://github.com/MDAnalysis/RotamerConvolveMD .

Our approach improves upon the existing method [Polyhach2011]_ by
increasing computational efficiency and implementing, via the
MDAnalysis library, analysis of ensembles of hundreds of structures,
which allowed us to estimate distance distributions for entire
simulation trajectories.

Examples of the application of this approach can be found in
[Stelzl2014]_.


Installation
============

Typically, install the package with ::

   python setup.py install --user

(see also the file ``INSTALL.rst``)

This will automatically install MDAnalysis and other dependencies. If
problems arise, try installing MDAnalysis first (see
http://www.mdanalysis.org for help).

Analysis is performed with the script ``convolve-mtss-rotamers.py``,
which will have been installed in your ``bin`` directory. You might
have to add the bin directory to your ``PATH``. Consult your Python
documentation for the details although often this will be
``~/.local/bin`` (for a ``--user`` installation or ``/usr/local/bin``
for site-wide installation).


Usage
=====

Analysis is performed with the script
``convolve-mtss-rotamers.py``. It takes as input

* a topology or structure file (PSF, GRO, PDB, ... any `topology
  format`_ recognized by MDAnalysis)
* a trajectory (DCD, XTC, TRR, ... any `trajectory format`_ that
  MDAnalysis can read)

A typical invocation::

   convolve-mtss-rotamers.py --resid 47 330  --histogramBins 0 80 1  --clashDistance 2.2  \
          --output "dat/peptso-xrd"  --dcdfilename "dcd/peptso-xrd-47-330" \
          peptso.gro peptso.xtc

It loads the MD trajectory from the topology ``peptso.gro`` and the
trajectory ``peptso.xtc``. The ``--resid`` pair is required and
denotes the residue numbers (in the topology) to which the MTSSL spin
labels would be attached. Rotamers that overlap with protein atoms as
measured by an atom-atom distance smaller than the ``--clashDistance``
will be discarded and not counted in the distance calculations. For
further explanations see the ``--help`` option.

For an example, see ``doc/example`` in the source distribution. The
example can also be run to test the installation as reference output
is provided.


Help
====

If you have questions or problems installing the package then ask on
the MDAnalysis user mailing list:
http://groups.google.com/group/mdnalysis-discussion


References
==========

.. Links
.. -----

.. _MDAnalysis: http://www.mdanalysis.org
.. _Gunnar Jeschke: http://www.epr.ethz.ch/
.. _topology format:
   https://pythonhosted.org/MDAnalysis/documentation_pages/topology/init.html#supported-topology-formats
.. _trajectory format:
   https://pythonhosted.org/MDAnalysis/documentation_pages/coordinates/init.html#id1

.. Articles
.. --------

.. [Michaud-Agrawal2011] N. Michaud-Agrawal, E. J. Denning,
   T. B. Woolf, and O. Beckstein. MDAnalysis: A toolkit for the
   analysis of molecular dynamics simulations. J Comp Chem,
   32:2319-2327, 2011. doi:`10.1002/jcc.21787`_. http://www.mdanalysis.org

.. _`10.1002/jcc.21787`: http://doi.org/10.1002/jcc.21787

.. [Polyhach2011] Y. Polyhach, E. Bordignon, and G. Jeschke. Rotamer
   libraries of spin labelled cysteines for protein
   studies. Phys. Chem. Chem. Phys., 13:2356-2366, 2011.
   doi: `10.1039/C0CP01865A`_.

.. _`10.1039/C0CP01865A`: http://dx.doi.org/10.1039/C0CP01865A

.. [Stelzl2014] L. S. Stelz, P. W. Fowler, M. S. P. Sansom, and
   O. Beckstein. Flexible gates generate occluded intermediates in the
   transport cycle of LacY. J Mol Biol, 426:735-751, 2013.
   doi: `10.1016/j.jmb.2013.10.024`_

.. _`10.1016/j.jmb.2013.10.024`: http://dx.doi.org/10.1016/j.jmb.2013.10.024



"""

VERSION = 1, 2, 0

def get_version_string():
    return ".".join([str(x) for x in VERSION])
