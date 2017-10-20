.. -*- mode: rst; coding: utf-8 -*-

======================================
 MTSL Rotamer fitting to trajectories
======================================

:Author:    Philip W Fowler, Oliver Beckstein, Katrin Reichel, and AUTHORS_
:Year:      2013
:Licence:   GNU Public Licence, version 2 (or higher)
:Copyright: © 2013 Philip W Fowler, Oliver Beckstein,
            © 2014–2017 AUTHORS_
:Citation:  LS Stelzl, PW Fowler, MSP Sansom, O Beckstein. J Mol Biol
            426 (2014), 735-751. doi: 10.1016/j.jmb.2013.10.024
:Documentation: |docs|
	    
.. _AUTHORS:
   https://raw.githubusercontent.com/MDAnalysis/RotamerConvolveMD/master/AUTHORS

This package contains the *MTSL rotamer library R1A_298K (2011)* and
an *updated MTSL rotamer library R1A_298K (2015)* provided by `Gunnar
Jeschke`_, which is also published under the GPL with his
permission. The updated rotamer library was sent by `Gunnar Jeschke`_
after personal discussion.

Summary
=======

This package analyses molecular dynamics trajectories or
conformational ensembles in terms of spin-label distances as probed in
double electron-electron resonance (DEER) experiments and spin-label
to amide protin distances as probed in paramagnetic relaxation
enhancement (PRE) experiments. The spin labels are fitted on
trajectories and the spin label mobility is taken into account using a
rotamer library.

For further details see the `RotamerConvolveMD documentation`_.



Background
==========

Site-directed spin labeling (SDSL) is a common technique to
investigate structure and dynamics of macromolecular
systems. Covalentry attached spin labels are introduced to the system
and induce electron spin resonance.  Double electron electron spin
resonance (DEER) is an EPR technique for measuring distances between
two spin labels that have been covalently attached to a protein. Two
cysteine residues are introduced into the protein and subsequently
labelled. Paramagnetic relaxation enhancement (PRE) is an NMR
technique for measuring distances between a spin label and the amide
protons of the protein backbone. One cysteine residue is introduced at
the position of the label The positions are chosen to report on the
expected conformational change. A commonly used spin label is
(1-oxyl-2,2,5,5-tetramethylpyrroline-3-methyl)-methanethiosulfonate
(MTSL). MTSL has a linker with five rotatable bonds and is therefore
very flexible. The distance distributions between the two spin labels
(DEER) or one spin label and the amide protons are measured by
experiments are typically broad and often multi-modal. The
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
simulation trajectories. In the case of PRE measurements, it enables 
the user to calculate back the transverse relaxation enhancement 
to compare raw data without calculating the distances based on the
experiment.

In the case of MTSL, the distances are determined by considering the position of the free electron
located between nitrogen (N1) and (O1).

Examples of the application of this approach can be found in
[Stelzl2014]_.


Installation
============

Typically, install the package with pip_ ::

   pip install RotamerConvolveMD

(see also the file ``INSTALL.rst``)

This will automatically install MDAnalysis and other dependencies. If
problems arise, try installing MDAnalysis first (see
http://www.mdanalysis.org for help).

Analysis is performed with the script ``convolve-mtss-rotamers.py`` 
(DEER) or ``convolve-mtss-rotamers_pre.py`` (PRE), 
which will have been installed in your ``bin`` directory. You might
have to add the bin directory to your ``PATH``. Consult your Python
documentation for the details although often this will be
``~/.local/bin`` (for a ``--user`` installation or ``/usr/local/bin``
for site-wide installation). 


Usage
=====

DEER
----

Analysis for standard DEER experiments with MTSL spin labels is
performed with the script ``convolve-mtss-rotamers.py``. It takes as
input

* a topology or structure file (psf, gro, pdb, ... any `topology
  format`_ recognized by mdanalysis)
* a trajectory (dcd, xtc, trr, ... any `trajectory format`_ that
  mdanalysis can read)

a typical invocation::

    convolve-mtss-rotamers.py \
        --resid 47 330  \
        --histogramBins 0 80 1  \
        --clashDistance 2.2  \
        --output "dat/peptso-xrd" \
        --plotname "dat/peptso-xrd.pdf" \
        --outputRawDistances "dat/peptso-xrd" \
        --dcdfilename "dcd/peptso-xrd" \
        --dcdfilenameNoClashes "dcd/peptso-xrd" \
        --useNOelectron \
        --libname "MTSSL 298K 2015" \
        peptso.gro peptso.xtc

It loads the MD trajectory from the topology ``peptso.gro`` and the
trajectory ``peptso.xtc``. The ``--resid`` pair is required and
denotes the residue numbers (in the topology) to which the MTSSL spin
labels would be attached. Rotamers that overlap with protein atoms as
measured by an atom-atom distance smaller than the ``--clashDistance``
will be discarded and not counted in the distance calculations. 
The user can decide to use either N1 ``--no-useNOelectron`` or the 
geometric midpointis N1 and O1  ``--useNOelectron``  to calculate 
the distances. Two libraries are offered: `MTSSL 298K 2011` or 
`MTSSL 298K 2015` (default) by defining the option ``--libname``. 
For further explanations see the ``--help`` option.

For an example, see ``doc/example`` in the source distribution. The
example can also be run to test the installation as reference output
is provided.


PRE
---

Analysis for standard PRE experiments with MTSL spin label is performed 
with the script ``convolve-mtss-rotamers_pre.py``. Similar to the 
analysis of DEER experiments, it takes as inputs:

* a topology or structure file (psf, gro, pdb, ... any `topology
  format`_ recognized by mdanalysis)
* a trajectory (dcd, xtc, trr, ... any `trajectory format`_ that
  mdanalysis can read)

a typical invocation::

    convolve-mtss-rotamers_pre.py \
        --resid 47  \
        --clashDistance 2.2  \
        --plotname "dat/peptso-xrd-47.pdf" \
        --outputRawDistances "dat/peptso-xrd" \
        --dcdfilenameAll "dcd/peptso-xrd" \
        --dcdfilenameNoClashes "dcd/peptso-xrd" \
        --useNOelectron \
        --libname "MTSSL 298K 2015" \
        peptso.gro peptso.xtc 

The ``--resid`` is required and denotes the residue number (in the topology) 
to which the MTSSL spin label would be attached. Rotamers that overlap 
with protein atoms as measured by an atom-atom distance smaller than 
the ``--clashDistance`` will be discarded and not counted in the distance 
calculations. The user can decide to use either N1 ``--no-useNOelectron`` 
or the geometric midpointis N1 and O1  ``--useNOelectron``  to calculate 
the distances. Two libraries are offered: `MTSSL 298K 2011` or 
`MTSSL 298K 2015` (default) by defining the option ``--libname``. 
For further explanations see the ``--help`` option.


Help
====

If you have questions or problems installing the package then ask on
the MDAnalysis user mailing list:
http://groups.google.com/group/mdnalysis-discussion

	
References
==========

.. Links
.. -----

.. _`RotamerConvolveMD documentation`:
   https://www.mdanalysis.org/RotamerConvolveMD
.. _MDAnalysis: http://www.mdanalysis.org
.. _Gunnar Jeschke: http://www.epr.ethz.ch/
.. _topology format: 
   https://pythonhosted.org/MDAnalysis/documentation_pages/topology/init.html#supported-topology-formats
.. _trajectory format:
   https://pythonhosted.org/MDAnalysis/documentation_pages/coordinates/init.html#id1
.. _pip: https://pip.pypa.io/

.. Badges
.. ------
.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
   :alt: Documentation (latest release)
   :target: `RotamerConvolveMD documentation`_

   
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

.. [Stelzl2014] L. S. Stelzl, P. W. Fowler, M. S. P. Sansom, and
   O. Beckstein. Flexible gates generate occluded intermediates in the
   transport cycle of LacY. J Mol Biol, 426:735-751, 2013. 
   doi: `10.1016/j.jmb.2013.10.024`_ 

.. _`10.1016/j.jmb.2013.10.024`: http://dx.doi.org/10.1016/j.jmb.2013.10.024


