.. -*- coding: utf-8 -*-

============
 Background
============

Site-directed spin labeling (SDSL) is a common technique to
investigate structure and dynamics of macromolecular
systems. Covalentry attached spin labels are introduced to the system
and induce electron spin resonance.

*Double electron electron spin resonance* (**DEER**) is an EPR_
technique for measuring distances between two spin labels that have
been covalently attached to a protein. Two cysteine residues are
introduced into the protein and subsequently labelled.

*Paramagnetic relaxation enhancement* (**PRE**) is an NMR_ technique for
measuring distances between a spin label and the amide protons of the
protein backbone. One cysteine residue is introduced at the position
of the label.

The positions are chosen to report on the expected conformational
change. A commonly used spin label is
(1-oxyl-2,2,5,5-tetramethylpyrroline-3-methyl)-methanethiosulfonate
(MTSL). MTSL has a linker with five rotatable bonds and is therefore
very flexible.


Method
======

The distance distributions between the two spin labels (DEER) or one
spin label and the amide protons are measured by experiments (PRE) are
typically broad and often multi-modal. The distributions are therefore
a convolution of the flexibility of the MTSL spin label and the
conformational spread of the proteins in the sample. In order to
compare a structural ensemble (e.g., from MD simulations or modelling)
to experimental DEER or PRE data, we developed a post-processing
method that

1. maps rotamer libraries of the MTSL spin label onto each position,

2. discards those rotamers that sterically clash with the protein
   (typically distances <2 Å) and

3. calculates all (weighted) distance pairs between the remaining
   rotamers and 

4. thereby estimates a distance distribution for that structure. 

The code was written in Python using the MDAnalysis_ library
[Michaud-Agrawal2011]_ [Gowers2016]_ and :ref:`rotamer-libraries` for
MTSL (2011, 2015) [Polyhach2011]_.

Our approach improves upon the existing method [Polyhach2011]_ by
increasing computational efficiency and implementing, via the
MDAnalysis library, analysis of ensembles of hundreds of structures,
which allowed us to estimate distance distributions for entire
simulation trajectories. In the case of PRE measurements, it enables
the user to calculate back the transverse relaxation enhancement to
compare raw data without calculating the distances based on the
experiment.

In the case of MTSL, the distances are determined by considering the
position of the free electron located between nitrogen (N1) and (O1).


Example applications
====================

Examples of the application of this approach for comparing
experimental DEER data to molecular dynamics (MD) simulations can be
found in [Stelzl2014]_ and [Fowler2015]_.


.. _MDAnalysis: https://www.mdanalysis.org
.. _EPR: https://en.wikipedia.org/wiki/Electron_paramagnetic_resonance
.. _NMR: https://en.wikipedia.org/wiki/Nuclear_magnetic_resonance
