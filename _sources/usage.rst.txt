.. -*- coding: utf-8 -*-

=========================
 Using RotamerConvolveMD
=========================

When installing the package, two scripts are installed in the ``bin``
directory:

* :program:`convolve-mtss-rotamers.py` analyzes standard DEER
  experiments with MTSL spin labels.
* :program:`convolve-mtss-rotamers_pre.py` analyzes PRE experiments
  with MTSL spin labels.

Scripts are run from the command line.  


DEER
----

Analysis for standard DEER experiments with MTSL spin labels is
performed with the script ``convolve-mtss-rotamers.py``. It takes as
input

* a topology or structure file (psf, gro, pdb, ... any `topology
  format`_ recognized by mdanalysis)
* a trajectory (dcd, xtc, trr, ... any `trajectory format`_ that
  mdanalysis can read)

A typical invocation::

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
        peptso.gro peptso.xtc

It loads the MD trajectory from the topology ``peptso.gro`` and the
trajectory ``peptso.xtc``. The ``--resid`` pair is required and
denotes the residue numbers (in the topology) to which the MTSSL spin
labels would be attached. Rotamers that overlap with protein atoms as
measured by an atom-atom distance smaller than the ``--clashDistance``
will be discarded and not counted in the distance calculations. 
The user can decide to use either N1 ``--useNOelectron`` or the 
geometric midpointis N1 and O1  ``--no-useNOelectron``  to calculate 
the distances. For further explanations see the ``--help`` option.

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

A typical invocation::

    convolve-mtss-rotamers_pre.py \
        --resid 47  \
        --clashDistance 2.2  \
        --plotname "dat/peptso-xrd-47.pdf" \
        --outputRawDistances "dat/peptso-xrd" \
        --dcdfilenameAll "dcd/peptso-xrd" \
        --dcdfilenameNoClashes "dcd/peptso-xrd" \
        --useNOelectron \
        peptso.gro peptso.xtc 

The ``--resid`` is required and denotes the residue number (in the topology) 
to which the MTSSL spin label would be attached. Rotamers that overlap 
with protein atoms as measured by an atom-atom distance smaller than 
the ``--clashDistance`` will be discarded and not counted in the distance 
calculations. The user can decide to use either N1 ``--useNOelectron`` 
or the geometric midpointis N1 and O1  ``--no-useNOelectron``  to calculate 
the distances. For further explanations see the ``--help`` option.


.. Links

.. _topology format: 
   https://www.mdanalysis.org/docs/documentation_pages/topology/init.html#supported-topology-formats
.. _trajectory format:
   https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#id1
   
