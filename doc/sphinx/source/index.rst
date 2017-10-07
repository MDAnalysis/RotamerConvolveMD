.. RotamerConvolve documentation master file, created by
   sphinx-quickstart on Thu Oct  5 15:44:48 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

RotamerConvolveMD documentation
===============================

This package analyses molecular dynamics trajectories or
conformational ensembles in terms of spin-label distances as probed in
double electron-electron resonance (DEER) experiments and spin-label
to amide protein distances as probed in paramagnetic relaxation
enhancement (PRE) experiments. The spin labels are fitted on
trajectories and the spin label mobility is taken into account using a
rotamer library.

When using this code please cite [Stelzl2014]_ and [Polyhach2011]_.

This package contains the *MTSL rotamer library R1A_298K (2011)* and an *updated MTSL
rotamer library R1A_298K (povided in 2017)* provided by `Gunnar Jeschke`_, which is also published under the GPL with his permission. The updated rotamer library was sent by `Gunnar Jeschke`_ after personal discussion.


If you have questions or problems installing the package then ask on
the MDAnalysis user mailing list:
http://groups.google.com/group/mdnalysis-discussion

Source code is available from
https://github.com/MDAnalysis/RotamerConvolveMD under the open source
GNU General Public License, version 2.

.. _Gunnar Jeschke: http://www.epr.ethz.ch/

Content
=======

.. toctree::
   :maxdepth: 1

   installation	      
   background
   usage

   
References
==========

.. [Michaud-Agrawal2011] N. Michaud-Agrawal, E. J. Denning,
   T. B. Woolf, and O. Beckstein. MDAnalysis: A toolkit for the
   analysis of molecular dynamics simulations. J Comp Chem,
   32:2319-2327, 2011. doi:`10.1002/jcc.21787`_.
   https://www.mdanalysis.org

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



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
