.. -*- coding: utf-8 -*-
.. RotamerConvolve documentation master file, created by
   sphinx-quickstart on Thu Oct  5 15:44:48 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================
 RotamerConvolveMD
===================

:Release: |release|
:Date: |today|

This package analyses molecular dynamics trajectories or
conformational ensembles in terms of spin-label distances as probed in
double electron-electron resonance (DEER) experiments and spin-label
to amide protein distances as probed in paramagnetic relaxation
enhancement (PRE) experiments. The spin labels are fitted on
trajectories and the spin label mobility is taken into account using a
rotamer library.

When using this code please cite [Stelzl2014]_ and [Polyhach2011]_.

This package contains MTSL :ref:`rotamer-libraries` provided by
`Gunnar Jeschke`_ (published under the GPL with his permission).


If you have questions or problems installing the package then ask on
the MDAnalysis user mailing list:
https://groups.google.com/group/mdnalysis-discussion

Source code is available from
https://github.com/MDAnalysis/RotamerConvolveMD under the open source
GNU General Public License, version 2.

.. _Gunnar Jeschke: http://www.epr.ethz.ch/

.. Content (in sidebar)
.. ====================

.. toctree::
   :maxdepth: 4
   :hidden:	      

   installation	      
   background
   libraries   
   usage
   
References
==========

.. [Michaud-Agrawal2011] N. Michaud-Agrawal, E. J. Denning,
   T. B. Woolf, and O. Beckstein. MDAnalysis: A toolkit for the
   analysis of molecular dynamics simulations. J Comp Chem,
   32:2319-2327, 2011. doi:`10.1002/jcc.21787`_.
   https://www.mdanalysis.org

.. _`10.1002/jcc.21787`: https://doi.org/10.1002/jcc.21787

.. [Gowers2016] R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L.
		Seyler, D. L. Dotson, J. Domański, S. Buchoux, I. M. Kenney,
                and O. Beckstein. `MDAnalysis: A Python package for the
                rapid analysis of molecular dynamics
                simulations`_. In S. Benthall and S. Rostrup, editors,
                Proceedings of the 15th Python in Science Conference,
                pages 102 – 109, Austin, TX, 2016. SciPy.

.. _`MDAnalysis: A Python package for the rapid analysis of molecular
     dynamics simulations`:
     http://conference.scipy.org/proceedings/scipy2016/oliver_beckstein.html

.. [Polyhach2011] Y. Polyhach, E. Bordignon, and G. Jeschke. Rotamer
   libraries of spin labelled cysteines for protein
   studies. Phys. Chem. Chem. Phys., 13:2356-2366, 2011. 
   doi: `10.1039/C0CP01865A`_.

.. _`10.1039/C0CP01865A`: https://doi.org/10.1039/C0CP01865A

.. [Stelzl2014] L. S. Stelz, P. W. Fowler, M. S. P. Sansom, and
   O. Beckstein. Flexible gates generate occluded intermediates in the
   transport cycle of LacY. J Mol Biol, 426:735-751, 2013. 
   doi: `10.1016/j.jmb.2013.10.024`_ 

.. _`10.1016/j.jmb.2013.10.024`: https://doi.org/10.1016/j.jmb.2013.10.024


.. [Fowler2015] P. W. Fowler, M. Orwick-Rydmark, S. Radestock, N. Solcan, P. M. Dijkman,
		J. A. Lyons, J. Kwok, M. Caffrey, A. Watts, L. R. Forrest,
                and S. Newstead. Gating topology of the proton-coupled
                oligopeptide symporters. Structure, 23:290–301,
                2015. doi:`10.1016/j.str.2014.12.012`_

.. _`10.1016/j.str.2014.12.012`: https://doi.org/10.1016/j.str.2014.12.012


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
