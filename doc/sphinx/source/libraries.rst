.. -*- coding: utf-8 -*-


.. _rotamer-libraries:

===================
 Rotamer libraries
===================

This package contains the *MTSL rotamer library R1A_298K (2011)*
[Polyhach2011]_ and an *updated MTSL rotamer library R1A_298K (2015)*
(which is the default library for MTSL in MMM_ since version 2015.1).

Both libraries were provided by `Gunnar Jeschke`_. The libraries are
published under the GPL with his permission.


.. _Gunnar Jeschke: http://www.epr.ethz.ch/
.. _MMM: http://www.epr.ethz.ch/software/mmm-older-versions.html


.. Table:: Rotamer Libraries
   :name: table-libraries 

   +-------------------+------------+---------+-----+-----------------+
   |name               | spin-label |ID       |year |citation         |
   +===================+============+=========+=====+=================+   
   | 'MTSSL 298K 2011' | MTSL       |R1A_298K |2011 |[Polyhach2011]_  |
   +-------------------+------------+---------+-----+-----------------+
   | 'MTSSL 298K 2015' | MTSL       |R1A_298K |2015 |[Polyhach2011]_  |
   +-------------------+------------+---------+-----+-----------------+
   
The included libraries are listed in Table :ref:`table-libraries` (see
also :data:`rotcon.libraries.LIBRARIES`). The *name* in the table is
used to select the library with the ``--libname NAME`` option of the
scripts.


