Installation instructions
=========================

From the sources::

     python setup.py install --user

From a tarball::

     pip install --user rotamers-1.0.tar.gz

or ::

     easy_install --user rotamers-1.0.tar.gz

Either form should automatically download required packages:

* numpy_
* MDAnalysis_

.. Note:: The package requires MDAnalysis 0.8 or a pre-release version
          as it makes use of some new features.

Please see http://www.mdanalysis.org for hints if you have problems
with the automatic installation of MDAnalysis_.

See README.rst for further notes.

.. _numpy: http://numpy.scipy.org/
.. _MDAnalysis: http://www.mdanalysis.org
