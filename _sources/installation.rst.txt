.. -*- coding: utf-8 -*-

==============
 Installation
==============

Latest release
==============

The two following methods (using :program:`conda` or :program:`pip`)
will both install all dependencies.

conda
-----

`Conda`_ packages are made available on conda-forge and can be
installed with ::

  conda config --add channels conda-forge
  conda install RotamerConvolveMD

Later updates are installed with ::

  conda update RotamerConvolveMD

.. _conda: https://conda.io/docs/  


pip (PyPi)
----------

The `PyPi RotamerConvolveMD`_ package can be installed with pip_::

  pip install --upgrade RotamerConvolveMD

or in a user directory  ::

  pip install --upgrade --user RotamerConvolveMD


.. _PyPi RotamerConvolveMD: https://pypi.python.org/pypi/RotamerConvolveMD
.. _pip: https://pip.pypa.io


From sources
============

Clone the GitHub repository and use the master branch (the default)::

  git clone https://github.com/MDAnalysis/RotamerConvolveMD.git
  pip install RotamerConvolveMD/

or with the classic ``python setup.py`` method ::

  cd RotamerConvolveMD
  python setup.py install

Either form should automatically download required packages:

* numpy_
* MDAnalysis_

.. Note:: The package requires at least MDAnalysis 0.16.2.

Please see https://www.mdanalysis.org for hints if you have problems
with the automatic installation of MDAnalysis_.

.. _numpy: http://numpy.scipy.org/
.. _MDAnalysis: https://www.mdanalysis.org
