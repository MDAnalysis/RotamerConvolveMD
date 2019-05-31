# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# Convolve MTSS rotamers with MD trajectory.
# Copyright (c) 2011-2017 Philip Fowler and AUTHORS
# Published under the GNU Public Licence, version 2 (or higher)
#
# Includes a rotamer library for MTSS at 298 K by Gunnar Jeschke,
# which is published under the same licence by permission.
"""\
Rotamer library handling
========================

:mod:`rotamers.library` contains the data (:data:`LIBRARIES`) to load
a rotamer library, represented by a :class:`RotamerLibrary`.

"""
from __future__ import absolute_import, division, print_function

import MDAnalysis

import logging
logger = logging.getLogger("MDAnalysis.app")

import numpy as np
import os.path
import pkg_resources

#: Name of the directory in the package that contains the library data.
LIBDIR = "data"

# This could be turned into a YAML file.
#: Registry of libraries, indexed by name.
LIBRARIES = {
    'MTSSL 298K 2011': {
        'topology': "rotamer1_R1A_298K_2011.pdb",
        'ensemble': "rotamer1_R1A_298K_2011.dcd",
        'populations': "R1A_298K_populations_2011.dat",
        'author': "Gunnar Jeschke",
        'licence': "GPL v2",
        'citation': "Polyhach Y, Bordignon E, Jeschke G. "
                    "Phys Chem Chem Phys. 2011; 13(6):2356-2366. doi: 10.1039/c0cp01865a",
    },
    'MTSSL 298K 2015': {
        'topology': "rotamer1_R1A_298K_2015.pdb",
        'ensemble': "rotamer1_R1A_298K_2015.dcd",
        'populations': "R1A_298K_populations_2015.dat",
        'author': "Gunnar Jeschke",
        'licence': "GPL v2",
        'citation': "Polyhach Y, Bordignon E, Jeschke G. "
                    "Phys Chem Chem Phys. 2011; 13(6):2356-2366. doi: 10.1039/c0cp01865a",
        'information': "updated version of the MTSSL rotamer library from 2015"
    },
}


def find_file(filename, pkglibdir=LIBDIR):
    """Return full path to file *filename*.

    1) If the path exists, return rooted canonical path.
    2) Try to find the file in the package.
    """
    if os.path.exists(filename):
        return MDAnalysis.core.util.realpath(filename)
    return pkg_resources.resource_filename(__name__, os.path.join(pkglibdir, filename))


class RotamerLibrary(object):
    """Rotamer library

    The library makes available the attributes :attr:`rotamers`, and :attr:`weights`.

    .. attribute:: rotamers
       :class:`MDAnalysis.core.AtomGroup.Universe` instance that
       records all rotamers as a trajectory

    .. attribute:: weights
       NumPy array containing the population of each rotomer.

    .. attribute:: name
       Name of the library.

    .. attribute:: lib
       Dictionary containing the file names and meta data for the library :attr:`name`.
    """

    def __init__(self, name):
        """RotamerLibrary(name)

        :Arguments:
           *name*
              name of the library (must exist in the registry of libraries, :data:`LIBRARIES`)
        """
        self.name = name
        self.lib = {}
        try:
            self.lib.update(LIBRARIES[name])  # make a copy
        except KeyError:
            raise ValueError("No rotamer library with name {0} known: must be one of {1}".format(name,
                                                                                                 LIBRARIES.keys()))
        logger.info("Using rotamer library '{0}' by {1[author]}".format(self.name, self.lib))
        logger.info("Please cite: {0[citation]}".format(self.lib))
        # adjust paths
        for k in 'ensemble', 'topology', 'populations':
            self.lib[k] = find_file(self.lib[k])
        logger.debug("[rotamers] ensemble = {0[ensemble]} with topology = {0[topology]}".format(self.lib))
        logger.debug("[rotamers] populations = {0[populations]}".format(self.lib))

        self.rotamers = MDAnalysis.Universe(self.lib['topology'], self.lib['ensemble'])
        self.weights = self.read_rotamer_weights(self.lib['populations'])

        if len(self.rotamers.trajectory) != len(self.weights):
            err_msg = "Discrepancy between number of rotamers ({0}) and weights ({1})".format(
                len(self.rotamers.trajectory), len(self.weights))
            logger.critical(err_msg)
            raise ValueError(err_msg)

    def read_rotamer_weights(self, filename):
        """read in the rotamer weights from *filename*

        There is one weight per conformer (frame) in the trajectory.
        """
        return np.loadtxt(filename)

    def __repr__(self):
        return "<RotamerLibrary '{0}' by {1} with {2} rotamers>".format(self.name, self.lib['author'],
                                                                        len(self.weights))
