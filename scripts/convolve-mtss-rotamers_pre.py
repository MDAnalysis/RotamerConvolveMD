#!/usr/bin/env python
# Convolve PRE spin label rotamers with MD trajectories.
# Copyright (c) 2011-2017 Philip Fowler and AUTHORS
# Published under the GNU Public Licence, version 2 (or higher)
#
# Includes a rotamer library for MTSSL at 298 K by Gunnar Jeschke,
# which is published under the same licence by permission
"""\
%prog [options] TOPOLOGY TRAJECTORY [TRAJECTORY ...]
       %prog [options] STRUCTURE

Calculate a distribution of PRE spin label distances from an MD trajectory
or an arbitrary ensemble of conformations by fitting Gunnar Jeschke's
rotamer library of MTSS (at 298 K).

Provide either a topology and trajectory data (multiple trajectories
are concatenated on the fly) or just a structure file such as a PDB or
GRO file.

For details on the method see Stelzl et al, J Mol Biol 426 (2014),
735-751, doi:10.1016/j.jmb.2013.10.024.
"""

import os.path

import MDAnalysis
import numpy as np
import time

import rotcon.convolve_pre

import logging
logger = logging.getLogger("MDAnalysis.app")


if __name__ == "__main__":
        from optparse import OptionParser

        parser = OptionParser(usage=__doc__)
        parser.add_option("--resid", type=int, dest="residue", default=None,
                          help="REQUIRED: residues to compute PRE distances")
        parser.add_option("--discard", dest="discardFrames", type=int, default=0,
                          help="discard the first N frames [%default]")
        parser.add_option("--clashDistance", dest="clashDistance", type=float, default=2.2,
                          help="the distance between heavy-atoms of the label and protein "
                          "within which they are assumed to be clashing [%default Angstroms]")
        parser.add_option("--outputRawDistances", dest="outputFileRawDistances",
                          default="outputRawDistances.dat", help="the path and name of "
                          "the output file with raw distances; the filename will "
                          "have resid 1 and resid 2 inserted before the extension [%default]")
        parser.add_option("--dcdfilenameAll", dest="dcdFilenameAll", metavar="FILENAME",
                          help="the path and stem of the DCD files of the fitted MTSL rotamers")
        parser.add_option("--dcdfilenameNoClashes", dest="dcdFilenameNoClashes", metavar="FILENAME",
                          help="the path and stem of the DCD files of the fitted MTSL rotamers "
                               "without clashes")
        parser.add_option("--libname", dest="libname", metavar="NAME", default="MTSSL 298K 2015",
                          help="name of the rotamer library [%default]")
        parser.add_option("--plotname", dest="plotname", metavar="FILENAME", default=None,
                          help="plot the histogram to FILENAME (the extensions determines the format) "
                               "By default <outputFile>.pdf.")
        parser.add_option("--no-plot", action="store_false", dest="with_plot", default=True,
                          help="suppress producing plots with --plotname")
        parser.add_option("--useNOelectron", action="store_true", dest="useNOelectron",
                          help="Set this flag, if the geometic midpoints of N1 and O1 atoms should be "
                          "used for distances measurements.")
        parser.add_option("--no-useNOelectron", action="store_false", dest="useNOelectron",
                          help="Set this flag, if N1 atoms should be used for distances measurements.")


        options, args = parser.parse_args()

        MDAnalysis.start_logging()
        logger.info(rotcon.get_info_string())
        logger.info(rotcon.get_license_string())
        logger.info(rotcon.get_citation_string())

        # load the reference protein structure
        try:
            proteinStructure = MDAnalysis.Universe(*args)
        except:
            logger.critical("protein structure and/or trajectory not correctly specified")
            raise
        if options.residue is None:
            raise ValueError("Provide residue id in --residue R1")

        logger.info("Loading trajectory data as Universe({0})".format(*args))

        if not options.dcdFilenameAll:
                options.dcdFilenameAll = options.outputFile + "-tmp"

        startTime = time.time()
        R = rotcon.convolve_pre.RotamerDistances(proteinStructure,
                                             options.residue,
                                             outputFileRawDistances=options.outputFileRawDistances,
                                             dcdFilenameAll=options.dcdFilenameAll,
                                             dcdFilenameNoClashes=options.dcdFilenameNoClashes,
                                             libname=options.libname,
                                             discardFrames=options.discardFrames,
                                             clashDistance=options.clashDistance,
                                             useNOelectron=options.useNOelectron)
        logger.info("DONE with analysis, elapsed time %6i s" % (int(time.time() - startTime)))

        if options.with_plot:
                if options.plotname is None:
                        root, ext = os.path.splitext(R.outputFile)
                        options.plotname = root + ".pdf"
                R.plot(filename=options.plotname, linewidth=2)

        MDAnalysis.stop_logging()
