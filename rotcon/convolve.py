# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# Convolve MTSS rotamers with MD trajectory.
# Copyright (c) 2011-2017 Philip Fowler and AUTHORS
# Published under the GNU Public Licence, version 2 (or higher)
#
# Includes a rotamer library for MTSS at 298 K by Gunnar Jeschke,
# which is published under the same licence by permission.
from __future__ import absolute_import, division, print_function

import MDAnalysis
import MDAnalysis.analysis.align
import MDAnalysis.lib.NeighborSearch as KDNS
import MDAnalysis.analysis.distances

import numpy as np
import os.path

from . import library

import logging
logger = logging.getLogger("MDAnalysis.app")

class RotamerDistancesBase(object):
    @staticmethod
    def fit_rotamers(rotamers, protein, site_resid, dcdfile):
        """Produce a temporary trajectory of the rotamers.

        The backbone of the rotamers is fitted to the backbone of the
        spin labelled residues.
        """
        # create an ordered list allowing the rotamer to be fitted onto the backbone of the protein
        fittingSelection = (["name C", "name CA", "name N"],
                            ["protein and name C and resid {0}".format(site_resid),
                             "protein and name CA and resid {0}".format(site_resid),
                             "protein and name N and resid {0}".format(site_resid)
                             ])
        # fit the rotamer library onto the protein
        MDAnalysis.analysis.align.AlignTraj(rotamers, protein,
                                            select=fittingSelection, weights="mass",
                                            filename=dcdfile,
                                            verbose=False).run()
        return dcdfile

    def find_clashing_rotamers(self, fitted_rotamers, protein, site_resid):
        """Detect any rotamer that clashes with the protein."""
        # make a KD tree of the protein neighbouring atoms
        proteinNotSite = protein.select_atoms("protein and not (name H* or name [123]H or type H) "
                  "and not (resid " + str(site_resid) +
                  " or (resid " + str(site_resid-1) + " and (name N or name CA or name C or name O)) "
                  "or (resid " + str(site_resid+1) + " and (name N or name CA or name C or name O)))")
        proteinNotSiteLookup = KDNS.AtomNeighborSearch(proteinNotSite)

        rotamerSel = fitted_rotamers.select_atoms("not (name H* or name [123]H or type H)")

        rotamer_clash = []
        for rotamer in fitted_rotamers.trajectory:
            bumps = proteinNotSiteLookup.search(rotamerSel, self.clashDistance)
            rotamer_clash.append(bool(bumps))
        return rotamer_clash, np.sum(rotamer_clash)


class RotamerDistances(RotamerDistancesBase):
    """Calculation of distance distributions between two spin labels."""
    def __init__(self, *args, **kwargs):
        """RotamerDistances(universe, residue_list, **kwargs)

        :Arguments:
           *universe*
              :class:`MDAnalysis.Universe`
           *residue_list*
              list of two residue numbers ``(r1, r2)`` that indicate
              the labelled sites

        :Keywords:
           *dcdFilename*
              name of the temporary files with rotamers fitted [``'trj'``]
           *dcdFilenameNoClashes*
              name of the temporary files with rotamers fitted [``'trj'``]
           *outputFile*
              stem of the name of the file containing the distance histogram
              (the final name will be ``<outputFile><resid_1>-<resid_2>.dat``
              [``'distances'``]
           *outputFileRawDistances*
              stem of the name of the file containing the raw distances
              (the final name will be ``<outputFile><resid_1>-<resid_2>_distances.dat``
              [``'distances'``]
           *libname*
              library name; the library is loaded with
              :class:`rotcon.library.RotamerLibrary` [``'MTSSL 298K'``]
           *discardFrames*
              skip initial frames < *discardFrames* [``0``]
           *clashDistance*
              discard rotamer if any distance between rotamer atoms
              and protein atoms is < *clashDistance*. Values down to
              1.5 Å are reasonable. The default is conservative. [``2.2`` Å]
           *histogramBins*
             tuple ``(dmin, dmax, delta)`` in Ångström [``(0.0, 100.0, 1.0)``]
           *useNOelectron*
            True = geometic midpoints of N1 and O1 atoms used for distance calculation
            False = N1 atoms used distance measurements,
        """
        proteinStructure = args[0]
        residues = args[1]
        if len(residues) != 2:
            raise ValueError("The residue_list must contain exactly 2 residue numbers: current "
                             "value {0}.".format(residues))

        outputFile, ext = os.path.splitext(kwargs.pop('outputFile', 'distances'))
        ext = ext or ".dat"
        self.outputFile = "{0}-{1[0]}-{1[1]}{2}".format(outputFile, residues, ext)

        outputFileRawDistances, ext = os.path.splitext(kwargs.pop('outputFileRawDistances', 'distances'))
        ext = ext or ".dat"
        self.outputFileRawDistances = "{0}-{1[0]}-{1[1]}-rawDistances{2}".format(outputFileRawDistances, residues, ext)


        dcdFilename, ext = os.path.splitext(kwargs.pop('dcdFilename', 'trj'))
        ext = ext or ".dcd"
        tmptrj = ["{0}-{1[0]}-{1[1]}-1{2}".format(dcdFilename, residues, ext),  # or make this temp files?
                  "{0}-{1[0]}-{1[1]}-2{2}".format(dcdFilename, residues, ext),  # or make this temp files?
                  ]

        dcdFilenameNoClashes, ext = os.path.splitext(kwargs.pop('dcdFilenameNoClashes', 'trj'))
        ext = ext or ".dcd"
        tmptrjNoClashes = ["{0}-{1[0]}-{1[1]}-noClashes-1{2}".format(dcdFilenameNoClashes, residues, ext),  # or make this temp files?
                           "{0}-{1[0]}-{1[1]}-noClashes-2{2}".format(dcdFilenameNoClashes, residues, ext),  # or make this temp files?
                          ]

        kwargs.setdefault('discardFrames', 0)
        self.clashDistance = kwargs.pop('clashDistance', 2.2)  # Ångström
        histogramBins = kwargs.pop('histogramBins', (0.0, 100.0, 1.0))
        useNOelectron = kwargs.pop('useNOelectron', True)

        self.lib = library.RotamerLibrary(kwargs.get('libname', 'MTSSL 298K'))

        # setup the main lists
        distances = []
        weights = []

        logger.info("Starting rotamer distance analysis of trajectory "
                    "{0}...".format(proteinStructure.trajectory.filename))
        logger.info("clashDistance = {0} A; rotamer library = '{1}'".format(self.clashDistance, self.lib.name))
        logger.debug("Temporary trajectories for rotamers 1 and 2 "
                     "(only last frame of MD trajectory): {0[0]} and {0[1]}".format(tmptrj))
        logger.debug("Results will be written to {0}.".format(self.outputFile))

        progressmeter = MDAnalysis.log.ProgressMeter(proteinStructure.trajectory.n_frames, interval=1)
        for protein in proteinStructure.trajectory:
            progressmeter.echo(protein.frame)
            if protein.frame < kwargs['discardFrames']:
                continue
            # define the atoms used to fit the rotamers. Note that an
            # ordered list has to be created as the ordering of C CA N is
            # different in both. Fit the rotamers onto the protein:
            self.fit_rotamers(self.lib.rotamers, proteinStructure, residues[0], tmptrj[0])
            rotamersSite1 = MDAnalysis.Universe(self.lib.rotamers.filename, tmptrj[0])
            (rotamer1_clash, rotamer1_clash_total) = self.find_clashing_rotamers(rotamersSite1,
                                                                                 proteinStructure, residues[0])

            self.fit_rotamers(self.lib.rotamers, proteinStructure, residues[1], tmptrj[1])
            rotamersSite2 = MDAnalysis.Universe(self.lib.rotamers.filename, tmptrj[1])
            (rotamer2_clash, rotamer2_clash_total) = self.find_clashing_rotamers(rotamersSite2, proteinStructure,
                                                                                 residues[1])

            # define the atoms to measure the distances between
            rotamer1nitrogen = rotamersSite1.select_atoms("name N1")
            rotamer2nitrogen = rotamersSite2.select_atoms("name N1")
            rotamer1oxygen = rotamersSite1.select_atoms("name O1")
            rotamer2oxygen = rotamersSite2.select_atoms("name O1")

            # define the atoms to measure the distances between
            rotamer1All = rotamersSite1.select_atoms("all")
            rotamer2All = rotamersSite2.select_atoms("all")

            with MDAnalysis.Writer("{}".format(tmptrjNoClashes[0]), rotamer1All.n_atoms) as S1:
                with MDAnalysis.Writer("{}".format(tmptrjNoClashes[1]), rotamer2All.n_atoms) as S2:
                    # loop over all the rotamers on the first site
                    for rotamer1 in rotamersSite1.trajectory:
                        if not rotamer1_clash[rotamer1.frame]:
                            # loop over all the rotamers on the second site
                            for rotamer2 in rotamersSite2.trajectory:
                                if not rotamer2_clash[rotamer2.frame]:
                                    S1.write(rotamersSite1.atoms)
                                    S2.write(rotamersSite2.atoms)
                                    # measure and record the distance
                                    (a, b, distance_nitrogen) = MDAnalysis.analysis.distances.dist(rotamer1nitrogen, rotamer2nitrogen)
                                    if useNOelectron == True:
                                        (a, b, distance_oxygen) = MDAnalysis.analysis.distances.dist(rotamer1oxygen, rotamer2oxygen)
                                        distance = np.mean([distance_nitrogen[0], distance_oxygen[0]])
                                    elif useNOelectron == False:
                                        distance = distance_nitrogen[0]
                                    distances.append(distance)
                                    # create the weights list
                                    try:
                                        weight = self.lib.weights[rotamer1.frame] * self.lib.weights[rotamer2.frame]
                                    except IndexError:
                                        logger.error("oppps: no weights for rotamer 1 #{0} - "
                                                    "rotamer 2 #{1}".format(rotamer1.frame, rotamer2.frame))
                                    weights.append(weight)

        # check that at least two distances have been measured
        if len(distances) < 2:
            logger.critical("no distances found between the spin pair!")
            raise RuntimeError("no distances found between the spin pair!")  # should this really be an exception?

        # calculate Nbins and min and max so that we cover at least
        # the requested lower and upper bounds with the given fixed
        # bin width
        bins = MDAnalysis.lib.util.fixedwidth_bins(histogramBins[2], histogramBins[0], histogramBins[1])
        # use numpy to histogram the distance data, weighted appropriately
        (a, b) = np.histogram(distances, weights=weights, density=True, bins=bins['Nbins'],
                              range=(bins['min'], bins['max']))

        with open(self.outputFile, 'w') as OUTPUT:
            for (i, j) in enumerate(a):
                OUTPUT.write("%6.2f %8.3e\n" % ((0.5*(b[i] + b[i+1])), j))
        logger.info("Distance distribution for residues {0[0]} - {0[1]} "
                    "was written to {1}".format(residues, self.outputFile))

        with open(self.outputFileRawDistances, 'w') as OUTPUT:
            for distance in distances:
                OUTPUT.write("%8.3e\n" % (distance))


    def plot(self, **kwargs):
        """Load data file and plot"""

        import matplotlib.pyplot as plt

        filename = kwargs.pop('filename', None)
        fig = kwargs.pop('fig', None)
        if fig is None:
            fig = plt.figure(figsize=(5, 5))
        ax = kwargs.pop('ax', None)
        if ax is None:
            ax = fig.add_subplot(111)

        dist, prob = np.loadtxt(self.outputFile, unpack=True)
        ax.plot(dist, prob, **kwargs)
        ax.set_xlabel(r"spin-label distance $d$ ($\AA$)")
        ax.set_ylabel("probability density")

        if filename:
            ax.figure.savefig(filename)
            logger.info("Plotted distance distribution to {0}".format(filename))

        return ax
