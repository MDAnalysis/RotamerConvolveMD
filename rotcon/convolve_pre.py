# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# Convolve MTSS rotamers with MD trajectory.
# Copyright (c) 2011-2013 Philip Fowler, Oliver Beckstein
# Published under the GNU Public Licence, version 2 (or higher)
#
# Includes a rotamer library for MTSS at 298 K by Gunnar Jeschke,
# which is published under the same licence by permission.

import MDAnalysis
import MDAnalysis.analysis.align
import MDAnalysis.lib.NeighborSearch as KDNS
import MDAnalysis.analysis.distances

import numpy as np
import os.path

import rotcon.library

import logging
logger = logging.getLogger("MDAnalysis.app")


def rms_fit_trj(*args, **kwargs):
    """Silenced :func:`MDAnalysis.analysis.align.rms_fit_trj`"""
    kwargs['quiet'] = True
    return MDAnalysis.analysis.align.rms_fit_trj(*args, **kwargs)


class RotamerDistances(object):
    """Calculation of distance distributions between two spin labels."""
    def __init__(self, *args, **kwargs):
        """RotamerDistances(universe, residue_list, **kwargs)

        :Arguments:
           *universe*
              :class:`MDAnalysis.Universe`
           *residue*
              residue number ``(r1)`` that indicate the labelled sites

        :Keywords:
           *dcdFilenameAll*
              name of the temporary files with rotamers fitted [``'trj'``]
           *dcdFilenameNoClashes*
              name of the temporary files with rotamers fitted [``'trj'``]
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
        """
        proteinStructure = args[0]
        residue = args[1]


        outputFileRawDistances, ext = os.path.splitext(kwargs.pop('outputFileRawDistances', 'distances'))
        ext = ext or ".dat"
        self.outputFileRawDistances = "{0}-{1}-rawDistances{2}".format(outputFileRawDistances, 
                                                                        residue, ext)
        
        
        dcdFilenameAll, ext = os.path.splitext(kwargs.pop('dcdFilenameAll', 'trj'))
        ext = ext or ".dcd"
        tmptrj = "{0}-{1}-all{2}".format(dcdFilenameAll, residue, ext)
        
        dcdFilenameNoClashes, ext = os.path.splitext(kwargs.pop('dcdFilenameNoClashes', 'trj'))
        ext = ext or ".dcd"
        tmptrjNoClashes = "{0}-{1}-noClashes{2}".format(dcdFilenameNoClashes, residue, ext)


        kwargs.setdefault('discardFrames', 0)
        self.clashDistance = kwargs.pop('clashDistance', 2.2)  # Ångström

        self.lib = rotcon.library.RotamerLibrary(kwargs.get('libname', 'MTSSL 298K'))

        # setup the main lists
        distances = []
        weights = []

        logger.info("Starting rotamer distance analysis of trajectory "
                    "{0}...".format(proteinStructure.trajectory.filename))
        logger.info("clashDistance = {0} A; rotamer library = '{1}'".format(self.clashDistance, 
                                                                            self.lib.name))
        logger.debug("Temporary trajectories for rotamers 1 and 2 "
                     "(only last frame of MD trajectory): {0[0]} and {0[1]}".format(tmptrj))

        progressmeter = MDAnalysis.log.ProgressMeter(proteinStructure.trajectory.n_frames, interval=1)
        for protein in proteinStructure.trajectory:
            progressmeter.echo(protein.frame)
            if protein.frame < kwargs['discardFrames']:
                continue
            # define the atoms used to fit the rotamers. Note that an
            # ordered list has to be created as the ordering of C CA N is
            # different in both. Fit the rotamers onto the protein:
            self.fit_rotamers(self.lib.rotamers, proteinStructure, residue, tmptrj)
            rotamersSite1 = MDAnalysis.Universe(self.lib.rotamers.filename, tmptrj)
            (rotamer1_clash, rotamer1_clash_total) = self.find_clashing_rotamers(rotamersSite1,
                                                                    proteinStructure, residue)
            proteinNH = proteinStructure.select_atoms("protein and name HN") # or HN
            
            # define the atoms to measure the distances between
            rotamer1nitrogen = rotamersSite1.select_atoms("name N1")

            # define the atoms to measure the distances between
            rotamer1All = rotamersSite1.select_atoms("all")
            
            with MDAnalysis.Writer("{}".format(tmptrjNoClashes), rotamer1All.n_atoms) as S1:
                # loop over all the rotamers on the first site
                for rotamer1 in rotamersSite1.trajectory:
                    if not rotamer1_clash[rotamer1.frame]:
                        S1.write(rotamersSite1.atoms)
                        for i in range(0,proteinNH.n_atoms):
                            nh = proteinNH.select_atoms("resid {}".format(i+1))
                            if nh.n_atoms == 1:
                                (a, b, distance) = MDAnalysis.analysis.distances.dist(rotamer1nitrogen, 
                                                                                      nh)
                                distances.append([i+1, distance[0]])



        # check that at least two distances have been measured
        if len(distances) < 2:
            logger.critical("no distances found between the spin pair!")
            raise RuntimeError("no distances found between the spin pair!")  
        # should this really be an exception?

        
        with open(self.outputFileRawDistances, 'w') as OUTPUT:
            for distance in distances:
                OUTPUT.write("{0[0]}\t{0[1]}\n".format(distance))

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

        data = np.loadtxt(self.outputFileRawDistances, unpack=True)

        dataResidues = dict()
        for j in range(0, len(data[0])):
            if int(data[0][j]) in dataResidues: 
                dataResidues[int(data[0][j])].append(data[1][j])
            else: 
                dataResidues[int(data[0][j])] = [data[1][j]]

        for data in dataResidues:
            ax.scatter(data, np.min(dataResidues[data]), color='blue')
            ax.scatter(data, np.max(dataResidues[data]), color='red')

        if filename:
            ax.figure.savefig(filename)
            logger.info("Plotted min and max distances to {0}".format(filename))

        return ax

    def fit_rotamers(self, rotamers, protein, site_resid, dcdfile):
        """Produce a temporary trajectory of the rotamers.

        The backbone of the rotamers is fitted to the backbone of the
        spin labelled residue.
        """
        # create an ordered list allowing the rotamer to be fitted onto the backbone of the protein
        fittingSelection = (["name C", "name CA", "name N"],
                            ["protein and name C and resid {0}".format(site_resid),
                             "protein and name CA and resid {0}".format(site_resid),
                             "protein and name N and resid {0}".format(site_resid)
                             ])
        # fit the rotamer library onto the protein
        rms_fit_trj(rotamers, protein, select=fittingSelection, mass_weighted=True, filename=dcdfile)
        return dcdfile

    def find_clashing_rotamers(self, fitted_rotamers, protein, site_resid):
        """Detect any rotamer that clashes with the protein."""
        # make a KD tree of the protein neighbouring atoms
        proteinNotSite = protein.selectAtoms("protein and not name H* and not (resid " + str(site_resid) +
                                             " or (resid " + str(site_resid-1) + " and (name C or name O)) "
                                                                                 "or (resid " + str(site_resid+1)
                                             + " and name N))")
        proteinNotSiteLookup = KDNS.AtomNeighborSearch(proteinNotSite)

        rotamerSel = fitted_rotamers.selectAtoms("not name H*")

        rotamer_clash = []
        for rotamer in fitted_rotamers.trajectory:
            bumps = proteinNotSiteLookup.search(rotamerSel, self.clashDistance)
            rotamer_clash.append(bool(bumps))
        return rotamer_clash, np.sum(rotamer_clash)
