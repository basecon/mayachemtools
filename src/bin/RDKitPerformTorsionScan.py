#!/usr/bin/env python
#
# File: RDKitPerformTorsionScan.py
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2020 Manish Sud. All rights reserved.
#
# The functionality available in this script is implemented using RDKit, an
# open source toolkit for cheminformatics developed by Greg Landrum.
#
# This file is part of MayaChemTools.
#
# MayaChemTools is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# MayaChemTools is distributed in the hope that it will be useful, but without
# any warranty; without even the implied warranty of merchantability of fitness
# for a particular purpose.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MayaChemTools; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation Inc., 59 Temple Place, Suite 330,
# Boston, MA, 02111-1307, USA.
#

from __future__ import print_function

# Add local python path to the global path and import standard library modules...
import os
import sys;  sys.path.insert(0, os.path.join(os.path.dirname(sys.argv[0]), "..", "lib", "Python"))
import time
import re
import glob
import multiprocessing as mp

import matplotlib.pyplot as plt
import seaborn as sns

# RDKit imports...
try:
    from rdkit import rdBase
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdFMCS
    from rdkit.Chem import rdMolTransforms
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import RDKit module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your RDKit environment and try again.\n\n")
    sys.exit(1)

# MayaChemTools imports...
try:
    from docopt import docopt
    import MiscUtil
    import RDKitUtil
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import MayaChemTools module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your MayaChemTools environment and try again.\n\n")
    sys.exit(1)

ScriptName = os.path.basename(sys.argv[0])
Options = {}
OptionsInfo = {}

def main():
    """Start execution of the script"""
    
    MiscUtil.PrintInfo("\n%s (RDK v%s; %s): Starting...\n" % (ScriptName, rdBase.rdkitVersion, time.asctime()))
    
    (WallClockTime, ProcessorTime) = MiscUtil.GetWallClockAndProcessorTime()
    
    # Retrieve command line arguments and options...
    RetrieveOptions()
    
    # Process and validate command line arguments and options...
    ProcessOptions()
    
    # Perform actions required by the script...
    PerformTorsionScan()

    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformTorsionScan():
    """Perform torsion scan."""
    
    # Setup a molecule reader for input file...
    MiscUtil.PrintInfo("\nProcessing file %s..." % OptionsInfo["Infile"])
    OptionsInfo["InfileParams"]["AllowEmptyMols"] = True
    Mols  = RDKitUtil.ReadMolecules(OptionsInfo["Infile"], **OptionsInfo["InfileParams"])
    
    PlotExt = OptionsInfo["OutPlotParams"]["OutExt"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Options["--outfile"])
    MiscUtil.PrintInfo("Generating output files %s_*.sdf, %s_*Torsion*Match*.sdf, %s_*Torsion*Match*Energies.csv, %s_*Torsion*Match*Plot.%s..." % (FileName, FileName, FileName, FileName, PlotExt))

    MolCount, ValidMolCount, MinimizationFailedCount, TorsionsMissingCount, TorsionsScanFailedCount = ProcessMolecules(Mols)
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of molecules failed during initial minimization: %d" % MinimizationFailedCount)
    MiscUtil.PrintInfo("Number of molecules without any matched torsions: %d" % TorsionsMissingCount)
    MiscUtil.PrintInfo("Number of molecules failed during torsion scan: %d" % TorsionsScanFailedCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount + TorsionsMissingCount + MinimizationFailedCount + TorsionsScanFailedCount))

def ProcessMolecules(Mols):
    """Process molecules to perform torsion scan. """

    if OptionsInfo["MPMode"]:
        return ProcessMoleculesUsingMultipleProcesses(Mols)
    else:
        return ProcessMoleculesUsingSingleProcess( Mols)

def ProcessMoleculesUsingSingleProcess(Mols):
    """Process molecules to perform torsion scan using a single process."""

    MolInfoText = "first molecule"
    if not OptionsInfo["FirstMolMode"]:
        MolInfoText = "all molecules"

    if OptionsInfo["TorsionMinimize"]:
        MiscUtil.PrintInfo("\nPeforming torsion scan on %s by generating conformation ensembles for specific torsion angles and constrained energy minimization of the ensembles..." % (MolInfoText))
    else:
        MiscUtil.PrintInfo("\nPeforming torsion scan on %s by skipping generation of conformation ensembles for specific torsion angles and constrained energy minimization of the ensembles..." % (MolInfoText))
    
    SetupTorsionsPatternsInfo()
    
    (MolCount, ValidMolCount, TorsionsMissingCount, MinimizationFailedCount, TorsionsScanFailedCount) = [0] * 5

    for Mol in Mols:
        MolCount += 1

        if OptionsInfo["FirstMolMode"] and MolCount > 1:
            MolCount -= 1
            break
        
        if Mol is None:
            continue
        
        if RDKitUtil.IsMolEmpty(Mol):
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolCount)
                MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
            continue
        ValidMolCount += 1

        Mol, MinimizationCalcStatus, TorsionsMatchStatus, TorsionsScanCalcStatus = PerformMinimizationAndTorsionScan(Mol, MolCount)
        
        if not MinimizationCalcStatus:
            MinimizationFailedCount += 1
            continue
        
        if not TorsionsMatchStatus:
            TorsionsMissingCount += 1
            continue

        if not TorsionsScanCalcStatus:
            TorsionsScanFailedCount += 1
            continue

    return (MolCount, ValidMolCount, MinimizationFailedCount, TorsionsMissingCount, TorsionsScanFailedCount)

def ProcessMoleculesUsingMultipleProcesses(Mols):
    """Process molecules to perform torsion scan using multiprocessing."""

    MolInfoText = "first molecule"
    if not OptionsInfo["FirstMolMode"]:
        MolInfoText = "all molecules"

    if OptionsInfo["TorsionMinimize"]:
        MiscUtil.PrintInfo("\nPeforming torsion scan on %s using multiprocessing by generating conformation ensembles for specific torsion angles and constrained energy minimization of the ensembles..." % (MolInfoText))
    else:
        MiscUtil.PrintInfo("\nPeforming torsion scan %s using multiprocessing by skipping generation of conformation ensembles for specific torsion angles and constrained energy minimization of the ensembles..." % (MolInfoText))
        
    
    MPParams = OptionsInfo["MPParams"]
    
    # Setup data for initializing a worker process...
    MiscUtil.PrintInfo("Encoding options info...")
    
    InitializeWorkerProcessArgs = (MiscUtil.ObjectToBase64EncodedString(Options), MiscUtil.ObjectToBase64EncodedString(OptionsInfo))

    if OptionsInfo["FirstMolMode"]:
        Mol = Mols[0]
        Mols = [Mol]

    # Setup a encoded mols data iterable for a worker process...
    WorkerProcessDataIterable = RDKitUtil.GenerateBase64EncodedMolStrings(Mols)

    # Setup process pool along with data initialization for each process...
    MiscUtil.PrintInfo("\nConfiguring multiprocessing using %s method..." % ("mp.Pool.imap()" if re.match("^Lazy$", MPParams["InputDataMode"], re.I) else "mp.Pool.map()"))
    MiscUtil.PrintInfo("NumProcesses: %s; InputDataMode: %s; ChunkSize: %s\n" % (MPParams["NumProcesses"], MPParams["InputDataMode"], ("automatic" if MPParams["ChunkSize"] is None else MPParams["ChunkSize"])))
    
    ProcessPool = mp.Pool(MPParams["NumProcesses"], InitializeWorkerProcess, InitializeWorkerProcessArgs)
    
    # Start processing...
    if re.match("^Lazy$", MPParams["InputDataMode"], re.I):
        Results = ProcessPool.imap(WorkerProcess, WorkerProcessDataIterable, MPParams["ChunkSize"])
    elif re.match("^InMemory$", MPParams["InputDataMode"], re.I):
        Results = ProcessPool.map(WorkerProcess, WorkerProcessDataIterable, MPParams["ChunkSize"])
    else:
        MiscUtil.PrintError("The value, %s, specified for \"--inputDataMode\" is not supported." % (MPParams["InputDataMode"]))

    (MolCount, ValidMolCount, TorsionsMissingCount, MinimizationFailedCount, TorsionsScanFailedCount) = [0] * 5
    
    for Result in Results:
        MolCount += 1
        
        MolIndex, EncodedMol, MinimizationCalcStatus, TorsionsMatchStatus, TorsionsScanCalcStatus = Result
        
        if EncodedMol is None:
            continue
        ValidMolCount += 1
    
        if not MinimizationCalcStatus:
            MinimizationFailedCount += 1
            continue
        
        if not TorsionsMatchStatus:
            TorsionsMissingCount += 1
            continue

        if not TorsionsScanCalcStatus:
            TorsionsScanFailedCount += 1
            continue
        
        Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
        
    return (MolCount, ValidMolCount, MinimizationFailedCount, TorsionsMissingCount, TorsionsScanFailedCount)
    
def InitializeWorkerProcess(*EncodedArgs):
    """Initialize data for a worker process."""
    
    global Options, OptionsInfo

    MiscUtil.PrintInfo("Starting process (PID: %s)..." % os.getpid())

    # Decode Options and OptionInfo...
    Options = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[0])
    OptionsInfo = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[1])

    # Initialize torsion patterns info...
    SetupTorsionsPatternsInfo()

def WorkerProcess(EncodedMolInfo):
    """Process data for a worker process."""

    MolIndex, EncodedMol = EncodedMolInfo

    (MinimizationCalcStatus, TorsionsMatchStatus, TorsionsScanCalcStatus) = [False] * 3
    
    if EncodedMol is None:
        return [MolIndex, None, MinimizationCalcStatus, TorsionsMatchStatus, TorsionsScanCalcStatus]
    
    Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
    if RDKitUtil.IsMolEmpty(Mol):
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, (MolIndex + 1))
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
        return [MolIndex, None, MinimizationCalcStatus, TorsionsMatchStatus, TorsionsScanCalcStatus]
    
    Mol, MinimizationCalcStatus, TorsionsMatchStatus, TorsionsScanCalcStatus = PerformMinimizationAndTorsionScan(Mol, (MolIndex + 1))
    
    return [MolIndex, RDKitUtil.MolToBase64EncodedMolString(Mol, PropertyPickleFlags = Chem.PropertyPickleOptions.MolProps | Chem.PropertyPickleOptions.PrivateProps), MinimizationCalcStatus, TorsionsMatchStatus, TorsionsScanCalcStatus]

def PerformMinimizationAndTorsionScan(Mol, MolNum = None):
    """Peform minimization and torsions scan."""
    
    if not OptionsInfo["Infile3D"]:
        Mol, MinimizationCalcStatus = MinimizeMolecule(Mol, MolNum)
        if not MinimizationCalcStatus:
            return (Mol, False, False, False)
        
    TorsionsMolInfo = SetupTorsionsMolInfo(Mol, MolNum)
    if TorsionsMolInfo["NumOfMatches"] == 0:
        return (Mol, True, False, False)
    
    Mol, ScanCalcStatus = ScanAllTorsionsInMol(Mol, TorsionsMolInfo, MolNum)
    if not ScanCalcStatus:
        return (Mol, True, True, False)
    
    return (Mol, True, True, True)

def ScanAllTorsionsInMol(Mol, TorsionsMolInfo,  MolNum = None):
    """Peform scans on all torsions in a molecule."""
    
    if TorsionsMolInfo["NumOfMatches"] == 0:
        return Mol, True

    Mol = AddHydrogens(Mol)
    MolName = RDKitUtil.GetMolName(Mol, MolNum)
    
    FirstTorsionMode = OptionsInfo["FirstTorsionMode"]
    TorsionsPatternsInfo = OptionsInfo["TorsionsPatternsInfo"]

    TorsionPatternCount, TorsionScanCount, TorsionMatchCount = [0] * 3
    TorsionMaxMatches = OptionsInfo["TorsionMaxMatches"]
    
    for TorsionID in TorsionsPatternsInfo["IDs"]:
        TorsionPatternCount +=  1
        TorsionPattern = TorsionsPatternsInfo["Pattern"][TorsionID]
        TorsionPatternMol = TorsionsPatternsInfo["Mol"][TorsionID]
        
        TorsionsMatches = TorsionsMolInfo["Matches"][TorsionID]

        if TorsionsMatches is None:
            continue
        
        if FirstTorsionMode and TorsionPatternCount > 1:
            if not OptionsInfo["QuietMode"]:
                MiscUtil.PrintWarning("Already scaned first torsion pattern, \"%s\" for molecule %s during \"%s\" value of \"--modeTorsions\" option . Abandoning torsion scan...\n" % (TorsionPattern, MolName, OptionsInfo["ModeTorsions"]))
            break

        for Index, TorsionMatches in enumerate(TorsionsMatches):
            TorsionMatchNum = Index + 1
            TorsionMatchCount +=  1

            if TorsionMatchCount > TorsionMaxMatches:
                if not OptionsInfo["QuietMode"]:
                    MiscUtil.PrintWarning("Already scaned a maximum of %s torsion matches for molecule %s specified by \"--torsionMaxMatches\" option. Abandoning torsion scan...\n" % (TorsionMaxMatches, MolName))
                break
            
            Mol, TorsionScanStatus,  TorsionMols, TorsionEnergies, TorsionAngles = ScanSingleTorsionInMol(Mol, TorsionID, TorsionPattern, TorsionPatternMol, TorsionMatches,  MolNum)
            if not TorsionScanStatus:
                continue
            
            TorsionScanCount +=  1
            GenerateOutputFiles(Mol, MolNum, TorsionID, TorsionMatchNum, TorsionMols, TorsionEnergies, TorsionAngles)
        
        if TorsionMatchCount > TorsionMaxMatches:
            break
    
    if  OptionsInfo["RemoveHydrogens"]:
        Mol = Chem.RemoveHs(Mol)

    if TorsionScanCount:
        GenerateStartingTorsionScanStructureOutfile(Mol, MolNum)

    Status = True if TorsionScanCount else False
    
    return (Mol, Status)
    
def ScanSingleTorsionInMol(Mol, TorsionID, TorsionPattern, TorsionPatternMol, TorsionMatches, MolNum = None):
    """Perform torsion scan for a molecule along with constrained energy minimization."""

    StartAngle = OptionsInfo["TorsionStart"]
    StopAngle = OptionsInfo["TorsionStop"]
    StepSize = OptionsInfo["TorsionStep"]
    TorsionMinimize = OptionsInfo["TorsionMinimize"]
    
    AtomIndex1, AtomIndex2, AtomIndex3, AtomIndex4 = TorsionMatches

    TorsionMols = []
    TorsionEnergies = []
    TorsionAngles = []
    
    Angles = [Angle for Angle in range(StartAngle, StopAngle, StepSize)]
    Angles.append(StopAngle)

    for Angle in Angles:
        TorsionMol = Chem.Mol(Mol)
        TorsionMolConf = TorsionMol.GetConformer(0)

        rdMolTransforms.SetDihedralDeg(TorsionMolConf, AtomIndex1, AtomIndex2, AtomIndex3, AtomIndex4, Angle)

        if TorsionMinimize:
            # Perform constrained minimization...
            TorsionMatchesMol = RDKitUtil.MolFromSubstructureMatch(TorsionMol, TorsionPatternMol, TorsionMatches)
            TorsionMol, CalcStatus, Energy = ConstrainAndMinimizeMolecule(TorsionMol, TorsionMatchesMol, MolNum)
            if not CalcStatus:
                if not OptionsInfo["QuietMode"]:
                    MolName = RDKitUtil.GetMolName(Mol, MolNum)
                    MiscUtil.PrintWarning("Failed to perform constrained minimization for molecule %s with torsion angle set to %s during torsion scan for torsion pattern %s. Abandoning torsion scan..." % (MolName, Angle, TorsionPattern))
                return (Mol, False, None, None, None)
        else:
            # Calculate energy...
            CalcStatus, Energy = GetEnergy(TorsionMol)
            if not CalcStatus:
                if not OptionsInfo["QuietMode"]:
                    MolName = RDKitUtil.GetMolName(Mol, MolNum)
                    MiscUtil.PrintWarning("Failed to retrieve calculated energy for molecule %s with torsion angle set to %s during torsion scan for torsion pattern %s. Abandoning torsion scan..." % (MolName, Angle, TorsionPattern))
                return (Mol, False, None, None, None)
        
        if  OptionsInfo["RemoveHydrogens"]:
            TorsionMol = Chem.RemoveHs(TorsionMol)
        
        TorsionMols.append(TorsionMol)
        TorsionEnergies.append(Energy)
        TorsionAngles.append(Angle)
    
    return (Mol, True, TorsionMols, TorsionEnergies, TorsionAngles)

def SetupTorsionsMolInfo(Mol, MolNum = None):
    """Setup torsions info for a molecule."""

    TorsionsPatternsInfo = OptionsInfo["TorsionsPatternsInfo"]
    
    # Initialize...
    TorsionsMolInfo = {}
    TorsionsMolInfo["IDs"] = []
    TorsionsMolInfo["NumOfMatches"] = 0
    TorsionsMolInfo["Matches"] = {}
    for TorsionID in TorsionsPatternsInfo["IDs"]:
        TorsionsMolInfo["IDs"].append(TorsionID)
        TorsionsMolInfo["Matches"][TorsionID] = None
    
    MolName = RDKitUtil.GetMolName(Mol, MolNum)
    UseChirality = OptionsInfo["UseChirality"]
    
    for TorsionID in TorsionsPatternsInfo["IDs"]:
        # Match torsions..
        TorsionPattern = TorsionsPatternsInfo["Pattern"][TorsionID]
        TorsionPatternMol = TorsionsPatternsInfo["Mol"][TorsionID]
        TorsionsMatches = RDKitUtil.FilterSubstructureMatchesByAtomMapNumbers(Mol, TorsionPatternMol, Mol.GetSubstructMatches(TorsionPatternMol, useChirality = UseChirality))
        
        # Validate tosion matches...
        ValidTorsionsMatches = []
        for Index, TorsionMatch in enumerate(TorsionsMatches):
            if len(TorsionMatch) != 4:
                if not OptionsInfo["QuietMode"]:
                    MiscUtil.PrintWarning("Ignoring invalid torsion match to atom indices, %s, for torsion pattern, %s, in molecule %s: It must match exactly 4 atoms." % (TorsionMatch, TorsionPattern, MolName))
                continue

            if not RDKitUtil.AreAtomIndicesSequentiallyConnected(Mol, TorsionMatch):
                if not OptionsInfo["QuietMode"]:
                    MiscUtil.PrintWarning("Ignoring invalid torsion match to atom indices, %s, for torsion pattern, %s, in molecule %s: Matched atom indices must be sequentially connected." % (TorsionMatch, TorsionPattern, MolName))
                continue

            Bond = Mol.GetBondBetweenAtoms(TorsionMatch[1], TorsionMatch[2])
            if Bond.IsInRing():
                if not OptionsInfo["QuietMode"]:
                    MiscUtil.PrintWarning("Ignoring invalid torsion match to atom indices, %s, for torsion pattern, %s, in molecule %s: Matched atom indices, %s and %s, are not allowed to be in a ring." % (TorsionMatch, TorsionPattern, MolName, TorsionMatch[1], TorsionMatch[2]))
                continue
            
            ValidTorsionsMatches.append(TorsionMatch)
        
        # Track valid matches...
        if len(ValidTorsionsMatches):
            TorsionsMolInfo["NumOfMatches"] += len(ValidTorsionsMatches)
            TorsionsMolInfo["Matches"][TorsionID] = ValidTorsionsMatches
        
    if TorsionsMolInfo["NumOfMatches"] == 0:
        if not OptionsInfo["QuietMode"]:
            MiscUtil.PrintWarning("Failed to match any torsions  in molecule %s" % (MolName))

    return TorsionsMolInfo

def SetupTorsionsPatternsInfo():
    """Setup torsions patterns info."""

    TorsionsPatternsInfo = {}
    TorsionsPatternsInfo["IDs"] = []
    TorsionsPatternsInfo["Pattern"] = {}
    TorsionsPatternsInfo["Mol"] = {}

    TorsionID = 0
    for TorsionPattern in OptionsInfo["TorsionPatternsList"]:
        TorsionID += 1
        
        TorsionMol = Chem.MolFromSmarts(TorsionPattern)
        if TorsionMol is None:
            MiscUtil.PrintError("Failed to create torsion pattern molecule. The torsion SMILES/SMARTS pattern, \"%s\", specified using \"-t, --torsions\" option is not valid." % (TorsionPattern))
        
        TorsionsPatternsInfo["IDs"].append(TorsionID)
        TorsionsPatternsInfo["Pattern"][TorsionID] = TorsionPattern
        TorsionsPatternsInfo["Mol"][TorsionID] = TorsionMol

    OptionsInfo["TorsionsPatternsInfo"] = TorsionsPatternsInfo
    
def MinimizeMolecule(Mol, MolNum = None):
    "Generate and minimize conformers for a molecule to get the lowest energy conformer."

    # Add hydrogens before  minimization. No need to remove them after minimization
    # as they will be used during torsion match and constrained minimization...
    #
    Mol = AddHydrogens(Mol)
    
    ConfIDs = EmbedMolecule(Mol, MolNum)
    if not len(ConfIDs):
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, MolNum)
            MiscUtil.PrintWarning("Minimization couldn't be performed for molecule %s: Embedding failed...\n" % MolName)
        return (Mol, False)

    CalcEnergyMap = {}
    for ConfID in ConfIDs:
        try:
            if OptionsInfo["UseUFF"]:
                Status = AllChem.UFFOptimizeMolecule(Mol, confId = ConfID, maxIters = OptionsInfo["MaxIters"])
            elif OptionsInfo["UseMMFF"]:
                Status = AllChem.MMFFOptimizeMolecule(Mol, confId = ConfID, maxIters = OptionsInfo["MaxIters"], mmffVariant = OptionsInfo["MMFFVariant"])
            else:
                MiscUtil.PrintError("Minimization couldn't be performed: Specified forcefield, %s, is not supported" % OptionsInfo["ForceField"])
        except (ValueError, RuntimeError, Chem.rdchem.KekulizeException) as ErrMsg:
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolNum)
                MiscUtil.PrintWarning("Minimization couldn't be performed for molecule %s:\n%s\n" % (MolName, ErrMsg))
            return (Mol, False)
        
        EnergyStatus, Energy = GetEnergy(Mol, ConfID)
        if not EnergyStatus:
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolNum)
                MiscUtil.PrintWarning("Failed to retrieve calculated energy for conformation number %d of molecule %s. Try again after removing any salts or cleaing up the molecule...\n" % (ConfID, MolName))
            return (Mol, False)
        
        if Status != 0:
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolNum)
                MiscUtil.PrintWarning("Minimization failed to converge for conformation number %d of molecule %s in %d steps. Try using higher value for \"--maxIters\" option...\n" % (ConfID, MolName, OptionsInfo["MaxIters"]))
            
        CalcEnergyMap[ConfID] = Energy
    
    SortedConfIDs = sorted(ConfIDs, key = lambda ConfID: CalcEnergyMap[ConfID])
    MinEnergyConfID = SortedConfIDs[0]
        
    for ConfID in [Conf.GetId() for Conf in Mol.GetConformers()]:
        if ConfID == MinEnergyConfID:
            continue
        Mol.RemoveConformer(ConfID)
    
    # Set ConfID to 0 for MinEnergyConf...
    Mol.GetConformer(MinEnergyConfID).SetId(0)

    return (Mol, True)

def ConstrainAndMinimizeMolecule(Mol, RefMolCore, MolNum = None):
    "Constrain and Minimize molecule."

    # Setup forcefield function to use for constrained minimization...
    ForceFieldFunction = None
    ForceFieldName = None
    if OptionsInfo["UseUFF"]:
        ForceFieldFunction = lambda mol, confId = -1 : AllChem.UFFGetMoleculeForceField(mol, confId = confId)
        ForeceFieldName = "UFF"
    else:
        ForceFieldFunction = lambda mol, confId = -1 : AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol, mmffVariant = OptionsInfo["MMFFVariant"]) , confId = confId)
        ForeceFieldName = "MMFF"

    if ForceFieldFunction is None:
        if not OptionsInfo["QuietMode"]:
            MiscUtil.PrintWarning("Failed to setup forcefield %s for molecule: %s\n" % (ForceFieldName, RDKitUtil.GetMolName(Mol, MolNum)))
        return (None, False, None)
        
    MaxConfs = OptionsInfo["MaxConfsTorsion"]
    EnforceChirality = OptionsInfo["EnforceChirality"]
    UseExpTorsionAnglePrefs = OptionsInfo["UseExpTorsionAnglePrefs"]
    UseBasicKnowledge = OptionsInfo["UseBasicKnowledge"]
    UseTethers = OptionsInfo["UseTethers"]

    CalcEnergyMap = {}
    MolConfsMap = {}
    ConfIDs = [ConfID for ConfID in range(0, MaxConfs)]

    for ConfID in ConfIDs:
        try:
            MolConf = Chem.Mol(Mol)
            AllChem.ConstrainedEmbed(MolConf, RefMolCore, useTethers = UseTethers, coreConfId = -1, randomseed = ConfID, getForceField = ForceFieldFunction, enforceChirality = EnforceChirality, useExpTorsionAnglePrefs = UseExpTorsionAnglePrefs, useBasicKnowledge = UseBasicKnowledge)
        except (ValueError, RuntimeError, Chem.rdchem.KekulizeException)  as ErrMsg:
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolNum)
                MiscUtil.PrintWarning("Constrained embedding coupldn't  be performed for molecule %s:\n%s\n" % (RDKitUtil.GetMolName(Mol, MolNum), ErrMsg))
            return (None, False, None)
        
        EnergyStatus, Energy = GetEnergy(MolConf)
        
        if not EnergyStatus:
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolNum)
                MiscUtil.PrintWarning("Failed to retrieve calculated energy for conformation number %d of molecule %s. Try again after removing any salts or cleaing up the molecule...\n" % (ConfID, MolName))
            return (None, False, None)
        
        CalcEnergyMap[ConfID] = Energy
        MolConfsMap[ConfID] = MolConf

    SortedConfIDs = sorted(ConfIDs, key = lambda ConfID: CalcEnergyMap[ConfID])
    MinEnergyConfID = SortedConfIDs[0]
    
    MinEnergy = CalcEnergyMap[MinEnergyConfID]
    MinEnergyMolConf = MolConfsMap[MinEnergyConfID]
    
    MinEnergyMolConf.ClearProp('EmbedRMS')
    
    return (MinEnergyMolConf, True, MinEnergy)

def GetEnergy(Mol, ConfID = None):
    "Calculate energy."

    Status = True
    Energy = None

    if ConfID is None:
        ConfID = -1
    
    if OptionsInfo["UseUFF"]:
        UFFMoleculeForcefield = AllChem.UFFGetMoleculeForceField(Mol, confId = ConfID)
        if UFFMoleculeForcefield is None:
            Status = False
        else:
            Energy = UFFMoleculeForcefield.CalcEnergy()
    elif OptionsInfo["UseMMFF"]:
        MMFFMoleculeProperties = AllChem.MMFFGetMoleculeProperties(Mol, mmffVariant = OptionsInfo["MMFFVariant"])
        MMFFMoleculeForcefield = AllChem.MMFFGetMoleculeForceField(Mol, MMFFMoleculeProperties, confId = ConfID)
        if MMFFMoleculeForcefield is None:
            Status = False
        else:
            Energy = MMFFMoleculeForcefield.CalcEnergy()
    else:
        MiscUtil.PrintError("Couldn't retrieve conformer energy: Specified forcefield, %s, is not supported" % OptionsInfo["ForceField"])
    
    return (Status, Energy)

def EmbedMolecule(Mol, MolNum = None):
    "Embed conformations"

    ConfIDs = []
    
    MaxConfs = OptionsInfo["MaxConfs"]
    RandomSeed = OptionsInfo["RandomSeed"]
    EnforceChirality = OptionsInfo["EnforceChirality"]
    UseExpTorsionAnglePrefs = OptionsInfo["UseExpTorsionAnglePrefs"]
    UseBasicKnowledge = OptionsInfo["UseBasicKnowledge"]
    
    try:
        ConfIDs = AllChem.EmbedMultipleConfs(Mol, numConfs = MaxConfs, randomSeed = RandomSeed, enforceChirality = EnforceChirality, useExpTorsionAnglePrefs = UseExpTorsionAnglePrefs, useBasicKnowledge = UseBasicKnowledge)
    except ValueError as ErrMsg:
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, MolNum)
            MiscUtil.PrintWarning("Embedding failed  for molecule %s:\n%s\n" % (MolName, ErrMsg))
        ConfIDs = []
    
    return ConfIDs

def GenerateStartingTorsionScanStructureOutfile(Mol, MolNum):
    """Write out the structure of molecule used for starting tosion scan. """
    
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Options["--outfile"])
    MolName = GetOutputFileMolName(Mol, MolNum)
    
    Outfile  = "%s_%s.%s" % (FileName, MolName, FileExt)
    
    # Set up a molecule writer...
    Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintWarning("Failed to setup a writer for output fie %s " % Outfile)
        return
    
    Writer.write(Mol)
    
    if Writer is not None:
        Writer.close()

def GenerateOutputFiles(Mol, MolNum, TorsionID, TorsionMatchNum, TorsionMols, TorsionEnergies, TorsionAngles):
    """Generate output files. """
    
    StructureOutfile, EnergyTextOutfile, PlotOutfile = SetupOutputFileNames(Mol, MolNum, TorsionID, TorsionMatchNum)
    
    GenerateScannedTorsionsStructureOutfile(StructureOutfile, Mol, MolNum, TorsionID, TorsionMatchNum, TorsionMols, TorsionEnergies, TorsionAngles)
    GenerateEnergyTextOutfile(EnergyTextOutfile, Mol, MolNum, TorsionID, TorsionMatchNum, TorsionMols, TorsionEnergies, TorsionAngles)
    GeneratePlotOutfile(PlotOutfile, Mol, MolNum, TorsionID, TorsionMatchNum, TorsionMols, TorsionEnergies, TorsionAngles)

def GenerateScannedTorsionsStructureOutfile(Outfile, Mol, MolNum, TorsionID, TorsionMatchNum, TorsionMols, TorsionEnergies, TorsionAngles):
    """Write out structures generated after torsion scan along with associated data."""

    # Set up a molecule writer...
    Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintWarning("Failed to setup a writer for output fie %s " % Outfile)
        return
    
    MolName = RDKitUtil.GetMolName(Mol, MolNum)
    
    for Index, TorsionMol in enumerate(TorsionMols):
        TorsionAngle = "%s" % TorsionAngles[Index]
        TorsionMol.SetProp("Torsion_Angle", TorsionAngle)
        
        TorsionEnergy = "%.2f" % TorsionEnergies[Index]
        TorsionMol.SetProp(OptionsInfo["EnergyLabel"], TorsionEnergy)

        TorsionMolName = "%s_Deg%s" % (MolName, TorsionAngle)
        TorsionMol.SetProp("_Name", TorsionMolName)
        
        Writer.write(TorsionMol)
        
    if Writer is not None:
        Writer.close()
    
def GenerateEnergyTextOutfile(Outfile, Mol, MolNum, TorsionID, TorsionMatchNum, TorsionMols, TorsionEnergies, TorsionAngles):
    """Write out torsion angles and energies."""

    # Setup a writer...
    Writer = open(Outfile, "w")
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
    
    # Write headers...
    Writer.write("TorsionAngle,%s\n" % (OptionsInfo["EnergyLabel"]))

    for Index, TorsionAngle in enumerate(TorsionAngles):
        Writer.write("%d,%.2f\n" % (TorsionAngle, TorsionEnergies[Index]))

    if Writer is not None:
        Writer.close()
    
def GeneratePlotOutfile(Outfile, Mol, MolNum, TorsionID, TorsionMatchNum, TorsionMols, TorsionEnergies, TorsionAngles):
    """Generate a plot corresponding to torsion angles and energies."""

    OutPlotParams = OptionsInfo["OutPlotParams"]

    # Initialize seaborn and matplotlib paramaters...
    if not OptionsInfo["OutPlotInitialized"]:
        OptionsInfo["OutPlotInitialized"] = True
        RCParams = {"figure.figsize":(OutPlotParams["Width"], OutPlotParams["Height"]),
                    "axes.titleweight": OutPlotParams["TitleWeight"],
                    "axes.labelweight": OutPlotParams["LabelWeight"]}
        sns.set(context = OutPlotParams["Context"], style = OutPlotParams["Style"], palette = OutPlotParams["Palette"], font = OutPlotParams["Font"], font_scale = OutPlotParams["FontScale"], rc = RCParams)

    # Create a new figure...
    plt.figure()

    # Draw plot...
    PlotType = OutPlotParams["Type"]
    if re.match("linepoint", PlotType, re.I):
        Axis = sns.lineplot(x = TorsionAngles, y = TorsionEnergies, marker = "o",  legend = False)
    elif re.match("scatter", PlotType, re.I):
        Axis = sns.scatterplot(x = TorsionAngles, y = TorsionEnergies, legend = False)
    elif re.match("line", PlotType, re.I):
        Axis = sns.lineplot(x = TorsionAngles, y = TorsionEnergies, legend = False)
    else:
        MiscUtil.PrintError("The value, %s, specified for \"type\" using option \"--outPlotParams\" is not supported. Valid plot types: linepoint, scatter or line" % (PlotType))

    # Setup title and labels...
    Title = OutPlotParams["Title"]
    if OptionsInfo["OutPlotTitleTorsionSpec"]:
        TorsionPattern = OptionsInfo["TorsionsPatternsInfo"]["Pattern"][TorsionID]
        Title = "%s: %s" % (OutPlotParams["Title"], TorsionPattern)

    # Set labels and title...
    Axis.set(xlabel = OutPlotParams["XLabel"], ylabel = OutPlotParams["YLabel"], title = Title)
    
    # Save figure...
    plt.savefig(Outfile)

    # Close the plot...
    plt.close()

def SetupOutputFileNames(Mol, MolNum, TorsionID, TorsionMatchNum):
    """Setup names of output files. """
    
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Options["--outfile"])
    MolName = GetOutputFileMolName(Mol, MolNum)
    
    OutfileRoot  = "%s_%s_Torsion%s_Match%s" % (FileName, MolName, TorsionID, TorsionMatchNum)
    
    StructureOutfile = "%s.%s" % (OutfileRoot, FileExt)
    EnergyTextOutfile = "%s_Energies.csv" % (OutfileRoot)
    PlotExt = OptionsInfo["OutPlotParams"]["OutExt"]
    PlotOutfile = "%s_Plot.%s" % (OutfileRoot, PlotExt)

    return (StructureOutfile, EnergyTextOutfile, PlotOutfile)

def GetOutputFileMolName(Mol, MolNum):
    """Get output file prefix. """
    
    MolName = "Mol%s" % MolNum
    if OptionsInfo["OutfileMolName"]:
        MolName = re.sub("[^a-zA-Z0-9]", "_", RDKitUtil.GetMolName(Mol, MolNum), re.I)

    return MolName

def AddHydrogens(Mol, AddCoords = True):
    """Check and add hydrogens. """
    
    if  not OptionsInfo["AddHydrogens"]:
        return Mol

    if OptionsInfo["HydrogensAdded"]:
        return Mol

    OptionsInfo["HydrogensAdded"] = True
    return Chem.AddHs(Mol, addCoords = AddCoords)

def ProcessOptions():
    """Process and validate command line arguments and options."""
    
    MiscUtil.PrintInfo("Processing options...")

    # Validate options...
    ValidateOptions()
    
    OptionsInfo["ModeMols"] = Options["--modeMols"]
    OptionsInfo["FirstMolMode"] = True if re.match("^First$", Options["--modeMols"], re.I) else False
    
    OptionsInfo["ModeTorsions"] = Options["--modeTorsions"]
    OptionsInfo["FirstTorsionMode"] = True if re.match("^First$", Options["--modeTorsions"], re.I) else False
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    OptionsInfo["Infile3D"] = True if re.match("^yes$", Options["--infile3D"], re.I) else False
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"])

    OptionsInfo["OutfileMolName"] = True if re.match("^yes$", Options["--outfileMolName"], re.I) else False
    
    OptionsInfo["OutPlotTitleTorsionSpec"] = True if re.match("^yes$", Options["--outPlotTitleTorsionSpec"], re.I) else False
    
    # The default width and height, 10.0 and 7.5, map to aspect raito of 16/9 (1.778)...
    DefaultValues = {'Type': 'linepoint', 'Width': 10.0, 'Height': 5.6, 'Title': 'Torsion Scan', 'XLabel': 'Torsion Angle (degrees)', 'YLabel': 'Energy (kcal/mol)'}
    OptionsInfo["OutPlotParams"] = MiscUtil.ProcessOptionSeabornPlotParameters("--outPlotParams", Options["--outPlotParams"], DefaultValues)
    if not re.match("^(linepoint|scatter|Line)$", OptionsInfo["OutPlotParams"]["Type"], re.I):
        MiscUtil.PrintError("The value, %s, specified for \"type\" using option \"--outPlotParams\" is not supported. Valid plot types: linepoint, scatter or line" % (OptionsInfo["OutPlotParams"]["Type"]))
    
    OptionsInfo["OutPlotInitialized"] = False
    
    # Procsss and validate specified SMILES/SMARTS torsion patterns...
    TorsionPatterns = Options["--torsions"]
    TorsionPatternsList = []
    for TorsionPattern in TorsionPatterns.split(","):
        TorsionPattern = TorsionPattern.strip()
        if not len(TorsionPattern):
            MiscUtil.PrintError("Empty value specified for SMILES/SMARTS pattern in  \"-t, --torsions\" option: %s" % TorsionPatterns)
        
        TorsionMol = Chem.MolFromSmarts(TorsionPattern)
        if TorsionMol is None:
            MiscUtil.PrintError("Failed to create torsion pattern molecule. The torsion SMILES/SMARTS pattern, \"%s\", specified using \"-t, --torsions\" option, \"%s\",  is not valid." % (TorsionPattern, TorsionPatterns))
        TorsionPatternsList.append(TorsionPattern)
    
    OptionsInfo["TorsionPatterns"] = TorsionPatterns
    OptionsInfo["TorsionPatternsList"] = TorsionPatternsList
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    OptionsInfo["AddHydrogens"] = True if re.match("^yes$", Options["--addHydrogens"], re.I) else False
    OptionsInfo["HydrogensAdded"] = False
    
    if re.match("^ETDG$", Options["--conformerGenerator"], re.I):
        ConformerGenerator = "ETDG"
        UseExpTorsionAnglePrefs = True
        UseBasicKnowledge = False
    elif re.match("^KDG$", Options["--conformerGenerator"], re.I):
        ConformerGenerator = "KDG"
        UseExpTorsionAnglePrefs = False
        UseBasicKnowledge = True
    elif re.match("^ETKDG$", Options["--conformerGenerator"], re.I):
        ConformerGenerator = "ETKDG"
        UseExpTorsionAnglePrefs = True
        UseBasicKnowledge = True
    elif re.match("^SDG$", Options["--conformerGenerator"], re.I):
        ConformerGenerator = "SDG"
        UseExpTorsionAnglePrefs = False
        UseBasicKnowledge = False
    else:
        MiscUtil.PrintError("The value, %s, specified for option \"-c, --conformerGenerator\" is not supported." % (Options["--conformerGenerator"]))
    
    OptionsInfo["ConformerGenerator"] = ConformerGenerator
    OptionsInfo["UseExpTorsionAnglePrefs"] = UseExpTorsionAnglePrefs
    OptionsInfo["UseBasicKnowledge"] = UseBasicKnowledge

    if re.match("^UFF$", Options["--forceField"], re.I):
        ForceField = "UFF"
        UseUFF = True
        UseMMFF = False
    elif re.match("^MMFF$", Options["--forceField"], re.I):
        ForceField = "MMFF"
        UseUFF = False
        UseMMFF = True
    else:
        MiscUtil.PrintError("The value, %s, specified for \"--forceField\" is not supported." % (Options["--forceField"],))
    
    MMFFVariant = "MMFF94" if re.match("^MMFF94$", Options["--forceFieldMMFFVariant"], re.I) else "MMFF94s"
    
    OptionsInfo["ForceField"] = ForceField
    OptionsInfo["MMFFVariant"] = MMFFVariant
    OptionsInfo["UseMMFF"] = UseMMFF
    OptionsInfo["UseUFF"] = UseUFF
    
    if UseMMFF:
        OptionsInfo["EnergyLabel"] = "%s_Energy" % MMFFVariant
    else:
        OptionsInfo["EnergyLabel"] = "%s_Energy" % ForceField
    
    OptionsInfo["EnforceChirality"] = True if re.match("^yes$", Options["--enforceChirality"], re.I) else False
    
    OptionsInfo["MaxConfs"] = int(Options["--maxConfs"])
    OptionsInfo["MaxConfsTorsion"] = int(Options["--maxConfsTorsion"])
    OptionsInfo["MaxIters"] = int(Options["--maxIters"])
    
    OptionsInfo["MPMode"] = True if re.match("^yes$", Options["--mp"], re.I) else False
    OptionsInfo["MPParams"] = MiscUtil.ProcessOptionMultiprocessingParameters("--mpParams", Options["--mpParams"])
    
    OptionsInfo["QuietMode"] = True if re.match("^yes$", Options["--quiet"], re.I) else False
    
    RandomSeed = -1
    if not re.match("^auto$", Options["--randomSeed"], re.I):
        RandomSeed = int(Options["--randomSeed"])
    OptionsInfo["RandomSeed"] = RandomSeed
    
    OptionsInfo["RemoveHydrogens"] = True if re.match("^yes$", Options["--removeHydrogens"], re.I) else False
    
    OptionsInfo["TorsionMaxMatches"] = int(Options["--torsionMaxMatches"])
    OptionsInfo["TorsionMinimize"] = True if re.match("^yes$", Options["--torsionMinimize"], re.I) else False

    TorsionRange = Options["--torsionRange"]
    TorsionRangeWords = TorsionRange.split(",")
    
    TorsionStart = int(TorsionRangeWords[0])
    TorsionStop = int(TorsionRangeWords[1])
    TorsionStep = int(TorsionRangeWords[2])
    
    if TorsionStart >= TorsionStop:
        MiscUtil.PrintError("The start value, %d, specified for option \"--torsionRange\" in string \"%s\" must be less than stop value, %s." % (TorsionStart, Options["--torsionRange"], TorsionStop))
    if TorsionStep == 0:
        MiscUtil.PrintError("The step value, %d, specified for option \"--torsonRange\" in string \"%s\" must be > 0." % (TorsionStep, Options["--torsionRange"]))
    if TorsionStep >= (TorsionStop - TorsionStart):
        MiscUtil.PrintError("The step value, %d, specified for option \"--torsonRange\" in string \"%s\" must be less than, %s." % (TorsionStep, Options["--torsionRange"], (TorsionStop - TorsionStart)))
    
    if TorsionStart < 0:
        if TorsionStart < -180:
            MiscUtil.PrintError("The start value, %d, specified for option \"--torsionRange\" in string \"%s\" must be  >= -180 to use scan range from -180 to 180." % (TorsionStart, Options["--torsionRange"]))
        if TorsionStop > 180:
            MiscUtil.PrintError("The stop value, %d, specified for option \"--torsionRange\" in string \"%s\" must be <= 180 to use scan range from -180 to 180." % (TorsionStop, Options["--torsionRange"]))
    else:
        if TorsionStop > 360:
            MiscUtil.PrintError("The stop value, %d, specified for option \"--torsionRange\" in string \"%s\" must be  <= 360 to use scan range from 0 to 360." % (TorsionStop, Options["--torsionRange"]))
    
    OptionsInfo["TorsionRange"] = TorsionRange
    OptionsInfo["TorsionStart"] = TorsionStart
    OptionsInfo["TorsionStop"] = TorsionStop
    OptionsInfo["TorsionStep"] = TorsionStep
    
    OptionsInfo["UseTethers"] = True if re.match("^yes$", Options["--useTethers"], re.I) else False
    OptionsInfo["UseChirality"] = True if re.match("^yes$", Options["--useTethers"], re.I) else False

def RetrieveOptions():
    """Retrieve command line arguments and options"""
    
    # Get options...
    global Options
    Options = docopt(_docoptUsage_)
    
    # Set current working directory to the specified directory...
    WorkingDir = Options["--workingdir"]
    if WorkingDir:
        os.chdir(WorkingDir)
    
    # Handle examples option...
    if "--examples" in Options and Options["--examples"]:
        MiscUtil.PrintInfo(MiscUtil.GetExamplesTextFromDocOptText(_docoptUsage_))
        sys.exit(0)

def ValidateOptions():
    """Validate option values"""

    MiscUtil.ValidateOptionTextValue("-a, --addHydrogens", Options["--addHydrogens"], "yes no")
    MiscUtil.ValidateOptionTextValue("-c, --conformerGenerator", Options["--conformerGenerator"], "SDG ETDG KDG ETKDG")
    
    MiscUtil.ValidateOptionTextValue("-f, --forceField", Options["--forceField"], "UFF MMFF")
    MiscUtil.ValidateOptionTextValue(" --forceFieldMMFFVariant", Options["--forceFieldMMFFVariant"], "MMFF94 MMFF94s")
    
    MiscUtil.ValidateOptionTextValue("--enforceChirality ", Options["--enforceChirality"], "yes no")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi txt csv tsv")
    MiscUtil.ValidateOptionTextValue("--infile3D", Options["--infile3D"], "yes no")

    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    if not Options["--overwrite"]:
        FileDir, FileName, FileExt = MiscUtil.ParseFileName(Options["--outfile"])
        FileNames = glob.glob("%s_*" % FileName)
        if len(FileNames):
            MiscUtil.PrintError("The outfile names, %s_*, generated from file specified, %s, for option \"-o, --outfile\" already exist. Use option \"--overwrite\" or \"--ov\"  and try again.\n" % (FileName, Options["--outfile"]))
    
    MiscUtil.ValidateOptionTextValue("--outPlotTitleTorsionSpec", Options["--outPlotTitleTorsionSpec"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--outfileMolName ", Options["--outfileMolName"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--modeMols", Options["--modeMols"], "First All")
    MiscUtil.ValidateOptionTextValue("--modeTorsions", Options["--modeTorsions"], "First All")

    MiscUtil.ValidateOptionIntegerValue("--maxConfs", Options["--maxConfs"], {">": 0})
    MiscUtil.ValidateOptionIntegerValue("--maxConfsTorsion", Options["--maxConfsTorsion"], {">": 0})
    MiscUtil.ValidateOptionIntegerValue("--maxIters", Options["--maxIters"], {">": 0})
    
    MiscUtil.ValidateOptionTextValue("--mp", Options["--mp"], "yes no")
    MiscUtil.ValidateOptionTextValue("-q, --quiet", Options["--quiet"], "yes no")
    
    if not re.match("^auto$", Options["--randomSeed"], re.I):
        MiscUtil.ValidateOptionIntegerValue("--randomSeed", Options["--randomSeed"], {})
    
    MiscUtil.ValidateOptionTextValue("-r, --removeHydrogens", Options["--removeHydrogens"], "yes no")
    
    MiscUtil.ValidateOptionIntegerValue("--torsionMaxMatches", Options["--torsionMaxMatches"], {">": 0})
    MiscUtil.ValidateOptionTextValue("--torsionMinimize", Options["--torsionMinimize"], "yes no")
    MiscUtil.ValidateOptionNumberValues("--torsionRange", Options["--torsionRange"], 3, ",", "integer", {})
    
    MiscUtil.ValidateOptionTextValue("--useChirality", Options["--useChirality"], "yes no")
    MiscUtil.ValidateOptionTextValue("--useTethers", Options["--useTethers"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitPerformTorsionScan.py - Perform torsion scan

Usage:
    RDKitPerformTorsionScan.py [--addHydrogens <yes or no>] [--conformerGenerator <SDG, ETDG, KDG, ETKDG>]
                               [--forceField <UFF, or MMFF>] [--forceFieldMMFFVariant <MMFF94 or MMFF94s>]
                               [--enforceChirality <yes or no>] [--infile3D <yes or no>] [--infileParams <Name,Value,...>]
                               [--modeMols  <First or All>] [--modeTorsions  <First or All> ] [--maxConfs <number>]
                               [--maxConfsTorsion <number>] [--maxIters <number>] [--mp <yes or no>] [--mpParams <Name.Value,...>]
                               [--outfileMolName  <yes or no>] [--outfileParams <Name,Value,...>] [--outPlotParams <Name,Value,...>]
                               [--outPlotTitleTorsionSpec <yes or no>] [--overwrite]  [--quiet <yes or no>] [--removeHydrogens <yes or no>]
                               [--randomSeed <number>] [--torsionMaxMatches <number>] [--torsionMinimize <yes or no>]
                               [--torsionRange <Start,Stop,Step>] [--useChirality <yes or no>] [--useTethers  <yes or no>]
                               [-w <dir>] -t <torsions> -i <infile>  -o <outfile> 
    RDKitPerformTorsionScan.py -h | --help | -e | --examples

Description:
    Perform torsion scan for molecules around torsion angles specified using
    SMILES/SMARTS patterns. A molecule is optionally minimized before performing
    a torsion scan. A set of initial 3D structures are generated for a molecule
    by scanning the torsion angle across the specified range and updating the 3D
    coordinates of the molecule. A conformation ensemble is optionally generated
    for each 3D structure representing a specific torsion angle. The conformation
    with the lowest energy is selected to represent the torsion angle. An option
    is available to skip the generation of the conformation ensemble and simply
    calculate the energy for the initial 3D structure for a specific torsion angle

    The torsions are specified using SMILES or SMARTS patterns. A substructure match
    is performed to select torsion atoms in a molecule. The SMILES pattern match must
    correspond to four torsion atoms. The SMARTS patterns containing atom indices may
    match  more than four atoms. The atoms indices, however, must match exactly four
    torsion atoms. For example: [s:1][c:2]([aX2,cH1])!@[CX3:3](O)=[O:4] for thiophene
    esters and carboxylates as specified in Torsion Library (TorLib) [Ref 146].

    A set of four output files is generated for each torsion match in each
    molecule. The names of the output files are generated using the root of
    the specified output file. They may either contain sequential molecule
    numbers or molecule names as shown below:
        
        <OutfileRoot>_Mol<Num>.sdf
        <OutfileRoot>_Mol<Num>_Torsion<Num>_Match<Num>.sdf
        <OutfileRoot>_Mol<Num>_Torsion<Num>_Match<Num>_Energies.csv
        <OutfileRoot>_Mol<Num>_Torsion<Num>_Match<Num>_Plot.<ImgExt>
        
        or
        
        <OutfileRoot>_<MolName>.sdf
        <OutfileRoot>_<MolName>_Torsion<Num>_Match<Num>.sdf
        <OutfileRoot>_<MolName>_Torsion<Num>_Match<Num>_Energies.csv
        <OutfileRoot>_<MolName>_Torsion<Num>_Match<Num>_Plot.<ImgExt>
        
    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), .csv, .tsv .txt)

    The supported output file formats are: SD (.sdf, .sd)

Options:
    -a, --addHydrogens <yes or no>  [default: yes]
        Add hydrogens before minimization.
    -c, --conformerGenerator <SDG, ETDG, KDG, ETKDG>  [default: ETKDG]
        Conformation generation methodology for generating initial 3D structure
        of a molecule and conformation ensemble representing a specific torsion
        angle. No conformation ensemble is generated for 'No' value of
        '--torsionMinimize' option.
        
        Possible values: Standard Distance Geometry, (SDG), Experimental Torsion-angle
        preference with Distance Geometry (ETDG), basic Knowledge-terms with Distance
        Geometry (KDG),  and Experimental Torsion-angle preference along with basic
        Knowledge-terms with Distance Geometry (ETKDG) [Ref 129] .
    -f, --forceField <UFF, MMFF>  [default: MMFF]
        Forcefield method to use for  energy minimization of initial 3D structure
        of a molecule and conformation ensemble representing a specific torsion.
        No conformation ensemble is generated during for 'No' value of '--torsionMinimze'
        option and constrained energy minimization is not performed. Possible values:
        Universal Force Field (UFF) [ Ref 81 ] or Merck Molecular Mechanics Force
        Field [ Ref 83-87 ] .
    --forceFieldMMFFVariant <MMFF94 or MMFF94s>  [default: MMFF94]
        Variant of MMFF forcefield to use for energy minimization.
    --enforceChirality <yes or no>  [default: Yes]
        Enforce chirality for defined chiral centers during generation of conformers.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infile <infile>
        Input file name.
    --infile3D <yes or no>  [default: no]
        Skip generation and minimization of initial 3D structures for molecules in
        input file containing 3D coordinates.
    --infileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for reading
        molecules from files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD, MOL: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                smilesTitleLine,auto,sanitize,yes
            
        Possible values for smilesDelimiter: space, comma or tab.
    --modeMols <First or All>  [default: First]
        Perform torsion scan for the first molecule or all molecules in input
        file.
    --modeTorsions <First or All>  [default: First]
        Perform torsion scan for the first or all specified torsion pattern in
        molecules up to a maximum number of matches for each torsion
        specification as indicated by '--torsionMaxMatches' option. 
    --maxConfs <number>  [default: 250]
        Maximum number of conformations to generate for initial 3D structure of a
        molecule. The lowest energy conformation is written to the output file.
    --maxConfsTorsion <number>  [default: 50]
        Maximum number of conformations to generate for conformation ensemble
        representing a specific torsion. A constrained minimization is performed
        using the coordinates of the specified torsion and the lowest energy
        conformation is written to the output file.
    --maxIters <number>  [default: 500]
        Maximum number of iterations to perform for a molecule during minimization
        to generation initial 3D structures. This option is ignored during 'yes' value
        of  '--infile3D' option.
    --mp <yes or no>  [default: no]
        Use multiprocessing.
         
        By default, input data is retrieved in a lazy manner via mp.Pool.imap()
        function employing lazy RDKit data iterable. This allows processing of
        arbitrary large data sets without any additional requirements memory.
        
        All input data may be optionally loaded into memory by mp.Pool.map()
        before starting worker processes in a process pool by setting the value
        of 'inputDataMode' to 'InMemory' in '--mpParams' option.
        
        A word to the wise: The default 'chunkSize' value of 1 during 'Lazy' input
        data mode may adversely impact the performance. The '--mpParams' section
        provides additional information to tune the value of 'chunkSize'.
    --mpParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for to
        configure multiprocessing.
        
        The supported parameter names along with their default and possible
        values are shown below:
        
            chunkSize, auto
            inputDataMode, Lazy   [ Possible values: InMemory or Lazy ]
            numProcesses, auto   [ Default: mp.cpu_count() ]
        
        These parameters are used by the following functions to configure and
        control the behavior of multiprocessing: mp.Pool(), mp.Pool.map(), and
        mp.Pool.imap().
        
        The chunkSize determines chunks of input data passed to each worker
        process in a process pool by mp.Pool.map() and mp.Pool.imap() functions.
        The default value of chunkSize is dependent on the value of 'inputDataMode'.
        
        The mp.Pool.map() function, invoked during 'InMemory' input data mode,
        automatically converts RDKit data iterable into a list, loads all data into
        memory, and calculates the default chunkSize using the following method
        as shown in its code:
        
            chunkSize, extra = divmod(len(dataIterable), len(numProcesses) * 4)
            if extra: chunkSize += 1
        
        For example, the default chunkSize will be 7 for a pool of 4 worker processes
        and 100 data items.
        
        The mp.Pool.imap() function, invoked during 'Lazy' input data mode, employs
        'lazy' RDKit data iterable to retrieve data as needed, without loading all the
        data into memory. Consequently, the size of input data is not known a priori.
        It's not possible to estimate an optimal value for the chunkSize. The default 
        chunkSize is set to 1.
        
        The default value for the chunkSize during 'Lazy' data mode may adversely
        impact the performance due to the overhead associated with exchanging
        small chunks of data. It is generally a good idea to explicitly set chunkSize to
        a larger value during 'Lazy' input data mode, based on the size of your input
        data and number of processes in the process pool.
        
        The mp.Pool.map() function waits for all worker processes to process all
        the data and return the results. The mp.Pool.imap() function, however,
        returns the the results obtained from worker processes as soon as the
        results become available for specified chunks of data.
        
        The order of data in the results returned by both mp.Pool.map() and 
        mp.Pool.imap() functions always corresponds to the input data.
    -o, --outfile <outfile>
        Output file name. The output file root is used for generating the names
        of the output files corresponding to structures, energies, and plots during
        the torsion scan.
    --outfileMolName <yes or no>  [default: no]
        Append molecule name to output file root during the generation of the names
        for output files. The default is to use <MolNum>. The non alphabetical
        characters in molecule names are replaced by underscores.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: kekulize,no
            
    --outPlotParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for generating
        plots using Seaborn module. The supported parameter names along with their
        default values are shown below:
            
            type,linepoint,outExt,svg,width,10,height,5.6,
            title,auto,xlabel,auto,ylabel,auto,titleWeight,bold,labelWeight,bold
            style,darkgrid,palette,deep,font,sans-serif,fontScale,1,
            context,notebook
            
        Possible values:
            
            type: linepoint, scatter, or line. Both points and lines are drawn
                for linepoint plot type.
            outExt: Any valid format supported by Python module Matplotlib.
                For example: PDF (.pdf), PNG (.png), PS (.ps), SVG (.svg)
            titleWeight, labelWeight: Font weight for title and axes labels.
                Any valid value.
            style: darkgrid, whitegrid, dark, white, ticks
            palette: deep, muted, pastel, dark, bright, colorblind
            font: Any valid font name
        
     --outPlotTitleTorsionSpec <yes or no>  [default: yes]
        Append torsion specification to the title of the torsion plot.
    --overwrite
        Overwrite existing files.
    -q, --quiet <yes or no>  [default: no]
        Use quiet mode. The warning and information messages will not be printed.
    --randomSeed <number>  [default: auto]
        Seed for the random number generator for generating initial 3D coordinates.
        Default is to use a random seed.
    --removeHydrogens <yes or no>  [default: Yes]
        Remove hydrogens after minimization.
    -t, --torsions <SMILES/SMARTS,...,...>
        SMILES/SMARTS patterns corresponding to torsion specifications. It's a 
        comma delimited list of valid SMILES/SMART patterns.
        
        A substructure match is performed to select torsion atoms in a molecule.
        The SMILES pattern match must correspond to four torsion atoms. The
        SMARTS patterns contain atom indices may match  more than four atoms.
        The atoms indices, however, must match exactly four torsion atoms. For example:
        [s:1][c:2]([aX2,cH1])!@[CX3:3](O)=[O:4] for thiophene esters and carboxylates
        as specified in Torsion Library (TorLib) [Ref 146].
    --torsionMaxMatches <number>  [default: 5]
        Maximum number of torsions to match for each torsion specification in a
        molecule.
    --torsionMinimize <yes or no>  [default: no]
        Perform constrained energy minimization on a conformation ensemble
        for  a specific torsion angle and select the lowest energy conformation
        representing the torsion angle.
    --torsionRange <Start,Stop,Step>  [default: 0,360,5]
        Start, stop, and step size angles in degrees for a torsion scan. In addition,
        you may specify values using start and stop angles from -180 to 180.
    --useChirality <yes or no>  [default: no]
        Use chirrality during substructure matches for identification of torsions.
     --useTethers <yes or no>  [default: yes]
        Use tethers to optimize the final conformation by applying a series of extra forces
        to align matching atoms to the positions of the core atoms. Otherwise, use simple
        distance constraints during the optimization.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To perform a torsion scan on first molecule in a SMILES file using a minimum
    energy structure of the molecule selected from an ensemble of conformations,
    skipping generation of conformation ensembles for specific torsion angles and
    constrained energy minimization of the ensemble, generate output files
    corresponding to structure, energy and torsion plot, type:
    
        % RDKitPerformTorsionScan.py  -t "O=CNC" -i SampleSeriesD3R.smi 
          -o SampleOut.sdf

    To run the previous example on all molecules in a SD file, type:
    
        % RDKitPerformTorsionScan.py  -t "O=CNC" --modeMols All
          -i SampleSeriesD3R.sdf -o SampleOut.sdf

    To perform a torsion scan on first molecule in a SMILES file using a minimum
    energy structure of the molecule selected from an ensemble of conformations,
    generation of conformation ensembles for specific torsion angles and constrained
    energy minimization of the ensemble, generate output files corresponding to
    structure, energy and torsion plot, type:
    
        % RDKitPerformTorsionScan.py  -t "O=CNC" --torsionMinimize Yes
           -i SampleSeriesD3R.smi -o SampleOut.sdf

    To run the previous example on all molecules in a SD file, type:
    
        % RDKitPerformTorsionScan.py  -t "O=CNC" --modeMols All
           --torsionMinimize Yes -i SampleSeriesD3R.sdf -o SampleOut.sdf

    To run the previous example in multiprocessing mode on all available CPUs
    without loading all data into memory and write out a SD file, type:

        % RDKitPerformTorsionScan.py  -t "O=CNC" -i SampleSeriesD3R.smi 
          -o SampleOut.sdf --modeMols All --torsionMinimize Yes --mp yes
    
    To run the previous example in multiprocessing mode on all available CPUs
    by loading all data into memory and write out a SD file, type:

        % RDKitPerformTorsionScan.py  -t "O=CNC" -i SampleSeriesD3R.smi 
          -o SampleOut.sdf --modeMols All --torsionMinimize Yes --mp yes
          --mpParams "inputDataMode,InMemory"
    
    To run the previous example in multiprocessing mode on specific number of
    CPUs and chunk size without loading all data into memory and write out a SD file,
    type:

        % RDKitPerformTorsionScan.py  -t "O=CNC" -i SampleSeriesD3R.smi 
          -o SampleOut.sdf --modeMols All --torsionMinimize Yes --mp yes
          --mpParams "inputDataMode,Lazy,numProcesses,4,chunkSize,8"

    To perform a torsion scan on first molecule in a SD file containing 3D coordinates,
    skipping generation of conformation ensembles for specific torsion angles and
    constrained energy minimization of the ensemble, generate output files
    corresponding to structure, energy and torsion plot, type:
    
        % RDKitPerformTorsionScan.py  -t "O=CNC"  --infile3D yes
          -i SampleSeriesD3R3D.sdf -o SampleOut.sdf

    To perform a torsion scan using multiple torsion specifications on all molecules in
    a SD file containing 3D coordinates, generation of conformation ensembles for specific
    torsion angles and constrained energy minimization of the ensemble, generate output files
    corresponding to structure, energy and torsion plot, type:
    
        % RDKitPerformTorsionScan.py  -t "O=CNC,[O:1]=[C:2](c)[N:3][C:4]"
          --infile3D yes --modeMols All  --modeTorsions All
          --torsionMinimize Yes -i SampleSeriesD3R3D.sdf -o SampleOut.sdf

    To run the previous example using a specific torsion scan range, type:

        % RDKitPerformTorsionScan.py  -t "O=CNC,[O:1]=[C:2](c)[N:3][C:4]"
          --infile3D yes --modeMols All --modeTorsions All --torsionMinimize
          Yes --torsionRange 0,360,10 -i SampleSeriesD3R.smi -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitCalculateRMSD.py, RDKitCalculateMolecularDescriptors.py, RDKitCompareMoleculeShapes.py,
    RDKitConvertFileFormat.py, RDKitPerformConstrainedMinimization.py

Copyright:
    Copyright (C) 2020 Manish Sud. All rights reserved.

    The functionality available in this script is implemented using RDKit, an
    open source toolkit for cheminformatics developed by Greg Landrum.

    This file is part of MayaChemTools.

    MayaChemTools is free software; you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option) any
    later version.

"""

if __name__ == "__main__":
    main()
