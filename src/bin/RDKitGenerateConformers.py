#!/usr/bin/env python
#
# File: RDKitGenerateConformers.py
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
import multiprocessing as mp

# RDKit imports...
try:
    from rdkit import rdBase
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Descriptors
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
    GenerateConformers()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def GenerateConformers():
    """Generate conformers."""
    
    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % OptionsInfo["Infile"])
    Mols  = RDKitUtil.ReadMolecules(OptionsInfo["Infile"], **OptionsInfo["InfileParams"])
    
    # Set up a molecule writer...
    Writer = RDKitUtil.MoleculesWriter(OptionsInfo["Outfile"], **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % OptionsInfo["Outfile"])
    MiscUtil.PrintInfo("Generating file %s..." % OptionsInfo["Outfile"])

    MolCount, ValidMolCount, ConfGenFailedCount = ProcessMolecules(Mols, Writer)

    if Writer is not None:
        Writer.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of molecules failed during conformation generation or minimization: %d" % ConfGenFailedCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount + ConfGenFailedCount))

def ProcessMolecules(Mols, Writer):
    """Process molecules to generate conformers."""
    
    if OptionsInfo["MPMode"]:
        return ProcessMoleculesUsingMultipleProcesses(Mols, Writer)
    else:
        return ProcessMoleculesUsingSingleProcess(Mols, Writer)

def ProcessMoleculesUsingSingleProcess(Mols, Writer):
    """Process molecules to generate conformers using a single process."""
    
    if OptionsInfo["SkipForceFieldMinimization"]:
        MiscUtil.PrintInfo("\nGenerating conformers without performing energy minimization...")
    else:
        MiscUtil.PrintInfo("\nGenerating conformers and performing energy minimization...")
    
    (MolCount, ValidMolCount, ConfGenFailedCount) = [0] * 3
    for Mol in Mols:
        MolCount += 1
        
        if Mol is None:
            continue
        
        if RDKitUtil.IsMolEmpty(Mol):
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolCount)
                MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
            continue
        ValidMolCount += 1

        ConformerMol, CalcStatus, ConfIDs, ConfEnergies = GenerateMolConformers(Mol, MolCount)
        
        if not CalcStatus:
            ConfGenFailedCount += 1
            continue
        
        WriteMolConformers(Writer, ConformerMol, MolCount, ConfIDs, ConfEnergies)
    
    return (MolCount, ValidMolCount, ConfGenFailedCount)

def ProcessMoleculesUsingMultipleProcesses(Mols, Writer):
    """Process and minimize molecules using multiprocessing."""
    
    if OptionsInfo["SkipForceFieldMinimization"]:
        MiscUtil.PrintInfo("\nGenerating conformers without performing energy minimization using multiprocessing...")
    else:
        MiscUtil.PrintInfo("\nGenerating conformers and performing energy minimization using multiprocessing...")
    
    MPParams = OptionsInfo["MPParams"]
    
    # Setup data for initializing a worker process...
    InitializeWorkerProcessArgs = (MiscUtil.ObjectToBase64EncodedString(Options), MiscUtil.ObjectToBase64EncodedString(OptionsInfo))
    
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
    
    (MolCount, ValidMolCount, ConfGenFailedCount) = [0] * 3
    for Result in Results:
        MolCount += 1
        MolIndex, EncodedMol, CalcStatus, ConfIDs, ConfEnergies = Result
        
        if EncodedMol is None:
            continue
        ValidMolCount += 1
        
        if not CalcStatus:
            ConfGenFailedCount += 1
            continue
            
        Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
        WriteMolConformers(Writer, Mol, MolCount, ConfIDs, ConfEnergies)
    
    return (MolCount, ValidMolCount, ConfGenFailedCount)

def InitializeWorkerProcess(*EncodedArgs):
    """Initialize data for a worker process."""
    
    global Options, OptionsInfo

    MiscUtil.PrintInfo("Starting process (PID: %s)..." % os.getpid())
    
    # Decode Options and OptionInfo...
    Options = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[0])
    OptionsInfo = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[1])
    
def WorkerProcess(EncodedMolInfo):
    """Process data for a worker process."""

    MolIndex, EncodedMol = EncodedMolInfo
    
    CalcStatus = False
    ConfIDs = None
    ConfEnergies = None
    
    if EncodedMol is None:
        return [MolIndex, None, CalcStatus, ConfIDs, ConfEnergies]

    Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
    if RDKitUtil.IsMolEmpty(Mol):
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, (MolIndex + 1))
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
        return [MolIndex, None, CalcStatus, ConfIDs, ConfEnergies]
    
    Mol, CalcStatus, ConfIDs, ConfEnergies = GenerateMolConformers(Mol, (MolIndex + 1))

    return [MolIndex, RDKitUtil.MolToBase64EncodedMolString(Mol, PropertyPickleFlags = Chem.PropertyPickleOptions.MolProps | Chem.PropertyPickleOptions.PrivateProps), CalcStatus, ConfIDs, ConfEnergies]

def GenerateMolConformers(Mol, MolNum = None):
    "Generate conformers for a molecule"
    
    if OptionsInfo["SkipForceFieldMinimization"]:
        return GenerateMolConformersWithoutMinimization(Mol, MolNum)
    else:
        return GenerateMolConformersWithMinimization(Mol, MolNum)

def GenerateMolConformersWithoutMinimization(Mol, MolNum = None):
    "Generate conformers for a molecule without performing minimization."
    
    ConfIDs = EmbedMolecule(Mol, MolNum)
    if not len(ConfIDs):
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, MolNum)
            MiscUtil.PrintWarning("Conformation generation couldn't be performed for molecule %s: Embedding failed...\n" % MolName)
        return [Mol, False, None, None]
    
    if OptionsInfo["AlignConformers"]:
        AllChem.AlignMolConformers(Mol)
    
    if not OptionsInfo["QuietMode"]:
        MolName = RDKitUtil.GetMolName(Mol, MolNum)
        MiscUtil.PrintInfo("\nNumber of conformations generated for %s: %d" % (MolName, len(ConfIDs)))

    # Convert ConfIDs into a list...
    ConfIDsList = [ConfID for ConfID in ConfIDs]

    # Setup conformation energies...
    ConfEnergies = None
    if OptionsInfo["EnergyOut"]:
        ConfEnergies = []
        for ConfID in ConfIDsList:
            EnergyStatus, Energy = GetConformerEnergy(Mol, ConfID)
            
            Energy = "%.2f" % Energy if EnergyStatus else "NotAvailable"
            ConfEnergies.append(Energy)
            
            if not EnergyStatus:
                if not OptionsInfo["QuietMode"]:
                    MolName = RDKitUtil.GetMolName(Mol, MolNum)
                    MiscUtil.PrintWarning("Failed to retrieve calculated energy for conformation number %d of molecule %s. Try again after removing any salts or cleaing up the molecule...\n" % (ConfID, MolName))
    
    return [Mol, True, ConfIDsList, ConfEnergies]

def GenerateMolConformersWithMinimization(Mol, MolNum):
    "Generate and mininize conformers for a molecule."
    
    if  OptionsInfo["AddHydrogens"]:
        Mol = Chem.AddHs(Mol)
    
    ConfIDs = EmbedMolecule(Mol, MolNum)
    if not len(ConfIDs):
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, MolNum)
            MiscUtil.PrintWarning("Conformation generation couldn't be performed for molecule %s: Embedding failed...\n" % MolName)
        return [Mol, False, None, None]

    CalcEnergyMap = {}
    for ConfID in ConfIDs:
        try:
            if OptionsInfo["UseUFF"]:
                Status = AllChem.UFFOptimizeMolecule(Mol, confId = ConfID, maxIters = OptionsInfo["MaxIters"])
            elif OptionsInfo["UseMMFF"]:
                Status = AllChem.MMFFOptimizeMolecule(Mol, confId = ConfID, maxIters = OptionsInfo["MaxIters"], mmffVariant = OptionsInfo["MMFFVariant"])
            else:
                MiscUtil.PrintError("Minimization couldn't be performed: Specified forcefield, %s, is not supported" % OptionsInfo["ForceField"])
        except (RuntimeError, Chem.rdchem.KekulizeException) as ErrMsg:
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolNum)
                MiscUtil.PrintWarning("Minimization couldn't be performed for molecule %s:\n%s\n" % (MolName, ErrMsg))
            return [Mol, False, None, None]
        
        EnergyStatus, Energy = GetConformerEnergy(Mol, ConfID)
        if not EnergyStatus:
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolNum)
                MiscUtil.PrintWarning("Failed to retrieve calculated energy for conformation number %d of molecule %s. Try again after removing any salts or cleaing up the molecule...\n" % (ConfID, MolName))
            return [Mol, False, None, None]
        
        if Status != 0:
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolNum)
                MiscUtil.PrintWarning("Minimization failed to converge for conformation number %d of molecule %s in %d steps. Try using higher value for \"--maxIters\" option...\n" % (ConfID, MolName, OptionsInfo["MaxIters"]))
            
        CalcEnergyMap[ConfID] = Energy

    if  OptionsInfo["RemoveHydrogens"]:
        Mol = Chem.RemoveHs(Mol)

    # Align molecules after minimization...
    if OptionsInfo["AlignConformers"]:
        AllChem.AlignMolConformers(Mol)
    
    SortedConfIDs = sorted(ConfIDs, key = lambda ConfID: CalcEnergyMap[ConfID])
    
    MinEnergyConfID = SortedConfIDs[0]
    MinConfEnergy = CalcEnergyMap[MinEnergyConfID]
    EnergyWindow = OptionsInfo["EnergyWindow"]

    EnergyRMSDCutoff = OptionsInfo["EnergyRMSDCutoff"]
    ApplyEnergyRMSDCutoff = False
    if EnergyRMSDCutoff > 0:
        ApplyEnergyRMSDCutoff = True

    # Calculate RMSD values for conformers...
    PreAligned = False
    if OptionsInfo["AlignConformers"]:
        PreAligned = True
    
    CalcRMSDMap = {}
    if ApplyEnergyRMSDCutoff:
        for ConfID in SortedConfIDs:
            RMSD = AllChem.GetConformerRMS(Mol, MinEnergyConfID, ConfID, prealigned=PreAligned)
            CalcRMSDMap[ConfID] = RMSD
    
    # Track conformers with in the specified energy window  from the lowest
    # energy conformation along with applying RMSD cutoff as needed...
    #
    SelectedConfIDs = []
    
    ConfCount = 0
    IgnoredByEnergyConfCount = 0
    IgnoredByRMSDConfCount = 0
    
    FirstConf = True
    
    for ConfID in SortedConfIDs:
        if FirstConf:
            FirstConf = False
            SelectedConfIDs.append(ConfID)
            continue
            
        ConfEnergyDiff = abs(CalcEnergyMap[ConfID] - MinConfEnergy )
        if  ConfEnergyDiff > EnergyWindow:
            IgnoredByEnergyConfCount += 1
            continue
        
        if ApplyEnergyRMSDCutoff:
            if CalcRMSDMap[ConfID] < EnergyRMSDCutoff:
                IgnoredByRMSDConfCount += 1
                continue

        ConfCount += 1
        SelectedConfIDs.append(ConfID)
    
    if not OptionsInfo["QuietMode"]:
        MolName = RDKitUtil.GetMolName(Mol, MolNum)
        MiscUtil.PrintInfo("\nTotal Number of conformations generated for %s: %d" % (MolName, ConfCount))
        MiscUtil.PrintInfo("Number of conformations ignored due to energy window cutoff: %d" % (IgnoredByEnergyConfCount))
        if ApplyEnergyRMSDCutoff:
            MiscUtil.PrintInfo("Number of conformations ignored due to energy RMSD cutoff:  %d" % (IgnoredByRMSDConfCount))

    SelectedConfEnergies = None
    if OptionsInfo["EnergyOut"]:
        SelectedConfEnergies = ["%.2f" % CalcEnergyMap[ConfID] for ConfID in SelectedConfIDs]
        
    return [Mol, True, SelectedConfIDs, SelectedConfEnergies]
    
def GetConformerEnergy(Mol, ConfID):
    "Calculate conformer energy"

    Status = True
    Energy = 0.0
    
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
    
    # Figure out the number of conformations to embded...
    if re.match("^Auto$", OptionsInfo["MaxConfs"], re.I):
        NumOfRotBonds = Descriptors.NumRotatableBonds(Mol)
        if NumOfRotBonds <= 5:
            MaxConfs = 100
        elif NumOfRotBonds >=6 and NumOfRotBonds <= 10:
            MaxConfs = 200
        else:
            MaxConfs = 300
    else:
        MaxConfs = int(OptionsInfo["MaxConfs"])
        
    RandomSeed = OptionsInfo["RandomSeed"]
    EnforceChirality = OptionsInfo["EnforceChirality"]
    UseExpTorsionAnglePrefs = OptionsInfo["UseExpTorsionAnglePrefs"]
    UseBasicKnowledge = OptionsInfo["UseBasicKnowledge"]
    EmbedRMSDCutoff = OptionsInfo["EmbedRMSDCutoff"]
    
    try:
        ConfIDs = AllChem.EmbedMultipleConfs(Mol, numConfs = MaxConfs, randomSeed = RandomSeed, pruneRmsThresh = EmbedRMSDCutoff, enforceChirality = EnforceChirality, useExpTorsionAnglePrefs = UseExpTorsionAnglePrefs, useBasicKnowledge = UseBasicKnowledge)
    except ValueError as ErrMsg:
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, MolNum)
            MiscUtil.PrintWarning("Embedding failed  for molecule %s:\n%s\n" % (MolName, ErrMsg))
        ConfIDs = []

    return ConfIDs

def WriteMolConformers(Writer, Mol, MolNum, ConfIDs, ConfEnergies = None):
    """Write molecule coformers. """

    if ConfIDs is None:
        return
    
    MolName = RDKitUtil.GetMolName(Mol, MolNum)

    for Index, ConfID in enumerate(ConfIDs):
        SetConfMolName(Mol, MolName, ConfID)
        
        if ConfEnergies is not None:
            Mol.SetProp(OptionsInfo["EnergyLabel"], ConfEnergies[Index])
            
        Writer.write(Mol, confId = ConfID)
    
def SetConfMolName(Mol, MolName, ConfCount):
    """Set conf mol name"""

    ConfName = "%s_Conf%d" % (MolName, ConfCount)
    Mol.SetProp("_Name", ConfName)

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    OptionsInfo["AddHydrogens"] = True
    if re.match("^no$", Options["--addHydrogens"], re.I):
        OptionsInfo["AddHydrogens"] = False

    OptionsInfo["AlignConformers"] = True
    if re.match("^no$", Options["--alignConformers"], re.I):
        OptionsInfo["AlignConformers"] = False

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
    
    OptionsInfo["ConformerGenerator"] = ConformerGenerator
    OptionsInfo["UseExpTorsionAnglePrefs"] = UseExpTorsionAnglePrefs
    OptionsInfo["UseBasicKnowledge"] = UseBasicKnowledge

    if re.match("^UFF$", Options["--forceField"], re.I):
        ForceField = "UFF"
        UseUFF = True
        UseMMFF = False
        SkipForceFieldMinimization = False
    elif re.match("^MMFF$", Options["--forceField"], re.I):
        ForceField = "MMFF"
        UseUFF = False
        UseMMFF = True
        SkipForceFieldMinimization = False
    else:
        ForceField = "None"
        UseUFF = False
        UseMMFF = False
        SkipForceFieldMinimization = True
        
    MMFFVariant = "MMFF94" if re.match("^MMFF94$", Options["--forceFieldMMFFVariant"], re.I) else "MMFF94s"
    
    OptionsInfo["SkipForceFieldMinimization"] = SkipForceFieldMinimization
    OptionsInfo["ForceField"] = ForceField
    OptionsInfo["MMFFVariant"] = MMFFVariant
    OptionsInfo["UseMMFF"] = UseMMFF
    OptionsInfo["UseUFF"] = UseUFF
        
    OptionsInfo["EnergyOut"] = True if re.match("^yes$", Options["--energyOut"], re.I) else False
    if UseUFF:
        EnergyLabel = "UFF_Energy"
    elif UseMMFF:
        EnergyLabel = "%s_Energy" % MMFFVariant
    else:
        EnergyLabel = "Energy"
    OptionsInfo["EnergyLabel"] = EnergyLabel
    
    OptionsInfo["EnforceChirality"] = True
    if re.match("^no$", Options["--enforceChirality"], re.I):
        OptionsInfo["EnforceChirality"] = False
    
    OptionsInfo["EnergyWindow"] = float(Options["--energyWindow"])
    
    EmbedRMSDCutoff = -1.0
    if not re.match("^none$", Options["--embedRMSDCutoff"], re.I):
        EmbedRMSDCutoff = float(Options["--embedRMSDCutoff"])
    OptionsInfo["EmbedRMSDCutoff"] = EmbedRMSDCutoff
    
    EnergyRMSDCutoff = -1.0
    if not re.match("^none$", Options["--energyRMSDCutoff"], re.I):
        EnergyRMSDCutoff = float(Options["--energyRMSDCutoff"])
    OptionsInfo["EnergyRMSDCutoff"] = EnergyRMSDCutoff
    
    OptionsInfo["MaxIters"] = int(Options["--maxIters"])
    OptionsInfo["MaxConfs"] = Options["--maxConfs"]
    
    OptionsInfo["MPMode"] = True if re.match("^yes$", Options["--mp"], re.I) else False
    OptionsInfo["MPParams"] = MiscUtil.ProcessOptionMultiprocessingParameters("--mpParams", Options["--mpParams"])
    
    OptionsInfo["QuietMode"] = True if re.match("^yes$", Options["--quiet"], re.I) else False
    
    RandomSeed = -1
    if not re.match("^auto$", Options["--randomSeed"], re.I):
        RandomSeed = int(Options["--randomSeed"])
    OptionsInfo["RandomSeed"] = RandomSeed
    
    OptionsInfo["RemoveHydrogens"] = True
    if re.match("^no$", Options["--removeHydrogens"], re.I):
        OptionsInfo["RemoveHydrogens"] = False
    
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
    MiscUtil.ValidateOptionTextValue("--alignConformers", Options["--alignConformers"], "yes no")
    MiscUtil.ValidateOptionTextValue("-c, --conformerGenerator", Options["--conformerGenerator"], "SDG ETDG KDG ETKDG")
    
    MiscUtil.ValidateOptionTextValue("-f, --forceField", Options["--forceField"], "UFF MMFF None")
    MiscUtil.ValidateOptionTextValue(" --forceFieldMMFFVariant", Options["--forceFieldMMFFVariant"], "MMFF94 MMFF94s")
    
    MiscUtil.ValidateOptionTextValue("--energyOut", Options["--energyOut"], "yes no")
    MiscUtil.ValidateOptionTextValue("--enforceChirality ", Options["--enforceChirality"], "yes no")
    MiscUtil.ValidateOptionFloatValue("--energyWindow", Options["--energyWindow"], {">": 0.0})
    
    if not re.match("^none$", Options["--embedRMSDCutoff"], re.I):
        MiscUtil.ValidateOptionFloatValue("--embedRMSDCutoff", Options["--embedRMSDCutoff"], {">": 0.0})
    
    if not re.match("^none$", Options["--energyRMSDCutoff"], re.I):
        MiscUtil.ValidateOptionFloatValue("--energyRMSDCutoff", Options["--energyRMSDCutoff"], {">": 0.0})
        # Make sure that the alignConformers option is being used...
        if not re.match("^yes$", Options["--alignConformers"], re.I):
            MiscUtil.PrintError("%s value of \"--alignConformers\" is not allowed for %s value of \"--energyRMSDCutoff\" option " % (Options["--alignConformers"], Options["--energyRMSDCutoff"]))
        
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi txt csv tsv")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
        
    if not re.match("^auto$", Options["--maxConfs"], re.I):
        MiscUtil.ValidateOptionIntegerValue("--maxConfs", Options["--maxConfs"], {">": 0})
    
    MiscUtil.ValidateOptionIntegerValue("--maxIters", Options["--maxIters"], {">": 0})
    
    MiscUtil.ValidateOptionTextValue("--mp", Options["--mp"], "yes no")
    MiscUtil.ValidateOptionTextValue("-q, --quiet", Options["--quiet"], "yes no")
    
    if not re.match("^auto$", Options["--randomSeed"], re.I):
        MiscUtil.ValidateOptionIntegerValue("--randomSeed", Options["--randomSeed"], {})
        RandomSeed = int(Options["--randomSeed"])
    
    MiscUtil.ValidateOptionTextValue("-r, --removeHydrogens", Options["--removeHydrogens"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitGenerateConformers.py - Generate molecular conformations

Usage:
    RDKitGenerateConformers.py [--alignConformers <yes or no>] [--addHydrogens <yes or no>]
                               [--conformerGenerator <SDG, ETDG, KDG, ETKDG>] [--embedRMSDCutoff <number>]
                               [--energyOut  <yes or no>] [--enforceChirality <yes or no>] [--energyRMSDCutoff <number>]
                               [--energyWindow <number> ]  [--forceField <UFF, MMFF, None>] [--forceFieldMMFFVariant <MMFF94 or MMFF94s>]
                               [--infileParams <Name,Value,...>] [--maxConfs <number>]
                               [--mp <yes or no>] [--mpParams <Name.Value,...>]
                               [--maxIters <number>]  [ --outfileParams <Name,Value,...> ]  [--overwrite]
                               [--quiet <yes or no>] [ --removeHydrogens <yes or no>] [--randomSeed <number>]
                               [-w <dir>] -i <infile> -o <outfile> 
    RDKitGenerateConformers.py -h | --help | -e | --examples

Description:
    Generate molecular conformations using a combination of distance geometry and
    forcefield minimization. The forcefield minimization may be skipped to only generate
    conformations by available distance geometry based methodologies.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi,
    .csv, .tsv .txt)

    The supported output file format are: SD (.sdf, .sd)

Options:
    -a, --addHydrogens <yes or no>  [default: yes]
        Add hydrogens before minimization.
    --alignConformers <yes or no>  [default: yes]
        Align conformers for each molecule.
    -c, --conformerGenerator <SDG, ETDG, KDG, ETKDG>  [default: ETKDG]
        Conformation generation methodology for generating initial 3D coordinates of a
        molecule. Possible values: Standard Distance Geometry, (SDG), Experimental
        Torsion-angle preference with Distance Geometry (ETDG)  [Ref 129] , basic Knowledge-terms
        with Distance Geometry (KDG), and Experimental Torsion-angle preference along
        with basic Knowledge-terms and Distance Geometry (ETKDG). 
    --embedRMSDCutoff <number>  [default: none]
        RMSD cutoff for retaining conformations after embedding and before energy minimization.
        All embedded conformations are kept by default. Otherwise, only those conformations
        which are different from each other by the specified RMSD cutoff are kept. The first
        embedded conformation is always retained.
    --energyOut <yes or no>  [default: No]
        Write out energy values.
    --enforceChirality <yes or no>  [default: Yes]
        Enforce chirality for defined chiral centers.
    --energyRMSDCutoff <number>  [default: none]
        RMSD cutoff for retaining conformations after energy minimization. By default,
        all minimized conformations with in the specified energy window from the lowest energy
        conformation are kept. Otherwise, only those conformations which are different from
        the lowest energy conformation by the specified RMSD cutoff and are with in the
        specified energy window are kept. The lowest energy conformation is always retained.
    --energyWindow <number>  [default: 20]
        Energy window in kcal/mol for selecting conformers. This option is ignored during
        'None' value of '-f, --forcefield' option.
    -e, --examples
        Print examples.
    -f, --forceField <UFF, MMFF, None>  [default: MMFF]
        Forcefield method to use for energy minimization. Possible values: Universal Force
        Field (UFF) [Ref 81],  Merck Molecular Mechanics Force Field (MMFF) [Ref 83-87] or
        None.
    --forceFieldMMFFVariant <MMFF94 or MMFF94s>  [default: MMFF94]
        Variant of MMFF forcefield to use for energy minimization.
    -h, --help
        Print this help message.
    -i, --infile <infile>
        Input file name.
    --infileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for reading
        molecules from files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD, MOL: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                smilesTitleLine,auto,sanitize,yes
            
        Possible values for smilesDelimiter: space, comma or tab.
    --maxConfs <number>  [default: auto]
        Maximum number of conformations to generate for each molecule by conformation
        generation methodology. The conformations are minimized using the specified
        forcefield as needed and written to the output file. The default value for maximum
        number of conformations is dependent on the number of rotatable bonds in molecules:
        RotBonds <= 5, maxConfs = 100; RotBonds >=6 and <= 10, MaxConfs = 200;
        RotBonds >= 11, maxConfs = 300
    --maxIters <number>  [default: 250]
        Maximum number of iterations to perform for each molecule during forcefield
        minimization.
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
        Output file name.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: kekulize,no
            
    --overwrite
        Overwrite existing files.
    -q, --quiet <yes or no>  [default: no]
        Use quiet mode. The warning and information messages will not be printed.
    -r, --removeHydrogens <yes or no>  [default: Yes]
        Remove hydrogens after minimization.
    --randomSeed <number>  [default: auto]
        Seed for the random number generator for reproducing 3D coordinates.
        Default is to use a random seed.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To generate conformers using Experimental Torsion-angle preference along
    with basic Knowledge-terms and Distance Geometry (ETKDG) followed by
    MMFF minimization with automatic determination of maximum number of
    conformers for each molecule and write out a SD file, type:

        % RDKitGenerateConformers.py  -i Sample.smi -o SampleOut.sdf

    To rerun the first example in a quiet mode and write out a SD file, type:

        % RDKitGenerateConformers.py -q yes -i Sample.smi -o SampleOut.sdf

    To rerun the first example in multiprocessing mode on all available CPUs
    without loading all data into memory and write out a SD file, type:

        % RDKitGenerateConformers.py --mp yes -i Sample.smi -o SampleOut.sdf

    To run the first example in multiprocessing mode on all available CPUs
    by loading all data into memory and write out a SD file, type:

        % RDKitGenerateConformers.py --mp yes --mpParams "inputDataMode,
          InMemory" -i Sample.smi -o SampleOut.sdf

    To rerun the first example in multiprocessing mode on specific number of
    CPUs and chunk size without loading all data into memory and write out a SD file,
    type:

        % RDKitGenerateConformers.py --mp yes --mpParams "inputDataMode,Lazy,
          numProcesses,4,chunkSize,8" -i Sample.smi -o SampleOut.sdf

    To generate up to 150 conformers for each molecule using ETKDG and UFF forcefield
    minimization along with conformers within 25 kcal/mol energy window and write out a
    SD file, type:

        % RDKitGenerateConformers.py  --energyWindow 25 -f UFF --maxConfs 150
          -i Sample.smi -o SampleOut.sdf

    To generate up to 50 conformers for each molecule using KDG without any forcefield
    minimization and alignment of conformers and write out a SD file, type:

        % RDKitGenerateConformers.py  -f none --maxConfs 50 --alignConformers no
          -i Sample.sdf -o SampleOut.sdf

    To generate up to 50 conformers using SDG without any forcefield minimization
    and alignment of conformers for molecules in a  CSV SMILES file, SMILES strings
    in column 1, name in column 2, and write out a SD file, type:

        % RDKitGenerateConformers.py  --maxConfs 50  --maxIters 50 -c SDG
          --alignConformers no -f none --infileParams "smilesDelimiter,comma,
          smilesTitleLine,yes, smilesColumn,1,smilesNameColumn,2"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitCalculateRMSD.py, RDKitCalculateMolecularDescriptors.py,
    RDKitCompareMoleculeShapes.py, RDKitConvertFileFormat.py,
    RDKitGenerateConstrainedConformers.py, RDKitPerformMinimization.py

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
