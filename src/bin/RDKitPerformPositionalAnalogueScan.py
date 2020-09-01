#!/usr/bin/env python
#
# File: RDKitPerformPositionalAnalogueScan.py
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
from itertools import combinations
from itertools import permutations

# RDKit imports...
try:
    from rdkit import rdBase
    from rdkit import Chem
    from rdkit.Chem import AllChem
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
    PerformPositionalAnalogueScan()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformPositionalAnalogueScan():
    """Perform positional analogue scan."""
    
    # Setup a molecule reader for input file...
    MiscUtil.PrintInfo("\nProcessing file %s..." % OptionsInfo["Infile"])
    OptionsInfo["InfileParams"]["AllowEmptyMols"] = True
    Mols  = RDKitUtil.ReadMolecules(OptionsInfo["Infile"], **OptionsInfo["InfileParams"])

    # Set up a molecule writer...
    Writer = RDKitUtil.MoleculesWriter(OptionsInfo["Outfile"], **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % OptionsInfo["Outfile"])
    MiscUtil.PrintInfo("Generating file %s..." % OptionsInfo["Outfile"])

    MolCount, ValidMolCount, SearchPatternMissingCount, AnaloguesGenerationFailedCount, AnaloguesCount = ProcessMolecules(Mols, Writer)

    if Writer is not None:
        Writer.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    if not OptionsInfo["RxnSMARTSMode"]:
        MiscUtil.PrintInfo("Number of molecules with missing search pattern: %d" % SearchPatternMissingCount)
    MiscUtil.PrintInfo("Number of molecules failed during generation of analogues: %d" % AnaloguesGenerationFailedCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount + SearchPatternMissingCount + AnaloguesGenerationFailedCount))
    
    MiscUtil.PrintInfo("\nTotal number of analogues for valid molecules: %d" % AnaloguesCount)
    MiscUtil.PrintInfo("Average number of analogues per valid molecule: %d" % (int(AnaloguesCount/(ValidMolCount))))

def ProcessMolecules(Mols, Writer):
    """Process molecules to generate analogues."""
    
    if OptionsInfo["MPMode"]:
        return ProcessMoleculesUsingMultipleProcesses( Mols, Writer)
    else:
        return ProcessMoleculesUsingSingleProcess(Mols, Writer)

def ProcessMoleculesUsingSingleProcess(Mols, Writer):
    """Process molecules to generate analogues using a single process."""

    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    SetSMILESMolProps = OptionsInfo["OutfileParams"]["SetSMILESMolProps"]

    (MolCount, ValidMolCount, SearchPatternMissingCount, AnaloguesGenerationFailedCount, AnaloguesCount) = [0] * 5

    # Initialize search and target patterns info...
    SetupSearchAndTargetPatternsInfo()

    FirstMol = True
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

        MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus = GenerateMolAnalogues(Mol, MolCount)
        
        if FirstMol:
            FirstMol = False
            if SetSMILESMolProps:
                if Writer is not None:
                    RDKitUtil.SetWriterMolProps(Writer, Mol)
        
        if not OptionsInfo["RxnSMARTSMode"]:
            if SearchPatternMissingStatus:
                SearchPatternMissingCount += 1
                WriteMolecule(Writer, Mol, Compute2DCoords)
                continue
        
        if not AnaloguesGenerationStatus:
            AnaloguesGenerationFailedCount += 1
            WriteMolecule(Writer, Mol, Compute2DCoords)
            continue
        
        # Write out molecues and analogues...
        WriteMolAnalogues(Writer, Mol, MolCount, MolAnalogues, Compute2DCoords)
        
        AnaloguesCount += len(MolAnalogues)

    return (MolCount, ValidMolCount, SearchPatternMissingCount, AnaloguesGenerationFailedCount, AnaloguesCount)

def ProcessMoleculesUsingMultipleProcesses(Mols, Writer):
    """Process molecules to generate analogues using multiprocessing."""

    MPParams = OptionsInfo["MPParams"]
    
    # Setup data for initializing a worker process...
    MiscUtil.PrintInfo("Encoding options info...")

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
    
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    SetSMILESMolProps = OptionsInfo["OutfileParams"]["SetSMILESMolProps"]

    (MolCount, ValidMolCount, SearchPatternMissingCount, AnaloguesGenerationFailedCount, AnaloguesCount) = [0] * 5
    FirstMol = True
    for Result in Results:
        MolCount += 1
        MolIndex, EncodedMol, EncodedMolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus  = Result
        
        if EncodedMol is None:
            continue
        ValidMolCount += 1

        Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
        
        if FirstMol:
            FirstMol = False
            if SetSMILESMolProps:
                if Writer is not None:
                    RDKitUtil.SetWriterMolProps(Writer, Mol)
        
        if not OptionsInfo["RxnSMARTSMode"]:
            if SearchPatternMissingStatus:
                SearchPatternMissingCount += 1
                WriteMolecule(Writer, Mol, Compute2DCoords)
                continue
        
        if not AnaloguesGenerationStatus:
            AnaloguesGenerationFailedCount += 1
            WriteMolecule(Writer, Mol, Compute2DCoords)
            continue
        
        MolAnalogues = [RDKitUtil.MolFromBase64EncodedMolString(EncodedMolAnalogue) for EncodedMolAnalogue in EncodedMolAnalogues]
        
        # Write out molecues and analogues...
        WriteMolAnalogues(Writer, Mol, MolCount, MolAnalogues, Compute2DCoords)
        
        AnaloguesCount += len(MolAnalogues)
    
    return (MolCount, ValidMolCount, SearchPatternMissingCount, AnaloguesGenerationFailedCount, AnaloguesCount)
    
def InitializeWorkerProcess(*EncodedArgs):
    """Initialize data for a worker process."""
    
    global Options, OptionsInfo

    MiscUtil.PrintInfo("Starting process (PID: %s)..." % os.getpid())

    # Decode Options and OptionInfo...
    Options = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[0])
    OptionsInfo = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[1])
    
    # Initialize search and target  patterns info...
    SetupSearchAndTargetPatternsInfo()

def WorkerProcess(EncodedMolInfo):
    """Process data for a worker process."""

    MolIndex, EncodedMol = EncodedMolInfo

    MolAnalogues = None
    SearchPatternMissingStatus = False
    AnaloguesGenerationStatus = True
    
    if EncodedMol is None:
        return [MolIndex, None, MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus]

    Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
    if RDKitUtil.IsMolEmpty(Mol):
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, (MolIndex + 1))
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
        return [MolIndex, None, MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus]

    MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus = GenerateMolAnalogues(Mol, (MolIndex +1))

    EncodedMolAnalogues = None
    if MolAnalogues is not None:
        EncodedMolAnalogues = [RDKitUtil.MolToBase64EncodedMolString(MolAnalogue, PropertyPickleFlags = Chem.PropertyPickleOptions.MolProps | Chem.PropertyPickleOptions.PrivateProps) for MolAnalogue in MolAnalogues]
    
    return [MolIndex, EncodedMol, EncodedMolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus]

def GenerateMolAnalogues(Mol, MolNum = None):
    """Generate analogues."""

    if OptionsInfo["AttachAtomsMode"] or OptionsInfo["ReplaceAtomsMode"] or OptionsInfo["AttachSMILESMode"]:
        return GenerateMolAnaloguesForMultipleScanModes(Mol, MolNum)
    elif OptionsInfo["RxnSMARTSMode"]:
        return GenerateMolAnaloguesForRxnSMARTSMode(Mol, MolNum)
    elif OptionsInfo["TandemScanMode"]:
        return GenerateMolAnaloguesForTandemScanMode(Mol, MolNum)
    else:
        MiscUtil.PrintError("Failed to generate analogues.  The mode value, %s, specified using \"-m, --mode\" is not supported." % (OptionsInfo["Mode"]))

def GenerateMolAnaloguesForMultipleScanModes(Mol, MolNum = None):
    """Generate analogues for the following scan modes: AttachAtom, ReplaceAtom
    or AttachSMILES."""

    MolAnalogues = None
    SearchPatternMissingStatus = False
    AnaloguesGenerationStatus = True

    # Retrieve matched atoms...
    MolMatchedAtomIndices = SetupMolMatchedAtomIndices(Mol, MolNum)
    if MolMatchedAtomIndices is None:
        SearchPatternMissingStatus = True
        return (MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus)

    AtomicNumber = OptionsInfo["TargetPatternAtomicNumber"]
    SMILESMol = OptionsInfo["TargetPatternSMILESMol"]
    
    MolAnalogues = []
    MolAnaloguesSet = set()
    for AtomIndicesCombination in combinations(MolMatchedAtomIndices, OptionsInfo["ComboSearchAtoms"]):
        MolAnalogue = Chem.RWMol(Mol)

        for AtomIndex in AtomIndicesCombination:
            if OptionsInfo["AttachAtomsMode"]:
                MolAnalogue = AttachAtom(MolAnalogue, AtomIndex, AtomicNumber, BondOrder = Chem.rdchem.BondType.SINGLE)
            elif OptionsInfo["ReplaceAtomsMode"]:
                MolAnalogue = ReplaceAtom(MolAnalogue, AtomIndex, AtomicNumber)
            elif OptionsInfo["AttachSMILESMode"]:
                AttachAtomIndex = MolAnalogue.GetNumAtoms()
                MolAnalogue = AttachMol(MolAnalogue, AtomIndex, SMILESMol, AttachAtomIndex, BondOrder = Chem.rdchem.BondType.SINGLE)
            
        if not SanitizeMolAnalogue(MolAnalogue, MolNum):
            continue

        TrackAnalogues(MolAnalogue, MolAnaloguesSet, MolAnalogues)
        
    if len(MolAnalogues) == 0:
        AnaloguesGenerationStatus = False
        
    return (MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus)

def GenerateMolAnaloguesForRxnSMARTSMode(Mol, MolNum = None):
    """Generate analogues for reaction SMARTS mode."""

    MolAnalogues = None
    SearchPatternMissingStatus = False
    AnaloguesGenerationStatus = True
    
    Rxn = OptionsInfo["TargetPatternRxn"]
    RxnProducts = Rxn.RunReactants((Mol,))

    if not len(RxnProducts):
        MolAnalogues = None
        AnaloguesGenerationStatus = False
        return (MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus)

    MolAnalogues = []
    MolAnaloguesSet = set()
    for RxnProduct in RxnProducts:
        MolAnalogue = RxnProduct[0]
        
        if not SanitizeMolAnalogue(MolAnalogue, MolNum):
            continue
        
        TrackAnalogues(MolAnalogue, MolAnaloguesSet, MolAnalogues)
    
    if len(MolAnalogues) == 0:
        AnaloguesGenerationStatus = False
        
    return (MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus)

def GenerateMolAnaloguesForTandemScanMode(Mol, MolNum = None):
    """Generate analogues for tandem scan mode."""

    MolAnalogues = None
    SearchPatternMissingStatus = False
    AnaloguesGenerationStatus = True

    # Retrieve matched atoms...
    MolMatchedAtomIndices = SetupMolMatchedAtomIndices(Mol, MolNum)
    if MolMatchedAtomIndices is None:
        SearchPatternMissingStatus = True
        return (MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus)

    TargetPatternTandemScanInfo = OptionsInfo["TargetPatternTandemScanInfo"]
    NumOfTandemOperations = len(TargetPatternTandemScanInfo["Names"])
    if  len(MolMatchedAtomIndices) < NumOfTandemOperations:
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, MolNum)
            MiscUtil.PrintWarning("Number of matched search pattern atoms, %s, must be >= number of specified operations, %s, during TandemScan mode for molecule %s..." % (len(MolMatchedAtomIndices), NumOfTandemOperations, MolName))
        SearchPatternMissingStatus = True
        return (MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus)
    
    MolAnalogues = []
    MolAnaloguesSet = set()
    for AtomIndicesCombination in combinations(MolMatchedAtomIndices, NumOfTandemOperations):
        AtomIndicesPermutations = []
        if OptionsInfo["PermuteTandemScan"]:
            for AtomIndicesPermutation in permutations(AtomIndicesCombination, NumOfTandemOperations):
                AtomIndicesPermutations.append(list(AtomIndicesPermutation))
        else:
            AtomIndicesPermutations.append(list(AtomIndicesCombination))
        
        for AtomIndicesPermutation in AtomIndicesPermutations:
            MolAnalogue = Chem.RWMol(Mol)
            
            for NameIndex, Name in enumerate(TargetPatternTandemScanInfo["Names"]):
                AtomIndex = AtomIndicesPermutation[NameIndex]
                
                if re.match("^AttachAtoms$", Name, re.I):
                    MolAnalogue = AttachAtom(MolAnalogue, AtomIndex, TargetPatternTandemScanInfo["PatternAtomicNumber"][Name], BondOrder = Chem.rdchem.BondType.SINGLE)
                elif re.match("^ReplaceAtoms$", Name, re.I):
                    MolAnalogue = ReplaceAtom(MolAnalogue, AtomIndex, TargetPatternTandemScanInfo["PatternAtomicNumber"][Name])
                elif re.match("^AttachSMILES$", Name, re.I):
                    AttachAtomIndex = MolAnalogue.GetNumAtoms()
                    MolAnalogue = AttachMol(MolAnalogue, AtomIndex, TargetPatternTandemScanInfo["PatternSMILESMol"][Name], AttachAtomIndex, BondOrder = Chem.rdchem.BondType.SINGLE)
                    MolAnalogue = Chem.RWMol(MolAnalogue)
            
            if not SanitizeMolAnalogue(MolAnalogue, MolNum):
                continue
            
            TrackAnalogues(MolAnalogue, MolAnaloguesSet, MolAnalogues)
    
    if len(MolAnalogues) == 0:
        AnaloguesGenerationStatus = False
        
    return (MolAnalogues, SearchPatternMissingStatus, AnaloguesGenerationStatus)

def TrackAnalogues(MolAnalogue, MolAnaloguesSet, MolAnalogues):
    """Track analogues. """
    
    if OptionsInfo["UniqueAnalogues"]:
        AnalogueSMILES = Chem.MolToSmiles(Chem.RemoveHs(MolAnalogue))
        if AnalogueSMILES not in MolAnaloguesSet:
            MolAnaloguesSet.add(AnalogueSMILES)
            MolAnalogues.append(MolAnalogue)
    else:
        MolAnalogues.append(MolAnalogue)
    
def ReplaceAtom(Mol, AtomIndex, AtomicNumber):
    """Replace an existing atom with a new atom in a molecule."""
    
    Atom = Mol.GetAtomWithIdx(AtomIndex)
    Atom.SetAtomicNum(AtomicNumber)

    return Mol

def AttachAtom(Mol, AtomIndex, AtomicNumber, BondOrder = Chem.rdchem.BondType.SINGLE):
    """Attach atom to an existing atom in a molecule."""
    
    NewAtomIndex = Mol.AddAtom(Chem.Atom(AtomicNumber))
    Mol.AddBond(AtomIndex, NewAtomIndex, BondOrder)

    return Mol

def AttachMol(Mol, AtomIndex, AttachMol, AttachAtomIndex, BondOrder = Chem.rdchem.BondType.SINGLE):
    """Attach molecule to an exisiting molecule."""
    
    NewMol = Chem.CombineMols(Mol, AttachMol)
    
    NewMol = Chem.EditableMol(NewMol)
    NewMol.AddBond(AtomIndex, AttachAtomIndex, order = BondOrder)
    NewMol = NewMol.GetMol()
    
    return NewMol

def SanitizeMolAnalogue(Mol, MolNum):
    """Sanitize analogue and catch any errors. """
    
    Status = True
    try:
        Chem.SanitizeMol(Mol)
    except Exception as ErrMsg:
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, MolNum)
            MiscUtil.PrintInfo("\n%s\nFailed to sanitize analogue of molecule %s...\nIgnoring analogue..." % (ErrMsg, MolName))
        Status = False

    return Status

def SetupMolMatchedAtomIndices(Mol, MolNum):
    """Setup a list of matched atom indices corresponding to first atom in search
    SMARTS pattern for a molecule."""

    MatchedAtoms = []
    
    # Retrieve first atom match for PAS...
    MatchedAtomIndices = [AtomIndices[0] for AtomIndices in Mol.GetSubstructMatches(OptionsInfo["SearchPatternMol"])]
    if len(MatchedAtomIndices) == 0:
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, MolNum)
            MiscUtil.PrintWarning("No search pattern match for molecule %s..." % MolName)
        return None
    
    return MatchedAtomIndices

def WriteMolAnalogues(Writer, Mol, MolNum, MolAnalogues, Compute2DCoords):
    """Write molecule analogues. """

    # Write out molecule...
    WriteMolecule(Writer, Mol, Compute2DCoords)
    
    if MolAnalogues is None:
        return

    # Write out analogues...
    for Index, MolAnalogue in enumerate(MolAnalogues):
        MolName = RDKitUtil.GetMolName(Mol, MolNum)
        SetMolAnalogueName(MolAnalogue, MolName, (Index + 1))

        WriteMolecule(Writer, MolAnalogue, Compute2DCoords)
    
def SetMolAnalogueName(MolAnalogue, MolName, AnalogueCount):
    """Set analogue name."""
    
    NameSuffix = OptionsInfo["NameSuffix"]
    AnalogueName = "%s_%s%d" % (MolName, NameSuffix, AnalogueCount)
    MolAnalogue.SetProp("_Name", AnalogueName)

def WriteMolecule(Writer, Mol, Compute2DCoords):
    """Write out molecule."""
    
    if Compute2DCoords:
        AllChem.Compute2DCoords(Mol)
    
    Writer.write(Mol)

def SetupSearchAndTargetPatternsInfo():
    """Setup search and target patterns info."""

    # Setup search pattern molecule...
    SearchPattern = OptionsInfo["SearchPattern"]
    SearchPatternMol = None
    if SearchPattern is not None:
        SearchPatternMol = Chem.MolFromSmarts(SearchPattern)
        if SearchPatternMol is None:
            MiscUtil.PrintError("Failed to create search pattern molecule. The torsion SMILES/SMARTS pattern, \"%s\", specified using \"-s, --searchPattern\" option is not valid." % (SearchPattern))
    OptionsInfo["SearchPatternMol"] = SearchPatternMol

    # Setup target pattern information...
    Mode = OptionsInfo["Mode"]
    TargetPattern = OptionsInfo["TargetPattern"]
    
    OptionsInfo["TargetPatternAtomicNumber"] = None
    OptionsInfo["TargetPatternSMILESMol"] = None
    OptionsInfo["TargetPatternRxn"] = None
    
    if re.match("^(AttachAtoms|ReplaceAtoms)$", Mode, re.I):
        OptionsInfo["TargetPatternAtomicNumber"] = Chem.GetPeriodicTable().GetAtomicNumber(TargetPattern)
    elif re.match("^AttachSMILES$", Mode, re.I):
        OptionsInfo["TargetPatternSMILESMol"] = Chem.MolFromSmiles(TargetPattern)
    elif re.match("^RxnSMARTS$", Mode, re.I):
        OptionsInfo["TargetPatternRxn"] = AllChem.ReactionFromSmarts(TargetPattern)
    elif re.match("^TandemScan$", Mode, re.I):
        # Setup atomic numbers and molecules...
        TargetPatternTandemScanInfo = OptionsInfo["TargetPatternTandemScanInfo"]
        
        TargetPatternTandemScanInfo["PatternAtomicNumber"] = {}
        TargetPatternTandemScanInfo["PatternSMILESMol"] = {}
        
        for Name in TargetPatternTandemScanInfo["Names"]:
            Pattern = TargetPatternTandemScanInfo["Pattern"][Name]
            if re.match("^(AttachAtoms|ReplaceAtoms)$", Name, re.I):
                TargetPatternTandemScanInfo["PatternAtomicNumber"][Name] = Chem.GetPeriodicTable().GetAtomicNumber(Pattern)
                TargetPatternTandemScanInfo["PatternSMILESMol"][Name] = None
            elif re.match("^AttachSMILES$", Name, re.I):
                TargetPatternTandemScanInfo["PatternSMILESMol"][Name] = Chem.MolFromSmiles(Pattern)
                TargetPatternTandemScanInfo["PatternAtomicNumber"][Name] = None
    
def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()

    OptionsInfo["ComboSearchAtoms"] = int(Options["--comboSearchAtoms"])

    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    OptionsInfo["Mode"] = Options["--mode"]
    (AttachAtomsMode, ReplaceAtomsMode, AttachSMILESMode, RxnSMARTSMode, TandemScanMode) = [False] * 5
    Mode = OptionsInfo["Mode"]
    if re.match("^AttachAtoms$", Mode, re.I):
        AttachAtomsMode = True
    elif re.match("^ReplaceAtoms$", Mode, re.I):
        ReplaceAtomsMode = True
    elif re.match("^AttachSMILES$", Mode, re.I):
        AttachSMILESMode = True
    elif re.match("^RxnSMARTS$", Mode, re.I):
        RxnSMARTSMode = True
    elif re.match("^TandemScan$", Mode, re.I):
        TandemScanMode = True
    else:
        MiscUtil.PrintError("The mode value , %s, specified using \"-m, --mode\" option is not valid" % (Mode))
    OptionsInfo["AttachAtomsMode"] = AttachAtomsMode
    OptionsInfo["ReplaceAtomsMode"] = ReplaceAtomsMode
    OptionsInfo["AttachSMILESMode"] = AttachSMILESMode
    OptionsInfo["RxnSMARTSMode"] = RxnSMARTSMode
    OptionsInfo["TandemScanMode"] = TandemScanMode

    OptionsInfo["MPMode"] = True if re.match("^yes$", Options["--mp"], re.I) else False
    OptionsInfo["MPParams"] = MiscUtil.ProcessOptionMultiprocessingParameters("--mpParams", Options["--mpParams"])
    
    OptionsInfo["NameSuffix"] = Options["--nameSuffix"]
    OptionsInfo["PermuteTandemScan"] = True if re.match("^yes$", Options["--permuteTandemScan"], re.I) else False
    OptionsInfo["QuietMode"] = True if re.match("^yes$", Options["--quiet"], re.I) else False
    
    OptionsInfo["UniqueAnalogues"] = True if re.match("^yes$", Options["--uniqueAnalogues"], re.I) else False
    
    # Process and validate search pattern...
    SearchPattern = re.sub(" ", "", Options["--searchPattern"])
    if not SearchPattern:
        MiscUtil.PrintError("No value specified using \"-s, --searchPattern\" option." % (SearchPattern))
    if re.match("^auto$", SearchPattern, re.I):
        if re.match("^(AttachAtoms|ReplaceAtoms|AttachSMILES|TandemScan)$", Mode, re.I):
            SearchPattern = "[cH]"
        else:
            SearchPattern = None
    else:
        if re.match("^RxnSMARTS$", Mode, re.I):
            MiscUtil.PrintError("The specification of SMARTS search pattern, %s,  using \"-s, --searchPattern\" option is not allowed during, %s, value of \"-m, --mode\"." % (SearchPattern, Mode))

    if SearchPattern is not None:
        SearchPatternMol = Chem.MolFromSmarts(SearchPattern)
        if SearchPatternMol is None:
            MiscUtil.PrintError("Failed to create search pattern molecule. The torsion SMILES/SMARTS pattern, \"%s\", specified using \"-s, --searchPattern\" option is not valid." % (SearchPattern))
    OptionsInfo["SearchPattern"] = SearchPattern
        
    # Process and validate target pattern...
    TargetPatternTandemScanInfo = None
    TargetPattern = re.sub(" ", "", Options["--targetPattern"])
    if not TargetPattern:
        MiscUtil.PrintError("No value specified using \"-t, --targetPattern\" option." % (TargetPattern))
        
    if re.match("^auto$", TargetPattern, re.I):
        # Setup default target pattern...
        if re.match("^AttachAtoms$", Mode, re.I):
            TargetPattern = "F"
        elif re.match("^ReplaceAtoms$", Mode, re.I):
            TargetPattern = "N"
        elif re.match("^AttachSMILES$", Mode, re.I):
            TargetPattern = "C(F)(F)(F)"
        elif re.match("^RxnSMARTS$", Mode, re.I):
            TargetPattern = "[cH:1]>>[N:1]"
        elif re.match("^TandemScan$", Mode, re.I):
            TargetPattern = "ReplaceAtoms,N,AttachAtoms,F"
        else:
            MiscUtil.PrintError("Failed to setup default target pattern for mode, %s, specified using \"-m, --mode\"" % (Mode))
    
    if re.match("^(AttachAtoms|ReplaceAtoms)$", Mode, re.I):
        if not RDKitUtil.IsValidElementSymbol(TargetPattern):
            MiscUtil.PrintError("The target pattern, \"%s\", specified using \"-t, --targetPattern\" must be a valid atomic symbol during, %s, value of \"-m, --mode\" option." % (TargetPattern, Mode))
    elif re.match("^AttachSMILES$", Mode, re.I):
        ValidateSMILESTargetPattern(TargetPattern, Mode)
    elif re.match("^RxnSMARTS$", Mode, re.I):
        try:
            Rxn = AllChem.ReactionFromSmarts(TargetPattern)
        except Exception as ErrMsg:
            MiscUtil.PrintError("The target pattern, \"%s\", specified using \"-t, --targetPattern\" must be a valid RxnSMARTS  during, %s, value of \"-m, --mode\" option.\nErrMsg: %s" % (TargetPattern, Mode, ErrMsg))
    elif re.match("^TandemScan$", Mode, re.I):
        TargetPatternTandemScanInfo = {}
        TargetPatternTandemScanInfo["Names"] = []
        TargetPatternTandemScanInfo["Pattern"] = {}
        
        TargetPatternWords = TargetPattern.split(",")
        if len(TargetPatternWords) % 2:
            MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified using \"-t --targetPattern\" option must be an even number during, %s, value of \"-m, --mode\" option." % (len(TargetPatternWords), Mode))
        
        ValidOperationNames = ["AttachAtoms", "ReplaceAtoms", "AttachSMILES"]
        CanonicalOperationNamesMap = {}
        for OperationName in ValidOperationNames:
            CanonicalOperationNamesMap[OperationName.lower()] = OperationName
        
        # Validate PAS operation name and value pairs...
        for Index in range(0, len(TargetPatternWords), 2):
            OperationName = TargetPatternWords[Index]
            OperationPattern = TargetPatternWords[Index + 1]
        
            CanonicalOperationName = OperationName.lower()
            if  not CanonicalOperationName in CanonicalOperationNamesMap:
                MiscUtil.PrintError("The operation name, %s, specified using \"-t, --targetPattern\" is not a valid during, %s, value of \"-m, --mode\" option. Supported operation  names: %s" % (OperationName, Mode, " ".join(ValidOperationNames)))
            
            if not OperationPattern:
                MiscUtil.PrintError("Empty target pattern specified using \"-t, --targetPattern\" must be a valid atomic symbol for operation, %s, during, %s, value of \"-m, --mode\" option." % (OperationName, Mode))
            
            if re.match("^(AttachAtoms|ReplaceAtoms)$", OperationName, re.I):
                if not RDKitUtil.IsValidElementSymbol(OperationPattern):
                    MiscUtil.PrintError("The target pattern, \"%s\", specified using \"-t, --targetPattern\" must be a valid atomic symbol for operation, %s, during, %s, value of \"-m, --mode\" option." % (OperationPattern, OperationName, Mode))
            elif re.match("^AttachSMILES$", OperationName, re.I):
                ValidateSMILESTargetPattern(OperationPattern, Mode, OperationName)
            
            # Track info...
            Name = CanonicalOperationNamesMap[CanonicalOperationName]
            TargetPatternTandemScanInfo["Names"].append(Name)
            TargetPatternTandemScanInfo["Pattern"][Name] = OperationPattern
    
    OptionsInfo["TargetPattern"] = TargetPattern
    OptionsInfo["TargetPatternTandemScanInfo"] = TargetPatternTandemScanInfo
    
def ValidateSMILESTargetPattern(TargetPattern, Mode, OperationName = None):
    """Validate SMILES pattern."""

    Mol = Chem.MolFromSmiles(TargetPattern)
    if Mol is None:
        if OperationName is None:
            MiscUtil.PrintError("The target pattern, \"%s\", specified using \"-t, --targetPattern\" must be a valid SMILES during, %s, value of \"-m, --mode\" option." % (TargetPattern, Mode))
        else:
            MiscUtil.PrintError("The target pattern, \"%s\", specified using \"-t, --targetPattern\" must be a valid SMILES for operation, %s, during, %s, value of \"-m, --mode\" option." % (TargetPattern, OperationName, Mode))
    
    if RDKitUtil.AreAtomMapNumbersPresentInMol(Mol):
        if OperationName is None:
            MiscUtil.PrintError("Atom map numbers are not allowed in the target pattern, \"%s\", specified using \"-t, --targetPattern\"  during, %s, value of \"-m, --mode\" option." % (TargetPattern, Mode))
        else:
            MiscUtil.PrintError("Atom map number are not allowed in  target pattern, \"%s\", specified using \"-t, --targetPattern\" must be a valid SMILES for operation, %s, during, %s, value of \"-m, --mode\" option." % (TargetPattern, OperationName, Mode))
            
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
    
    MiscUtil.ValidateOptionIntegerValue("-c, --comboSearchAtoms", Options["--comboSearchAtoms"], {">": 0})
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi txt csv tsv")

    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi txt csv tsv")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "ReplaceAtoms AttachAtoms AttachSMILES RxnSMARTS TandemScan")
    
    MiscUtil.ValidateOptionTextValue("--mp", Options["--mp"], "yes no")
    MiscUtil.ValidateOptionTextValue("-p, --permuteTandemScan", Options["--permuteTandemScan"], "yes no")
    MiscUtil.ValidateOptionTextValue("-q, --quiet", Options["--quiet"], "yes no")
    
    MiscUtil.ValidateOptionTextValue(" -u, --uniqueAnalogues", Options["--uniqueAnalogues"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitPerformPositionalAnalogueScan.py - Positional analogue scanning.

Usage:
    RDKitPerformPositionalAnalogueScan.py [--comboSearchAtoms <number>] [--infileParams <Name,Value,...>]
                                          [--mode <ReplaceAtoms, AttachAtoms, AttachSMILES, RxnSMARTS, TandemScan>]
                                          [--mp <yes or no>] [--mpParams <Name,Value,...>] [ --nameSuffix <text>] [ --outfileParams <Name,Value,...> ]
                                          [--overwrite] [--permuteTandemScan <yes or no>] [--quiet <yes or no>] [--searchPattern <SMARTSPattern>]
                                          [--targetPattern <ElementSymbol, SMILES, RxnSMARTS,...>] [--uniqueAnalogues <yes or no> ]
                                          [-w <dir>] -i <infile>  -o <outfile> 
    RDKitPerformPositionalAnalogueScan.py -h | --help | -e | --examples

Description:
    Perform Positional Analogue Scanning (PAS) to generate analogues of molecules
    by applying chemical transformations to molecules [ Ref 147-148 ]. The chemical
    transformations are defined using SMARTS, element symbols, SMILES, and 
    RxnSMARTS. Four different types of chemical transformations are available for
    for generating analogues of molecules: replace atoms, attach atoms, attach SMILES,
    and RxnSMARTS. Tandem positional analogue scanning may be performed by the
    concurrent application of multiple chemical transformations.
    
    A SMARTS search pattern identifies atoms in molecules for attachment or
    replacement points during positional analogue scanning. It may retrieve multiple
    substructure matches in a molecule. The first matched atom in each substructure
    match comprises a set of attachment or replacement points.
    
    A target pattern encompasses information regarding element symbol, SMILES, and
    reaction SMARTS for replacing and attaching atoms, attaching SMILES, and applying
    reaction SMARTS. In addition, multiple concurrent chemical transformations may
    be specified during tandem positional analogue scanning. 
    
    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi,
    .csv, .tsv .txt)

    The supported output file formats are: SD (.sdf, .sd), SMILES (.smi, .csv,
    .tsv .txt)

Options:
    -c, --comboSearchAtoms <number>  [default: 1]
        Number of concurrent search pattern match atoms to use as attachment or
        replacement points during positional analogue scanning in 'AttachAtoms', 
        'AttachSMILES', and 'ReplaceAtoms' modes. This value is ignored during 
        'RxnSMARTS' and 'TandemScan' values of  '-m, --mode' option.
    -e, --examples
        Print examples.
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
    -m, --mode <ReplaceAtoms, AttachAtoms,...>  [default: ReplaceAtoms]
        Type of operations to perform during positional analogue scanning. The
        supported values, along with a brief description, are shown below:
            
            Mode           Description
            ReplaceAtoms   Replace atoms with new atoms
            AttachAtoms    Attach new atoms to atoms
            AttachSMILES   Attach SMILES to atoms
            RxnSMARTS      Run reaction SMARTS
            TandemScan     Perform tandem scan by combining ReplaceAtoms,
                AttachAtoms, and AttachSMILES                    
             
        The chemical transformations of input molecules is dependent on the
        values of '-s, --searchPattern' and  '-t, --targetPattern' options. For
        example, nitrogen-walk or nitrogen scan is performed by '[cH]' and 'N'
        values for '-s, --searchPattern' and  '-t, --targetPattern' options.
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
    -n, --nameSuffix <text>  [default: Analogue]
        Name suffix for generating molecule names of analogues. Format of analogue
        names: <MolName>_<NameSuffix>_<MolNum>
    -o, --outfile <outfile>
        Output file name.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: compute2DCoords,auto,kekulize,no
            SMILES: kekulize,no,smilesDelimiter,space, smilesIsomeric,yes,
                smilesTitleLine,yes,smilesMolName,yes,smilesMolProps,no
            
    --overwrite
        Overwrite existing files.
    -p, --permuteTandemScan <yes or no>  [default: yes]
        Permute atom positions matched by SMARTS search pattern in a molecule
        to generate all possible analogues during tandem positional analogue
        scanning.
        
        This option is only valid for 'TandemScan' value of '-m, --mode' option.
    -q, --quiet <yes or no>  [default: no]
        Use quiet mode. The warning and information messages will not be printed.
    -s, --searchPattern <SMARTSPattern>  [default: auto]
        SMARTS search pattern identifying atoms in molecules for attachment or
        replacement points during positional analogue scanning. The SMARTS search
        pattern may retrieve multiple substructure matches in a molecule. The first
        matched atom in each substructure match comprises a set of attachment or
        replacement points.
        
        The default values, dependent on the value of '-m, --mode' option, are
        shown below:
            
            Mode            Default   Description
            ReplaceAtoms    [cH]      Aromatic carbon  
            AttachAtoms     [cH]      Aromatic carbon
            AttachSMILES    [cH]      Aromatic carbon
            RxnSMARTS       None      Not applicable
            TandemScan      [cH]      Aromatic carbon
            
        This option is ignored during 'RxnSMARTS' value of '-m, --mode' option.
    -t, --targetPattern <ElementSymbol, SMILES, RxnSMARTS...>  [default: auto]
        Target pattern for performing chemical transformations during positional
        analogue scanning. These values are used in conjunction with the value of
        '-s, --searchPattern' to generate appropriate analogues.
        
        The default values, dependent on the values of '-m, --mode' and 
        '-s, --searchPattern' options, are shown below:
            
            Mode            Default        Description
            ReplaceAtoms    N              Element symbol for nitrogen
            AttachAtoms     F              Element symbol for fluorine
            AttachSMILES    C(F)(F)(F)     SMILES for CF3
            RxnSMARTS       [cH:1]>>[N:1]  Replace aromatic carbon by nitrogen
            TandemScan      ReplaceAtoms,N,AttachAtoms,F  Replace and attach
            
        Multiple concurrent chemical transformations are allowed during 'TandemScan'. The
        target pattern specification for 'TandemScan' is a comma delimited list of operation
        type and target pattern. Format: OperationType,TargetPattern,...
        
        The supported operation types and target pattern are shown below:
      
            ReplaceAtoms,<ElementSymbol>
            AttachAtoms,<ElementSymbol>
            AttachSMILES,<SMILES> 
            
        For example:
      
            ReplaceAtoms,N,AttachAtoms,F
            ReplaceAtoms,N,AttachAtoms,F,AttachSMILES,C(F)(F)(F) 
             
        The number of chemical transformations  in 'TandemScan' must be less than or
        equal to the total number atoms matched by SMARTS search pattern in a molecule.
        Otherwise, it is not possible to perform a 'TandemScan'. The matched atom positions
        may be optionally permuted to generate all possible analogues during  positional
        analogue scanning using '-p, --permuteTandemScan' option.
    -u, --uniqueAnalogues <yes or no>  [default: yes]
        Keep only unique analogues of a molecule corresponding to unique SMILES
        strings. The duplicate SMILES string may be generated during PAS due to
        symmetric replacement or attachment points in molecules.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To perform a nitrogen-walk or nitrogen scan by replacing aromatic carbons by
    nitrogens in molecules in a SMILES file and write out a SMILES file, type:

        % RDKitPerformPositionalAnalogueScan.py  -i Sample.smi -o SampleOut.smi

    To run the first example by explicity specifying search and target patterns for
    replacing aromatic carbons by nitogens in molecules in a SD file and write out
    a SD file, type:

        % RDKitPerformPositionalAnalogueScan.py  -m ReplaceAtoms -s "[cH]"
          -t "N" -i Sample.sdf -o SampleOut.sdf

    To run the first example in multiprocessing mode on all available CPUs
    without loading all data into memory and write out a SD file, type:

        % RDKitPerformPositionalAnalogueScan.py  -i Sample.smi -o SampleOut.sdf
          --mp yes
     
    To run the previous example in multiprocessing mode on all available CPUs
    by loading all data into memory and write out a SD file, type:

        % RDKitPerformPositionalAnalogueScan.py  -i Sample.smi -o SampleOut.sdf
          --mp yes --mpParams "inputDataMode,InMemory"
     
    To run the previous example in multiprocessing mode on specific number of
    CPUs and chunk size without loading all data into memory and write out a SD file,
    type:

        % RDKitPerformPositionalAnalogueScan.py  -i Sample.smi -o SampleOut.sdf
          --mpParams "inputDataMode,Lazy,numProcesses,4,chunkSize,8"
     
    To perform positional analogue scanning by simultaneously attaching fluorines
    to two aromatic carbons in molecules in a SMILES file and write out a SD file,
    type:

        % RDKitPerformPositionalAnalogueScan.py  -m AttachAtoms -s "[cH]"
          -t "F" -c 2 -i Sample.smi -o SampleOut.sdf

    To perform positional analogue scanning by attaching SMILES for CF3 to aromatic
    carbons in molecules in a SMILES file and write out a SD file, type:

        % RDKitPerformPositionalAnalogueScan.py  -m AttachSMILES -s "[cH]"
          -t "C(F)(F)(F)" -i Sample.smi -o SampleOut.sdf

    To perform a nitrogen-walk or nitrogen scan by using reaction SMARTS to replace
    aromatic carbons by nitrogens in molecules in a SMILES file and write out a  SMILES
    file,  type:

        % RDKitPerformPositionalAnalogueScan.py  -m RxnSMARTS
          -t "[cH:1]>>[N:1]" -i Sample.smi -o SampleOut.smi

    To perform a tandem positional analogue scan by concurrently applying multiple
    chemical transformations to aromatic carbons, permute all matched search
    atom positions during analogue generation, and write out a SD file, type:

        % RDKitPerformPositionalAnalogueScan.py  -m TandemScan -s "[cH]"
          -t "ReplaceAtoms,N,AttachAtoms,F,AttachSMILES,OC"
          -p yes  -i Sample.smi -o SampleOut.smi

    To perform a nitrogen-walk or nitrogen scan by replacing aromatic carbons by
    nitrogens in molecules in a SMILES CSV fileS, MILES strings in column 1, name
    in column 2, and write out a SD file, type:

        % RDKitPerformPositionalAnalogueScan.py  -m ReplaceAtoms -s "[cH]"
          -t "N" --infileParams "smilesDelimiter,comma, smilesTitleLine,yes,
          smilesColumn,1,smilesNameColumn,2"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitEnumerateCompoundLibrary.py,
    RDKitPerformTorsionScan.py

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
