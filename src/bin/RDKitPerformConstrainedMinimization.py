#!/usr/bin/env python
#
# File: RDKitPerformConstrainedMinimization.py
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
    from rdkit.Chem import rdFMCS
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
    PerformConstrainedMinimization()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformConstrainedMinimization():
    """Perform constrained minimization."""
    
    # Read and validate reference molecule...
    RefMol = RetrieveReferenceMolecule()
    
    # Setup a molecule reader for input file...
    MiscUtil.PrintInfo("\nProcessing file %s..." % OptionsInfo["Infile"])
    OptionsInfo["InfileParams"]["AllowEmptyMols"] = True
    Mols  = RDKitUtil.ReadMolecules(OptionsInfo["Infile"], **OptionsInfo["InfileParams"])

    # Set up a molecule writer...
    Writer = RDKitUtil.MoleculesWriter(OptionsInfo["Outfile"], **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % OptionsInfo["Outfile"])
    MiscUtil.PrintInfo("Generating file %s..." % OptionsInfo["Outfile"])

    MolCount, ValidMolCount, CoreScaffoldMissingCount, MinimizationFailedCount = ProcessMolecules(RefMol, Mols, Writer)

    if Writer is not None:
        Writer.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of molecules with missing core scaffold: %d" % CoreScaffoldMissingCount)
    MiscUtil.PrintInfo("Number of molecules failed during conformation generation or minimization: %d" % MinimizationFailedCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount + CoreScaffoldMissingCount + MinimizationFailedCount))

def ProcessMolecules(RefMol, Mols, Writer):
    """Process and minimize molecules. """
    
    if OptionsInfo["MPMode"]:
        return ProcessMoleculesUsingMultipleProcesses(RefMol, Mols, Writer)
    else:
        return ProcessMoleculesUsingSingleProcess(RefMol, Mols, Writer)

def ProcessMoleculesUsingSingleProcess(RefMol, Mols, Writer):
    """Process and minimize molecules using a single process."""

    (MolCount, ValidMolCount, CoreScaffoldMissingCount, MinimizationFailedCount) = [0] * 4
    
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

        # Setup a reference molecule core containing common scaffold atoms...
        RefMolCore = SetupCoreScaffold(RefMol, Mol, MolCount)
        if RefMolCore is None:
            CoreScaffoldMissingCount += 1
            continue
            
        Mol, CalcStatus, Energy, ScaffoldEmbedRMSD = ConstrainAndMinimizeMolecule(Mol, RefMolCore, MolCount)
        
        if not CalcStatus:
            MinimizationFailedCount += 1
            continue
        
        WriteMolecule(Writer, Mol, Energy, ScaffoldEmbedRMSD)

    return (MolCount, ValidMolCount, CoreScaffoldMissingCount, MinimizationFailedCount)

def ProcessMoleculesUsingMultipleProcesses(RefMol, Mols, Writer):
    """Process and minimize molecules using multiprocessing."""

    MPParams = OptionsInfo["MPParams"]
    
    # Setup data for initializing a worker process...
    MiscUtil.PrintInfo("Encoding options info and reference molecule...")
    
    OptionsInfo["EncodedRefMol"] = RDKitUtil.MolToBase64EncodedMolString(RefMol)
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
    
    (MolCount, ValidMolCount, CoreScaffoldMissingCount, MinimizationFailedCount) = [0] * 4
    for Result in Results:
        MolCount += 1
        MolIndex, EncodedMol, CoreScaffoldMissingStatus, CalcStatus, Energy, ScaffoldEmbedRMSD  = Result
        
        if EncodedMol is None:
            continue
        ValidMolCount += 1

        if CoreScaffoldMissingStatus:
            CoreScaffoldMissingStatus += 1
            continue
        
        if not CalcStatus:
            MinimizationFailedCount += 1
            continue
            
        Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
        WriteMolecule(Writer, Mol, Energy, ScaffoldEmbedRMSD)
    
    return (MolCount, ValidMolCount, CoreScaffoldMissingCount, MinimizationFailedCount)
    
def InitializeWorkerProcess(*EncodedArgs):
    """Initialize data for a worker process."""
    
    global Options, OptionsInfo

    MiscUtil.PrintInfo("Starting process (PID: %s)..." % os.getpid())

    # Decode Options and OptionInfo...
    Options = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[0])
    OptionsInfo = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[1])

    # Decode RefMol...
    OptionsInfo["RefMol"] = RDKitUtil.MolFromBase64EncodedMolString(OptionsInfo["EncodedRefMol"])
    
def WorkerProcess(EncodedMolInfo):
    """Process data for a worker process."""

    MolIndex, EncodedMol = EncodedMolInfo

    CoreScaffoldMissingStatus = False
    CalcStatus = False
    Energy = None
    ScaffoldEmbedRMSD = None
    
    if EncodedMol is None:
        return [MolIndex, None, CoreScaffoldMissingStatus, CalcStatus, Energy, ScaffoldEmbedRMSD]

    RefMol = OptionsInfo["RefMol"]
    
    Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
    if RDKitUtil.IsMolEmpty(Mol):
        if not OptionsInfo["QuietMode"]:
            MolName = RDKitUtil.GetMolName(Mol, (MolIndex + 1))
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
        return [MolIndex, None, CoreScaffoldMissingStatus, CalcStatus, Energy, ScaffoldEmbedRMSD]
    
    # Setup a reference molecule core containing common scaffold atoms...
    RefMolCore = SetupCoreScaffold(RefMol, Mol, (MolIndex + 1))
    if RefMolCore is None:
        CoreScaffoldMissingStatus = True
        return [MolIndex, None, CalcStatus, CoreScaffoldMissingStatus, Energy, ScaffoldEmbedRMSD]
    
    Mol, CalcStatus, Energy, ScaffoldEmbedRMSD = ConstrainAndMinimizeMolecule(Mol, RefMolCore, (MolIndex + 1))

    return [MolIndex, RDKitUtil.MolToBase64EncodedMolString(Mol, PropertyPickleFlags = Chem.PropertyPickleOptions.MolProps | Chem.PropertyPickleOptions.PrivateProps), CoreScaffoldMissingStatus, CalcStatus, Energy, ScaffoldEmbedRMSD]

def RetrieveReferenceMolecule():
    """Retrieve and validate reference molecule """
    
    RefFile = OptionsInfo["RefFile"]
    
    MiscUtil.PrintInfo("\nProcessing file %s..." % (RefFile))
    OptionsInfo["InfileParams"]["AllowEmptyMols"] = False
    ValidRefMols, RefMolCount, ValidRefMolCount  = RDKitUtil.ReadAndValidateMolecules(RefFile, **OptionsInfo["InfileParams"])
    
    if ValidRefMolCount == 0:
        MiscUtil.PrintError("The reference file, %s, contains no valid molecules." % RefFile)
    elif ValidRefMolCount > 1:
        MiscUtil.PrintWarning("The reference file, %s, contains, %d, valid molecules. Using first molecule as the reference molecule..." % (RefFile, ValidRefMolCount))
    
    RefMol = ValidRefMols[0]

    if OptionsInfo["UseScaffoldSMARTS"]:
        ScaffoldPatternMol = Chem.MolFromSmarts(OptionsInfo["ScaffoldSMARTS"])
        if ScaffoldPatternMol is None:
            MiscUtil.PrintError("Failed to create scaffold pattern molecule. The scaffold SMARTS pattern, %s, specified using \"-s, --scaffold\" option is not valid." % (OptionsInfo["ScaffoldSMARTS"]))
        
        if not RefMol.HasSubstructMatch(ScaffoldPatternMol):
            MiscUtil.PrintError("The scaffold SMARTS pattern, %s, specified using \"-s, --scaffold\" option, is missing in the first valid reference molecule." % (OptionsInfo["ScaffoldSMARTS"]))
            
    return RefMol

def SetupCoreScaffold(RefMol, Mol, MolCount):
    """Setup a reference molecule core containing common scaffold atoms between
    a pair of molecules."""

    if OptionsInfo["UseScaffoldMCS"]:
        return SetupCoreScaffoldByMCS(RefMol, Mol, MolCount)
    elif OptionsInfo["UseScaffoldSMARTS"]:
        return SetupCoreScaffoldBySMARTS(RefMol, Mol, MolCount)
    else:
        MiscUtil.PrintError("The  value, %s, specified for  \"-s, --scaffold\" option is not supported." % (OptionsInfo["Scaffold"]))
        
def SetupCoreScaffoldByMCS(RefMol, Mol, MolCount):
    """Setup a reference molecule core containing common scaffold atoms between
    a pair of molecules using MCS."""
    
    MCSParams = OptionsInfo["MCSParams"]
    Mols = [RefMol, Mol]

    MCSResultObject = rdFMCS.FindMCS(Mols, maximizeBonds = MCSParams["MaximizeBonds"], threshold = MCSParams["Threshold"], timeout = MCSParams["TimeOut"], verbose = MCSParams["Verbose"], matchValences = MCSParams["MatchValences"], ringMatchesRingOnly = MCSParams["RingMatchesRingOnly"], completeRingsOnly = MCSParams["CompleteRingsOnly"], matchChiralTag = MCSParams["MatchChiralTag"], atomCompare = MCSParams["AtomCompare"], bondCompare = MCSParams["BondCompare"], seedSmarts = MCSParams["SeedSMARTS"]) 
    
    if MCSResultObject.canceled:
        if not OptionsInfo["QuietMode"]:
            MiscUtil.PrintWarning("MCS failed to identify a common core scaffold between reference moecule and input molecule %s. Specify a different set of parameters using \"-m, --mcsParams\" option and try again." % (RDKitUtil.GetMolName(Mol, MolCount)))
        return None
    
    CoreNumAtoms = MCSResultObject.numAtoms
    CoreNumBonds = MCSResultObject.numBonds
    
    SMARTSCore = MCSResultObject.smartsString
    
    if not len(SMARTSCore):
        if not OptionsInfo["QuietMode"]:
            MiscUtil.PrintWarning("MCS failed to identify a common core scaffold between reference moecule and input molecule %s. Specify a different set of parameters using \"-m, --mcsParams\" option and try again." % (RDKitUtil.GetMolName(Mol, MolCount)))
        return None
        
    if CoreNumAtoms < MCSParams["MinNumAtoms"]:
        if not OptionsInfo["QuietMode"]:
            MiscUtil.PrintWarning("Number of atoms, %d, in core scaffold identified by MCS is less than, %d, as specified by \"minNumAtoms\" parameter in  \"-m, --mcsParams\" option." % (CoreNumAtoms, MCSParams["MinNumAtoms"]))
        return None
    
    if CoreNumBonds < MCSParams["MinNumBonds"]:
        if not OptionsInfo["QuietMode"]:
            MiscUtil.PrintWarning("Number of bonds, %d, in core scaffold identified by MCS is less than, %d, as specified by \"minNumBonds\" parameter in  \"-m, --mcsParams\" option." % (CoreNumBonds, MCSParams["MinNumBonds"]))
        return None

    return GenerateCoreMol(RefMol, SMARTSCore)
    
def SetupCoreScaffoldBySMARTS(RefMol, Mol, MolCount):
    """Setup a reference molecule core containing common scaffold atoms between
    a pair of molecules using specified SMARTS."""
    
    if OptionsInfo["ScaffoldPatternMol"] is None:
        OptionsInfo["ScaffoldPatternMol"] = Chem.MolFromSmarts(OptionsInfo["ScaffoldSMARTS"])
        
    if not Mol.HasSubstructMatch(OptionsInfo["ScaffoldPatternMol"]):
        if not OptionsInfo["QuietMode"]:
            MiscUtil.PrintWarning("The scaffold SMARTS pattern, %s, specified using \"-s, --scaffold\" option is missing in input molecule,  %s." % (OptionsInfo["ScaffoldSMARTS"], RDKitUtil.GetMolName(Mol, MolCount)))
        return None

    return GenerateCoreMol(RefMol, OptionsInfo["ScaffoldSMARTS"])

def GenerateCoreMol(RefMol, SMARTSCore):
    """Generate core molecule for embedding. """

    # Create a molecule corresponding to core atoms...
    SMARTSCoreMol = Chem.MolFromSmarts(SMARTSCore)

    # Setup a ref molecule containing core atoms with dummy atoms as
    # attachment points for atoms around the core atoms...
    Core = AllChem.ReplaceSidechains(Chem.RemoveHs(RefMol), SMARTSCoreMol)

    # Delete any substructures containing dummy atoms..
    RefMolCore = AllChem.DeleteSubstructs(Core, Chem.MolFromSmiles('*'))
    RefMolCore.UpdatePropertyCache()
    
    return RefMolCore

def ConstrainAndMinimizeMolecule(Mol, RefMolCore, MolNum = None):
    "Constrain and Minimize molecule."

    if  OptionsInfo["AddHydrogens"]:
        Mol = Chem.AddHs(Mol, addCoords = True)

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
        return (None, False, None, None)
        
    MaxConfs = OptionsInfo["MaxConfs"]
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
            return (None, False, None, None)
        
        EnergyStatus, Energy = GetEnergy(MolConf)
        
        if not EnergyStatus:
            if not OptionsInfo["QuietMode"]:
                MolName = RDKitUtil.GetMolName(Mol, MolNum)
                MiscUtil.PrintWarning("Failed to retrieve calculated energy for conformation number %d of molecule %s. Try again after removing any salts or cleaing up the molecule...\n" % (ConfID, MolName))
            return (None, False, None, None)
        
        CalcEnergyMap[ConfID] = Energy
        MolConfsMap[ConfID] = MolConf

    SortedConfIDs = sorted(ConfIDs, key = lambda ConfID: CalcEnergyMap[ConfID])
    MinEnergyConfID = SortedConfIDs[0]
    
    MinEnergy = "%.2f" % CalcEnergyMap[MinEnergyConfID]  if OptionsInfo["EnergyOut"] else None
    MinEnergyMolConf = MolConfsMap[MinEnergyConfID]
    
    ScaffoldEmbedRMSD = "%.4f" % float(MinEnergyMolConf.GetProp('EmbedRMS')) if OptionsInfo["ScaffoldRMSDOut"] else None
    MinEnergyMolConf.ClearProp('EmbedRMS')
    
    if  OptionsInfo["RemoveHydrogens"]:
        MinEnergyMolConf = Chem.RemoveHs(MinEnergyMolConf)
        
    return (MinEnergyMolConf, True, MinEnergy, ScaffoldEmbedRMSD)

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
    
def WriteMolecule(Writer, Mol, Energy = None, ScaffoldEmbedRMSD = None, ConfID = None,):
    """Write molecule. """

    if ScaffoldEmbedRMSD is not None:
        Mol.SetProp("CoreScaffoldEmbedRMSD", ScaffoldEmbedRMSD)
            
    if Energy is not None:
        Mol.SetProp(OptionsInfo["EnergyLabel"], Energy)
            
    if ConfID is None:
        Writer.write(Mol)
    else:
        Writer.write(Mol, confId = ConfID)
    
def ProcessMCSParameters():
    """Set up and process MCS parameters."""

    SetupMCSParameters()
    ProcessSpecifiedMCSParameters()

def SetupMCSParameters():
    """Set up default MCS parameters."""
    
    OptionsInfo["MCSParams"] = {"MaximizeBonds": True, "Threshold": 0.9, "TimeOut": 3600, "Verbose": False, "MatchValences": True, "MatchChiralTag": False, "RingMatchesRingOnly": True, "CompleteRingsOnly": True, "AtomCompare": rdFMCS.AtomCompare.CompareElements, "BondCompare": rdFMCS.BondCompare.CompareOrder, "SeedSMARTS": "", "MinNumAtoms": 1, "MinNumBonds": 0}
    
def ProcessSpecifiedMCSParameters():
    """Process specified MCS parameters."""

    if re.match("^auto$", OptionsInfo["SpecifiedMCSParams"], re.I):
        # Nothing to process...
        return
    
    # Parse specified parameters...
    MCSParams = re.sub(" ", "", OptionsInfo["SpecifiedMCSParams"])
    if not MCSParams:
        MiscUtil.PrintError("No valid parameter name and value pairs specified using \"-m, --mcsParams\" option.")

    MCSParamsWords = MCSParams.split(",")
    if len(MCSParamsWords) % 2:
        MiscUtil.PrintError("The number of comma delimited paramater names and values, %d, specified using \"-m, --mcsParams\" option must be an even number." % (len(MCSParamsWords)))
    
    # Setup  canonical parameter names...
    ValidParamNames = []
    CanonicalParamNamesMap = {}
    for ParamName in sorted(OptionsInfo["MCSParams"]):
        ValidParamNames.append(ParamName)
        CanonicalParamNamesMap[ParamName.lower()] = ParamName

    # Validate and set paramater names and value...
    for Index in range(0, len(MCSParamsWords), 2):
        Name = MCSParamsWords[Index]
        Value = MCSParamsWords[Index + 1]

        CanonicalName = Name.lower()
        if  not CanonicalName in CanonicalParamNamesMap:
            MiscUtil.PrintError("The parameter name, %s, specified using \"-m, --mcsParams\" option is not a valid name. Supported parameter names: %s" % (Name,  " ".join(ValidParamNames)))

        ParamName = CanonicalParamNamesMap[CanonicalName]
        if re.match("^Threshold$", ParamName, re.I):
            Value = float(Value)
            if Value <= 0.0 or Value > 1.0 :
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: > 0 and <= 1.0" % (Value, Name))
            ParamValue = Value
        elif re.match("^Timeout$", ParamName, re.I):
            Value = int(Value)
            if Value <= 0:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: > 0" % (Value, Name))
            ParamValue = Value
        elif re.match("^MinNumAtoms$", ParamName, re.I):
            Value = int(Value)
            if Value < 1:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: >= 1" % (Value, Name))
            ParamValue = Value
        elif re.match("^MinNumBonds$", ParamName, re.I):
            Value = int(Value)
            if Value < 0:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: >=0 " % (Value, Name))
            ParamValue = Value
        elif re.match("^AtomCompare$", ParamName, re.I):
            if re.match("^CompareAny$", Value, re.I):
                ParamValue = rdFMCS.AtomCompare.CompareAny
            elif re.match("^CompareElements$", Value, re.I):
                ParamValue = Chem.rdFMCS.AtomCompare.CompareElements
            elif re.match("^CompareIsotopes$", Value, re.I):
                ParamValue = Chem.rdFMCS.AtomCompare.CompareIsotopes
            else:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: CompareAny CompareElements CompareIsotopes" % (Value, Name))
        elif re.match("^BondCompare$", ParamName, re.I):
            if re.match("^CompareAny$", Value, re.I):
                ParamValue = Chem.rdFMCS.BondCompare.CompareAny
            elif re.match("^CompareOrder$", Value, re.I):
                ParamValue = rdFMCS.BondCompare.CompareOrder
            elif re.match("^CompareOrderExact$", Value, re.I):
                ParamValue = rdFMCS.BondCompare.CompareOrderExact
            else:
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: CompareAny CompareOrder CompareOrderExact" % (Value, Name))
        elif re.match("^SeedSMARTS$", ParamName, re.I):
            if not len(Value):
                MiscUtil.PrintError("The parameter value specified using \"-m, --mcsParams\" option  for parameter, %s, is empty. " % (Name))
            ParamValue = Value
        else:
            if not re.match("^(Yes|No|True|False)$", Value, re.I):
                MiscUtil.PrintError("The parameter value, %s, specified using \"-m, --mcsParams\" option  for parameter, %s, is not a valid value. Supported values: Yes No True False" % (Value, Name))
            ParamValue = False
            if re.match("^(Yes|True)$", Value, re.I):
                ParamValue = True
        
        # Set value...
        OptionsInfo["MCSParams"][ParamName] = ParamValue

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["RefFile"] = Options["--reffile"]

    OptionsInfo["Scaffold"] = Options["--scaffold"]
    if re.match("^auto$", Options["--scaffold"], re.I):
        UseScaffoldMCS = True
        UseScaffoldSMARTS = False
        ScaffoldSMARTS = None
    else:
        UseScaffoldMCS = False
        UseScaffoldSMARTS = True
        ScaffoldSMARTS = OptionsInfo["Scaffold"]
    
    OptionsInfo["UseScaffoldMCS"] = UseScaffoldMCS
    OptionsInfo["UseScaffoldSMARTS"] = UseScaffoldSMARTS
    OptionsInfo["ScaffoldSMARTS"] = ScaffoldSMARTS
    OptionsInfo["ScaffoldPatternMol"] = None

    OptionsInfo["SpecifiedMCSParams"] = Options["--mcsParams"]
    ProcessMCSParameters()
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    OptionsInfo["AddHydrogens"] = True if re.match("^yes$", Options["--addHydrogens"], re.I) else False
    
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
    
    OptionsInfo["ScaffoldRMSDOut"] = True if re.match("^yes$", Options["--scaffoldRMSDOut"], re.I) else False
    
    OptionsInfo["EnergyOut"] = True if re.match("^yes$", Options["--energyOut"], re.I) else False
    if UseMMFF:
        OptionsInfo["EnergyLabel"] = "%s_Energy" % MMFFVariant
    else:
        OptionsInfo["EnergyLabel"] = "%s_Energy" % ForceField
    
    OptionsInfo["EnforceChirality"] = True if re.match("^yes$", Options["--enforceChirality"], re.I) else False
    
    OptionsInfo["MaxConfs"] = int(Options["--maxConfs"])
    
    OptionsInfo["MPMode"] = True if re.match("^yes$", Options["--mp"], re.I) else False
    OptionsInfo["MPParams"] = MiscUtil.ProcessOptionMultiprocessingParameters("--mpParams", Options["--mpParams"])
    
    OptionsInfo["QuietMode"] = True if re.match("^yes$", Options["--quiet"], re.I) else False
    
    OptionsInfo["RemoveHydrogens"] = True if re.match("^yes$", Options["--removeHydrogens"], re.I) else False
    OptionsInfo["UseTethers"] = True if re.match("^yes$", Options["--useTethers"], re.I) else False

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
    
    MiscUtil.ValidateOptionTextValue("--scaffoldRMSDOut", Options["--scaffoldRMSDOut"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--energyOut", Options["--energyOut"], "yes no")
    MiscUtil.ValidateOptionTextValue("--enforceChirality ", Options["--enforceChirality"], "yes no")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi txt csv tsv")

    MiscUtil.ValidateOptionFilePath("-r, --reffile", Options["--reffile"])
    MiscUtil.ValidateOptionFileExt("-r, --reffile", Options["--reffile"], "sdf sd mol")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
        
    MiscUtil.ValidateOptionIntegerValue("--maxConfs", Options["--maxConfs"], {">": 0})
    
    MiscUtil.ValidateOptionTextValue("--mp", Options["--mp"], "yes no")
    MiscUtil.ValidateOptionTextValue("-q, --quiet", Options["--quiet"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("-r, --removeHydrogens", Options["--removeHydrogens"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("-u, --useTethers", Options["--useTethers"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitPerformConstrainedMinimization.py - Perform constrained minimization

Usage:
    RDKitPerformConstrainedMinimization.py [--addHydrogens <yes or no>] [--conformerGenerator <SDG, ETDG, KDG, ETKDG>]
                                           [--forceField <UFF, or MMFF>] [--forceFieldMMFFVariant <MMFF94 or MMFF94s>]
                                           [--energyOut  <yes or no>] [--enforceChirality <yes or no>] [--infileParams <Name,Value,...>]
                                           [--maxConfs <number>]  [--mcsParams <Name,Value,...>] [--mp <yes or no>] [--mpParams <Name.Value,...>]
                                           [ --outfileParams <Name,Value,...> ] [--overwrite] [--quiet <yes or no>] [ --removeHydrogens <yes or no>]
                                           [--scaffold <auto or SMARTS>]  [--scaffoldRMSDOut  <yes or no>] [--useTethers  <yes or no>] 
                                           [-w <dir>] -i <infile> -r <reffile> -o <outfile> 
    RDKitPerformConstrainedMinimization.py -h | --help | -e | --examples

Description:
    Generate 3D structures for molecules by performing a constrained energy minimization
    against a reference molecule. An initial set of 3D conformers are generated for the
    input molecules using distance geometry. A common core scaffold, corresponding to
    a Maximum Common Substructure (MCS) or an explicit SMARTS pattern,  is identified
    between a pair of input and reference molecules. The core scaffold atoms in input
    molecules are aligned against the same atoms in the reference molecule. The energy
    of aligned structures are minimized using the forcefield to generate the final 3D structures.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd)
    .csv, .tsv .txt)

    The supported output file formats are: SD (.sdf, .sd)

Options:
    -a, --addHydrogens <yes or no>  [default: yes]
        Add hydrogens before minimization.
    -c, --conformerGenerator <SDG, ETDG, KDG, ETKDG>  [default: ETKDG]
        Conformation generation methodology for generating initial 3D coordinates
        for molecules in input file. A common core scaffold is identified between a
        a pair of input and reference molecules. The atoms in common core scaffold 
        of input molecules are aligned against the reference molecule followed by
        energy minimization to generate final 3D structure.
        
        Possible values: Standard Distance Geometry, (SDG), Experimental Torsion-angle
        preference with Distance Geometry (ETDG), basic Knowledge-terms with Distance
        Geometry (KDG),  and Experimental Torsion-angle preference along with basic
        Knowledge-terms with Distance Geometry (ETKDG) [Ref 129] .
    -f, --forceField <UFF, MMFF>  [default: MMFF]
        Forcefield method to use for  constrained energy minimization. Possible values:
        Universal Force Field (UFF) [ Ref 81 ] or Merck Molecular Mechanics Force
        Field [ Ref 83-87 ] .
    --forceFieldMMFFVariant <MMFF94 or MMFF94s>  [default: MMFF94]
        Variant of MMFF forcefield to use for energy minimization.
    --energyOut <yes or no>  [default: No]
        Write out energy values.
    --enforceChirality <yes or no>  [default: Yes]
        Enforce chirality for defined chiral centers.
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
    --maxConfs <number>  [default: 250]
        Maximum number of conformations to generate for each molecule by conformation
        generation methodology for initial 3D coordinates. A constrained minimization is
        performed using the specified forcefield and the lowest energy conformation is written
        to the output file.
    --mcsParams <Name,Value,...>  [default: auto]
        Parameter values to use for identifying a maximum common substructure
        (MCS) in between a pair of reference and input molecules.In general, it is a
        comma delimited list of parameter name and value pairs. The supported
        parameter names along with their default values are shown below:
            
            atomCompare,CompareElements,bondCompare,CompareOrder,
            maximizeBonds,yes,matchValences,yes,matchChiralTag,no,
            minNumAtoms,1,minNumBonds,0,ringMatchesRingOnly,yes,
            completeRingsOnly,yes,threshold,1.0,timeOut,3600,seedSMARTS,none
            
        Possible values for atomCompare: CompareAny, CompareElements,
        CompareIsotopes. Possible values for bondCompare: CompareAny,
        CompareOrder, CompareOrderExact.
        
        A brief description of MCS parameters taken from RDKit documentation is
        as follows:
            
            atomCompare - Controls match between two atoms
            bondCompare - Controls match between two bonds
            maximizeBonds - Maximize number of bonds instead of atoms
            matchValences - Include atom valences in the MCS match
            matchChiralTag - Include atom chirality in the MCS match
            minNumAtoms - Minimum number of atoms in the MCS match
            minNumBonds - Minimum number of bonds in the MCS match
            ringMatchesRingOnly - Ring bonds only match other ring bonds
            completeRingsOnly - Partial rings not allowed during the match
            threshold - Fraction of the dataset that must contain the MCS
            seedSMARTS - SMARTS string as the seed of the MCS
            timeout - Timeout for the MCS calculation in seconds
            
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
    -r, --reffile <reffile>
        Reference input file name containing a 3D reference molecule. A common
        core scaffold must be present in a pair of an input and reference molecules.
        Otherwise, no constrained minimization is performed on the input molecule.
    --removeHydrogens <yes or no>  [default: Yes]
        Remove hydrogens after minimization.
    -s, --scaffold <auto or SMARTS>  [default: auto]
        Common core scaffold between a pair of input and reference molecules used for
        constrained minimization of molecules in input file. Possible values: Auto or a
        valid SMARTS pattern. The common core scaffold is automatically detected
        corresponding to the Maximum Common Substructure (MCS) between a pair of
        reference and input molecules. A valid SMARTS pattern may be optionally specified
        for the common core scaffold.
    --scaffoldRMSDOut <yes or no>  [default: No]
        Write out RMSD value for common core alignment between a pair of input and
        reference molecules.
    -u, --useTethers <yes or no>  [default: yes]
        Use tethers to optimize the final conformation by applying a series of extra forces
        to align matching atoms to the positions of the core atoms. Otherwise, use simple
        distance constraints during the optimization.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To perform constrained energy minimization for molecules in a SMILES file against
    a reference 3D molecule in a SD file using a common core scaffold between pairs of
    input and reference molecules identified using MCS, generating up to 250 conformations
    using ETKDG methodology followed by MMFF forcefield minimization, and write out
    a SD file containing minimum energy structure corresponding to each constrained
    molecule, type:

        % RDKitPerformConstrainedMinimization.py  -i SampleSeriesD3R.smi
          -r SampleSeriesRef3D.sdf  -o SampleOut.sdf

    To rerun the first example in a quiet mode and write out a SD file, type:

        % RDKitPerformConstrainedMinimization.py  -q yes -i SampleSeriesD3R.smi
          -r SampleSeriesRef3D.sdf  -o SampleOut.sdf

    To run the first example in multiprocessing mode on all available CPUs
    without loading all data into memory and write out a SD file, type:

        % RDKitPerformConstrainedMinimization.py  --mp yes
          -i SampleSeriesD3R.smi -r SampleSeriesRef3D.sdf  -o SampleOut.sdf

    To rerun the first example in multiprocessing mode on all available CPUs
    by loading all data into memory and write out a SD file, type:

        % RDKitPerformConstrainedMinimization.py  --mp yes --mpParams
          "inputDataMode,InMemory" -i SampleSeriesD3R.smi
          -r SampleSeriesRef3D.sdf  -o SampleOut.sdf

    To rerun the first example using an explicit SMARTS string for a common core
    scaffold and write out a SD file, type:

        % RDKitPerformConstrainedMinimization.py  --scaffold
          "c1c(C(N(C(c2cc(-c3nc(N)ncc3)cn2))))cccc1" -i SampleSeriesD3R.smi -r
          SampleSeriesRef3D.sdf -o SampleOut.sdf 

    To rerun the first example using molecules in a CSV SMILES file, SMILES
    strings in column 1, name in column2, and write out a SD file, type:

        % RDKitPerformConstrainedMinimization.py  --infileParams "smilesDelimiter,
          comma,smilesTitleLine,yes,smilesColumn,1,smilesNameColumn,2"
          -i SampleSeriesD3R.csv -r SampleSeriesRef3D.sdf  -o SampleOut.sdf

    To perform constrained energy minimization for molecules in a SD file against
    a reference 3D molecule in a SD file using a common core scaffold between pairs of
    input and reference molecules identified using MCS, generating up to 50 conformations
    using SDG methodology followed by UFF forcefield minimization, and write out
    a SD file containing minimum energy structure along with energy and embed RMS values
    corresponding to each constrained molecule, type:

        % RDKitPerformConstrainedMinimization.py  --maxConfs 50  -c SDG -f UFF
          --scaffoldRMSDOut yes --energyOut yes -i SampleSeriesD3R.sdf
          -r SampleSeriesRef3D.sdf  -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitCalculateRMSD.py, RDKitCalculateMolecularDescriptors.py, RDKitCompareMoleculeShapes.py,
    RDKitConvertFileFormat.py, RDKitGenerateConstrainedConformers.py, RDKitPerformMinimization.py

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
