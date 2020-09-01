#!/bin/env python
#
# File: RDKitSearchFunctionalGroups.py
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
    from rdkit.Chem import FunctionalGroups
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

FunctionalGroupsMap = {}

def main():
    """Start execution of the script"""
    
    MiscUtil.PrintInfo("\n%s (RDK v%s; %s): Starting...\n" % (ScriptName, rdBase.rdkitVersion, time.asctime()))
    
    (WallClockTime, ProcessorTime) = MiscUtil.GetWallClockAndProcessorTime()
    
    # Retrieve command line arguments and options...
    RetrieveOptions()
    
    # Process and validate command line arguments and options...
    ProcessOptions()
    
    # Perform actions required by the script...
    PerformFunctionalGroupsSearch()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformFunctionalGroupsSearch():
    """Retrieve functional groups information and perform search."""

    # Process functional groups info...
    ProcessFunctionalGroupsInfo()
    
    # Setup pattern mols for functional group SMARTS...
    GroupsPatternMols = SetupFunctionalGroupsSMARTSPatterns()
    
    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % OptionsInfo["Infile"])
    Mols  = RDKitUtil.ReadMolecules(OptionsInfo["Infile"], **OptionsInfo["InfileParams"])
    
    # Set up  molecule writers...
    Writer,  GroupOutfilesWriters = SetupMoleculeWriters()

    MolCount, ValidMolCount, RemainingMolCount,  GroupsPatternMatchCountList = ProcessMolecules(Mols, GroupsPatternMols, Writer, GroupOutfilesWriters)
    
    if Writer is not None:
        Writer.close()
    for GroupOutfileWriter in GroupOutfilesWriters:
        GroupOutfileWriter.close()
        
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

    MiscUtil.PrintInfo("\nTotal number of molecules matched against specified match criteria: %d" % RemainingMolCount)
    
    MiscUtil.PrintInfo("\nNumber of molecuels matched against individual functional groups:")
    MiscUtil.PrintInfo("FunctionalGroupName,MatchCount")
    
    for GroupIndex in range(0, len(OptionsInfo["SpecifiedFunctionalGroups"])):
        GroupName = OptionsInfo["SpecifiedFunctionalGroups"][GroupIndex]
        if OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"][GroupIndex]:
            GroupName = '!' + GroupName
        GroupMatchCount = GroupsPatternMatchCountList[GroupIndex]
        MiscUtil.PrintInfo("%s,%d" % (GroupName, GroupMatchCount))

def ProcessMolecules(Mols, GroupsPatternMols, Writer, GroupOutfilesWriters):
    """Process and search molecules. """

    if OptionsInfo["MPMode"]:
        return ProcessMoleculesUsingMultipleProcesses(Mols, GroupsPatternMols, Writer, GroupOutfilesWriters)
    else:
        return ProcessMoleculesUsingSingleProcess(Mols, GroupsPatternMols, Writer, GroupOutfilesWriters)
    
def ProcessMoleculesUsingSingleProcess(Mols, GroupsPatternMols, Writer, GroupOutfilesWriters):
    """Process and search molecules using a single process."""

    MiscUtil.PrintInfo("\nSearching functional groups...")
    
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    CombineMatchResults = OptionsInfo["CombineMatchResults"]
    SetSMILESMolProps = OptionsInfo["OutfileParams"]["SetSMILESMolProps"]
    
    GroupsPatternsMatchCountList = [0] * len(OptionsInfo["SpecifiedFunctionalGroups"])
    (MolCount, ValidMolCount, RemainingMolCount) = [0] * 3
    
    FirstMol = True
    for Mol in Mols:
        MolCount += 1
        
        if Mol is None:
            continue
        
        if RDKitUtil.IsMolEmpty(Mol):
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
            continue
        
        ValidMolCount += 1
        if FirstMol:
            FirstMol = False
            if SetSMILESMolProps:
                if Writer is not None:
                    RDKitUtil.SetWriterMolProps(Writer, Mol)
                for GroupOutfileWriter in GroupOutfilesWriters:
                    if GroupOutfileWriter is not None:
                        RDKitUtil.SetWriterMolProps(GroupOutfileWriter, Mol)
        
        # Match molecule against functional group patterns...
        MolMatched, GroupsPatternMatchStatusList = MatchMolecule(Mol, GroupsPatternMols)

        # Update functional group match count...
        for GroupIndex, MatchStatus in enumerate(GroupsPatternMatchStatusList):
            if MatchStatus:
                GroupsPatternsMatchCountList[GroupIndex] += 1
        
        if not MolMatched:
            continue
        
        RemainingMolCount += 1
        WriteMolecule(Writer, GroupOutfilesWriters, Mol,  Compute2DCoords, CombineMatchResults, GroupsPatternMatchStatusList)
    
    return (MolCount, ValidMolCount, RemainingMolCount,  GroupsPatternsMatchCountList)

def ProcessMoleculesUsingMultipleProcesses(Mols, GroupsPatternMols, Writer, GroupOutfilesWriters):
    """Process and search molecules using multiprocessing."""

    MiscUtil.PrintInfo("\nSearching functional groups  using multiprocessing...")
    
    MPParams = OptionsInfo["MPParams"]
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    CombineMatchResults = OptionsInfo["CombineMatchResults"]
    SetSMILESMolProps = OptionsInfo["OutfileParams"]["SetSMILESMolProps"]
    
    # Setup data for initializing a worker process...
    MiscUtil.PrintInfo("Encoding options info and functional groups pattern molecules...")
    OptionsInfo["EncodedGroupPatternMols"] = [RDKitUtil.MolToBase64EncodedMolString(PatternMol) for PatternMol in GroupsPatternMols]
    InitializeWorkerProcessArgs = (MiscUtil.ObjectToBase64EncodedString(Options), MiscUtil.ObjectToBase64EncodedString(OptionsInfo), MiscUtil.ObjectToBase64EncodedString(FunctionalGroupsMap))

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
    
    GroupsPatternsMatchCountList = [0] * len(OptionsInfo["SpecifiedFunctionalGroups"])
    
    (MolCount, ValidMolCount, RemainingMolCount) = [0] * 3
    
    FirstMol = True
    for Result in Results:
        MolCount += 1
        MolIndex, EncodedMol, MolMatched, GroupsPatternMatchStatusList = Result
        
        if EncodedMol is None:
            continue
        ValidMolCount += 1
        
        Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
        
        if FirstMol:
            FirstMol = False
            if SetSMILESMolProps:
                if Writer is not None:
                    RDKitUtil.SetWriterMolProps(Writer, Mol)
                for GroupOutfileWriter in GroupOutfilesWriters:
                    if GroupOutfileWriter is not None:
                        RDKitUtil.SetWriterMolProps(GroupOutfileWriter, Mol)
        
        # Update functional group match count...
        for GroupIndex, MatchStatus in enumerate(GroupsPatternMatchStatusList):
            if MatchStatus:
                GroupsPatternsMatchCountList[GroupIndex] += 1
        
        if not MolMatched:
            continue
        
        RemainingMolCount += 1
        WriteMolecule(Writer, GroupOutfilesWriters, Mol,  Compute2DCoords, CombineMatchResults, GroupsPatternMatchStatusList)
    
    return (MolCount, ValidMolCount, RemainingMolCount,  GroupsPatternsMatchCountList)

def InitializeWorkerProcess(*EncodedArgs):
    """Initialize data for a worker process."""
    
    global Options, OptionsInfo, FunctionalGroupsMap

    MiscUtil.PrintInfo("Starting process (PID: %s)..." % os.getpid())
    
    # Decode Options and OptionInfo...
    Options = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[0])
    OptionsInfo = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[1])
    FunctionalGroupsMap = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[2])

    # Decode ChEMBLPatternMols...
    OptionsInfo["GroupPatternMols"] = [RDKitUtil.MolFromBase64EncodedMolString(EncodedMol) for EncodedMol in OptionsInfo["EncodedGroupPatternMols"]]
    
def WorkerProcess(EncodedMolInfo):
    """Process data for a worker process."""

    MolIndex, EncodedMol = EncodedMolInfo

    if EncodedMol is None:
        return [MolIndex, None, False, None]
    
    Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
    if RDKitUtil.IsMolEmpty(Mol):
        MolName = RDKitUtil.GetMolName(Mol, (MolIndex + 1))
        MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
        return [MolIndex, None, False, None]
        
    # Match molecule against functional group patterns...
    MolMatched, GroupsPatternMatchStatusList = MatchMolecule(Mol, OptionsInfo["GroupPatternMols"])
    
    return [MolIndex, EncodedMol, MolMatched, GroupsPatternMatchStatusList]

def MatchMolecule(Mol, GroupsPatternMols):
    """Search for functional groups in a molecule."""

    GroupsPatternMatchStatusList = []
    
    # Match pattern mols...
    for GroupIndex in range(0, len(OptionsInfo["SpecifiedFunctionalGroups"])):
        Status = DoesPatternMolMatch(GroupsPatternMols[GroupIndex], Mol, OptionsInfo["UseChirality"], OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"][GroupIndex])
        GroupsPatternMatchStatusList.append(Status)
        
    # Match mol against all specified criteria...
    MolMatched = DoesMolMeetSpecifiedMatchCriteria(GroupsPatternMatchStatusList, OptionsInfo["CombineMatchResults"], OptionsInfo["AndCombineOperatorMode"])

    return (MolMatched, GroupsPatternMatchStatusList)

def DoesMolMeetSpecifiedMatchCriteria(GroupsPatternMolsMatchStatus,  CombineMatchResults, AndCombineOperatorMode):
    """Match molecule using specified match criteia."""

    if CombineMatchResults and AndCombineOperatorMode:
        # Must match all specified SMARTS
        Status = True
        for MatchStatus in GroupsPatternMolsMatchStatus:
            if not MatchStatus:
                Status = False
                break
    else:
        # One match is enough...
        Status = False
        for MatchStatus in GroupsPatternMolsMatchStatus:
            if MatchStatus:
                Status = True
                break
    
    return Status
    
def WriteMolecule(Writer, GroupOutfilesWriters, Mol, Compute2DCoords, CombineMatchResults, GroupsPatternMatchStatusList):
    """Write out molecule."""
    
    if OptionsInfo["CountMode"]:
        return
    
    if Compute2DCoords:
        AllChem.Compute2DCoords(Mol)
    
    if CombineMatchResults:
        Writer.write(Mol)
    else:
        for GroupIndex in range(0, len(GroupsPatternMatchStatusList)):
            if GroupsPatternMatchStatusList[GroupIndex]:
                GroupOutfilesWriters[GroupIndex].write(Mol)
    
def SetupMoleculeWriters():
    """Set up molecule writers for output files."""

    Writer = None
    GroupOutfilesWriters = []
    
    if OptionsInfo["CountMode"]:
        return (Writer, GroupOutfilesWriters)
    
    Outfile = OptionsInfo["Outfile"]
    CombineMatchResults = OptionsInfo["CombineMatchResults"]
    GroupsOutfiles = OptionsInfo["SpecifiedFunctionalGroupsOutfiles"]
    
    if CombineMatchResults:
        Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
        if Writer is None:
            MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
        MiscUtil.PrintInfo("Generating file %s..." % Outfile)
    else:
        for GroupOutfile in GroupsOutfiles:
            GroupOutfileWriter = RDKitUtil.MoleculesWriter(GroupOutfile, **OptionsInfo["OutfileParams"])
            if GroupOutfileWriter is None:
                MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Writer)
            GroupOutfilesWriters.append(GroupOutfileWriter)
        
        GroupsCount = len(GroupsOutfiles)
        if GroupsCount > 4:
            MiscUtil.PrintInfo("Generating %d output files with the following file name format: %s<GroupName>.%s" % (GroupsCount, OptionsInfo["OutfileBasename"], OptionsInfo["OutfileExt"]))
        else:
            Delmiter = ', '
            OutfileNames = Delmiter.join(GroupsOutfiles)
            MiscUtil.PrintInfo("Generating %d output files: %s..." % (GroupsCount, OutfileNames))

    return (Writer, GroupOutfilesWriters)
    
def DoesPatternMolMatch(PatternMol, Mol, UseChirality, NegateMatch):
    """Perform a substructure match for the presence of pattern molecule in a molecule."""

    MolMatched = Mol.HasSubstructMatch(PatternMol, useChirality = UseChirality)
    if NegateMatch:
        if MolMatched:
            MolMatched = False
        else:
            MolMatched = True
    
    return MolMatched
    
def ProcessFunctionalGroupsInfo():
    """Process functional groups information."""
    
    RetrieveFunctionalGroupsInfo()
    ProcessSpecifiedFunctionalGroups()
    
    SetupFunctionalGroupsOutputFileNames()

def ProcessSpecifiedFunctionalGroups():
    """Process and validate specified functional groups"""
    
    OptionsInfo["SpecifiedFunctionalGroups"] = []
    OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"] = []
    
    if re.match("^All$", OptionsInfo["FunctionalGroups"], re.I):
        OptionsInfo["SpecifiedFunctionalGroups"] = FunctionalGroupsMap['Names']
        OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"] = [False] * len(OptionsInfo["SpecifiedFunctionalGroups"])
        return
    
    # Set up a map of valid group names for checking specified group names...
    CanonicalGroupNameMap = {}
    for GroupName in FunctionalGroupsMap['Names']:
        CanonicalGroupNameMap[GroupName.lower()] = GroupName
        
    # Parse and validate specified names...
    GroupNames = re.sub(" ", "", OptionsInfo["FunctionalGroups"])
    if not GroupNames:
        MiscUtil.PrintError("No functional group name specified for \"-f, --functionalGroups\" option")
    
    SpecifiedFunctionalGroups = []
    SpecifiedNegateMatchStatus = []
    
    for GroupName in GroupNames.split(","):
        CanonicalGroupName = GroupName.lower()
        NegateMatchStatus = False
        if re.match("^!", CanonicalGroupName, re.I):
            NegateMatchStatus = True
            CanonicalGroupName = re.sub("^!", "", CanonicalGroupName)
        if CanonicalGroupName in CanonicalGroupNameMap:
            SpecifiedFunctionalGroups.append(CanonicalGroupNameMap[CanonicalGroupName])
            SpecifiedNegateMatchStatus.append(NegateMatchStatus)
        else:
            MiscUtil.PrintWarning("The functional group name, %s, specified using \"-f, --functionalGroups\" option is not a valid name." % (GroupName))

    if not len(SpecifiedFunctionalGroups):
        MiscUtil.PrintError("No valid functional group names specified for \"-f, --functionalGroups\" option")
        
    OptionsInfo["SpecifiedFunctionalGroups"] = SpecifiedFunctionalGroups
    OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"] = SpecifiedNegateMatchStatus

def SetupFunctionalGroupsSMARTSPatterns():
    """Setup SMARTS patterns for specified functional groups."""

    OptionsInfo["SpecifiedFunctionalGroupsSMARTSPatterns"] = []
    FunctionalGroupsPatternMols = []
    
    for Name in OptionsInfo["SpecifiedFunctionalGroups"]:
        SMARTSPattern = FunctionalGroupsMap['SMARTSPattern'][Name]
        PatternMol = Chem.MolFromSmarts(SMARTSPattern)
        if PatternMol is None:
            MiscUtil.PrintError("Failed to parse SMARTS pattern, %s, for function group, %s" % (SMARTSPattern, Name))
        
        OptionsInfo["SpecifiedFunctionalGroupsSMARTSPatterns"].append(SMARTSPattern)
        FunctionalGroupsPatternMols.append(PatternMol)
    
    return FunctionalGroupsPatternMols

def SetupFunctionalGroupsOutputFileNames():
    """Setup output file names for specified functional group names."""

    OptionsInfo["SpecifiedFunctionalGroupsOutfiles"] = []
    
    if OptionsInfo["CountMode"]:
        # No need of any output file...
        return
    
    if OptionsInfo["CombineMatchResults"]:
        # No need of output files for specified functional groups...
        return
    
    OutfileBasename = OptionsInfo["OutfileBasename"]
    OutfileExt = OptionsInfo["OutfileExt"]
    SpecifiedFunctionalGroupsOutfiles = []

    GroupsCount = len(OptionsInfo["SpecifiedFunctionalGroups"])
    for GroupIndex in range(0, GroupsCount):
        GroupName = OptionsInfo["SpecifiedFunctionalGroups"][GroupIndex]
        if OptionsInfo["SpecifiedFunctionalGroupsNegateMatch"][GroupIndex]:
            GroupName = "Not" + GroupName
        GroupName = re.sub("\.", "", GroupName)
        
        GroupOutfile = "%s%s.%s" % (OutfileBasename, GroupName, OutfileExt)
        SpecifiedFunctionalGroupsOutfiles.append(GroupOutfile)

    OptionsInfo["SpecifiedFunctionalGroupsOutfiles"] = SpecifiedFunctionalGroupsOutfiles
    
def RetrieveFunctionalGroupsInfo():
    """Retrieve functional groups information"""

    MiscUtil.PrintInfo("\nRetrieving data from default RDKit functional groups hierarchy file Functional_Group_Hierarchy.txt...")
    
    FunctionalGroupNamesFile = OptionsInfo["GroupNamesFile"]
    FunctionalGroupsNodes = FunctionalGroups.BuildFuncGroupHierarchy(FunctionalGroupNamesFile)

    FunctionalGroupsMap['Names'] = []
    FunctionalGroupsMap['SMARTSPattern'] = {}
    
    RetrieveDataFromFunctionalGroupsHierarchy(FunctionalGroupsNodes)

    if not len(FunctionalGroupsMap['Names']):
        MiscUtil.PrintError("Failed to retrieve any functional group names and SMARTS patterns...")
        
    MiscUtil.PrintInfo("Total number of functional groups present functional group hierarchy: %d" % (len(FunctionalGroupsMap['Names'])))
    
def RetrieveDataFromFunctionalGroupsHierarchy(FGNodes):
    """Retrieve functional groups data by recursively visiting functional group nodes."""
    
    for FGNode in FGNodes:
        Name = FGNode.label
        SMARTSPattern = FGNode.smarts

        if Name in FunctionalGroupsMap['SMARTSPattern']:
            MiscUtil.PrintWarning("Ignoring duplicate functional group name: %s..." % Name)
        else:
            FunctionalGroupsMap['Names'].append(Name)
            FunctionalGroupsMap['SMARTSPattern'][Name] = SMARTSPattern

        RetrieveDataFromFunctionalGroupsHierarchy(FGNode.children)

def ListFunctionalGroupsInfo():
    """List functional groups information"""

    MiscUtil.PrintInfo("\nListing available functional groups names and SMARTS patterns...")
    MiscUtil.PrintInfo("\nFunctionalGroupName\tSMARTSPattern")
    
    for Name in sorted(FunctionalGroupsMap['Names']):
        SMARTSPattern = FunctionalGroupsMap['SMARTSPattern'][Name]
        MiscUtil.PrintInfo("%s\t%s" % (Name, SMARTSPattern))
    
    MiscUtil.PrintInfo("")

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["CombineMatches"] = Options["--combineMatches"]
    
    OptionsInfo["CombineMatchResults"] = True
    if re.match("^No$", Options["--combineMatches"], re.I):
        OptionsInfo["CombineMatchResults"] = False
        if Options["--outfile"]:
            FileDir, FileName, FileExt = MiscUtil.ParseFileName(Options["--outfile"])
            OptionsInfo["OutfileBasename"] = FileName
            OptionsInfo["OutfileExt"] = FileExt
    
    OptionsInfo["CombineOperator"] = Options["--combineOperator"]
    OptionsInfo["AndCombineOperatorMode"] = True
    if re.match("^or$", Options["--combineOperator"], re.I):
        OptionsInfo["AndCombineOperatorMode"] = False
    
    OptionsInfo["GroupNamesFile"] = None
    if not re.match("^auto$", Options["--groupNamesFile"], re.I):
        OptionsInfo["GroupNamesFile"] = Options["--groupNamesFile"]
        
    OptionsInfo["FunctionalGroups"] = Options["--functionalGroups"]
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    OptionsInfo["CountMode"] = False
    if re.match("^count$", Options["--mode"], re.I):
        OptionsInfo["CountMode"] = True
    
    OptionsInfo["MPMode"] = True if re.match("^yes$", Options["--mp"], re.I) else False
    OptionsInfo["MPParams"] = MiscUtil.ProcessOptionMultiprocessingParameters("--mpParams", Options["--mpParams"])
    
    OptionsInfo["UseChirality"] = False
    if re.match("^yes$", Options["--useChirality"], re.I):
        OptionsInfo["UseChirality"] = True

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
    
    # Handle listing of functional group information...
    if  Options and Options["--list"]:
        ProcessListFunctionalGroupsOption()
        sys.exit(0)

def ProcessListFunctionalGroupsOption():
    """Process list functional groups information."""

    # Validate and process dataFile option for listing functional groups information...
    OptionsInfo["GroupNamesFile"] = None
    if not re.match("^auto$", Options["--groupNamesFile"], re.I):
        MiscUtil.ValidateOptionFilePath("-g, --groupNamesFile", Options["--groupNamesFile"])
        OptionsInfo["GroupNamesFile"] = Options["--groupNamesFile"]
    
    RetrieveFunctionalGroupsInfo()
    ListFunctionalGroupsInfo()
    
def ValidateOptions():
    """Validate option values"""
    
    MiscUtil.ValidateOptionTextValue("-c, --combineMatches", Options["--combineMatches"], "yes no")
    MiscUtil.ValidateOptionTextValue("--combineOperator", Options["--combineOperator"], "and or")
    
    if not re.match("^auto$", Options["--groupNamesFile"], re.I):
        MiscUtil.ValidateOptionFilePath("-g, groupNamesFile", Options["--groupNamesFile"])
        
    if re.match("^none$", Options["--functionalGroups"], re.I):
        MiscUtil.PrintError("The name(s) of functional groups must be specified using \"-f, --functionalGroups\" option")
        
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd smi txt csv tsv")
    if Options["--outfile"]:
        MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi")
        MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
        MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "retrieve count")
    if re.match("^retrieve$", Options["--mode"], re.I):
        if not Options["--outfile"]:
            MiscUtil.PrintError("The outfile must be specified using \"-o, --outfile\" during \"retrieve\" value of \"-m, --mode\" option")
        
    MiscUtil.ValidateOptionTextValue("--mp", Options["--mp"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--useChirality", Options["--useChirality"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitSearchFunctionalGroups.py - Search for functional groups using SMARTS patterns

Usage:
    RDKitSearchFunctionalGroups.py  [--combineMatches <yes or no>] [--combineOperator <and or or>]
                                           [--groupNamesFile <FileName or auto>] [--infileParams <Name,Value,...>]
                                           [--mode <retrieve or count>] [--mp <yes or no>] [--mpParams <Name.Value,...>]
                                           [--negate <yes or no>] [--outfileParams <Name,Value,...>] [--overwrite]
                                           [--useChirality <yes or no>] [-w <dir>] [-o <outfile>] -i <infile> -f <Name1,Name2,Name3... or All>
    RDKitSearchFunctionalGroups.py [--groupNamesFile <FileName or auto>] -l | --list
    RDKitSearchFunctionalGroups.py -h | --help | -e | --examples

Description:
    Perform a substructure search in an input file using SMARTS patterns for functional
    groups and write out the matched molecules to an output file or simply count the
    number of matches.

    The SMARTS patterns for specified functional group(s) are retrieved from file,
    Functional_Group_Hierarchy.txt, available in RDKit data directory.

    The names of valid functional groups and hierarchies  are dynamically retrieved from the
    functional groups hierarchy file and are shown below:

        AcidChloride, AcidChloride.Aromatic, AcidChloride.Aliphatic
        Alcohol, Alcohol.Aromatic, Alcohol.Aliphatic
        Aldehyde, Aldehyde.Aromatic, Aldehyde.Aliphatic
        Amine, Amine.Primary, Amine.Primary.Aromatic, Amine.Primary.Aliphatic,
        Amine.Secondary, Amine.Secondary.Aromatic, Amine.Secondary.Aliphatic
        Amine.Tertiary, Amine.Tertiary.Aromatic, Amine.Tertiary.Aliphatic
        Amine.Aromatic, Amine.Aliphatic, Amine.Cyclic
        Azide, Azide.Aromatic, Azide.Aliphatic
        BoronicAcid, BoronicAcid.Aromatic, BoronicAcid.Aliphatic
        CarboxylicAcid, CarboxylicAcid.Aromatic, CarboxylicAcid.Aliphatic,
        CarboxylicAcid.AlphaAmino
        Halogen, Halogen.Aromatic, Halogen.Aliphatic
        Halogen.NotFluorine, Halogen.NotFluorine.Aliphatic,
        Halogen.NotFluorine.Aromatic
        Halogen.Bromine, Halogen.Bromine.Aliphatic, Halogen.Bromine.Aromatic,
        Halogen.Bromine.BromoKetone
        Isocyanate, Isocyanate.Aromatic, Isocyanate.Aliphatic
        Nitro, Nitro.Aromatic, Nitro.Aliphatic,
        SulfonylChloride, SulfonylChloride.Aromatic, SulfonylChloride.Aliphatic
        TerminalAlkyne

    The supported input file formats are: SD (.sdf, .sd), SMILES (.smi, .csv, .tsv, .txt)

    The supported output file formats are: SD (.sdf, .sd), SMILES (.smi)

Options:
    -c, --combineMatches <yes or no>  [default: yes]
        Combine search results for matching SMARTS patterns of specified functional groups
        against a molecule. Possible values: yes or no.
        
        The matched molecules are written to a single output file for "yes" value. Otherwise,
        multiple output files are generated, one for each functional group. The names of  
        these files correspond to a combination of the basename of the specified output file
        and the name of the functional group.
        
        No output files are generated during "count" value of "-m, --mode" option.
    --combineOperator <and or or>  [default: and]
        Logical operator to use for combining match results corresponding to specified
        functional group names before writing out a single file. This option is ignored
        during "No" value of  "-c, --combineMatches" option.
    -e, --examples
        Print examples.
    -g, --groupNamesFile <FileName or auto>  [default: auto]
        Specify a file name containing data for functional groups hierarchy or use functional
        group hierarchy file, Functional_Group_Hierarchy.txt, available in RDKit data directory.
        
        RDKit data format: Name<tab>Smarts<tab>Label<tab>RemovalReaction (optional)
        
        The format of data in local functional group hierarchy must match format of the
        data in functional group file available in RDKit data directory.
    -f, --functionalGroups <Name1,Name2,Name3... or All>  [default: none]
        Functional group names for performing substructure SMARTS search. Possible values:
        Comma delimited list of valid functional group names or All. The current set of valid
        functional group names are listed in the description section.
        
        The match results for multiple functional group names are combined using 'and'
        operator before writing them out to single file. No merging of match results takes
        place during generation of individual result files corresponding to fictional group
        names.
        
        The functional group name may be started with an exclamation mark to negate
        the match result for that fictional group.
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
    -l, --list
        List functional groups information without performing any search.
    -m, --mode <retrieve or count>  [default: retrieve]
        Specify whether to retrieve and write out matched molecules to an output
        file or simply count the number of matches.
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
            
            SD: compute2DCoords,auto,kekulize,no
            SMILES: kekulize,no,smilesDelimiter,space, smilesIsomeric,yes,
                smilesTitleLine,yes,smilesMolName,yes,smilesMolProps,no
            
        Default value for compute2DCoords: yes for SMILES input file; no for all other
        file types.
    --overwrite
        Overwrite existing files.
    -u, --useChirality <yes or no>  [default: no]
        Use stereochemistry information for SMARTS search.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To list names of all available functional groups along with their SMARTS
    patterns, type:

        % RDKitSearchFunctionalGroups.py -l

    To retrieve molecules containing amine functional group and write out a
    SMILES file, type: 

        % RDKitSearchFunctionalGroups.py -f Amine -i Sample.smi -o SampleOut.smi

    To retrieve molecules containing amine functional group, perform search in
    multiprocessing mode on all  available CPUs without loading all data into
    memory, and write out a SMILES file, type: 

        % RDKitSearchFunctionalGroups.py --mp yes -f Amine -i Sample.smi
          -o SampleOut.smi

    To retrieve molecules containing amine functional group, perform search in
    multiprocessing mode on all  available CPUs by loading all data into memory,
    and write out a SMILES file, type: 

        % RDKitSearchFunctionalGroups.py --mp yes --mpParams "inputDataMode,
          InMemory" -f Amine -i Sample.smi -o SampleOut.smi

    To retrieve molecules containing amine functional group, perform search in
    multiprocessing mode on specific number of CPUs and chunksize without loading
    all data into memory, and write out a SMILES file, type: 

        % RDKitSearchFunctionalGroups.py --mp yes --mpParams "inputDataMode,
          lazy,numProcesses,4,chunkSize,8" -f Amine -i Sample.smi -o
          SampleOut.smi

    To retrieve molecules containing amine functional group but not halogens and carboxylic
    acid functional groups and write out a SMILES file, type: 

        % RDKitSearchFunctionalGroups.py -f 'Amine,!Halogen,!CarboxylicAcid'
          -i Sample.smi -o SampleOut.smi

    To retrieve molecules containing amine, halogens or carboxylic  acid functional groups
    and write out a SMILES file, type: 

        % RDKitSearchFunctionalGroups.py -f 'Amine,Halogen,CarboxylicAcid'
          --combineOperator or -i Sample.smi -o SampleOut.smi

    To retrieve molecules containing amine and carboxylic acid functional groups defined in
    a local functional groups hierarchy file and write out individual SD files for each
    funcitonal group, type: 

        % RDKitSearchFunctionalGroups.py -f 'Amine,CarboxylicAcid' -i Sample.sdf 
          -g Custom_Functional_Group_Hierarchy.txt --combineMatches No -o SampleOut.sdf

    To count number of all functional groups in molecules without writing out an output
    files, type:

        % RDKitSearchFunctionalGroups.py -m count -f All --combineMatches no -i Sample.smi

    To retrieve molecule not containing aromatic alcohol and aromatic halogen functional
    group along with the use of chirality during substructure search and write out individual
    SMILES files for each functional group, type: 

        % RDKitSearchFunctionalGroups.py --combineMatches no -u yes
           -f '!Alcohol.Aromatic,!Halogen.Aromatic' -i Sample.smi -o SampleOut.smi

    To retrieve molecule containing amine functional group from a CSV SMILES file,
    SMILES strings in column 1, name in column 2, and write out a SD file, type: 

        % RDKitSearchFunctionalGroups.py -f Amine --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitFilterPAINS.py, RDKitSearchSMARTS.py

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
