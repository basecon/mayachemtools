#!/usr/bin/env python
#
# File: RDKitCalculateMolecularDescriptors.py
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
    from rdkit.Chem import rdMolDescriptors
    from rdkit.Chem import Descriptors
    from rdkit.Chem import Descriptors3D
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import RDKit module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your RDKit environment and try again.\n\n")
    sys.exit(1)

# RDKit dependency imports...
import numpy

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

DescriptorNamesMap = {}

def main():
    """Start execution of the script"""
    
    MiscUtil.PrintInfo("\n%s (RDK v%s; %s): Starting...\n" % (ScriptName, rdBase.rdkitVersion, time.asctime()))
    
    (WallClockTime, ProcessorTime) = MiscUtil.GetWallClockAndProcessorTime()
    
    # Retrieve command line arguments and options...
    RetrieveOptions()
    
    # Process and validate command line arguments and options...
    ProcessOptions()
    
    # Perform actions required by the script...
    CalculateMolecularDescriptors()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def CalculateMolecularDescriptors():
    """Calculate molecular descriptors."""

    ProcessMolecularDescriptorsInfo()
    PerformCalculations()

def ProcessMolecularDescriptorsInfo():
    """Process descriptors information."""

    RetrieveMolecularDescriptorsInfo()
    ProcessSpecifiedDescriptorNames()

def PerformCalculations():
    """Calculate descriptors for a specified list of descriptors."""

    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % OptionsInfo["Infile"])
    Mols = RDKitUtil.ReadMolecules(OptionsInfo["Infile"], **OptionsInfo["InfileParams"])
    
    # Setup a writer...
    Writer = SetupMoleculeWriter()
    
    # Process molecules...
    MolCount, ValidMolCount = ProcessMolecules(Mols, Writer)

    if Writer is not None:
        Writer.close()
        
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

def ProcessMolecules(Mols, Writer):
    """Process and filter molecules. """
    
    if OptionsInfo["MPMode"]:
        return ProcessMoleculesUsingMultipleProcesses(Mols, Writer)
    else:
        return ProcessMoleculesUsingSingleProcess(Mols, Writer)

def ProcessMoleculesUsingSingleProcess(Mols, Writer):
    """Process molecules and calculate descriptors using a single process."""

    DescriptorsCount = len(OptionsInfo["SpecifiedDescriptorNames"])
    MiscUtil.PrintInfo("\nCalculating %d molecular %s for each molecule..." % (DescriptorsCount, ("descroptors" if DescriptorsCount > 1 else "descriptor")))
    
    (MolCount, ValidMolCount) = [0] * 2
    for MolIndex, Mol in enumerate(Mols):
        MolCount += 1
        if Mol is None:
            continue
        
        if RDKitUtil.IsMolEmpty(Mol):
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
            continue
        ValidMolCount += 1

        # Calculate and write descriptor values...
        CalculatedValues = CalculateDescriptorValues(MolIndex, Mol)
        WriteDescriptorValues(Mol, MolCount, Writer, CalculatedValues)
    
    return (MolCount, ValidMolCount)
    
def ProcessMoleculesUsingMultipleProcesses(Mols, Writer):
    """Process molecules and calculate descriptors using multiprocessing."""

    DescriptorsCount = len(OptionsInfo["SpecifiedDescriptorNames"])
    MiscUtil.PrintInfo("\nCalculating %d molecular %s for each molecule using multiprocessing......" % (DescriptorsCount, ("descroptors" if DescriptorsCount > 1 else "descriptor")))

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

    (MolCount, ValidMolCount) = [0] * 2
    for Result in Results:
        MolCount += 1
        MolIndex, EncodedMol, CalculatedValues = Result
        
        if EncodedMol is None:
            continue
        ValidMolCount += 1
        
        Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
        
        # Write descriptor values...
        WriteDescriptorValues(Mol, MolCount, Writer, CalculatedValues)
    
    return (MolCount, ValidMolCount)

def InitializeWorkerProcess(*EncodedArgs):
    """Initialize data for a worker process."""
    
    global Options, OptionsInfo

    MiscUtil.PrintInfo("Starting process (PID: %s)..." % os.getpid())
    
    # Decode Options and OptionInfo...
    Options = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[0])
    OptionsInfo = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[1])
    
    RetrieveMolecularDescriptorsInfo(PrintInfo = False)
    
def WorkerProcess(EncodedMolInfo):
    """Process data for a worker process."""

    MolIndex, EncodedMol = EncodedMolInfo
    
    CalculatedValues = []
    if EncodedMol is None:
        return [MolIndex, None, CalculatedValues]

    Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
    if RDKitUtil.IsMolEmpty(Mol):
        MolName = RDKitUtil.GetMolName(Mol, (MolIndex + 1))
        MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
        return [MolIndex, None, CalculatedValues]

    return [MolIndex, EncodedMol, CalculateDescriptorValues(MolIndex, Mol)]

def CalculateDescriptorValues(MolIndex, Mol):
    """Calculate descriptor values. """

    return [FormatCalculatedValue(CalculateDescriptorValue(MolIndex, Mol, Name), OptionsInfo["Precision"]) for Name in OptionsInfo["SpecifiedDescriptorNames"]]

def CalculateDescriptorValue(MolIndex, Mol, Name):
    """Calculate value for a specific descriptor along with handling any calculation failure."""

    try:
        Value = DescriptorNamesMap["ComputeFunction"][Name](Mol)
    except ValueError as ErrMsg:
        MiscUtil.PrintWarning("Failed to calculate descriptor %s for molecule %s:\n%s\n" % (Name, (MolIndex + 1), ErrMsg))
        Value = "NA"
    
    return Value
    
def WriteDescriptorValues(Mol, MolNum, Writer, CalculatedValues):
    """Write out calculated descriptor values."""

    if OptionsInfo["TextOutFileMode"]:
        LineWords = []
        if OptionsInfo["SMILESOut"]:
            SMILES = Chem.MolToSmiles(Mol, isomericSmiles = True, canonical = True)
            LineWords.append(SMILES)
            
        MolName = RDKitUtil.GetMolName(Mol, MolNum)
        LineWords.append(MolName)
        LineWords.extend(CalculatedValues)
        Line = OptionsInfo["TextOutFileDelim"].join(LineWords)
        Writer.write("%s\n" % Line)
    else:
        for Index, Name in enumerate(OptionsInfo["SpecifiedDescriptorNames"]):
            Mol.SetProp(Name, CalculatedValues[Index])

        if OptionsInfo["OutfileParams"]["Compute2DCoords"]:
            AllChem.Compute2DCoords(Mol)
        
        Writer.write(Mol)

def SetupMoleculeWriter():
    """Setup a molecule writer. """
    
    if OptionsInfo["TextOutFileMode"]:
        Writer = open(OptionsInfo["Outfile"], "w")
    else:
        Writer = RDKitUtil.MoleculesWriter(OptionsInfo["Outfile"], **OptionsInfo["OutfileParams"])
        
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % OptionsInfo["Outfile"])

    MiscUtil.PrintInfo("Generating file %s..." % OptionsInfo["Outfile"])

    # Wite out headers for a text file...
    if OptionsInfo["TextOutFileMode"]:
        LineWords = []
        if OptionsInfo["SMILESOut"]:
            LineWords.append("SMILES")
        LineWords.append("MolID")
        LineWords.extend(OptionsInfo["SpecifiedDescriptorNames"])
        Line = OptionsInfo["TextOutFileDelim"].join(LineWords)
        Writer.write("%s\n" % Line)
    
    return Writer

def FormatCalculatedValue(Value, Precision):
    """Format calculated value of descriptor based on its type."""

    if (type(Value) is float) or (type(Value) is numpy.float64):
        FormattedValue = "%.*f" % (Precision, Value)
        if not re.search("[1-9]", FormattedValue):
            FormattedValue = "0.0"
    elif type(Value) is list:
        FormattedValue = "%s" % Value
        FormattedValue = re.sub('\[|\]|,', '', FormattedValue)
    else:
        FormattedValue = "%s" % Value

    return FormattedValue

def ProcessSpecifiedDescriptorNames():
    """Process and validate specified decriptor names."""

    OptionsInfo["SpecifiedDescriptorNames"] = []

    if not re.match("^(2D|3D|All|FragmentCountOnly|Specify)$", OptionsInfo["Mode"], re.I):
        MiscUtil.PrintError("Mode value, %s, using \"-m, --mode\" option is not a valid value." % OptionsInfo["Mode"])
    
    if re.match("^2D$", OptionsInfo["Mode"], re.I):
        OptionsInfo["SpecifiedDescriptorNames"] = DescriptorNamesMap["2D"]["Names"]
        if OptionsInfo["FragmentCount"]:
            OptionsInfo["SpecifiedDescriptorNames"].extend(DescriptorNamesMap["FragmentCount"]["Names"])
        return
    elif re.match("^3D$", OptionsInfo["Mode"], re.I):
        OptionsInfo["SpecifiedDescriptorNames"] = DescriptorNamesMap["3D"]["Names"]
        return
    elif re.match("^All$", OptionsInfo["Mode"], re.I):
        OptionsInfo["SpecifiedDescriptorNames"] = DescriptorNamesMap["2D"]["Names"]
        if OptionsInfo["FragmentCount"]:
            OptionsInfo["SpecifiedDescriptorNames"].extend(DescriptorNamesMap["FragmentCount"]["Names"])
        OptionsInfo["SpecifiedDescriptorNames"].extend(DescriptorNamesMap["3D"]["Names"])
        return
    elif re.match("^FragmentCountOnly$", OptionsInfo["Mode"], re.I):
        OptionsInfo["SpecifiedDescriptorNames"] = DescriptorNamesMap["FragmentCount"]["Names"]
        return

    # Set up a canonical descriptor names map for checking specified names...
    CanonicalNameMap = {}
    for Name in  DescriptorNamesMap["ComputeFunction"]:
        CanonicalNameMap[Name.lower()] = Name
    
    # Parse and validate specified names...
    DescriptorNames = re.sub(" ", "", OptionsInfo["DescriptorNames"])
    if not DescriptorNames:
        MiscUtil.PrintError("No descriptor names specified for \"-d, --descriptorNames\" option")

    SMILESInfile = MiscUtil.CheckFileExt(Options["--infile"], "smi")
    Canonical3DNameMap = {}
    if SMILESInfile:
        for Name in DescriptorNamesMap["3D"]["Names"]:
            Canonical3DNameMap[Name.lower()] = Name
            
    SpecifiedDescriptorNames = []
    for Name in DescriptorNames.split(","):
        CanonicalName = Name.lower()
        if CanonicalName in CanonicalNameMap:
            SpecifiedDescriptorNames.append(CanonicalNameMap[CanonicalName])
        else:
            MiscUtil.PrintError("The descriptor name, %s, specified using \"-d, --descriptorNames\" option is not a valid name." % (Name))
        if SMILESInfile:
            if CanonicalName in Canonical3DNameMap:
                MiscUtil.PrintError("The 3D descriptor name, %s, specified using \"-d, --descriptorNames\" option is not a valid for SMILES input file." % (Name))
                
    if not len(SpecifiedDescriptorNames):
        MiscUtil.PrintError("No valid descriptor name specified for \"-d, --descriptorNames\" option")
    
    OptionsInfo["SpecifiedDescriptorNames"] = SpecifiedDescriptorNames

def RetrieveMolecularDescriptorsInfo(PrintInfo = True):
    """Retrieve descriptors information."""

    if PrintInfo:
        MiscUtil.PrintInfo("\nRetrieving information for avalible molecular descriptors...")
    
    # Initialze data for 2D, FragmentCount and 3D descriptors...
    DescriptorNamesMap["Types"] = ["2D", "FragmentCount", "3D"]
    DescriptorNamesMap["ComputeFunction"] = {}

    Autocorr2DExclude = OptionsInfo["Autocorr2DExclude"]
    
    for Type in DescriptorNamesMap["Types"]:
        DescriptorNamesMap[Type] = {}
        DescriptorNamesMap[Type]["Names"] = []
    
    # Setup data for 2D and FragmentCount...
    DescriptorsInfo = Descriptors.descList
    for DescriptorInfo in DescriptorsInfo:
        Name = DescriptorInfo[0]
        ComputeFunction = DescriptorInfo[1]

        Type = "2D"
        if re.match("^fr_", Name, re.I):
            Type = "FragmentCount"
        elif re.match("^Autocorr2D$", Name, re.I):
            if Autocorr2DExclude:
                continue
        
        if Name in DescriptorNamesMap["ComputeFunction"]:
            if PrintInfo:
                MiscUtil.PrintWarning("Ignoring duplicate descriptor name: %s..." % Name)
        else:
            DescriptorNamesMap[Type]["Names"].append(Name)
            DescriptorNamesMap["ComputeFunction"][Name] = ComputeFunction

    # Add new 2D decriptor name to the list...
    Type = "2D"
    if not Autocorr2DExclude:
        try:
            Name = "Autocorr2D"
            ComputeFunction =  rdMolDescriptors.CalcAUTOCORR2D
            if not Name in DescriptorNamesMap["ComputeFunction"]:
                DescriptorNamesMap[Type]["Names"].append(Name)
                DescriptorNamesMap["ComputeFunction"][Name] = ComputeFunction
        except (AttributeError):
            if PrintInfo:
                MiscUtil.PrintInfo("2D descriptor, %s, is not available in your current version of RDKit." % Name)
        
    # Set data for 3D descriptors...
    Type = "3D"
    NameToComputeFunctionMap = {'PMI1' : Descriptors3D.PMI1, 'PMI2' : Descriptors3D.PMI2, 'PMI3' : Descriptors3D.PMI3, 'NPR1' : Descriptors3D.NPR1, 'NPR2' : Descriptors3D.NPR2, 'RadiusOfGyration' : Descriptors3D.RadiusOfGyration, 'InertialShapeFactor' :  Descriptors3D.InertialShapeFactor, 'Eccentricity' : Descriptors3D.Eccentricity, 'Asphericity' : Descriptors3D.Asphericity, 'SpherocityIndex' : Descriptors3D.SpherocityIndex}
    
    for Name in NameToComputeFunctionMap:
        ComputeFunction = NameToComputeFunctionMap[Name]
        if Name in DescriptorNamesMap["ComputeFunction"]:
            if PrintInfo:
                MiscUtil.PrintWarning("Ignoring duplicate descriptor name: %s..." % Name)
        else:
            DescriptorNamesMap[Type]["Names"].append(Name)
            DescriptorNamesMap["ComputeFunction"][Name] = ComputeFunction

    # Check and add new 3D descriptors not directly available through Descriptors3D module...
    Type = "3D"
    AvailableName3DMap = {}
    for Name in ['Autocorr3D', 'RDF', 'MORSE', 'WHIM', 'GETAWAY']:
        ComputeFunction = None
        try:
            if re.match("^Autocorr3D$", Name, re.I):
                ComputeFunction =  rdMolDescriptors.CalcAUTOCORR3D
            elif re.match("^RDF$", Name, re.I):
                ComputeFunction =  rdMolDescriptors.CalcRDF
            elif re.match("^MORSE$", Name, re.I):
                ComputeFunction =  rdMolDescriptors.CalcMORSE
            elif re.match("^WHIM$", Name, re.I):
                ComputeFunction =  rdMolDescriptors.CalcWHIM
            elif re.match("^GETAWAY$", Name, re.I):
                ComputeFunction =  rdMolDescriptors.CalcGETAWAY
            else:
                ComputeFunction = None
        except (AttributeError):
            if PrintInfo:
                MiscUtil.PrintWarning("3D descriptor, %s, is not available in your current version of RDKit" % Name)
        
        if ComputeFunction is not None:
            AvailableName3DMap[Name] = ComputeFunction

    for Name in AvailableName3DMap:
        ComputeFunction = AvailableName3DMap[Name]
        if not Name in DescriptorNamesMap["ComputeFunction"]:
            DescriptorNamesMap[Type]["Names"].append(Name)
            DescriptorNamesMap["ComputeFunction"][Name] = ComputeFunction

    Count = 0
    TypesCount = []
    for Type in DescriptorNamesMap["Types"]:
        TypeCount = len(DescriptorNamesMap[Type]["Names"])
        TypesCount.append(TypeCount)
        Count += TypeCount

    if not Count:
        if PrintInfo:
            MiscUtil.PrintError("Failed to retrieve any  molecular descriptors...")

    # Sort descriptor names...
    for Type in DescriptorNamesMap["Types"]:
        DescriptorNamesMap[Type]["Names"] = sorted(DescriptorNamesMap[Type]["Names"])
    
    if PrintInfo:
        MiscUtil.PrintInfo("\nTotal number of availble molecular descriptors: %d" % Count)
    for Index in range(0, len(DescriptorNamesMap["Types"])):
        Type = DescriptorNamesMap["Types"][Index]
        TypeCount = TypesCount[Index]
        if PrintInfo:
            MiscUtil.PrintInfo("Number of %s molecular descriptors: %d" % (Type, TypeCount))
    
def ListMolecularDescriptorsInfo():
    """List descriptors information."""

    MiscUtil.PrintInfo("\nListing information for avalible molecular descriptors...")

    Delimiter = ", "
    for Type in DescriptorNamesMap["Types"]:
        Names = DescriptorNamesMap[Type]["Names"]
        MiscUtil.PrintInfo("\n%s descriptors: %s" % (Type, Delimiter.join(Names)))

    MiscUtil.PrintInfo("")

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()

    OptionsInfo["Autocorr2DExclude"] = True
    if not re.match("^Yes$", Options["--autocorr2DExclude"], re.I):
        OptionsInfo["Autocorr2DExclude"] = False
    
    OptionsInfo["FragmentCount"] = True
    if not re.match("^Yes$", Options["--fragmentCount"], re.I):
        OptionsInfo["FragmentCount"] = False
    
    OptionsInfo["DescriptorNames"] = Options["--descriptorNames"]
    OptionsInfo["Mode"] = Options["--mode"]
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    TextOutFileMode = False
    TextOutFileDelim = ""
    if MiscUtil.CheckFileExt(Options["--outfile"], "csv"):
        TextOutFileMode = True
        TextOutFileDelim = ","
    elif MiscUtil.CheckFileExt(Options["--outfile"], "tsv txt"):
        TextOutFileMode = True
        TextOutFileDelim = "\t"
    OptionsInfo["TextOutFileMode"] = TextOutFileMode
    OptionsInfo["TextOutFileDelim"] = TextOutFileDelim
    
    OptionsInfo["MPMode"] = True if re.match("^yes$", Options["--mp"], re.I) else False
    OptionsInfo["MPParams"] = MiscUtil.ProcessOptionMultiprocessingParameters("--mpParams", Options["--mpParams"])
    
    OptionsInfo["Precision"] = int(Options["--precision"])
    
    OptionsInfo["SMILESOut"] = False
    if re.match("^Yes$", Options["--smilesOut"], re.I):
        OptionsInfo["SMILESOut"] = True

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
        
    # Handle listing of descriptor information...
    if  Options["--list"]:
        ProcessListMolecularDescriptorsOption()
        sys.exit(0)

def ProcessListMolecularDescriptorsOption():
    """Process list descriptors information."""

    RetrieveMolecularDescriptorsInfo()
    ListMolecularDescriptorsInfo()

def ValidateOptions():
    """Validate option values"""

    MiscUtil.ValidateOptionTextValue("-a, --autocorr2DExclude", Options["--autocorr2DExclude"], "yes no")
    MiscUtil.ValidateOptionTextValue("-f, --fragmentCount", Options["--fragmentCount"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "2D 3D All FragmentCountOnly Specify")
    
    if re.match("^Specify$", Options["--mode"], re.I):
        if re.match("^none$", Options["--descriptorNames"], re.I):
            MiscUtil.PrintError("The name(s) of molecular descriptors must be specified using \"-d, --descriptorNames\" option during \"Specify\" value of \"-m, --mode\" option.")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd mol smi csv tsv txt")
    
    if re.match("^3D|All$", Options["--mode"], re.I):
        if MiscUtil.CheckFileExt(Options["--infile"], "smi"):
            MiscUtil.PrintError("The input SMILES file, %s, is not valid for  \"3D or All\" value of \"-m, --mode\" option.")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd csv tsv txt")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionTextValue("--mp", Options["--mp"], "yes no")
    
    MiscUtil.ValidateOptionIntegerValue("-p, --precision", Options["--precision"], {">": 0})
    MiscUtil.ValidateOptionTextValue("-s, --smilesOut", Options["--smilesOut"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitCalculateMolecularDescriptors.py - Calculate 2D/3D molecular descriptors

Usage:
    RDKitCalculateMolecularDescriptors.py [--autocorr2DExclude <yes or no>] [--fragmentCount <yes or no>]
                                          [--descriptorNames <Name1,Name2,...>] [--infileParams <Name,Value,...>]
                                          [--mode <2D, 3D, All...>] [--mp <yes or no>] [--mpParams <Name.Value,...>]
                                          [--outfileParams <Name,Value,...>] [--overwrite] [--precision <number>]
                                          [--smilesOut <yes or no>] [-w <dir>] -i <infile> -o <outfile> 
    RDKitCalculateMolecularDescriptors.py -l | --list
    RDKitCalculateMolecularDescriptors.py -h | --help | -e | --examples

Description:
    Calculate 2D/3D molecular descriptors for molecules and write them out to a SD or
    CSV/TSV text file.

    The complete list of currently available molecular descriptors may be obtained by
    using '-l, --list' option. The names of valid 2D, fragment count, and 3D molecular
    descriptors are shown below:

    2D descriptors: Autocorr2D, BalabanJ, BertzCT, Chi0, Chi1, Chi0n - Chi4n, Chi0v - Chi4v,
    EState_VSA1 - EState_VSA11, ExactMolWt, FpDensityMorgan1, FpDensityMorgan2, FpDensityMorgan3,
    FractionCSP3, HallKierAlpha, HeavyAtomCount, HeavyAtomMolWt, Ipc, Kappa1 - Kappa3,
    LabuteASA, MaxAbsEStateIndex, MaxAbsPartialCharge, MaxEStateIndex, MaxPartialCharge,
    MinAbsEStateIndex, MinAbsPartialCharge, MinEStateIndex, MinPartialCharge, MolLogP,
    MolMR, MolWt, NHOHCount, NOCount, NumAliphaticCarbocycles, NumAliphaticHeterocycles,
    NumAliphaticRings, NumAromaticCarbocycles, NumAromaticHeterocycles, NumAromaticRings,
    NumHAcceptors, NumHDonors, NumHeteroatoms, NumRadicalElectrons, NumRotatableBonds,
    NumSaturatedCarbocycles, NumSaturatedHeterocycles, NumSaturatedRings, NumValenceElectrons,
    PEOE_VSA1 - PEOE_VSA14,  RingCount, SMR_VSA1 - SMR_VSA10, SlogP_VSA1 - SlogP_VSA12,
    TPSA, VSA_EState1 - VSA_EState10, qed

    FragmentCount 2D descriptors: fr_Al_COO, fr_Al_OH, fr_Al_OH_noTert, fr_ArN, fr_Ar_COO,
    fr_Ar_N, fr_Ar_NH, fr_Ar_OH, fr_COO, fr_COO2, fr_C_O, fr_C_O_noCOO, fr_C_S, fr_HOCCN,
    fr_Imine, fr_NH0, fr_NH1, fr_NH2, fr_N_O, fr_Ndealkylation1, fr_Ndealkylation2, fr_Nhpyrrole,
    fr_SH, fr_aldehyde, fr_alkyl_carbamate, fr_alkyl_halide, fr_allylic_oxid, fr_amide, fr_amidine,
    fr_aniline, fr_aryl_methyl, fr_azide, fr_azo, fr_barbitur, fr_benzene, fr_benzodiazepine,
    fr_bicyclic, fr_diazo, fr_dihydropyridine, fr_epoxide, fr_ester, fr_ether, fr_furan, fr_guanido,
    fr_halogen, fr_hdrzine, fr_hdrzone, fr_imidazole, fr_imide, fr_isocyan, fr_isothiocyan, fr_ketone,
    fr_ketone_Topliss, fr_lactam, fr_lactone, fr_methoxy, fr_morpholine, fr_nitrile, fr_nitro,
    fr_nitro_arom, fr_nitro_arom_nonortho, fr_nitroso, fr_oxazole, fr_oxime, fr_para_hydroxylation,
    fr_phenol, fr_phenol_noOrthoHbond, fr_phos_acid, fr_phos_ester, fr_piperdine, fr_piperzine,
    fr_priamide, fr_prisulfonamd, fr_pyridine, fr_quatN, fr_sulfide, fr_sulfonamd, fr_sulfone,
    fr_term_acetylene, fr_tetrazole, fr_thiazole, fr_thiocyan, fr_thiophene, fr_unbrch_alkane, fr_urea

    3D descriptors: Asphericity, Autocorr3D, Eccentricity, GETAWAY, InertialShapeFactor, MORSE,
    NPR1, NPR2, PMI1, PMI2, PMI3, RDF, RadiusOfGyration, SpherocityIndex, WHIM

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd), SMILES (.smi,
    .txt, .csv, .tsv)

    The supported output file formats are: SD File (.sdf, .sd), CSV/TSV (.csv, .tsv, .txt)

Options:
    -a, --autocorr2DExclude <yes or no>  [default: yes]
        Exclude Autocorr2D descriptor from the calculation of 2D descriptors. 
    -f, --fragmentCount <yes or no>  [default: yes]
        Include 2D fragment count descriptors during the calculation. These descriptors are
        counted using SMARTS patterns specified in FragmentDescriptors.csv file distributed
        with RDKit. This option is only used during '2D' or 'All' value of '-m, --mode' option.
    -d, --descriptorNames <Name1,Name2,...>  [default: none]
        A comma delimited list of supported molecular descriptor names to calculate.
        This option is only used during 'Specify' value of '-m, --mode' option.
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
    -l, --list
        List molecular descriptors without performing any calculations.
    -m, --mode <2D, 3D, All, FragmentCountOnly, or Specify>  [default: 2D]
        Type of molecular descriptors to calculate. Possible values: 2D, 3D,
        All or Specify. The name of molecular descriptors must be specified using
        '-d, --descriptorNames' for 'Specify'. 2D descriptors also include 1D descriptors.
        The structure  of molecules must contain 3D coordinates for the  calculation
        of 3D descriptors.
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
            
        Default value for compute2DCoords: yes for SMILES input file; no for all other
        file types.
    -p, --precision <number>  [default: 3]
        Floating point precision for writing the calculated descriptor values.
    -s, --smilesOut <yes or no>  [default: no]
        Write out SMILES string to CSV/TSV text output file.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To compute all available 2D descriptors except Autocorr2D descriptor and
    write out a CSV file, type:

        % RDKitCalculateMolecularDescriptors.py  -i Sample.smi -o SampleOut.csv

    To compute all available 2D descriptors except Autocorr2D descriptor in
    multiprocessing mode on all available CPUs without loading all data into
    memory, and write out a CSV file, type:

        % RDKitCalculateMolecularDescriptors.py  --mp yes -i Sample.smi
          -o SampleOut.csv

    To compute all available 2D descriptors except Autocorr2D descriptor in
    multiprocessing mode on all available CPUs by loading all data into memory,
    and write out a CSV file, type:

        % RDKitCalculateMolecularDescriptors.py  --mp yes --mpParams
          "inputDataMode,InMemory" -i Sample.smi -o SampleOut.csv

    To compute all available 2D descriptors except Autocorr2D descriptor in
    multiprocessing mode on specific number of CPUs and chunk size without
    loading all data into memory, and write out a SDF file, type:

        % RDKitCalculateMolecularDescriptors.py  --mp yes --mpParams
          "inputDataMode,Lazy,numProcesses,4,chunkSize,8" -i Sample.smi
          -o SampleOut.sdf

    To compute all available 2D descriptors including Autocorr2D descriptor and
    excluding fragment count descriptors, and write out a TSV file, type:

        % RDKitCalculateMolecularDescriptors.py  -m 2D -a no -f no
          -i Sample.smi -o SampleOut.tsv

    To compute all available 3D descriptors and write out a SD file, type:

        % RDKitCalculateMolecularDescriptors.py  -m 3D -i Sample3D.sdf
          -o Sample3DOut.sdf

    To compute only fragment count 2D descriptors and write out a SD
    file file, type:

        % RDKitCalculateMolecularDescriptors.py  -m FragmentCountOnly
          -i Sample.sdf -o SampleOut.sdf

    To compute all available 2D and 3D descriptors including fragment count and
    Autocorr2D and write out a CSV file, type:

        % RDKitCalculateMolecularDescriptors.py  -m All -a no -i Sample.sdf
          -o SampleOut.csv

    To compute a specific set of 2D and 3D descriptors and write out a
    write out a TSV file, type:

        % RDKitCalculateMolecularDescriptors.py  -m specify
          -d 'MolWt,MolLogP,NHOHCount, NOCount,RadiusOfGyration'
          -i Sample3D.sdf -o SampleOut.csv

    To compute all available 2D descriptors except Autocorr2D descriptor for 
    molecules in a CSV SMILES file, SMILES strings in column 1, name in
    column 2, and write out a SD file without calculation of 2D coordinates, type:

        % RDKitCalculateMolecularDescriptors.py --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,no"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitCalculateRMSD.py, RDKitCompareMoleculeShapes.py, RDKitConvertFileFormat.py,
    RDKitGenerateConformers.py, RDKitPerformMinimization.py

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
