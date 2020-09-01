#!/usr/bin/env python
#
# File: RDKitFilterPAINS.py
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
    PerformFiltering()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformFiltering():
    """Filter molecules using SMARTS specified in PAINS filter file."""

    # Setup PAINS patterns and pattern mols...
    MiscUtil.PrintInfo("\nSetting up PAINS pattern molecules for performing substructure search...")
    PAINSPatterns = RetrievePAINSPatterns()
    PAINSPatternMols = SetupPAINSPatternMols(PAINSPatterns)
    
    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % OptionsInfo["Infile"])
    Mols  = RDKitUtil.ReadMolecules(OptionsInfo["Infile"], **OptionsInfo["InfileParams"])
    
    # Set up molecule writers...
    Writer, WriterFiltered = SetupMoleculeWriters()
    
    MolCount, ValidMolCount, RemainingMolCount = ProcessMolecules(Mols, PAINSPatternMols, Writer, WriterFiltered)
    
    if Writer is not None:
        Writer.close()
    if WriterFiltered is not None:
        WriterFiltered.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))

    MiscUtil.PrintInfo("\nNumber of remaining molecules: %d" % RemainingMolCount)
    MiscUtil.PrintInfo("Number of filtered molecules: %d" % (ValidMolCount - RemainingMolCount))

def ProcessMolecules(Mols, PAINSPatternMols, Writer, WriterFiltered):
    """Process and filter molecules. """
    
    if OptionsInfo["MPMode"]:
        return ProcessMoleculesUsingMultipleProcesses(Mols, PAINSPatternMols, Writer, WriterFiltered)
    else:
        return ProcessMoleculesUsingSingleProcess(Mols, PAINSPatternMols, Writer, WriterFiltered)

def ProcessMoleculesUsingSingleProcess(Mols, PAINSPatternMols, Writer, WriterFiltered):
    """Process and filter molecules using a single process."""
    
    NegateMatch = OptionsInfo["NegateMatch"]
    OutfileFilteredMode = OptionsInfo["OutfileFilteredMode"]
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    SetSMILESMolProps = OptionsInfo["OutfileParams"]["SetSMILESMolProps"]
    
    MiscUtil.PrintInfo("\nFiltering molecules...")
    
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
                if WriterFiltered is not None:
                    RDKitUtil.SetWriterMolProps(WriterFiltered, Mol)
        
        MolMatched = DoesMoleculeContainsPAINSPattern(Mol, PAINSPatternMols)
        if MolMatched == NegateMatch:
            RemainingMolCount += 1
            WriteMolecule(Writer, Mol, Compute2DCoords)
        else:
            if OutfileFilteredMode:
                WriteMolecule(WriterFiltered, Mol, Compute2DCoords)
    
    return (MolCount, ValidMolCount, RemainingMolCount)
    
def ProcessMoleculesUsingMultipleProcesses(Mols, PAINSPatternMols, Writer, WriterFiltered):
    """Process and filter molecules using multiprocessing."""
    
    MiscUtil.PrintInfo("\nFiltering molecules using multiprocessing...")
    
    MPParams = OptionsInfo["MPParams"]
    NegateMatch = OptionsInfo["NegateMatch"]
    OutfileFilteredMode = OptionsInfo["OutfileFilteredMode"]
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    SetSMILESMolProps = OptionsInfo["OutfileParams"]["SetSMILESMolProps"]
    
    # Setup data for initializing a worker process...
    MiscUtil.PrintInfo("Encoding options info and PAINS pattern molecules...")
    OptionsInfo["EncodedPAINSPatternMols"] = [RDKitUtil.MolToBase64EncodedMolString(PatternMol) for PatternMol in PAINSPatternMols]
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
    
    (MolCount, ValidMolCount, RemainingMolCount) = [0] * 3
    FirstMol = True
    for Result in Results:
        MolCount += 1
        MolIndex, EncodedMol, MolMatched = Result
        
        if EncodedMol is None:
            continue
        ValidMolCount += 1
        
        Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
        
        if FirstMol:
            FirstMol = False
            if SetSMILESMolProps:
                if Writer is not None:
                    RDKitUtil.SetWriterMolProps(Writer, Mol)
                if WriterFiltered is not None:
                    RDKitUtil.SetWriterMolProps(WriterFiltered, Mol)
        
        if MolMatched == NegateMatch:
            RemainingMolCount += 1
            WriteMolecule(Writer, Mol, Compute2DCoords)
        else:
            if OutfileFilteredMode:
                WriteMolecule(WriterFiltered, Mol, Compute2DCoords)
    
    return (MolCount, ValidMolCount, RemainingMolCount)

def InitializeWorkerProcess(*EncodedArgs):
    """Initialize data for a worker process."""

    global Options, OptionsInfo
    
    MiscUtil.PrintInfo("Starting process (PID: %s)..." % os.getpid())

    # Decode Options and OptionInfo...
    Options = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[0])
    OptionsInfo = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[1])

    # Decode PAINSPatternMols...
    OptionsInfo["PAINSPatternMols"] = [RDKitUtil.MolFromBase64EncodedMolString(EncodedMol) for EncodedMol in OptionsInfo["EncodedPAINSPatternMols"]]
    
def WorkerProcess(EncodedMolInfo):
    """Process data for a worker process."""
    
    MolIndex, EncodedMol = EncodedMolInfo
    
    if EncodedMol is None:
        return [MolIndex, None, False]
        
    Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
    if RDKitUtil.IsMolEmpty(Mol):
        MolName = RDKitUtil.GetMolName(Mol, (MolIndex + 1))
        MiscUtil.PrintWarning("Ignoring empty molecule: %s" % MolName)
        return [MolIndex, None, False]
        
    MolMatched = DoesMoleculeContainsPAINSPattern(Mol, OptionsInfo["PAINSPatternMols"])

    return [MolIndex, EncodedMol, MolMatched]
    
def WriteMolecule(Writer, Mol, Compute2DCoords):
    """Write out molecule."""
    
    if OptionsInfo["CountMode"]:
        return
    
    if Compute2DCoords:
        AllChem.Compute2DCoords(Mol)
    
    Writer.write(Mol)

def SetupMoleculeWriters():
    """Setup molecule writers."""
    
    Writer = None
    WriterFiltered = None

    if OptionsInfo["CountMode"]:
        return (Writer, WriterFiltered)

    Writer = RDKitUtil.MoleculesWriter(OptionsInfo["Outfile"], **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % OptionsInfo["Outfile"])
    MiscUtil.PrintInfo("Generating file %s..." % OptionsInfo["Outfile"])
    
    if OptionsInfo["OutfileFilteredMode"]:
        WriterFiltered = RDKitUtil.MoleculesWriter(OptionsInfo["OutfileFiltered"], **OptionsInfo["OutfileParams"])
        if WriterFiltered is None:
            MiscUtil.PrintError("Failed to setup a writer for output fie %s " % OptionsInfo["OutfileFiltered"])
        MiscUtil.PrintInfo("Generating file %s..." % OptionsInfo["OutfileFiltered"])
    
    return (Writer, WriterFiltered)

def DoesMoleculeContainsPAINSPattern(Mol, PAINSPatternMols):
    """Check presence of PAINS pattern in the molecule"""

    MolMatched = False
    
    for PatternMol in PAINSPatternMols:
        if Mol.HasSubstructMatch(PatternMol, useChirality = True):
            MolMatched = True
            break
        
    return MolMatched
    
def RetrievePAINSPatterns():
    """Retrieve PAINS patterns for specified PAINS mode"""

    SMARTSPatterns = []
    for FilterType in OptionsInfo["SpecifiedFilterTypes"]:
        SMARTSPatterns.extend(OptionsInfo["PAINSFiltersMap"]["SMARTS"][FilterType])

    return SMARTSPatterns

def SetupPAINSPatternMols(PAINSPatterns):
    """Set up PAINS pattern mols for substructure search"""

    PatternMols = []
    for Pattern in PAINSPatterns:
        PatternMol = Chem.MolFromSmarts(Pattern)
        if PatternMol is None:
            MiscUtil.PrintWarning("Failed to convert PAINS pattern, %s, into a molecule..." % Pattern)
            continue
        PatternMols.append(PatternMol)
        
    return PatternMols    

def ProcessPAINSMode():
    """Process specified PAINS mode. """
    
    # Retrieve filetrs information...
    RetrievePAINSFiltersInfo()
    
    # Process PAINS mode...
    OptionsInfo["SpecifiedFilterTypes"] = OptionsInfo["PAINSFiltersMap"]["FilterTypes"]
    if re.match("^All$", OptionsInfo["PAINSMode"], re.I):
        return
    
    PAINSMode = re.sub(" ", "", OptionsInfo["PAINSMode"])
    if not len(PAINSMode):
        MiscUtil.PrintError("The PAINSMode mode specified using \"-p, --painsMode\" option are empty.")

    CanonicalFilterTypesMap = {}
    for FilterType in OptionsInfo["PAINSFiltersMap"]["FilterTypes"]:
        CanonicalFilterTypesMap[FilterType.lower()] = FilterType

    SpecifiedFilterTypes = []
    for FilterType in PAINSMode.split(","):
        CanonicalFilterType = FilterType.lower()
        if not CanonicalFilterType in CanonicalFilterTypesMap:
            MiscUtil.PrintError("The PAINS mode, %s, specified using \"-p, --PAINSMode\" is not valid. Supported PAINS modes: %s" % (FilterType, ", ".join(OptionsInfo["PAINSFiltersMap"]["FilterTypes"])))

        SpecifiedFilterTypes.append(CanonicalFilterTypesMap[CanonicalFilterType])

    OptionsInfo["SpecifiedFilterTypes"] = SpecifiedFilterTypes

def RetrievePAINSFiltersInfo():
    """Retrieve information for PAINS filters."""
    
    MayaChemToolsDataDir = MiscUtil.GetMayaChemToolsLibDataPath()
    PAINSFiltersFilePath = os.path.join(MayaChemToolsDataDir, "PAINSFilters.csv")
    
    MiscUtil.PrintInfo("\nRetrieving PAINS SMARTS patterns from file %s" % (PAINSFiltersFilePath))

    Delimiter = ','
    QuoteChar = '"'
    IgnoreHeaderLine = True
    FilterLinesWords = MiscUtil.GetTextLinesWords(PAINSFiltersFilePath, Delimiter, QuoteChar, IgnoreHeaderLine)

    PAINSFiltersMap = {}
    PAINSFiltersMap["FilterTypes"] = []
    PAINSFiltersMap["ID"] = {}
    PAINSFiltersMap["SMARTS"] = {}

    for LineWords in FilterLinesWords:
        FilterType = LineWords[0]
        ID = LineWords[1]
        SMARTS = LineWords[2]

        if not FilterType in PAINSFiltersMap["FilterTypes"]:
            PAINSFiltersMap["FilterTypes"].append(FilterType)
            PAINSFiltersMap["ID"][FilterType] = []
            PAINSFiltersMap["SMARTS"][FilterType] = []

        PAINSFiltersMap["ID"][FilterType].append(ID)
        PAINSFiltersMap["SMARTS"][FilterType].append(SMARTS)

    OptionsInfo["PAINSFiltersMap"] = PAINSFiltersMap
    
    MiscUtil.PrintInfo("\nTotal number filters: %d" % len(FilterLinesWords))
    MiscUtil.PrintInfo("Number of filter family types: %d\nFilter familty types: %s\n" % (len(PAINSFiltersMap["FilterTypes"]), ", ".join(PAINSFiltersMap["FilterTypes"])))

    for FilterType in PAINSFiltersMap["FilterTypes"]:
        MiscUtil.PrintInfo("Filter family type: %s; Number of filters: %d" % (FilterType, len(PAINSFiltersMap["ID"][FilterType])))

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"], Options["--infile"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"], Options["--infile"], Options["--outfile"])
    
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Options["--outfile"])
    OutfileFiltered = "%s_Filtered.%s" % (FileName, FileExt)
    OptionsInfo["OutfileFiltered"] = OutfileFiltered
    OptionsInfo["OutfileFilteredMode"] = True if re.match("^yes$", Options["--outfileFiltered"], re.I) else False
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    OptionsInfo["CountMode"] = True if re.match("^count$", Options["--mode"], re.I) else False
    OptionsInfo["NegateMatch"] = True if re.match("^yes$", Options["--negate"], re.I) else False

    OptionsInfo["MPMode"] = True if re.match("^yes$", Options["--mp"], re.I) else False
    OptionsInfo["MPParams"] = MiscUtil.ProcessOptionMultiprocessingParameters("--mpParams", Options["--mpParams"])

    OptionsInfo["PAINSMode"] = Options["--painsMode"]
    ProcessPAINSMode()
    
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
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "sdf sd smi txt csv tsv")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])

    MiscUtil.ValidateOptionTextValue("--outfileFiltered", Options["--outfileFiltered"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "filter count")
    if re.match("^filter$", Options["--mode"], re.I):
        if not Options["--outfile"]:
            MiscUtil.PrintError("The outfile must be specified using \"-o, --outfile\" during \"filter\" value of \"-m, --mode\" option")
        
    MiscUtil.ValidateOptionTextValue("--mp", Options["--mp"], "yes no")
    MiscUtil.ValidateOptionTextValue("-n, --negate", Options["--negate"], "yes no")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitFilterPAINS.py - Filter PAINS molecules

Usage:
    RDKitFilterPAINS.py  [--infileParams <Name,Value,...>] [--mode <filter or count>]
                         [--mp <yes or no>] [--mpParams <Name.Value,...>]
                         [--outfileFiltered <yes or no>] [ --outfileParams <Name,Value,...> ]
                         [--painsMode <All or A, B, C>] [--negate <yes or no>]
                         [--overwrite] [-w <dir>] -i <infile> -o <outfile>
    RDKitFilterPAINS.py -h | --help | -e | --examples

Description:
    Filter Pan-assay Interference molecules (PAINS) [ Ref 130 - 131 ] from an input
    file by performing a substructure search using SMARTS pattern specified in
    MAYACHEMTOOLS/lib/data/PAINSFilters.csv file and write out appropriate
    molecules to an output file or simply count the number of filtered molecules.

    The supported input file formats are: SD (.sdf, .sd), SMILES (.smi, .csv,
    .tsv, .txt)

    The supported output file formats are: SD (.sdf, .sd), SMILES (.smi)

Options:
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
            
            SD: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            SMILES: smilesColumn,1,smilesNameColumn,2,smilesDelimiter,space,
                smilesTitleLine,auto,sanitize,yes
            
        Possible values for smilesDelimiter: space, comma or tab.
    -m, --mode <filter or count>  [default: filter]
        Specify whether to filter the matched molecules and write out the rest of the 
        molecules to an outfile or simply count the number of matched molecules
        marked for filtering.
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
    -n, --negate <yes or no>  [default: no]
        Specify whether to filter molecules not matching the PAINS filters specified by
        SMARTS patterns.
    -o, --outfile <outfile>
        Output file name.
    --outfileFiltered <yes or no>  [default: no]
        Write out a file containing filtered molecules. Its name is automatically
        generated from the specified output file. Default: <OutfileRoot>_
        Filtered.<OutfileExt>.
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
    -p, --painsMode <All or A, B, or C>  [default: All]
        All or a comma delimited list of PAINS filter family type to used for
        filtering molecules. 
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To count the number of molecules not containing any substructure corresponding to
    PAINS SMARTS patterns and write out a SMILES file, type: 

        % RDKitFilterPAINS.py -i Sample.smi -o SampleOut.smi

    To count the number of molecules not containing any substructure corresponding to
    PAINS SMARTS patterns, perform filtering in multiprocessing mode on all available
    CPUs without loading all data into memory, and write out a SMILES file, type: 

        % RDKitFilterPAINS.py --mp yes -i Sample.smi -o SampleOut.smi

    To count the number of molecules not containing any substructure corresponding to
    PAINS SMARTS patterns, perform filtering in multiprocessing mode on all available
    CPUs by loading all data into memory, and write out a SMILES file, type: 

        % RDKitFilterPAINS.py --mp yes --mpParams "inputDataMode,InMemory"
          -i Sample.smi -o SampleOut.smi

    To count the number of molecules not containing any substructure corresponding to
    PAINS SMARTS patterns, perform filtering in multiprocessing mode on specific
    number of CPUs and chunk size without loading all data into memory, and
    write out a SMILES file, type: 

        % RDKitFilterPAINS.py --mp yes --mpParams "inputDataMode,Lazy,
          numProcesses,4,chunkSize,8" -i Sample.smi -o SampleOut.smi

    To count the number of molecules not containing any substructure corresponding to
    PAINS SMARTS patterns and write out a SMILES file containing these and filtered
    molecules, type: 

        % RDKitFilterPAINS.py --outfileFiltered yes -i Sample.smi
          -o SampleOut.smi

    To only count the number of molecules not containing any substructure corresponding
    to PAINS SMARTS patterns without writing out any file, type: 

        % RDKitFilterPAINS.py -m count -i Sample.sdf -o SampleOut.smi

    To count the number of molecules containing any substructure corresponding to
    PAINS SMARTS patterns and write out a SD file with computed 2D coordinates,
    type: 

        % RDKitFilterPAINS.py -n yes -i Sample.smi -o SampleOut.sdf

    To count the number of molecules not containing any substructure corresponding to
    PAINS SMARTS patterns family of Type A in a CSV SMILES file and write out a SD file, type: 

        % RDKitFilterPAINS.py --painsMode A --infileParams
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitFilterChEMBLAlerts.py, RDKitConvertFileFormat.py, RDKitSearchSMARTS.py

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
