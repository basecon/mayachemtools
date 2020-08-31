#!/bin/env python
#
# File: RDKitRemoveSalts.py
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
    from rdkit.Chem.SaltRemover import SaltRemover
    from rdkit.Chem.SaltRemover import InputFormat
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
    RemoveSalts()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def RemoveSalts():
    """Identify and remove salts from molecules"""
    
    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % OptionsInfo["Infile"])
    Mols  = RDKitUtil.ReadMolecules(OptionsInfo["Infile"], **OptionsInfo["InfileParams"])
    
    # Set up a molecule writer...
    Writer = SetupMoleculeWriter()

    MolCount, ValidMolCount, SaltsMolCount = ProcessMolecules(Mols, Writer)

    if Writer is not None:
        Writer.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of ignored molecules: %d" % (MolCount - ValidMolCount))
    
    MiscUtil.PrintInfo("\nNumber of molecules coontaining salts: %d" % (SaltsMolCount))

def ProcessMolecules(Mols, Writer):
    """Process and remove salts from molecules. """

    if OptionsInfo["MPMode"]:
        return ProcessMoleculesUsingMultipleProcesses(Mols, Writer)
    else:
        return ProcessMoleculesUsingSingleProcess(Mols, Writer)

def ProcessMoleculesUsingSingleProcess(Mols,  Writer):
    """Process and remove salts from molecules using a single process. """
    
    MiscUtil.PrintInfo("\nRemoving salts...")
    
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    SetSMILESMolProps = OptionsInfo["OutfileParams"]["SetSMILESMolProps"]
    
    # Set up a salt remover...
    Remover = SetupSaltRemover()
    
    (MolCount, ValidMolCount, SaltsMolCount) = [0] * 3
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
                RDKitUtil.SetWriterMolProps(Writer, Mol)
        
        UnsaltedMol, SaltyStatus = RemoveMolSalts(Mol, Remover, MolCount)
        
        if SaltyStatus:
            SaltsMolCount += 1

        WriteMolecule(Writer, UnsaltedMol, Compute2DCoords)
    
    return (MolCount, ValidMolCount, SaltsMolCount)
    
def ProcessMoleculesUsingMultipleProcesses(Mols, Writer):
    """Process and remove salts from molecules using  multiprocessing."""
    
    MiscUtil.PrintInfo("\nRemoving salts using multiprocessing...")
    
    MPParams = OptionsInfo["MPParams"]
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    
    # Setup data for initializing a worker process...
    InitializeWorkerProcessArgs = (MiscUtil.ObjectToBase64EncodedString(Options), MiscUtil.ObjectToBase64EncodedString(OptionsInfo))

    # Setup a encoded mols data iterable for a worker process by pickling only public
    # and private molecule properties...
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
    
    SetSMILESMolProps = OptionsInfo["OutfileParams"]["SetSMILESMolProps"]
    
    (MolCount, ValidMolCount, SaltsMolCount) = [0] * 3
    FirstMol = True
    for Result in Results:
        MolCount += 1
        MolIndex, EncodedMol, SaltyStatus = Result
        
        if EncodedMol is None:
            continue
        ValidMolCount += 1
        
        Mol = RDKitUtil.MolFromBase64EncodedMolString(EncodedMol)
        
        if FirstMol:
            FirstMol = False
            if SetSMILESMolProps:
                RDKitUtil.SetWriterMolProps(Writer, Mol)
        
        if SaltyStatus:
            SaltsMolCount += 1

        WriteMolecule(Writer, Mol, Compute2DCoords)
    
    return (MolCount, ValidMolCount, SaltsMolCount)

def InitializeWorkerProcess(*EncodedArgs):
    """Initialize data for a worker process."""

    global Options, OptionsInfo
    
    MiscUtil.PrintInfo("Starting process (PID: %s)..." % os.getpid())

    # Decode Options and OptionInfo...
    Options = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[0])
    OptionsInfo = MiscUtil.ObjectFromBase64EncodedString(EncodedArgs[1])

    # Set up salt remover...
    OptionsInfo["SaltRemover"] = SetupSaltRemover()

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
        
    Mol, SaltyStatus = RemoveMolSalts(Mol, OptionsInfo["SaltRemover"], (MolIndex + 1))
    EncodedMol = RDKitUtil.MolToBase64EncodedMolString(Mol, PropertyPickleFlags = Chem.PropertyPickleOptions.MolProps | Chem.PropertyPickleOptions.PrivateProps)

    return [MolIndex, EncodedMol, SaltyStatus]
    
def RemoveMolSalts(Mol, Remover, MolNum):
    """Remove salts from mol and return unsalted mol along with mol salty status."""

    UnsaltedMol = Mol
    SaltyStatus = False
    
    if Remover is not None:
        KeptMol, DeletedMols = Remover.StripMolWithDeleted(Mol, dontRemoveEverything = False)
        if len(DeletedMols) >= 1:
            SaltyStatus = True
        if RDKitUtil.IsMolEmpty(KeptMol):
            if len(DeletedMols) >= 1:
                # Take the larged fragment from DeletedMols
                UnsaltedMol = GetLargestMol(DeletedMols)
    else:
        # Use largest fragment as unsalted molecule...
        MolFrags = Chem.GetMolFrags(Mol, asMols = True)
        if len(MolFrags) > 1:
            # Keep the largest fragment as unsalted molecule...
            SaltyStatus = True
            UnsaltedMol = GetLargestMol(MolFrags)

    if SaltyStatus:
        Chem.SanitizeMol(UnsaltedMol)
        MolName = RDKitUtil.GetMolName(Mol, MolNum)
        if len(MolName):
            UnsaltedMol.SetProp("_Name", MolName)
    
    return (UnsaltedMol, SaltyStatus)

def GetLargestMol(Mols):
    """Get largest mol from list of mols"""

    LargestMol = None
    LargestMolSize = -1
    for Mol in Mols:
        Size = Mol.GetNumAtoms()
        if Size > LargestMolSize:
            LargestMol = Mol
            LargestMolSize = Size

    return LargestMol

def SetupSaltRemover():
    """Setup a salt removerr."""
    
    Remover = None
    if OptionsInfo["SaltsByComponentsMode"]:
        return Remover

    return SaltRemover(defnFilename = OptionsInfo["SaltsFile"], defnData = OptionsInfo["SaltsSMARTS"], defnFormat = InputFormat.SMARTS)

def WriteMolecule(Writer, Mol, Compute2DCoords):
    """Write out molecule."""
    
    if OptionsInfo["CountMode"]:
        return
    
    if Compute2DCoords:
        AllChem.Compute2DCoords(Mol)
    
    Writer.write(Mol)

def SetupMoleculeWriter():
    """Setup a molecule writer."""
    
    Writer = None
    if OptionsInfo["CountMode"]:
        return Writer

    Writer = RDKitUtil.MoleculesWriter(OptionsInfo["Outfile"], **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % OptionsInfo["Outfile"])
    MiscUtil.PrintInfo("Generating file %s..." % OptionsInfo["Outfile"])
    
    return Writer

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
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

    SaltsByComponentsMode = False
    SaltsBySMARTSFileMode = False
    SaltsBySMARTSMode = False
    if re.match("^ByComponent$", Options["--saltsMode"], re.I):
        SaltsByComponentsMode = True
    elif re.match("^BySMARTSFile$", Options["--saltsMode"], re.I):
        SaltsBySMARTSFileMode = False
    elif re.match("^BySMARTS$", Options["--saltsMode"], re.I):
        SaltsBySMARTSMode = True
    else:
        MiscUtil.PrintError("The salts mode specified, %s, using \"--saltsMode\" option is not valid." % Options["--saltsMode"])
    OptionsInfo["SaltsByComponentsMode"]  = SaltsByComponentsMode
    OptionsInfo["SaltsBySMARTSFileMode"]  = SaltsBySMARTSFileMode
    OptionsInfo["SaltsBySMARTSMode"]  = SaltsBySMARTSMode

    SaltsFile = None
    if re.match("^BySMARTSFile$", Options["--saltsMode"], re.I):
        if not re.match("^auto$", Options["--saltsFile"], re.I):
            SaltsFile = Options["--saltsFile"]
    OptionsInfo["SaltsFile"] = SaltsFile
    
    SaltsSMARTS = None
    if re.match("^BySMARTS$", Options["--saltsMode"], re.I):
        if not Options["--saltsSMARTS"]:
            MiscUtil.PrintError("No salts SMARTS pattern specified using \"--saltsSMARTS\" option during \"BySMARTS\" value of \"-s, --saltsMode\" option")
        SaltsSMARTS = Options["--saltsSMARTS"].strip(" ")
        if not len(SaltsSMARTS):
            MiscUtil.PrintError("Empty SMARTS pattern specified using \"--saltsSMARTS\" option during \"BySMARTS\" value of \"-s, --saltsMode\" option")
        if re.search(" ", SaltsSMARTS):
            SaltsSMARTS = re.sub('[ ]+', '\n', SaltsSMARTS)
        
    OptionsInfo["SaltsSMARTS"] = SaltsSMARTS
    
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
    
    if Options["--outfile"]:
        MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sdf sd smi")
        MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
        MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])

    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "remove count")
    if re.match("^remove$", Options["--mode"], re.I):
        if not Options["--outfile"]:
            MiscUtil.PrintError("The outfile must be specified using \"-o, --outfile\" during \"remove\" value of \"-m, --mode\" option")
    
    MiscUtil.ValidateOptionTextValue("--mp", Options["--mp"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--saltsMode", Options["--saltsMode"], "ByComponent BySMARTSFile BySMARTS")
    
    if re.match("^BySMARTSFile$", Options["--saltsMode"], re.I):
        if not re.match("^auto$", Options["--saltsFile"], re.I):
            MiscUtil.ValidateOptionFilePath("--saltsFile", Options["--saltsFile"])

# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitRemoveSalts.py - Remove salts

Usage:
    RDKitRemoveSalts.py  [--infileParams <Name,Value,...>] [--mode <remove or count>]
                         [--mp <yes or no>] [--mpParams <Name.Value,...>] [--outfileParams <Name,Value,...> ]
                         [--overwrite] [--saltsMode <ByComponent, BySMARTSFile, BySMARTS>]
                         [--saltsFile <FileName or auto>] [--saltsSMARTS <SMARTS>]
                         [-w <dir>] [-o <outfile>]  -i <infile>
    RDKitRemoveSalts.py -h | --help | -e | --examples

Description:
    Remove salts from molecules or simply count the number of molecules containing
    salts. Salts are identified and removed based on either SMARTS strings or by selecting
    the largest disconnected components in molecules as non-salt portion of molecules.

    The supported input file formats are: SD (.sdf, .sd), SMILES (.smi., csv, .tsv, .txt)

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
    -m, --mode <remove or count>  [default: remove]
        Specify whether to remove salts from molecules and write out molecules
        or or simply count the number of molecules containing salts.
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
    -s, --saltsMode <ByComponent, BySMARTSFile, BySMARTS>  [default: ByComponent]
        Specify whether to identify and remove salts based on SMARTS strings or
        by selecting the largest disconnected component as non-salt portion of a
        molecule. Possible values: ByComponent, BySMARTSFile or BySMARTS.
    --saltsFile <FileName or auto>  [default: auto]
        Specify a file name containing specification for SMARTS corresponding to salts or
        use default salts file, Salts.txt, available in RDKit data directory. This option is only
        used during 'BySMARTSFile' value of '-s, --saltsMode' option.
        
        RDKit data format: Smarts<tab>Name(optional)
        
        For example:
            
            [Cl,Br,I]
            [N](=O)(O)O
            [CH3]C(=O)O	  Acetic acid
            
    --saltsSMARTS <SMARTS text>
        Space delimited SMARTS specifications to use for salts identification instead
        their specifications in '--saltsFile'. This option is only used during 'BySMARTS'
        value of '-s, --saltsMode' option.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To remove salts from molecules in a SMILES file by keeping largest disconnected
    components as non-salt portion of molecules and write out a SMILES file, type:

        % RDKitRemoveSalts.py -i Sample.smi -o SampleOut.smi

    To remove salts from molecules in a SMILES file by keeping largest disconnected
    components as non-salt portion of molecules, perform salt removal in multiprocessing
    mode on all available CPUs without loading all data into memory, and write out a
    SMILES file, type:

        % RDKitRemoveSalts.py --mp yes -i Sample.smi -o SampleOut.smi

    To remove salts from molecules in a SMILES file by keeping largest disconnected
    components as non-salt portion of molecules, perform salt removal in multiprocessing
    mode on all available CPUs by loading all data into memory, and write out a
    SMILES file, type:

        % RDKitRemoveSalts.py --mp yes --mpParams "inputDataMode,InMemory"
          -i Sample.smi -o SampleOut.smi

    To remove salts from molecules in a SMILES file by keeping largest disconnected
    components as non-salt portion of molecules, perform salt removal in multiprocessing
    mode on specific number of CPUs and chunk size without loading all data into memory,
    and write out a SMILES file, type:

        % RDKitRemoveSalts.py --mp yes --mpParams "inputDataMode,Lazy,
          numProcesses,4,chunkSize,8" -i Sample.smi -o SampleOut.smi

    To count number of molecule containing salts from in a SD file, using largest
    components as non-salt portion of molecules, without generating any output
    file, type:

        % RDKitRemoveSalts.py -m count -i Sample.sdf

    To remove salts from molecules in a SMILES file using SMARTS strings in default
    Salts.txt distributed with RDKit to identify salts and write out a SMILES file, type:

        % RDKitRemoveSalts.py -m remove -s BySMARTSFile -i Sample.smi
          -o SampleOut.smi

    To remove salts from molecules in a SD file using SMARTS strings in a local
    CustomSalts.txt to identify salts and write out a SMILES file, type:

        % RDKitRemoveSalts.py -m remove -s BySMARTSFile --saltsFile
          CustomSalts.txt -i Sample.sdf -o SampleOut.smi

    To remove salts from molecules in a SD file using specified SMARTS to identify
    salts and write out a SD file, type:

        % RDKitRemoveSalts.py -m remove -s BySMARTS  --saltsSMARTS
          '[Cl,Br,I]  [N](=O)(O)O [N](=O)(O)O'
          -i Sample.sdf -o SampleOut.smi

    To remove salts form  molecules from a CSV SMILES file, SMILES strings in column 1,
    name in column 2, and generate output SD file, type:

        % RDKitRemoveSalts.py --infileParams 
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitRemoveDuplicateMolecules.py,
    RDKitRemoveInvalidMolecules.py, RDKitSearchFunctionalGroups.py,
    RDKitSearchSMARTS.py

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
