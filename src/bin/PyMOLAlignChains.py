#!/usr/bin/env python
#
# File: PyMOLAlignChains.py
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2020 Manish Sud. All rights reserved.
#
# The functionality available in this script is implemented using PyMOL, a
# molecular visualization system on an open source foundation originally
# developed by Warren DeLano.
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

# PyMOL imports...
try:
    import pymol
    # Finish launching PyMOL in  a command line mode for batch processing (-c)
    # along with the following options:  disable loading of pymolrc and plugins (-k);
    # suppress start up messages (-q)
    pymol.finish_launching(['pymol', '-ckq'])
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import PyMOL module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your PyMOL environment and try again.\n\n")
    sys.exit(1)

# MayaChemTools imports...
try:
    from docopt import docopt
    import MiscUtil
    import PyMOLUtil
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import MayaChemTools module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your MayaChemTools environment and try again.\n\n")
    sys.exit(1)

ScriptName = os.path.basename(sys.argv[0])
Options = {}
OptionsInfo = {}

def main():
    """Start execution of the script"""
    
    MiscUtil.PrintInfo("\n%s (PyMOL v%s; %s) Starting...\n" % (ScriptName, pymol.cmd.get_version()[0], time.asctime()))
    
    (WallClockTime, ProcessorTime) = MiscUtil.GetWallClockAndProcessorTime()
    
    # Retrieve command line arguments and options...
    RetrieveOptions()
    
    # Process and validate command line arguments and options...
    ProcessOptions()

    # Perform actions required by the script...
    PerformChainAlignment()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformChainAlignment():
    """Align chains and write out new files."""

    MiscUtil.PrintInfo("\nGenerating output files...")

    # Load reffile for alignment..
    SetupAlignReference()

    # Perform alignment for each input file and write it out...
    for FileIndex in range(0, len(OptionsInfo["InfilesInfo"]["InfilesNames"])):
        SetupInputObject(FileIndex)
        AlignInputObject(FileIndex)
        WriteAlignedInputObject(FileIndex)
        DeleteInputObject(FileIndex)

    # Delete reference object...
    pymol.cmd.delete(OptionsInfo["RefFileInfo"]["PyMOLObjectName"])
    
def SetupAlignReference():
    """Setup object for alignment reference."""

    RefFile = OptionsInfo["RefFileInfo"]["RefFileName"]
    RefName = OptionsInfo["RefFileInfo"]["PyMOLObjectName"]

    pymol.cmd.load(RefFile, RefName)

def SetupInputObject(FileIndex):
    """Setup a PyMOL object for input file."""

    InputFile = OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex]
    InputName = OptionsInfo["InfilesInfo"]["PyMOLObjectNames"][FileIndex]
    
    pymol.cmd.load(InputFile, InputName)

def AlignInputObject(FileIndex):
    """Align input object to reference object."""
    
    RefName = OptionsInfo["RefFileInfo"]["PyMOLObjectName"]
    FitName = OptionsInfo["InfilesInfo"]["PyMOLObjectNames"][FileIndex]
    
    MiscUtil.PrintInfo("\nAligning %s to %s..." % (FitName, RefName))

    if re.match("^FirstChain$", OptionsInfo["AlignMode"], re.I):
        RefFirstChainID = OptionsInfo["RefFileInfo"]["ChainIDs"][0]
        RefSelection = "(%s and chain %s)" % (RefName, RefFirstChainID)
        
        FitFirstChainID = RetrieveFirstChainID(FileIndex)
        FitSelection = "(%s and chain %s)" % (FitName, FitFirstChainID)
    else:
        RefSelection = RefName
        FitSelection = FitName

    if re.match("^align$", OptionsInfo["AlignMethod"], re.I):
        pymol.cmd.align(FitSelection, RefSelection)
    elif re.match("^cealign$", OptionsInfo["AlignMethod"], re.I):
        pymol.cmd.cealign(RefSelection, FitSelection)
    elif re.match("^super$", OptionsInfo["AlignMethod"], re.I):
        pymol.cmd.super(FitSelection, RefSelection)
    else:
        MiscUtil.PrintError("Invalid alignment method: %s" % OptionsInfo["AlignMethod"])

def WriteAlignedInputObject(FileIndex):
    """Write out aligned input object"""
    
    Outfile = OptionsInfo["InfilesInfo"]["OutfilesNames"][FileIndex]
    InputName = OptionsInfo["InfilesInfo"]["PyMOLObjectNames"][FileIndex]
    
    MiscUtil.PrintInfo("Generating aligned output file %s..." % Outfile)
    
    pymol.cmd.save(Outfile, InputName)
    
    if not os.path.exists(Outfile):
        MiscUtil.PrintWarning("Failed to generate aligned output file, %s..." % (Outfile))

def DeleteInputObject(FileIndex):
    """Delete aligned input object."""
    
    InputName = OptionsInfo["InfilesInfo"]["PyMOLObjectNames"][FileIndex]
    pymol.cmd.delete(InputName)

def RetrieveInfilesInfo():
    """Retrieve information for input files."""

    InfilesInfo = {}
    
    InfilesInfo["InfilesNames"] = []
    InfilesInfo["InfilesRoots"] = []
    InfilesInfo["ChainIDs"] = []
    InfilesInfo["PyMOLObjectNames"] = []
    
    InfilesInfo["OutfilesNames"] = []

    OutSuffix = OptionsInfo["OutSuffix"]
    
    for Infile in OptionsInfo["InfilesNames"]:
        MiscUtil.PrintInfo("\nRetrieving chains information for input file %s..." % Infile)
        
        FileDir, FileName, FileExt = MiscUtil.ParseFileName(Infile)
        InfileRoot = FileName
        
        ChainIDs = RetrieveChainIDs(Infile, InfileRoot)
        if not len(ChainIDs):
            if re.match("^FirstChain$", OptionsInfo["AlignMode"], re.I):
                MiscUtil.PrintError("The align mode, %s, can't be used for aligning chains: No non-empty chain IDs found in input file." % (OptionsInfo["AlignMode"]))
        
        InfilesInfo["InfilesNames"].append(Infile)
        InfilesInfo["InfilesRoots"].append(InfileRoot)
        InfilesInfo["ChainIDs"].append(ChainIDs)
        
        InfilesInfo["PyMOLObjectNames"].append(InfileRoot)

        # Setup outfile name...
        Outfile = "%s%s.pdb" % (InfileRoot, OutSuffix)
        InfilesInfo["OutfilesNames"].append(Outfile)
        if os.path.exists(Outfile):
            if not OptionsInfo["Overwrite"]:
                MiscUtil.PrintError("The output file, %s, already exist. Use option \"--ov\" or \"--overwrite\" and try again.\n" % (Outfile))
    
    OptionsInfo["InfilesInfo"] = InfilesInfo

def RetrieveRefFileInfo():
    """Retrieve information for ref file."""

    RefFileInfo = {}

    RefFile = OptionsInfo["RefFileName"]
    
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(RefFile)
    RefFileRoot = FileName
    
    if re.match("^FirstInputFile$", OptionsInfo["AlignRefFile"], re.I):
        ChainIDs = OptionsInfo["InfilesInfo"]["ChainIDs"][0]
    else:
        MiscUtil.PrintInfo("\nRetrieving chains information for alignment reference file %s..." % RefFile)
        ChainIDs = RetrieveChainIDs(RefFile, RefFileRoot)
        if not len(ChainIDs):
            if re.match("^FirstChain$", OptionsInfo["AlignMode"], re.I):
                MiscUtil.PrintError("The align mode, %s, can't be used for aligning chains: No non-empty chain IDs found in input file." % (OptionsInfo["AlignMode"]))

    RefFileInfo["RefFileName"] = RefFile
    RefFileInfo["RefFileRoot"] = RefFileRoot
    RefFileInfo["PyMOLObjectName"] = "AlignRef_%s" % RefFileRoot
    RefFileInfo["ChainIDs"] = ChainIDs
    
    OptionsInfo["RefFileInfo"] = RefFileInfo

def RetrieveChainIDs(Infile, InfileRoot):
    """Retrieve chains IDs for an input file."""

    pymol.cmd.reinitialize()
    
    MolName = InfileRoot
    pymol.cmd.load(Infile, MolName)

    ChainIDs = PyMOLUtil.GetChains(MolName)
    pymol.cmd.delete(MolName)

    if ChainIDs is None:
        ChainIDs = []

    # Print out chain and ligand IDs...
    ChainInfo = ", ".join(ChainIDs) if len(ChainIDs) else "None"
    MiscUtil.PrintInfo("Chain IDs: %s" % ChainInfo)
                         
    return ChainIDs

def RetrieveFirstChainID(FileIndex):
    """Get first chain ID."""
    
    ChainIDs = OptionsInfo["InfilesInfo"]["ChainIDs"][FileIndex]
    
    FirstChainID = None
    if len(ChainIDs):
        FirstChainID = ChainIDs[0]
    
    return FirstChainID

def ProcessOptions():
    """Process and validate command line arguments and options"""

    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["AlignMethod"] = Options["--alignMethod"].lower()
    OptionsInfo["AlignMode"] = Options["--alignMode"]
    
    OptionsInfo["Infiles"] = Options["--infiles"]
    OptionsInfo["InfilesNames"] =  Options["--infileNames"]

    OptionsInfo["AlignRefFile"] = Options["--alignRefFile"]
    if re.match("^FirstInputFile$", Options["--alignRefFile"], re.I):
        OptionsInfo["RefFileName"] = OptionsInfo["InfilesNames"][0]
    else:
        OptionsInfo["RefFileName"] = Options["--alignRefFile"]
    
    OptionsInfo["OutSuffix"] = Options["--outSuffix"]
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]

    RetrieveInfilesInfo()
    RetrieveRefFileInfo()

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
    
    MiscUtil.ValidateOptionTextValue("--alignMethod", Options["--alignMethod"], "align cealign super")
    MiscUtil.ValidateOptionTextValue("--alignMode", Options["--alignMode"], "FirstChain Complex")
    
    # Expand infiles to handle presence of multiple input files...
    InfileNames = MiscUtil.ExpandFileNames(Options["--infiles"], ",")
    if len(InfileNames) < 2:
        MiscUtil.PrintError("Number of input files specified for \"-i, --infiles\" option, %d, must be greater than 2..." % (len(InfileNames)))

    # Validate file extensions...
    for Infile in InfileNames:
        MiscUtil.ValidateOptionFilePath("-i, --infiles", Infile)
        MiscUtil.ValidateOptionFileExt("-i, --infiles", Infile, "pdb cif")
    Options["--infileNames"] = InfileNames

    if not re.match("^FirstInputFile$", Options["--alignRefFile"], re.I):
        AlignRefFile = Options["--alignRefFile"]
        MiscUtil.ValidateOptionFilePath("--alignRefFile", AlignRefFile)
        MiscUtil.ValidateOptionFileExt("--alignRefFile", AlignRefFile, "pdb cif")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLAlignChains.py - Align chains

Usage:
    PyMOLAlignChains.py [--alignMethod <align, cealign, super>]
                        [--alignMode <FirstChain or Complex>] [--alignRefFile <filename>]
                        [--outSuffix <text>] [--overwrite] [-w <dir>] -i <infile1,infile2,infile3...>
    PyMOLAlignChains.py -h | --help | -e | --examples

Description:
    Align chains in input files to a reference file and write out aligned files.

    The supported input and output file format are: PDB (.pdb), CIF(.cif)

    The names of the aligned output files are automatically generated from the
    names of input as shown below:
    
        <InfileRoot><OutSuffux>.pdb
        Default: <InfileRoot>_Aligned.pdb
    
Options:
    -a, --alignMethod <align, cealign, super>  [default: super]
        Alignment methodology to use for aligning input files to a
        reference file.
    --alignMode <FirstChain or Complex>  [default: FirstChain]
        Portion of input and reference files to use for spatial alignment of
        input files against reference file.  Possible values: FirstChain or
        Complex.
        
        The FirstChain mode allows alignment of the first chain in each input
        file to the first chain in the reference file along with moving the rest
        of the complex to coordinate space of the reference file. The complete
        complex in each input file is aligned to the complete complex in reference
        file for the Complex mode.
    --alignRefFile <filename>  [default: FirstInputFile]
        Reference input file name. The default is to use the first input file
        name specified using '-i, --infiles' option.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infiles <infile1,infile2,...>
        A comma delimited list of input files. The wildcards are also allowed
        in file names.
    --outSuffix <text>  [default: _Aligned]
        Suffix to append to input file root for generating name of output file.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To align first chain in all input files to the first chain in first input file
    and write out aligned output files, type:

        % PyMOLAlignChains.py -i "Sample3.pdb,Sample4.pdb,Sample5.pdb"

    To align first chain in all input files to the first chain in specific reference
    file and write out aligned output files, type:

        % PyMOLAlignChains.py --alignRefFile Sample5.pdb
          -i "Sample3.pdb,Sample4.pdb,Sample5.pdb"

    To align first chain in all input files to the first chain in first input file
    using a specific alignment method and write out aligned output files
    with specific suffix in names, type:

        % PyMOLAlignChains.py --alignMethod cealign --outSuffix "_aligned"
          -i "Sample3.pdb,Sample4.pdb,Sample5.pdb"

    To align all chains in each input files to all chains in first input file and
    write out aligned output files, type:

        % PyMOLAlignChains.py --alignMode Complex
          -i "Sample3.pdb,Sample4.pdb,Sample5.pdb"

Author:
    Manish Sud(msud@san.rr.com)

See also:
    PyMOLCalculateRMSD.py, PyMOLSplitChainsAndLigands.py,
    PyMOLVisualizeMacromolecules.py

Copyright:
    Copyright (C) 2020 Manish Sud. All rights reserved.

    The functionality available in this script is implemented using PyMOL, a
    molecular visualization system on an open source foundation originally
    developed by Warren DeLano.

    This file is part of MayaChemTools.

    MayaChemTools is free software; you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option) any
    later version.

"""

if __name__ == "__main__":
    main()
