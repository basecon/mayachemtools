#!/bin/env python
#
# File: PyMOLCalculateRMSD.py
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
    CalculateRMSDValues()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def CalculateRMSDValues():
    """Calculate RMSD between reference and probe files."""

    Outfile = OptionsInfo["Outfile"]
    OutDelim = OptionsInfo["OutDelim"]
    
    MiscUtil.PrintInfo("\nGenerating file %s..." % Outfile)
    OutFH = open(Outfile, "w")
    if OutFH is None:
        MiscUtil.PrintError("Couldn't open output file: %s.\n" % (Outfile))

    WriteColumnLabels(OutFH, OutDelim)
    
    pymol.cmd.reinitialize()
    if re.match("^OneToOne$", OptionsInfo["Mode"], re.I):
        CalculateOneToOneRMSDValues(OutFH, OutDelim)
    elif re.match("^AllToAll$", OptionsInfo["Mode"], re.I):
        CalculateAllToAllRMSDValues(OutFH, OutDelim)
    elif re.match("^FirstToAll$", OptionsInfo["Mode"], re.I):
        CalculateFirstToAllRMSDValues(OutFH, OutDelim)
    else:
        MiscUtil.PrintError("RMSD couldn't be calculated: Specified mode, %s, is not supported" % OptionsInfo["Mode"])

    OutFH.close()

def CalculateOneToOneRMSDValues(OutFH, OutDelim):
    """Calculate pairwise RMSD values."""
    
    RefFilesCount = len(OptionsInfo["RefFilesNames"])
    ProbeFilesCount = len(OptionsInfo["ProbeFilesNames"])

    FilesCount = ProbeFilesCount if RefFilesCount > ProbeFilesCount else RefFilesCount
    
    if RefFilesCount != ProbeFilesCount:
        MiscUtil.PrintWarning("Number of reference files, %d,  is not equal to number of probe files, %d .\n" % (RefFilesCount, ProbeFilesCount))
        MiscUtil.PrintWarning("Pairwise RMSD will be calculated only for first %s files.\n" % (FilesCount))

    # Process files...
    for FileIndex in range(0, FilesCount):
        RefFileIndex = FileIndex
        ProbeFileIndex = FileIndex
        
        LoadRefFile(RefFileIndex)
        LoadProbeFile(ProbeFileIndex)
        
        RMSD = CalculateRMSDValue(RefFileIndex, ProbeFileIndex)
        
        RefID = OptionsInfo["RefFilesInfo"]["FilesRoots"][RefFileIndex]
        ProbeID = OptionsInfo["ProbeFilesInfo"]["FilesRoots"][ProbeFileIndex]
        Line = "%s%s%s%s%s\n" % (RefID, OutDelim, ProbeID, OutDelim, RMSD)
        OutFH.write(Line)
        
        DeleteRefObject(RefFileIndex)
        DeleteProbeObject(ProbeFileIndex)
        
def CalculateAllToAllRMSDValues(OutFH, OutDelim):
    """Calculate RMSD values between all pairs of files."""
    
    RefFilesCount = len(OptionsInfo["RefFilesNames"])
    ProbeFilesCount = len(OptionsInfo["ProbeFilesNames"])
    OutMatrix = OptionsInfo["OutMatrix"]
    
    for RefFileIndex in range(0, RefFilesCount):
        LoadRefFile(RefFileIndex)
        RefID = OptionsInfo["RefFilesInfo"]["FilesRoots"][RefFileIndex]

        LineWords = []
        if OutMatrix:
            LineWords.append(RefID)
            
        for ProbeFileIndex in range(0, ProbeFilesCount):
            LoadProbeFile(ProbeFileIndex)
            RMSD = CalculateRMSDValue(RefFileIndex, ProbeFileIndex)
            DeleteProbeObject(ProbeFileIndex)
            
            if OutMatrix:
                LineWords.append(RMSD)
            else:
                ProbeID = OptionsInfo["ProbeFilesInfo"]["FilesRoots"][ProbeFileIndex]
                Line = "%s%s%s%s%s\n" % (RefID, OutDelim, ProbeID, OutDelim, RMSD)
                OutFH.write(Line)
        
        DeleteRefObject(RefFileIndex)
        
        if OutMatrix:
            Line = OutDelim.join(LineWords)
            OutFH.write("%s\n" % Line)

def CalculateFirstToAllRMSDValues(OutFH, OutDelim):
    """Calculate RMSD values between first reference file and all probe files. """

    # Setup reference...
    RefFileIndex = 0
    RefID = OptionsInfo["RefFilesInfo"]["FilesRoots"][RefFileIndex]
    LoadRefFile(RefFileIndex)

    # Go over probe files...
    for ProbeFileIndex in range(0, len(OptionsInfo["ProbeFilesNames"])):
        LoadProbeFile(ProbeFileIndex)
        
        RMSD = CalculateRMSDValue(RefFileIndex, ProbeFileIndex)
        
        ProbeID = OptionsInfo["ProbeFilesInfo"]["FilesRoots"][ProbeFileIndex]
        Line = "%s%s%s%s%s\n" % (RefID, OutDelim, ProbeID, OutDelim, RMSD)
        OutFH.write(Line)
        
        DeleteProbeObject(ProbeFileIndex)

    DeleteRefObject(RefFileIndex)

def WriteColumnLabels(OutFH, OutDelim):
    """Write out column labels. """

    ColLabels = []
    
    if re.match("^AllToAll$", OptionsInfo["Mode"], re.I) and OptionsInfo["OutMatrix"]:
        ColLabels.append("")
        ColLabels.extend(OptionsInfo["ProbeFilesInfo"]["FilesRoots"])
    else:
        ColLabels = ["RefFileID", "ProbeFileID", "RMSD"]
    
    Line = OutDelim.join(ColLabels)
    OutFH.write("%s\n" % Line)
    
def LoadRefFile(RefFileIndex):
    """Load reference file. """
    
    RefFile = OptionsInfo["RefFilesNames"][RefFileIndex]
    RefName = OptionsInfo["RefFilesInfo"]["PyMOLObjectNames"][RefFileIndex]
    LoadFile(RefFile, RefName)
    
def LoadProbeFile(ProbeFileIndex):
    """Load probe file. """
    
    ProbeFile = OptionsInfo["ProbeFilesNames"][ProbeFileIndex]
    ProbeName = OptionsInfo["ProbeFilesInfo"]["PyMOLObjectNames"][ProbeFileIndex]
    LoadFile(ProbeFile, ProbeName)

def LoadFile(FileName, ObjectName):
    """Load a file. """
    
    pymol.cmd.load(FileName, ObjectName)
    
def DeleteRefObject(RefFileIndex):
    """Delete reference object. """
    
    RefName = OptionsInfo["RefFilesInfo"]["PyMOLObjectNames"][RefFileIndex]
    DeleteObject(RefName)
    
def DeleteProbeObject(ProbeFileIndex):
    """Delete probe object. """
    
    ProbeName = OptionsInfo["ProbeFilesInfo"]["PyMOLObjectNames"][ProbeFileIndex]
    DeleteObject(ProbeName)

def DeleteObject(Name):
    """Delete PyMOL object. """
    
    pymol.cmd.delete(Name)
    
def CalculateRMSDValue(RefFileIndex, ProbeFileIndex):
    """Calculate RMSD value between referece and probe objects. """

    RefName = OptionsInfo["RefFilesInfo"]["PyMOLObjectNames"][RefFileIndex]
    ProbeName = OptionsInfo["ProbeFilesInfo"]["PyMOLObjectNames"][ProbeFileIndex]
    
    if re.match("^FirstChain$", OptionsInfo["AlignMode"], re.I):
        RefFirstChainID = OptionsInfo["RefFilesInfo"]["ChainIDs"][RefFileIndex][0]
        RefSelection = "(%s and chain %s)" % (RefName, RefFirstChainID)
        
        ProbeFirstChainID = OptionsInfo["ProbeFilesInfo"]["ChainIDs"][ProbeFileIndex][0]
        ProbeSelection = "(%s and chain %s)" % (ProbeName, ProbeFirstChainID)
    else:
        RefSelection = RefName
        ProbeSelection = ProbeName
        
    RMSD = CalculateRMSD(RefSelection, ProbeSelection, OptionsInfo["AlignMethod"])

    return RMSD

def CalculateRMSD(RefSelectionName, ProbeSelectionName, AlignMethod):
    """Calculate RMSD between two selections after aligning the selections."""
    
    if re.match("^super$", AlignMethod, re.I):
        Results = pymol.cmd.super(ProbeSelectionName, RefSelectionName)
        RMSD = Results[0]
    elif re.match("^cealign$", AlignMethod, re.I):
        Results = pymol.cmd.super(RefSelectionName, ProbeSelectionName)
        RMSD = Results["RMSD"]
    elif re.match("^align$", AlignMethod, re.I):
        Results = pymol.cmd.super(ProbeSelectionName, RefSelectionName)
        RMSD = Results[0]
    else:
        RMSD = None
        MiscUtil.PrintWarning("Failed to calculate RMSD. Unknown alignment method: %s" % AlignMethod)
    
    if RMSD is not None:
        RMSD = "%.2f" % RMSD
    
    return RMSD

def RetrieveProbeFilesInfo():
    """Retrieve information for probe input files."""

    RetrieveInfilesInfo("ProbeFiles")
    
def RetrieveRefFilesInfo():
    """Retrieve information for reference input files."""

    RetrieveInfilesInfo("RefFiles")
    
def RetrieveInfilesInfo(InfilesMode):
    """Retrieve information for input files."""

    if re.match("^ProbeFiles$", InfilesMode, re.I):
        MiscUtil.PrintInfo("Retrieving information for probe files...")
        InfilesNames = OptionsInfo["ProbeFilesNames"]
        NameSuffix = "_Probe"
    elif re.match("^RefFiles$", InfilesMode, re.I):
        MiscUtil.PrintInfo("Retrieving information for reference files...")
        InfilesNames = OptionsInfo["RefFilesNames"]
        NameSuffix = "_Ref"
    else:
        MiscUtil.PrintError("Internal Error: Unknown infiles mode: %s" % InfilesMode)
    
    InfilesInfo = {}
    
    InfilesInfo["FilesNames"] = []
    InfilesInfo["FilesRoots"] = []
    InfilesInfo["ChainIDs"] = []
    InfilesInfo["PyMOLObjectNames"] = []

    for Infile in InfilesNames:
        MiscUtil.PrintInfo("\nRetrieving chains information for input file %s..." % Infile)
        
        FileDir, FileName, FileExt = MiscUtil.ParseFileName(Infile)
        InfileRoot = FileName
        
        ChainIDs = RetrieveChainIDs(Infile, InfileRoot)
        if not len(ChainIDs):
            if re.match("^FirstChain$", OptionsInfo["AlignMode"], re.I):
                MiscUtil.PrintError("The align mode, %s, can't be used for calculating RMSD: No non-empty chain IDs found in input file." % (OptionsInfo["AlignMode"]))
        
        InfilesInfo["FilesNames"].append(Infile)
        InfilesInfo["FilesRoots"].append(InfileRoot)
        InfilesInfo["ChainIDs"].append(ChainIDs)

        Name = "%s%s" % (InfileRoot, NameSuffix)
        InfilesInfo["PyMOLObjectNames"].append(Name)

    if re.match("^ProbeFiles$", InfilesMode, re.I):
        OptionsInfo["ProbeFilesInfo"] = InfilesInfo
    elif re.match("^RefFiles$", InfilesMode, re.I):
        OptionsInfo["RefFilesInfo"] = InfilesInfo
        
def RetrieveChainIDs(Infile, InfileRoot):
    """Retrieve chains IDs for an input file."""

    pymol.cmd.reinitialize()
    
    MolName = InfileRoot
    pymol.cmd.load(Infile, MolName)

    ChainIDs = PyMOLUtil.GetChains(MolName, RemoveEmpty = True)
    pymol.cmd.delete(MolName)

    if ChainIDs is None:
        ChainIDs = []

    # Print out chain and ligand IDs...
    ChainInfo = ", ".join(ChainIDs) if len(ChainIDs) else "None"
    MiscUtil.PrintInfo("Chain IDs: %s" % ChainInfo)
                         
    return ChainIDs

def ProcessOptions():
    """Process and validate command line arguments and options"""

    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["AlignMethod"] = Options["--alignMethod"].lower()
    OptionsInfo["AlignMode"] = Options["--alignMode"]

    OptionsInfo["Mode"] = Options["--mode"]
    
    OptionsInfo["ProbeFiles"] = Options["--probefiles"]
    OptionsInfo["ProbeFilesNames"] =  Options["--probeFilesNames"]

    OptionsInfo["RefFiles"] = Options["--reffiles"]
    OptionsInfo["RefFilesNames"] =  Options["--refFilesNames"]
    
    RetrieveProbeFilesInfo()
    RetrieveRefFilesInfo()
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutMatrix"] = True if re.match("^Yes$", Options["--outMatrix"], re.I) else False
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    OptionsInfo["OutDelim"] = " "
    if MiscUtil.CheckFileExt(OptionsInfo["Outfile"], "csv"):
        OptionsInfo["OutDelim"] = ","
    elif MiscUtil.CheckFileExt(OptionsInfo["Outfile"], "tsv txt"):
        OptionsInfo["OutDelim"] = "\t"
    else:
        MiscUtil.PrintError("The file name specified , %s, for option \"--outfile\" is not valid. Supported file formats: csv tsv txt\n" % (OptionsInfo["Outfile"]))

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
    
    MiscUtil.ValidateOptionTextValue("-a, --alignMethod", Options["--alignMethod"], "align cealign super")
    MiscUtil.ValidateOptionTextValue("--alignMode", Options["--alignMode"], "FirstChain Complex")

    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "OneToOne AllToAll FirstToAll")
    
    # Expand reffiles to handle presence of multiple input files...
    RefFilesNames = MiscUtil.ExpandFileNames(Options["--reffiles"], ",")

    # Validate file extensions...
    for RefFile in RefFilesNames:
        MiscUtil.ValidateOptionFilePath("-r, --reffiles", RefFile)
        MiscUtil.ValidateOptionFileExt("-r, --reffiles", RefFile, "pdb cif")
    Options["--refFilesNames"] = RefFilesNames

    # Expand probefiles to handle presence of multiple input files...
    ProbeFilesNames = MiscUtil.ExpandFileNames(Options["--probefiles"], ",")

    # Validate file extensions...
    for ProbeFile in ProbeFilesNames:
        MiscUtil.ValidateOptionFilePath("-p, --probefiles", ProbeFile)
        MiscUtil.ValidateOptionFileExt("-p, --probefiles", ProbeFile, "pdb cif")
    Options["--probeFilesNames"] = ProbeFilesNames

    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "csv tsv txt")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    
    MiscUtil.ValidateOptionTextValue("--outMatrix", Options["--outMatrix"], "Yes No")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLCalculateRMSD.py - Calculate RMSD between macromolecules

Usage:
    PyMOLCalculateRMSD.py [--alignMethod <align, cealign, super>]
                          [--alignMode <FirstChain or Complex>] [--mode <OneToOne, AllToAll, FirstToAll>]
                          [--outMatrix <yes or no>] [--overwrite]
                          [-w <dir>] -p <probefile1,probefile2,probefile3...> -r <reffile1,reffile2,reffile3...> -o <outfile>
    PyMOLCalculateRMSD.py -h | --help | -e | --examples

Description:
    Calculate Root Mean Square Distance (RMSD) between a set of similar
    macromolecules in reference and probe input files. The probe and reference
    files are spatially aligned before the the calculation of RMSD values.

    The supported input  file format are: PDB (.pdb), mmCIF (.cif)

    The supported output file formats are:  CSV (.csv), TSV (.tsv, .txt)

Options:
    -a, --alignMethod <align, cealign, super>  [default: super]
        Alignment methodology to use for aligning probe input files to
        reference files.
    --alignMode <FirstChain or Complex>  [default: FirstChain]
        Portion of probe and reference files to use for spatial alignment of
        probe files against reference files.  Possible values: FirstChain or
        Complex.
        
        The FirstChain mode allows alignment of the first chain in probe files
        to the first chain in reference files along with moving the rest of the
        complex to coordinate space of the reference files. The complete
        complex in probe files is aligned to the complete complex in reference
        files for the Complex mode.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -m, --mode <OneToOne, AllToAll, FirstToAll>  [default: OneToOne]
        Specify how reference and probe input files are handled during the calculation
        of RMSD between reference and probe files.  Possible values: OneToOne,
        AllToAll and AllToFirst. For OneToOne mode, the number of reference input
        files must be equal to the number of probe input files. The RMSD is
        calculated for each pair of reference and probe file and written to the
        output file. For AllToAll mode, the RMSD is calculated for each reference
        input file against all probe input files. For FirstToAll mode, however, the RMSD
        is only calculated for the first reference input file against all probe files.
    -p, --probefiles <probefile1,probefile2,probelfile3...>
        A comma delimited list of probe input files. The wildcards are also allowed
        in file names.
    -r, --reffiles <reffile1,reffile2,reffile3...>
        A comma delimited list of reference input files. The wildcards are also allowed
        in file names.
    -o, --outfile <outfile>
        Output file name for writing out RMSD values. Supported text file extensions:
        csv, tsv or txt.
    --outMatrix <yes or no>  [default: yes]
        Output file in a matrix format during 'AllToAll' value for '-m, --mode' option.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To calculate RMSD between pair of macromolecules in reference and probe files
    using only first chain in each file and write out a CSV file containing calculated RMSD
    values along with IDs, type:

        % PyMOLCalculateRMSD.py  -r "Sample3.pdb,Sample4.pdb,Sample5.pdb"
          -p "Sample3.pdb,Sample4.pdb,Sample5.pdb" -o SampleOut.csv 

    To calculate RMSD between all macromolecules in reference and probe files using
    complete complex and write out a CSV matrix file, type:

        % PyMOLCalculateRMSD.py  -m AllToAll --alignMode Complex
           --outMatrix Yes -r "Sample3.pdb,Sample4.pdb,Sample5.pdb"
          -p "Sample3.pdb,Sample4.pdb" -o SampleOut.csv 

    To calculate RMSD between macromolecule in first reference against all probe files
    using only first chain in each file and write out a TSV file containing calculated RMSD
    values along with IDs, type:

        % PyMOLCalculateRMSD.py  -m FirstToAll
          -r "Sample3.pdb,Sample4.pdb,Sample5.pdb"
          -p "Sample3.pdb,Sample4.pdb,Sample5.pdb" -o SampleOut.tsv 

    To calculate RMSD between pair of macromolecules in reference and probe files
    using only first chain in each file along with a specific alignment method and write
    out a CSV file containing calculated RMSD values, type:

        % PyMOLCalculateRMSD.py  --alignMethod align
          -r "Sample3.pdb,Sample4.pdb,Sample5.pdb"
          -p "Sample3.pdb,Sample4.pdb,Sample5.pdb" -o SampleOut.csv 

Author:
    Manish Sud(msud@san.rr.com)

See also:
    PyMOLAlignChains.py, PyMOLSplitChainsAndLigands.py,
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
