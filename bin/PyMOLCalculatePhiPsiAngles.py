#!/bin/env python
#
# File: PyMOLCalculatePhiPsiAngles.py
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
    CalculatePhiPsiAngles()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def CalculatePhiPsiAngles():
    """Calculate phi and psi angles for macromolecules containing amino acids."""
    
    SetupOutputFiles()
    WriteColumnLabels()
    
    Infile = OptionsInfo["Infile"]
    MolName = OptionsInfo["InfileRoot"]
    
    MiscUtil.PrintInfo("\nCalculating phi and psi torsion angles for input file %s..." % Infile)
    
    # Load infile
    pymol.cmd.load(Infile, MolName)

    OutDelim = OptionsInfo["OutDelim"]
    Precision = OptionsInfo["Precision"]

    # Go over specified chain IDs..
    for ChainID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainIDs"]:
        # Write out information for combined file...
        PhiPsiInfo = PyMOLUtil.GetPhiPsiResiduesInfo(MolName, ChainID, Categorize = True)
        OptionsInfo["OutfileResCount"] += len(PhiPsiInfo["ResNums"])
        WritePhiPsiInfo(OptionsInfo["OutFH"], MolName, ChainID, PhiPsiInfo, OutDelim, Precision)

        # Write out information for category fies...
        if OptionsInfo["MultipleOutFiles"]:
            PhiPsiInfoList = []
            GeneralPhiPsiInfo, GlycinePhiPsiInfo, ProlinePhiPsiInfo, PreProlinePhiPsiInfo = PyMOLUtil.GetPhiPsiCategoriesResiduesInfo(MolName, ChainID)
            PhiPsiInfoList.extend([GeneralPhiPsiInfo, GlycinePhiPsiInfo, ProlinePhiPsiInfo, PreProlinePhiPsiInfo])
            
            for Index, Category in enumerate(OptionsInfo["Categories"]):
                OptionsInfo["CategoriesResCount"][Category] += len(PhiPsiInfoList[Index]["ResNums"])
                WritePhiPsiInfo(OptionsInfo["CategoriesOutFHs"][Category], MolName, ChainID, PhiPsiInfoList[Index], OutDelim, Precision)
    
    # Delete MolName object
    pymol.cmd.delete(MolName)

    # Close all files...
    CloseOutputFiles()

    # List number of phi and psi angles in output files...
    MiscUtil.PrintInfo("\nNumber of phi and psi angles in output file %s: %d" % (OptionsInfo["Outfile"],  OptionsInfo["OutfileResCount"]))
    if OptionsInfo["MultipleOutFiles"]:
        MiscUtil.PrintInfo("")
        for Index, Category in enumerate(OptionsInfo["Categories"]):
            MiscUtil.PrintInfo("Number of phi and psi angles in output file %s: %d" % (OptionsInfo["CategoriesOutfiles"][Category], OptionsInfo["CategoriesResCount"][Category]))

def WritePhiPsiInfo(OutFH, MolName, ChainID, PhiPsiResiduesInfo, OutDelim, Precision):
    """Write out phi and psi information."""

    for ResNum in PhiPsiResiduesInfo["ResNums"]:
        ResName = PhiPsiResiduesInfo["ResName"][ResNum]
        Phi = "%.*f" % (Precision, PhiPsiResiduesInfo["Phi"][ResNum])
        Psi = "%.*f" % (Precision, PhiPsiResiduesInfo["Psi"][ResNum])
        Category = PhiPsiResiduesInfo["Category"][ResNum]

        LineWords = []
        if OptionsInfo["OutChainID"]:
            LineWords.append(ChainID)
        
        LineWords.extend([ResNum, ResName, Phi, Psi])
        
        if OptionsInfo["OutCategory"]:
            LineWords.append(Category)
            
        Line = OutDelim.join(LineWords)
        
        OutFH.write("%s\n" % Line)

def WriteColumnLabels():
    """Write out column labels. """
    
    OutDelim = OptionsInfo["OutDelim"]

    ColLabels = []
    if OptionsInfo["OutChainID"]:
        ColLabels.append("ChainID")
    
    ColLabels.extend(["ResNum", "ResName", "Phi", "Psi"])
    
    if OptionsInfo["OutCategory"]:
        ColLabels.append("Category")
    
    Line = OutDelim.join(ColLabels)

    # Write column labels for combined file...
    OutFH = OptionsInfo["OutFH"]
    OutFH.write("%s\n" % Line)

    if not OptionsInfo["MultipleOutFiles"]:
        return
    
    # Write columns labels for category files...
    for Category in OptionsInfo["Categories"]:
        CategoryOutFH = OptionsInfo["CategoriesOutFHs"][Category]
        CategoryOutFH.write("%s\n" % Line)

def SetupOutputFiles():
    """Open output files."""

    if OptionsInfo["MultipleOutFiles"]:
        MiscUtil.PrintInfo("\nGenerating output files: %s" % (", ".join(OptionsInfo["OutfilesList"])))
    else:
        MiscUtil.PrintInfo("\nGenerating output file %s..." % (OptionsInfo["Outfile"]))
        
    # Open combined output file...
    Outfile = OptionsInfo["Outfile"]
    OutFH = open(Outfile, "w")
    if OutFH is None:
        MiscUtil.PrintError("Couldn't open output file: %s.\n" % (Outfile))
    OptionsInfo["OutFH"] = OutFH
    OptionsInfo["OutfileResCount"] = 0

    if not OptionsInfo["MultipleOutFiles"]:
        return
    
    # Open output files for different categories...
    OptionsInfo["CategoriesOutFHs"] = {}
    OptionsInfo["CategoriesResCount"] = {}
    for Category in OptionsInfo["Categories"]:
        CategoryOutfile = OptionsInfo["CategoriesOutfiles"][Category]
        CategoryOutFH = open(CategoryOutfile, "w")
        if CategoryOutfile is None:
            MiscUtil.PrintError("Couldn't open output file: %s.\n" % (CategoryOutfile))
        
        OptionsInfo["CategoriesOutFHs"][Category] = CategoryOutFH
        OptionsInfo["CategoriesResCount"][Category] = 0
    
def CloseOutputFiles():
    """Close output files."""
    
    OptionsInfo["OutFH"].close()
    
    if not OptionsInfo["MultipleOutFiles"]:
        return
    
    for Category in OptionsInfo["Categories"]:
        CategoryOutFH = OptionsInfo["CategoriesOutFHs"][Category]
        CategoryOutFH.close()

def RetrieveInfileInfo():
    """Retrieve information for input file."""
    
    Infile = OptionsInfo["Infile"]
    InfileRoot = OptionsInfo["InfileRoot"]
    
    ChainsAndLigandsInfo = PyMOLUtil.GetChainsAndLigandsInfo(Infile, InfileRoot)
    OptionsInfo["ChainsAndLigandsInfo"] = ChainsAndLigandsInfo

def ProcessChainIDs():
    """Process specified chain IDs for infile."""
    
    MiscUtil.PrintInfo("\nProcessing specified chain IDs for input file %s..." % OptionsInfo["Infile"])        
    ChainsAndLigandsInfo = OptionsInfo["ChainsAndLigandsInfo"]
    SpecifiedChainsAndLigandsInfo = PyMOLUtil.ProcessChainsAndLigandsOptionsInfo(ChainsAndLigandsInfo, "-c, --chainIDs", OptionsInfo["ChainIDs"], None, None)
    
    OptionsInfo["SpecifiedChainsAndLigandsInfo"] = SpecifiedChainsAndLigandsInfo
    
    MiscUtil.PrintInfo("Specified chain IDs: %s" % (", ".join(SpecifiedChainsAndLigandsInfo["ChainIDs"])))

def SetupCategoryOutfiles():
    """Setup output file names for different categories of phi and psi angles. """

    # Initialize...
    OptionsInfo["OutfilesList"] = []
    OptionsInfo["OutfilesList"].append(OptionsInfo["Outfile"])
    
    OptionsInfo["Categories"] = ["General", "Glycine", "Proline", "PreProline"]
    OptionsInfo["CategoriesOutfiles"] = {}
    for Category in OptionsInfo["Categories"]:
        OptionsInfo["CategoriesOutfiles"][Category] = None

    if not OptionsInfo["MultipleOutFiles"]:
        return

    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Outfile"])
    OutfileRoot = FileName
    OutfileExt = FileExt
    
    for Category in OptionsInfo["Categories"]:
        CategoryOutfile = "%s_%s.%s" % (OutfileRoot, Category, OutfileExt)
        if os.path.exists(CategoryOutfile):
            if not OptionsInfo["Overwrite"]:
                MiscUtil.PrintError("The category output file, %s, already exist. Use option \"--ov\" or \"--overwrite\" and try again.\n" % (CategoryOutfile))
                
        OptionsInfo["CategoriesOutfiles"][Category] = CategoryOutfile
        OptionsInfo["OutfilesList"].append(CategoryOutfile)
        
def ProcessOptions():
    """Process and validate command line arguments and options"""

    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()

    OptionsInfo["Infile"] = Options["--infile"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Infile"])
    OptionsInfo["InfileRoot"] = FileName
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Outfile"])
    OptionsInfo["OutfileRoot"] = FileName
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    OptionsInfo["OutDelim"] = " "
    if MiscUtil.CheckFileExt(OptionsInfo["Outfile"], "csv"):
        OptionsInfo["OutDelim"] = ","
    elif MiscUtil.CheckFileExt(OptionsInfo["Outfile"], "tsv txt"):
        OptionsInfo["OutDelim"] = "\t"
    else:
        MiscUtil.PrintError("The file name specified , %s, for option \"--outfile\" is not valid. Supported file formats: csv tsv txt\n" % (OptionsInfo["Outfile"]))
    
    OptionsInfo["OutMode"] = Options["--outMode"]
    OptionsInfo["MultipleOutFiles"] = True if re.match("^MultipleFiles$", OptionsInfo["OutMode"], re.I) else False
    
    OptionsInfo["OutChainID"] = True if re.match("^Yes$", Options["--outChainID"], re.I) else False
    OptionsInfo["OutCategory"] = True if re.match("^Yes$", Options["--outCategory"], re.I) else False
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    OptionsInfo["Precision"] = int(Options["--precision"])
    
    RetrieveInfileInfo()
    
    OptionsInfo["ChainIDs"] = Options["--chainIDs"]
    ProcessChainIDs()
    
    SetupCategoryOutfiles()

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
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "pdb cif")

    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "csv tsv txt")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    
    MiscUtil.ValidateOptionTextValue("--outMode", Options["--outMode"], "SingleFile MultipleFiles")
    MiscUtil.ValidateOptionTextValue("--outChainID", Options["--outChainID"], "yes no")
    MiscUtil.ValidateOptionTextValue("--outCategory", Options["--outCategory"], "yes no")
    MiscUtil.ValidateOptionIntegerValue("-p, --precision", Options["--precision"], {">": 0})
    
# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLCalculatePhiPsiAngles.py - Calculate phi and psi torsion angles

Usage:
    PyMOLCalculatePhiPsiAngles.py [--chainIDs <First, All or ID1,ID2...>]
                                  [--outMode <SingleFile or MultipleFies>] [--outChainID <yes or no>]
                                  [--outCategory <yes or no>] [--overwrite] [--precision <number>]
                                  [-w <dir>] -i <infile> -o <outfile>
    PyMOLCalculatePhiPsiAngles.py -h | --help | -e | --examples

Description:
    Calculate phi and psi torsion angels for amino acid residues present
    in macromolecules.

    The phi and psi angles are categorized into the following groups
    corresponding to four types of Ramachandran plots:
    
        General: All residues except glycine, proline, or pre-proline
        Glycine: Only glycine residues
        Proline: Only proline residues
        Pre-Proline: Only residues before proline not including glycine or
            proline
    
    The supported input  file format are: PDB (.pdb), mmCIF (.cif)

    The supported output file formats are:  CSV (.csv), TSV (.tsv, .txt)

Options:
    -c, --chainIDs <First, All or ID1,ID2...>  [default: All]
        List of chain IDs to use for calculating phi and psi angles for residues
        in chains. Possible values: First, All, or a comma delimited list of chain
        IDs. The default is to use all chain IDs in input file.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infile <infile>
        Input file name.
    -o, --outfile <outfile>
        Output file name for writing out calculated values. Supported text file
        extensions: csv, tsv or txt.
        
        In addition to the specified outfile containing phi and psi angles for all
        residues, a set of additional output files is generated for 'MultipleFiles'
        value of '--outMode' option. The names of these output files are
        automatically generated from the the name of the specified output
        file as shown below:
        
            General: <OutfileRoot>_General.<OutfileExt>
            Glycine: <OutfileRoot>_Glycine.<OutfileExt>
            Proline: <OutfileRoot>_Proline.<OutfileExt>
            Pre-Proline: <OutfileRoot>_PreProline.<OutfileExt>
        
    --outMode <SingleFile or MultipleFiles>  [default: SingleFile]
        A single output file containing phi and psi angles for all residues or
        multiple output files corresponding to different categories of angles.
        
        The phi and psi angles are categorized into the following groups
        corresponding to four types of Ramachandran plots:
        
            General: All residues except glycine, proline, or pre-proline
            Glycine: Only glycine residues
            Proline: Only proline residues
            Pre-Proline: Only residues before proline not including glycine or
                proline
        
        The output files contain the following information:
        
            ChainID ResNum ResName Phi Psi Category
        
    --outChainID <yes or no>  [default: yes]
        Write chain IDs to output file.
    --outCategory <yes or no>  [default: yes]
        Write phi and psi category to output file.
    --overwrite
        Overwrite existing files.
    -p, --precision <number>  [default: 2]
        Floating point precision for writing the calculated phi and psi angles.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To calculate phi and psi angles for all residues across all chains in input
    file and write out a single CSV file containing calculated values along with
    chain IDs, residue names and numbers, and category of angles corresponding
    to Ramachandran plots, type:

        % PyMOLCalculatePhiPsiAngles.py -i Sample3.pdb -o Sample3Out.csv

    To calculate phi and psi angles for all residues across all chains in input
    file and write out a multiple CSV files corresponding to categories of angles
    for Ramachandran plots along with other relevant information, type:

        % PyMOLCalculatePhiPsiAngles.py --outMode MultipleFiles -i Sample3.pdb
          -o Sample3Out.csv

    To calculate phi and psi angles for all residues in a specific chain in input
    file and write out a single TSV file containing calculated values along with
    other relevant information, type:

        % PyMOLCalculatePhiPsiAngles.py -c E  -i Sample3.pdb -o Sample3Out.csv

    To calculate phi and psi angles for all residues in a specific chain in input
    file and write out a multiple TSV files containing calculated values at a specific
    precision along with other relevant information, type:

        % PyMOLCalculatePhiPsiAngles.py --outMode MultipleFiles --chainIDs I
          -i Sample3.pdb -o Sample3Out.csv

Author:
    Manish Sud(msud@san.rr.com)

See also:
    DownloadPDBFiles.pl, PyMOLCalculateRMSD.py, PyMOLCalculateProperties.py,
    PyMOLGenerateRamachandranPlots.py

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
