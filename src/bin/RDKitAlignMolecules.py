#!/usr/bin/env python
#
# File: RDKitAlignMolecules.py
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

# RDKit imports...
try:
    from rdkit import rdBase
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
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
    AlignMolecules()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def AlignMolecules():
    """Align molecules."""
    
    if not re.match("^(OneToOne|FirstToAll)$", OptionsInfo["Mode"], re.I):
        MiscUtil.PrintError("Alignment couldn't be performed: Specified mode, %s, is not supported" % OptionsInfo["Mode"])
        
    RefFile = OptionsInfo["RefFile"]
    ProbeFile = OptionsInfo["ProbeFile"]
    
    Outfile = OptionsInfo["Outfile"]

    # Read reference and probe molecules...
    OptionsInfo["InfileParams"]["AllowEmptyMols"] = False
    
    MiscUtil.PrintInfo("\nProcessing file %s..." % (RefFile))
    ValidRefMols, RefMolCount, ValidRefMolCount  = RDKitUtil.ReadAndValidateMolecules(RefFile, **OptionsInfo["InfileParams"])
    
    MiscUtil.PrintInfo("Processing file %s..." % (ProbeFile))
    ValidProbeMols, ProbeMolCount, ValidProbeMolCount  = RDKitUtil.ReadAndValidateMolecules(ProbeFile, **OptionsInfo["InfileParams"])

    # Set up a molecule writer...
    Writer = RDKitUtil.MoleculesWriter(OptionsInfo["Outfile"], **OptionsInfo["OutfileParams"])
    if Writer is None:
        MiscUtil.PrintError("Failed to setup a writer for output fie %s " % OptionsInfo["Outfile"])
    MiscUtil.PrintInfo("Generating file %s..." % OptionsInfo["Outfile"])

    AlignmentFailedCount = 0
    if re.match("^OneToOne$", OptionsInfo["Mode"], re.I):
        AlignmentFailedCount = PerformOneToOneAlignment(ValidRefMols, ValidProbeMols, Writer)
    elif re.match("^FirstToAll$", OptionsInfo["Mode"], re.I):
         AlignmentFailedCount = PerformFirstToAllAlignment(ValidRefMols, ValidProbeMols, Writer)
    else:
        MiscUtil.PrintError("Alignment couldn't be performed: Specified mode, %s, is not supported" % OptionsInfo["Mode"])

    if Writer is not None:
        Writer.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: Reference - %d; Probe - %d" % (RefMolCount, ProbeMolCount))
    MiscUtil.PrintInfo("Number of valid molecules: Reference - %d; Probe - %d" % (ValidRefMolCount, ValidProbeMolCount))
    MiscUtil.PrintInfo("Number of probe molecules failed during alignment: %d" % AlignmentFailedCount)
    MiscUtil.PrintInfo("Number of ignored molecules:  Reference - %d; Probe - %d" % ((RefMolCount - ValidRefMolCount), (ProbeMolCount - ValidProbeMolCount + AlignmentFailedCount)))
    
def PerformOneToOneAlignment(ValidRefMols, ValidProbeMols, Writer):
    """Perform pairwise alignment"""
    
    ValidRefMolCount = len(ValidRefMols)
    ValidProbeMolCount = len(ValidProbeMols)
    
    MolCount = ValidRefMolCount
    if ValidRefMolCount > ValidProbeMolCount:
        MolCount = ValidProbeMolCount
    
    if ValidRefMolCount != ValidProbeMolCount:
        MiscUtil.PrintWarning("Number of valid reference molecules, %d,  is not equal to number of valid probe molecules, %d .\n" % (ValidRefMolCount, ValidProbeMolCount))
        MiscUtil.PrintWarning("Pairwise alignment will be performed only for first %s molecules.\n" % (MolCount))

    # Process molecules...
    AlignmentFailedCount = 0
    for MolIndex in range(0, MolCount):
        RefMol = ValidRefMols[MolIndex]
        ProbeMol = ValidProbeMols[MolIndex]

        RefMolName = RDKitUtil.GetMolName(RefMol, (MolIndex + 1))
        ProbeMolName = RDKitUtil.GetMolName(ProbeMol, (MolIndex + 1))

        Status = PerformAlignmentAndWrieOutput(RefMol, ProbeMol, RefMolName, ProbeMolName, Writer)
        if not Status:
            AlignmentFailedCount += 1
    
    return AlignmentFailedCount

def PerformFirstToAllAlignment(ValidRefMols, ValidProbeMols, Writer):
    """Perform alignment between first reference molecues and all probe molecules. """
    
    # Process molecules...
    RefMol = ValidRefMols[0]
    RefMolCount = 1
    RefMolName = RDKitUtil.GetMolName(RefMol, RefMolCount)

    ProbeMolCount = 0
    AlignmentFailedCount = 0
    for ProbeMol in ValidProbeMols:
        ProbeMolCount += 1
        ProbeMolName = RDKitUtil.GetMolName(ProbeMol, ProbeMolCount)
            
        Status = PerformAlignmentAndWrieOutput(RefMol, ProbeMol, RefMolName, ProbeMolName, Writer)
        if not Status:
            AlignmentFailedCount += 1

    return AlignmentFailedCount

def PerformAlignmentAndWrieOutput(RefMol, ProbeMol, RefMolName, ProbeMolName, Writer):
    """Perform alignment and write to output file."""

    Status = True
    try:
        if OptionsInfo["UseRMSD"]:
            RMSD = rdMolAlign.AlignMol(ProbeMol, RefMol, maxIters = OptionsInfo["MaxIters"])
        elif OptionsInfo["UseBestRMSD"]:
            RMSD = AllChem.GetBestRMS(RefMol, ProbeMol)
        elif OptionsInfo["UseOpen3A"]:
            O3A = rdMolAlign.GetO3A(ProbeMol, RefMol)
            Score = O3A.Align()
        elif OptionsInfo["UseCrippenOpen3A"]:
            CrippenO3A = rdMolAlign.GetCrippenO3A(ProbeMol, RefMol)
            Score = CrippenO3A.Align()
        else:
            MiscUtil.PrintError("Alignment couldn't be performed: Specified alignment value, %s, is not supported" % OptionsInfo["Alignment"])
    except (RuntimeError, ValueError):
        Status = False
        MiscUtil.PrintWarning("Alignment failed between reference molecule, %s, and probe molecule, %s.\nWriting unaligned probe molecule...\n" % (RefMolName, ProbeMolName))

    # Write out aligned probe molecule...
    Writer.write(ProbeMol)
    
    return Status

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Alignment"] = Options["--alignment"]
    
    OptionsInfo["UseRMSD"] = False
    OptionsInfo["UseBestRMSD"] = False
    
    OptionsInfo["UseOpen3A"] = False
    OptionsInfo["UseCrippenOpen3A"] = False

    # Process alignment
    if re.match("^RMSD$", OptionsInfo["Alignment"], re.I):
        OptionsInfo["UseRMSD"] = True
    elif re.match("^BestRMSD$", OptionsInfo["Alignment"], re.I):
        OptionsInfo["UseBestRMSD"] = True
    elif re.match("^Open3A$", OptionsInfo["Alignment"], re.I):
        OptionsInfo["UseOpen3A"] = True
    elif re.match("^CrippenOpen3A$", OptionsInfo["Alignment"], re.I):
        OptionsInfo["UseCrippenOpen3A"] = True
    
    OptionsInfo["MaxIters"] = int(Options["--maxIters"])
    
    OptionsInfo["Mode"] = Options["--mode"]
    
    OptionsInfo["RefFile"] = Options["--reffile"]
    OptionsInfo["ProbeFile"] = Options["--probefile"]

    # No need for any RDKit specific --outfileParams....
    OptionsInfo["InfileParams"] = MiscUtil.ProcessOptionInfileParameters("--infileParams", Options["--infileParams"])
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["OutfileParams"] = MiscUtil.ProcessOptionOutfileParameters("--outfileParams", Options["--outfileParams"])
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
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
    
    MiscUtil.ValidateOptionTextValue("--alignment", Options["--alignment"], "Open3A CrippenOpen3A RMSD BestRMSD")
    
    MiscUtil.ValidateOptionIntegerValue("--maxIters", Options["--maxIters"], {">": 0})
    MiscUtil.ValidateOptionTextValue("--mode", Options["--mode"], "OneToOne  FirstToAll")
    
    MiscUtil.ValidateOptionFilePath("-r, --reffile", Options["--reffile"])
    MiscUtil.ValidateOptionFileExt("-r, --reffile", Options["--reffile"], "sdf sd mol")
    
    MiscUtil.ValidateOptionFilePath("-p, --probefile", Options["--probefile"])
    MiscUtil.ValidateOptionFileExt("-p, --probefile", Options["--probefile"], "sdf sd mol")
        
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "sd sdf")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-r, --reffile", Options["--reffile"], "-o, --outfile", Options["--outfile"])
    MiscUtil.ValidateOptionsDistinctFileNames("-p, --probefile", Options["--probefile"], "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionsDistinctFileNames("-r, --reffile", Options["--reffile"], "-o, --outfile", Options["--outfile"])
    MiscUtil.ValidateOptionsDistinctFileNames("-p, --probefile", Options["--probefile"], "-o, --outfile", Options["--outfile"])
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitAlignMolecules.py - Align molecules by RMSD or shape

Usage:
    RDKitAlignMolecules.py [--alignment <Open3A, CrippenOpen3A, RMSD, BestRMSD>]
                           [--infileParams <Name,Value,...>] [--maxIters <number>]
                           [--mode <OneToOne, FirstToAll>] [ --outfileParams <Name,Value,...> ] 
                           [--overwrite] [-w <dir>] -r <reffile> -p <probefile> -o <outfile> 
    RDKitAlignMolecules.py -h | --help | -e | --examples

Description:
    Perform alignment between a set of similar molecules in reference and probe
    input files. The molecules are aligned either by Root Mean Square Distance (RMSD)
    between molecules or overlying their shapes (Open3A or CrippenOpen3A).
    The RDKit function fails to calculate RMSD values for dissimilar molecules. Consequently,
    unaligned probe molecules are written to the output file for dissimilar molecule pairs.

    The supported input file formats are: Mol (.mol), SD (.sdf, .sd)

    The supported output file formats are:  SD (.sdf, .sd)

Options:
    -a, --alignment <Open3A, CrippenOpen3A, RMSD, BestRMSD>  [default: Open3A]
        Alignment methodology to use for aligning molecules. Possible values: Open3A,
        CrippenOpen3A, RMSD, BestRMSD.
        
        The Open3A and CrippenOpen3A allow alignment of molecules using their shapes
        Open 3DAlign (Open3A) [ Ref 132 ] overlays molecules based on MMFF atom types
        and charges. Crippen Open 3DAlign (CrippenOpen3A) uses Crippen logP contributions
        to overlay molecules.
        
        During BestRMSMode mode, the RDKit 'function AllChem.GetBestRMS' is used to
        align and calculate RMSD. This function calculates optimal RMSD for aligning two
        molecules, taking symmetry into account. Otherwise, the RMSD value is calculated
        using 'AllChem.AlignMol function' without changing the atom order. A word to the
        wise from RDKit documentation: The AllChem.GetBestRMS function will attempt to
        align all permutations of matching atom orders in both molecules, for some molecules
        it will lead to 'combinatorial explosion'.
    --infileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for reading
        molecules from files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD, MOL: removeHydrogens,yes,sanitize,yes,strictParsing,yes
            
    --maxIters <number>  [default: 50]
        Maximum number of iterations to perform for each molecule pair during minimization
        of RMSD values. This option is ignored during BestRMSD mode.
    -m, --mode <OneToOne, FirstToAll>  [default: OneToOne]
        Specify how molecules are handled in reference and probe input files during
        alignment of molecules between reference and probe molecules.  Possible values:
        OneToOne and  FirstToAll. For OneToOne mode, the alignment is performed
        for each pair of molecules in the reference and probe file, and the aligned probe
        molecule is written the output file. For FirstToAll mode, the alignment is only
        performed between the first reference molecule against all probe molecules.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -o, --outfile <outfile>
        Output file name for writing out aligned probe molecules values. Supported
        file extensions: sdf or sd.
    --outfileParams <Name,Value,...>  [default: auto]
        A comma delimited list of parameter name and value pairs for writing
        molecules to files. The supported parameter names for different file
        formats, along with their default values, are shown below:
            
            SD: kekulize,no
            
    -p, --probefile <probefile>
        Probe input file name.
    -r, --reffile <reffile>
        Reference input file name.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To perform shape alignment using Open3A methodology between paris of molecules in
    reference and probe input 3D SD files and write out a SD file containing aligned
    molecules, type:

        % RDKitAlignMolecules.py  -r Sample3DRef.sdf -p Sample3DProb.sdf
          -o SampleOut.sdf

    To perform alignment using RMSD methodology between paris of molecules in
    reference and probe input 3D SD files and write out a SD file containing aligned
    molecules, type:

        % RDKitAlignMolecules.py  -a RMSD -r Sample3DRef.sdf -p Sample3DProb.sdf
          -o SampleOut.sdf

    To perform alignment using Open3A methodology  between first reference molecule
    against all probe molecules in 3D SD files without removing hydrogens , and write out
    a SD file containing aligned molecules, type:

        % RDKitAlignMolecules.py -m FirstToAll -a Open3A
          --infileParams "removeHydrogens,no" -r Sample3DRef.sdf
          -p Sample3DProb.sdf -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitCalculateMolecularDescriptors.py, RDKitCompareMoleculeShapes.py, RDKitCalculateRMSD.py,
    RDKitConvertFileFormat.py, RDKitGenerateConformers.py, RDKitPerformMinimization.py

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
