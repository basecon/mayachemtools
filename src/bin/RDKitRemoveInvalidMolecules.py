#!/usr/bin/env python
#
# File: RDKitRemoveInvalidMolecules.py
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
    RemoveInvalidMolecules()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def RemoveInvalidMolecules():
    """Identify and remove invalid molecules."""
    
    Infile = OptionsInfo["Infile"]
    Outfile = OptionsInfo["Outfile"]
    InvalidMolsOutfile = OptionsInfo["InvalidMolsOutfile"]
    
    CountMode = OptionsInfo["CountMode"]
    
    # Setup a molecule reader...
    MiscUtil.PrintInfo("\nProcessing file %s..." % Infile)
    Mols  = RDKitUtil.ReadMolecules(Infile, **OptionsInfo["InfileParams"])
    
    Writer = None
    MolNumWriter = None
    if not CountMode:
        # Set up a molecule writer...
        Writer = RDKitUtil.MoleculesWriter(Outfile, **OptionsInfo["OutfileParams"])
        if Writer is None:
            MiscUtil.PrintError("Failed to setup a writer for output fie %s " % Outfile)
    
        # Set up a invalid molecule number writer...
        InvalidMolsWriter = open(InvalidMolsOutfile, "w")
        if InvalidMolsWriter is None:
            MiscUtil.PrintError("Failed to open output fie %s " % InvalidMolsOutfile)
        InvalidMolsWriter.write("MolName\n")
    
        MiscUtil.PrintInfo("Generating files %s and %s..." % (Outfile, InvalidMolsOutfile))

    # Process molecules...
    MolCount = 0
    ValidMolCount = 0
    
    Compute2DCoords = OptionsInfo["OutfileParams"]["Compute2DCoords"]
    SetSMILESMolProps = OptionsInfo["OutfileParams"]["SetSMILESMolProps"]

    FirstMol = True
    for Mol in Mols:
        MolCount += 1
        
        if Mol is None:
            MolName = "Mol%s" % MolCount
            MiscUtil.PrintWarning("Ignoring invalid molecule:  %s" % MolName)
            if not CountMode:
                InvalidMolsWriter.write("%s\n" % MolName)
            continue
        
        if RDKitUtil.IsMolEmpty(Mol):
            MolName = RDKitUtil.GetMolName(Mol, MolCount)
            MiscUtil.PrintWarning("Ignoring invalid empty molecule: %s" % MolName)
            if not CountMode:
                InvalidMolsWriter.write("%s\n" % MolName)
            continue
        ValidMolCount += 1
        
        if FirstMol:
            FirstMol = False
            if not CountMode:
                if SetSMILESMolProps:
                    RDKitUtil.SetWriterMolProps(Writer, Mol)
        
        if Compute2DCoords:
            if not CountMode:
                AllChem.Compute2DCoords(Mol)
        
        if not CountMode:
            Writer.write(Mol)
    
    if Writer is not None:
        Writer.close()
    
    MiscUtil.PrintInfo("\nTotal number of molecules: %d" % MolCount)
    MiscUtil.PrintInfo("Number of valid molecules: %d" % ValidMolCount)
    MiscUtil.PrintInfo("Number of invalid molecules: %d" % (MolCount - ValidMolCount))

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

    OptionsInfo["CountMode"] = True if re.match("^count$", Options["--mode"], re.I) else False
    
    # Setup outfile for writing out molecule number for invalid molecules...
    OptionsInfo["InvalidMolsOutfile"] = ""
    if not OptionsInfo["CountMode"] :
        FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Outfile"])
        OptionsInfo["InvalidMolsOutfile"] = "%sInvalidMols.csv" % (FileName)

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
    
# Setup a usage string for docopt...
_docoptUsage_ = """
RDKitRemoveInvalidMolecules.py - Remove invalid molecules

Usage:
    RDKitRemoveInvalidMolecules.py  [--infileParams <Name,Value,...>]
                              [--mode <remove or count>] [ --outfileParams <Name,Value,...> ] 
                              [--overwrite] [-w <dir>] [-o <outfile>]  -i <infile>
    RDKitRemoveInvalidMolecules.py -h | --help | -e | --examples

Description:
    Identify and remove invalid molecules based on success or failure of RDKit molecule
    readers or simply count the number of invalid molecules.

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
        Specify whether to remove invalid molecules and write out filtered molecules
        to output file or or simply count the number of invalid molecules.
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
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To remove invalid molecules and generate an output file SMILES file
    containing valid molecules, type:

        % RDKitRemoveInvalidMolecules.py -i Sample.smi -o SampleOut.smi

    To count number of valid and invaid molecules without generating any
    output file, type:

        % RDKitRemoveInvalidMolecules.py -m count -i Sample.sdf

    To remove invalid molecules from a CSV SMILES file, SMILES strings in
    column 1, name in column 2, and generate output SD file containing valid
    molecules, type:

        % RDKitRemoveInvalidMolecules.py --infileParams 
          "smilesDelimiter,comma,smilesTitleLine,yes,smilesColumn,1,
          smilesNameColumn,2" --outfileParams "compute2DCoords,yes"
          -i SampleSMILES.csv -o SampleOut.sdf

Author:
    Manish Sud(msud@san.rr.com)

See also:
    RDKitConvertFileFormat.py, RDKitRemoveDuplicateMolecules.py,
    RDKitRemoveSalts, RDKitSearchFunctionalGroups.py,
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
