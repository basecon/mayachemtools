#!/usr/bin/env python
#
# File: PyMOLSplitChainsAndLigands.py
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
    SplitChainsAndLigands()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def SplitChainsAndLigands():
    """Split input file into output files corresponding to chains and ligands."""

    MiscUtil.PrintInfo("\nGenerating output files...")

    # Load macromolecule from input file...
    MolName = OptionsInfo["InfileRoot"]
    pymol.cmd.load(OptionsInfo["Infile"], MolName)
    
    for ChainID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainIDs"]:
        ChainFile = OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainOutfiles"][ChainID]
        WriteChainFile(MolName, ChainID, ChainFile)
        
        for LigandID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["LigandIDs"][ChainID]:
            LigandFile = OptionsInfo["SpecifiedChainsAndLigandsInfo"]["LigandOutfiles"][ChainID][LigandID]
            WriteLigandFile(MolName, ChainID, LigandID, LigandFile)
    
    # Delete macromolecule...
    pymol.cmd.delete(MolName)

def WriteChainFile(MolName, ChainID, ChainFile):
    """Write chain file."""

    MiscUtil.PrintInfo("\nGenerating output file %s..." % ChainFile)
    
    ChainName = "%s_Chain%s" % (MolName, ChainID)

    ChainSelection = "%s and (chain %s)" % (MolName, ChainID)
    if not OptionsInfo["ChainsMode"]:
        ChainSelection += " and (not organic)"
    
    if not OptionsInfo["KeepSolvents"]:
        ChainSelection += " and (not solvent)"
        
    if not OptionsInfo["KeepInorganics"]:
        ChainSelection += " and (not inorganic)"

    ChainSelection = "(%s)" % ChainSelection
    MiscUtil.PrintInfo("Chain selection: %s" % ChainSelection)
    
    pymol.cmd.create(ChainName, ChainSelection)
    pymol.cmd.save(ChainFile, ChainSelection)
    pymol.cmd.delete(ChainName)
    
    if not os.path.exists(ChainFile):
        MiscUtil.PrintWarning("Failed to generate Chain file, %s..." % (ChainFile))

def WriteLigandFile(MolName, ChainID, LigandID, LigandFile):
    """Write ligand file."""

    MiscUtil.PrintInfo("\nGenerating output file %s..." % LigandFile)
    
    LigandName = "%s_Chain%s_%s" % (MolName, ChainID, LigandID)
    LigandSelection = "(%s and (chain %s) and organic and (resn %s))" % (MolName, ChainID, LigandID)
    MiscUtil.PrintInfo("Ligand selection: %s" % LigandSelection)

    pymol.cmd.create(LigandName, LigandSelection)
    pymol.cmd.save(LigandFile, LigandSelection)
    pymol.cmd.delete(LigandName)
    
    if not os.path.exists(LigandFile):
        MiscUtil.PrintWarning("Failed to generate ligand file, %s..." % (LigandFile))
    
def ProcessChainAndLigandIDs():
    """Process chain and ligand IDs"""
    
    MolName = OptionsInfo["InfileRoot"]
    ChainsAndLigandsInfo = PyMOLUtil.GetChainsAndLigandsInfo(OptionsInfo["Infile"], MolName)
    OptionsInfo["ChainsAndLigandsInfo"] = ChainsAndLigandsInfo
    
    MiscUtil.PrintInfo("\nProcessing specified chain and ligand IDs for input file %s..." % OptionsInfo["Infile"])
    
    SpecifiedChainsAndLigandsInfo = PyMOLUtil.ProcessChainsAndLigandsOptionsInfo(ChainsAndLigandsInfo, "-c, --chainIDs", OptionsInfo["ChainIDs"], "-l, --ligandIDs", OptionsInfo["LigandIDs"])
    OptionsInfo["SpecifiedChainsAndLigandsInfo"] = SpecifiedChainsAndLigandsInfo
    
    CheckPresenceOfValidLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo)
    
def CheckPresenceOfValidLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo):
    """Check presence of valid ligand IDs."""

    MiscUtil.PrintInfo("\nSpecified chain IDs: %s" % (", ".join(SpecifiedChainsAndLigandsInfo["ChainIDs"])))
    
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
        if len (SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID]):
            MiscUtil.PrintInfo("Chain ID: %s; Specified LigandIDs: %s" % (ChainID, ", ".join(SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID])))
        else:
            MiscUtil.PrintInfo("Chain IDs: %s; Specified LigandIDs: None" % (ChainID))
            MiscUtil.PrintWarning("No valid ligand IDs found for chain ID, %s." % (ChainID))

def SetupChainAndLigandOutfiles():
    """Setup output file names for chains and ligands."""

    OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainOutfiles"] = {}
    OptionsInfo["SpecifiedChainsAndLigandsInfo"]["LigandOutfiles"] = {}

    InfileRoot = OptionsInfo["InfileRoot"]
    LigandFileExt = OptionsInfo["LigandFileExt"]
    
    for ChainID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainIDs"]:
        ChainOutfileRoot = "%s_Chain%s" % (InfileRoot, ChainID)
        ChainOutfile = "%s.pdb" % (ChainOutfileRoot)
        OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainOutfiles"][ChainID] = ChainOutfile
        if os.path.exists(ChainOutfile):
            if not OptionsInfo["Overwrite"]:
                MiscUtil.PrintError("The chain output file, %s, already exist. Use option \"--ov\" or \"--overwrite\" and try again.\n" % (ChainOutfile))
        
        OptionsInfo["SpecifiedChainsAndLigandsInfo"]["LigandOutfiles"][ChainID] = {}
        for LigandID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["LigandIDs"][ChainID]:
            LigandOutfile = "%s_%s.%s" % (ChainOutfileRoot, LigandID, LigandFileExt)
            OptionsInfo["SpecifiedChainsAndLigandsInfo"]["LigandOutfiles"][ChainID][LigandID] = LigandOutfile
            if os.path.exists(LigandOutfile):
                if not OptionsInfo["Overwrite"]:
                    MiscUtil.PrintError("The ligand output file, %s, already exist. Use option \"--ov\" or \"--overwrite\" and try again.\n" % (LigandOutfile))
    
def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()

    OptionsInfo["Mode"] = Options["--mode"]
    OptionsInfo["ChainsMode"] = False
    if re.match("^Chains$", OptionsInfo["Mode"], re.I):
        OptionsInfo["ChainsMode"] = True
    
    OptionsInfo["LigandFileFormat"] = Options["--ligandFileFormat"]
    LigandFileExt = "mol"
    if re.match("^PDB$", OptionsInfo["LigandFileFormat"], re.I):
        LigandFileExt = "pdb"
    elif re.match("^(SD|SDF)$", OptionsInfo["LigandFileFormat"], re.I):
        LigandFileExt = "sdf"
    elif re.match("^MOL$", OptionsInfo["LigandFileFormat"], re.I):
        LigandFileExt = "mol"
    OptionsInfo["LigandFileExt"] = LigandFileExt
    
    OptionsInfo["KeepInorganics"] = True if re.match("^Yes$", Options["--keepInorganics"], re.I) else False
    OptionsInfo["KeepSolvents"] = True if re.match("^Yes$", Options["--keepSolvents"], re.I) else False
    
    OptionsInfo["Infile"] = Options["--infile"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Infile"])
    OptionsInfo["InfileRoot"] = FileName

    OptionsInfo["Overwrite"] = Options["--overwrite"]

    OptionsInfo["ChainIDs"] = Options["--chainIDs"]
    OptionsInfo["LigandIDs"] = Options["--ligandIDs"]
    ProcessChainAndLigandIDs()

    SetupChainAndLigandOutfiles()

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
    
    MiscUtil.ValidateOptionTextValue("--ligandFileFormat", Options["--ligandFileFormat"], "PDB SDF SD MDLMOL")
    
    MiscUtil.ValidateOptionTextValue("-m, --mode", Options["--mode"], "Chains ChainsLigands")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "pdb cif")

    MiscUtil.ValidateOptionTextValue("--keepInorganics", Options["--keepInorganics"], "yes no")
    MiscUtil.ValidateOptionTextValue("--keepSolvents", Options["--keepSolvents"], "yes no")

# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLSplitChainsAndLigands.py - Split macromolecule into chains and ligands

Usage:
    PyMOLSplitChainsAndLigands.py [--chainIDs <First, All or ID1,ID2...>]
                                  [--ligandIDs <Largest, All or ID1,ID2...>] [--ligandFileFormat <PDB, SDF, MDLMOL>]
                                  [--mode <Chains or ChainsLigands>] [--keepInorganics <yes or no>]
                                  [--keepSolvents <yes or no>] [--overwrite] [-w <dir>] -i <infile>
    PyMOLSplitChainsAndLigands.py -h | --help | -e | --examples

Description:
    Spit a macromolecule into chains and ligands, and write them out to different
    files. The solvents and inorganic molecules may be optionally removed from
    chains. You may also skip the generation of ligand files and write out a chain
    along with associated ligands into the same chain file.
 
    The supported input file format is:  PDB (.pdb), CIF (.cif)
 
    The supported output file formats are: Chains - PDB (.pdb); Ligands: PDB (.pdb),
    SD file (.sdf, .sd), MDL MOL (.mol)

    The names of the output files are automatically generated from the name of
    input file as shown below:
    
        Chains: <InfileRoot>_<ChainID>.pdb
        Ligands: <InfileRoot>_<ChainID>.{pdb,sdf,sd,mol}
    
Options:
    -c, --chainIDs <First, All or ID1,ID2...>  [default: All]
        List of chain IDs for splitting input file. Possible values: First, All,
        or a comma delimited list of chain IDs. The default is to use
        all chain IDs in input file.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infile <infile>
        Input file name.
    -l, --ligandIDs <Largest, All or ID1,ID2...>  [default: Largest]
        List of ligand IDs present in chains for splitting input file. Possible
        values: Largest, All, or a comma delimited list of ligand IDs. The default
        is to use the largest ligand present in all or specified chains in input file.
        This option is ignored during 'Chains' value of '--mode' option.
        
        Ligands are identified using organic selection operator available in PyMOL.
        It'll also  identify buffer molecules as ligands. The largest ligand contains
        the highest number of heavy atoms.
    --ligandFileFormat <PDB, SDF, MDLMOL>  [default: SDF]
        Ligand file format.
    -m, --mode <Chains or ChainsLigands>  [default: ChainsLigands]
        Split input file into chains or chains and ligands. The ligands are kept
        together chains in the output files for 'Chains' mode. Separate files are
        generated for ligands during 'ChainsAndLigands' mode.
    --keepInorganics <yes or no>  [default: yes]
        Keep inorganic molecules during splitting of input file and write them to
        output files. The inorganic molecules are identified using inorganic selection
        operator available in PyMOL.
    --keepSolvents <yes or no>  [default: yes]
        Keep solvent molecules during splitting of input file and write them to
        output files. The solvent molecules are identified using solvent selection
        operator available in PyMOL.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To split a macromolecule into the first chain and the largest ligand in the
    first chain along with solvent and inorganic molecules, and write chain PDB
    and ligand SDF files, type:

        % PyMOLSplitChainsAndLigands.py -i Sample3.pdb

    To split a macromolecule into all chains and all ligands across all chains
    along with solvent and inorganic molecules, and write out corresponding
    chain and ligand files, type:

        % PyMOLSplitChainsAndLigands.py -i Sample3.pdb -c All -l All

    To split a macromolecule into all chains along with any associated ligands
    without any solvent and inorganic molecules, and write corresponding
    PDB files for chains and skipping generation of any ligand files, type:

        % PyMOLSplitChainsAndLigands.py -c all -m Chains --keepSolvents no
          --keepInorganics no -i Sample3.pdb

    To split a macromolecule into a specific chain and a specific ligand in the
    chain along with solvent and inorganic molecules, and write chain PDB
    and ligand MDLMOL files, type:

        % PyMOLSplitChainsAndLigands.py -c E -l ADP --ligandFileFormat MDLMOL
          -i Sample3.pdb 

Author:
    Manish Sud(msud@san.rr.com)

See also:
    PyMOLAlignChains.py, PyMOLVisualizeMacromolecules.py

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
