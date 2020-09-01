#!/usr/bin/env python
#
# File: PyMOLMutateAminoAcids.py
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
    PerformMutagenesis()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def PerformMutagenesis():
    """Mutate specified residues across chains and generate an output file."""

    MiscUtil.PrintInfo("\nApplying mutations...")

    # Load macromolecule from input file...
    MolName = OptionsInfo["InfileRoot"]
    LoadMolecule(OptionsInfo["Infile"], MolName)
    
    # Apply mutations...
    for Mutation, ChainID, ResName, ResNum, NewResName in OptionsInfo["SpecifiedMutationsInfo"] :
        ApplyMutation(Mutation, MolName, ChainID, ResName, ResNum, NewResName)

    #  Generate output file...
    Outfile = OptionsInfo["Outfile"]
    MiscUtil.PrintInfo("\nGenerating output file %s..." % Outfile)
    pymol.cmd.save(Outfile, MolName)
    
    # Delete macromolecule...
    DeleteMolecule(MolName)

def ApplyMutation(Mutation, MolName, ChainID, ResName, ResNum, NewResName):
    """Apply mutatation. """

    MiscUtil.PrintInfo("\nApplying mutation %s" % Mutation)
    
    # Setup wizard for amino acid mutagenesis...
    pymol.cmd.wizard("mutagenesis")
    
    pymol.cmd.refresh_wizard()

    # Setup residue to be mutated...
    ResSelection = "/%s//%s/%s" % (MolName, ChainID, ResNum)
    pymol.cmd.get_wizard().do_select(ResSelection)

    # Setup new mutated residue...
    pymol.cmd.get_wizard().set_mode("%s" % NewResName)
    
    # Mutate...
    pymol.cmd.get_wizard().apply()
    
    # Quit wizard...
    pymol.cmd.set_wizard()
    
def RetrieveChainsIDs():
    """Retrieve chain IDs. """

    MolName = OptionsInfo["InfileRoot"]
    Infile = OptionsInfo["Infile"]
    
    MiscUtil.PrintInfo("\nRetrieving chains information for input file %s..." % Infile)

    LoadMolecule(Infile, MolName)

    ChainIDs = PyMOLUtil.GetChains(MolName)
    
    DeleteMolecule(MolName)

    if ChainIDs is None:
        ChainIDs = []

    # Print out chain and ligand IDs...
    ChainInfo = ", ".join(ChainIDs) if len(ChainIDs) else "None"
    MiscUtil.PrintInfo("Chain IDs: %s" % ChainInfo)
                         
    OptionsInfo["ChainIDs"] = ChainIDs
    
def ProcessSpecifiedMutations():
    """Process specified mutations"""
    
    MiscUtil.PrintInfo("\nProcessing specified mutations...")

    SpecifiedMutationsInfo = []

    Mutations = re.sub(" ", "", OptionsInfo["Mutations"])
    MutationsWords = Mutations.split(",")
    if not len(MutationsWords):
        MiscUtil.PrintError("The number of comma delimited mutations specified using \"-m, --mutations\" option, \"%s\",  must be > 0." % (OptionsInfo["Mutations"]))

    # Load macromolecule from input file...
    MolName = OptionsInfo["InfileRoot"]
    LoadMolecule(OptionsInfo["Infile"], MolName)
    
    FirstMutation = True
    CurrentChainID = None
    CanonicalMutationMap = {}
    MutationsCount, ValidMutationsCount = [0] * 2
    
    for Mutation in MutationsWords:
        MutationsCount += 1
        if not len(Mutation):
            MiscUtil.PrintWarning("The mutation, \"%s\", specified using \"-m, --mutations\" option is empty.\nIgnoring mutation..." % (Mutation))
            continue

        CanonicalMutation = Mutation.lower()
        if CanonicalMutation in CanonicalMutationMap:
            MiscUtil.PrintWarning("The mutation, \"%s\", specified using \"-m, --mutations\" option already exist.\nIgnoring mutation..." % (Mutation))
            continue
        CanonicalMutationMap[CanonicalMutation] = Mutation
            
        # Match with a chain ID...
        MatchedResults = re.match(r"^([a-z0-9]+):([a-z]+)([0-9]+)([a-z]+)$", Mutation, re.I)
        if not MatchedResults:
            # Match without a chain ID...
            MatchedResults = re.match(r"^([a-z]+)([0-9]+)([a-z]+)$", Mutation, re.I)
            
        if not MatchedResults:
            MiscUtil.PrintWarning("The format of mutation, \"%s\", specified using \"-m, --mutations\" option is not valid. Supported format: <ChainID>:<ResName><ResNum><ResName> or <ResName><ResNum><ResName>\nIgnoring mutation..." % (Mutation))
            continue

        NumOfMatchedGroups =  len(MatchedResults.groups())
        if NumOfMatchedGroups == 3:
            ResName, ResNum, NewResName = MatchedResults.groups()
        elif NumOfMatchedGroups == 4:
            CurrentChainID, ResName, ResNum, NewResName = MatchedResults.groups()
        else:
            MiscUtil.PrintWarning("The format of mutation, \"%s\", specified using \"-m, --mutations\" option is not valid. Supported format: <ChainID>:<ResName><ResNum><ResName> or <ResName><ResNum><ResName>\nIgnoring mutation..." % (Mutation))
            continue
        
        if FirstMutation:
            FirstMutation = False
            if  CurrentChainID is None:
                MiscUtil.PrintError("The first mutation, \"%s\", specified using \"-m, --mutations\" option must be colon delimited and contain only two values, the first value corresponding to chain ID" % (Mutation))
        
        ResName = ResName.upper()
        NewResName = NewResName.upper()
        
        # Is ResNum and ResName present in input file?
        SelectionCmd = "%s and chain %s and resi %s and resn %s" % (MolName, CurrentChainID, ResNum, ResName)
        ResiduesInfo = PyMOLUtil.GetSelectionResiduesInfo(SelectionCmd)
        if (ResiduesInfo is None) or (not len(ResiduesInfo["ResNames"])):
            MiscUtil.PrintWarning("The residue name, %s, and residue number, %s, in mutation, \"%s\", specified using \"-m, --mutations\" option appears to be missing in input file.\nIgnoring mutation..." % (ResName, ResNum, Mutation))
            continue

        ValidMutationsCount += 1
        
        # Track mutation information...
        SpecifiedMutationsInfo.append([Mutation, CurrentChainID, ResName, ResNum, NewResName])
        
    # Delete macromolecule...
    DeleteMolecule(MolName)

    MiscUtil.PrintInfo("\nTotal number of mutations: %d" % MutationsCount)
    MiscUtil.PrintInfo("Number of valid mutations: %d" % ValidMutationsCount)
    MiscUtil.PrintInfo("Number of ignored mutations: %d" % (MutationsCount - ValidMutationsCount))
    
    if not len(SpecifiedMutationsInfo):
        MiscUtil.PrintError("No valid mutations, \"%s\" specified using \"-m, --mutations\" option." % (OptionsInfo["Mutations"]))
        
    OptionsInfo["SpecifiedMutationsInfo"] = SpecifiedMutationsInfo
    
def LoadMolecule(Infile, MolName):
    """Load input file. """
    
    pymol.cmd.reinitialize()
    pymol.cmd.load(Infile, MolName)
    
def DeleteMolecule(MolName):
    """Delete molecule."""
    
    pymol.cmd.delete(MolName)
    
def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()

    OptionsInfo["Infile"] = Options["--infile"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Infile"])
    OptionsInfo["InfileRoot"] = FileName

    OptionsInfo["Overwrite"] = Options["--overwrite"]
    OptionsInfo["Outfile"] = Options["--outfile"]

    RetrieveChainsIDs()
    
    Mutations = Options["--mutations"]
    if re.match("^None$", Mutations, re.I):
        MiscUtil.PrintError("No mutations specified using \"-m, --mutations\" option.")
    
    OptionsInfo["Mutations"] = Options["--mutations"]
    ProcessSpecifiedMutations()

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
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "pdb")
    
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "pdb")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])

# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLMutateAminoAcids.py - Mutate amino acids

Usage:
    PyMOLMutateAminoAcids.py [--mutations <Spec1,Spec2,...>]
                            [--overwrite] [-w <dir>] -i <infile> -o <outfile>
    PyMOLMutateAminoAcids.py -h | --help | -e | --examples

Description:
    Mutate amino acids in macromolecules. The mutations are performed using
    protein mutagenesis wizard available in PyMOL.
 
    The supported input and output file format is: PDB (.pdb)
 
Options:
    -m, --mutations <Spec1,Spec2,...>  [default: None]
        Comma delimited list of specifications for mutating amino acid residues
        in proteins.
        
        The format of mutation specification is as follows:
        
            <ChainID>:<ResName><ResNum><ResName>,...
        
        A chain ID in the first specification of a mutation is required. It may be
        skipped in subsequent specifications. The most recent chain ID is used
        for the missing chain ID. The first reside name corresponds to the residue
        to be mutated. The second residue name represents the new residue.
        The residue number corresponds to the first residue name and must be
        present in the current chain.
        
        Examples:
        
            E:LEU49CYS, E:SER53TYR
            E:LEU49CYS, SER53TYR
            E:LEU49CYS, SER53TYR, I:TYR7SER, ILE11VAL
              
        The residue names must be valid amino acid names. No validation is
        performed before mutating residues via protein mutagenesis wizard
        available in PyMOL.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infile <infile>
        Input file name.
    -o, --outfile <outfile>
        Output file name.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To mutate a single residue in a specific chain and write a PDB file, type:

        % PyMOLMutateAminoAcids.py -m "I:TYR7SER" -i Sample3.pdb
          -o Sample3Out.pdb

    To mutate multiple residues in a single chain and write a PDB file, type:

        % PyMOLMutateAminoAcids.py -m "I:TYR7SER, ILE11VAL" -i Sample3.pdb
          -o Sample3Out.pdb

    To mutate multiple residues across multiple chains and write a PDB file, type:

        % PyMOLMutateAminoAcids.py -m "E:LEU49CYS,SER53TYR,I:TYR7SER,ILE11VAL"
          -i Sample3.pdb -o Sample3Out.pdb

Author:
    Manish Sud(msud@san.rr.com)

See also:
    DownloadPDBFiles.pl, PyMOLMutateNucleicAcids.py,
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
