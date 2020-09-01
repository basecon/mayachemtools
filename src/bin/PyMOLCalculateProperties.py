#!/bin/env python
#
# File: PyMOLCalculateProperties.py
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
    CalculatePhysicochemicalProperties()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def CalculatePhysicochemicalProperties():
    """Calculate physicochemical properties for macromolecules."""

    Outfile = OptionsInfo["Outfile"]
    OutDelim = OptionsInfo["OutDelim"]
    Precision = OptionsInfo["Precision"]
    
    MiscUtil.PrintInfo("\nGenerating file %s...\n" % Outfile)
    OutFH = open(Outfile, "w")
    if OutFH is None:
        MiscUtil.PrintError("Couldn't open output file: %s.\n" % (Outfile))

    WriteColumnLabels(OutFH, OutDelim)
    
    InfilesInfo = OptionsInfo["InfilesInfo"]
    for FileIndex in range(0, len(InfilesInfo["InfilesNames"])):
        MiscUtil.PrintInfo("Calculating properties for input file %s..." % OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex])
        
        LoadInfile(FileIndex)
        AddHydrogens(FileIndex)

        CalculatedValues = CalculatePropertyValues(FileIndex)
        WriteCalculatedValues(FileIndex, OutFH, OutDelim, CalculatedValues)

        # Delete MolName object
        DeleteInfileObject(FileIndex)
    
    OutFH.close()

def CalculatePropertyValues(FileIndex):
    """Calculate property values."""

    MolName = OptionsInfo["InfilesInfo"]["InfilesRoots"][FileIndex]
    Precision = OptionsInfo["Precision"]
    
    CalculatedValues = []
    for PropertyName in OptionsInfo["SpecifiedPropertyNames"]:
        PropertyValue = GetFormattedPropertyValue(MolName, PropertyName)
        CalculatedValues.append(PropertyValue)

    return CalculatedValues

def GetFormattedPropertyValue(Selection, Name):
    """Calculate and return a formatted property value. """

    Quiet = OptionsInfo["Quiet"]
    Precision = OptionsInfo["Precision"]
    
    Value = None
    if re.match("^CenterOfMass$", Name, re.I):
        Value = PyMOLUtil.CalculateCenterOfMass(Selection, Quiet)
    elif re.match("^MolecularWeight$", Name, re.I):
        Value = pymol.util.compute_mass(Selection, implicit = False, quiet = Quiet)
    elif re.match("^MolecularSurfaceArea$", Name, re.I):
        Value = pymol.util.get_area(Selection, -1, 0, quiet = Quiet)
    elif re.match("^SumOfFormalCharges$", Name, re.I):
        Value = pymol.util.sum_formal_charges(Selection, quiet = Quiet)
    elif re.match("^SumOfPartialCharges$", Name, re.I):
        Value = pymol.util.sum_partial_charges(Selection, quiet = Quiet)
    elif re.match("^SolventAccessibleSurfaceArea$", Name, re.I):
        Value = pymol.util.get_sasa(Selection, quiet = Quiet)
    else:
        MiscUtil.PrintError("The property name specified, %s, using \"-m, --mode\" option is not a valid name." % Name)

    if Value is None:
        FormattedValue = "NA"
    else:
        if type(Value) is list:
            FormattedValues = []
            for ListElement in Value:
                FormattedListElement = "%.*f" % (Precision, ListElement)
                FormattedValues.append(FormattedListElement)
            FormattedValue = " ".join(FormattedValues)
        else:
            FormattedValue = "%.*f" % (Precision, Value)
        
    return FormattedValue

def WriteCalculatedValues(FileIndex, OutFH, OutDelim, CalculatedValues):
    """Write out calculated values. """

    PDBID = OptionsInfo["InfilesInfo"]["InfilesRoots"][FileIndex]
    LineWords = [PDBID]

    LineWords.extend(CalculatedValues)
    
    Line = OutDelim.join(LineWords)
    OutFH.write("%s\n" % Line)
    
def WriteColumnLabels(OutFH, OutDelim):
    """Write out column labels. """

    ColLabels = []

    ColLabels = ["PDBID"]
    ColLabels.extend(OptionsInfo["SpecifiedPropertyNames"])

    Line = OutDelim.join(ColLabels)
    OutFH.write("%s\n" % Line)
    
def LoadInfile(FileIndex):
    """Load a file. """

    Infile = OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex]
    MolName = OptionsInfo["InfilesInfo"]["InfilesRoots"][FileIndex]

    ChainSelections = OptionsInfo["InfilesInfo"]["ChainSelections"][FileIndex]
    NonChainSelections = OptionsInfo["NonChainSelections"]
    
    if ChainSelections is None and NonChainSelections is None:
        pymol.cmd.load(Infile, MolName)
        return

    TmpMolName = "Tmp%s" % MolName
    pymol.cmd.load(Infile, TmpMolName)

    MolSelections = []
    MolSelections.append(TmpMolName)
    if ChainSelections is not None:
        MolSelections.append(ChainSelections)
    if NonChainSelections is not None:
        MolSelections.append(NonChainSelections)
    
    MolSelection = " and ".join(MolSelections)
    MolSelection = "(%s)" % MolSelection
    pymol.cmd.create(MolName, MolSelection)

    pymol.cmd.delete(TmpMolName)
    
def DeleteInfileObject(FileIndex):
    """Delete PyMOL object. """
    
    MolName = OptionsInfo["InfilesInfo"]["InfilesRoots"][FileIndex]
    
    pymol.cmd.delete(MolName)

def AddHydrogens(FileIndex):
    """Add hydrogens."""

    if not OptionsInfo["Addhydrogens"]:
        return
    
    MolName = OptionsInfo["InfilesInfo"]["InfilesRoots"][FileIndex]

    pymol.cmd.h_add(MolName)
    pymol.cmd.sort("%s extend 1" % MolName)
    
def ProcessSpecifiedPropertyNames():
    """Process specified property names. """

    PropertyNames = RetrievePropertyNames()
    
    OptionsInfo["SpecifiedPropertyNames"] = []
    
    SpecifiedNames = re.sub(" ", "", OptionsInfo["Mode"])
    if not SpecifiedNames:
        MiscUtil.PrintError("No valid property names specifed  using \"-m, --mode\" option")
    
    if re.match("^All$", SpecifiedNames, re.I):
        OptionsInfo["SpecifiedPropertyNames"] = PropertyNames
        return

    # Validate propery names...
    CanonicalPropertyNamesMap = {}
    for Name in PropertyNames:
        CanonicalPropertyNamesMap[Name.lower()] = Name

    SpecifiedNamesWords = SpecifiedNames.split(",")
    for Name in SpecifiedNamesWords:
        CanonicalName = Name.lower()
        if CanonicalName not in CanonicalPropertyNamesMap:
            MiscUtil.PrintError("The property name specified, %s, using \"-m, --mode\" option is not a valid name." % Name)
        
        PropertyName = CanonicalPropertyNamesMap[CanonicalName]
        OptionsInfo["SpecifiedPropertyNames"].append(PropertyName)

def ProcessListPropertyNames():
    """List available property names. """
    
    PropertyNames = RetrievePropertyNames()
    
    MiscUtil.PrintInfo("\nListing available property names...")
    Delimiter = ", "
    MiscUtil.PrintInfo("\n%s" % (Delimiter.join(PropertyNames)))
    
    MiscUtil.PrintInfo("")

def RetrievePropertyNames():
    """Retrieve available property names. """
    
    PropertyNames = ["CenterOfMass", "MolecularWeight", "MolecularSurfaceArea", "SumOfFormalCharges", "SumOfPartialCharges", "SolventAccessibleSurfaceArea"]
    
    return PropertyNames

def RetrieveInfilesInfo():
    """Retrieve information for input files."""

    InfilesInfo = {}
    
    InfilesInfo["InfilesNames"] = []
    InfilesInfo["InfilesRoots"] = []
    InfilesInfo["ChainsAndLigandsInfo"] = []
    
    for Infile in OptionsInfo["InfilesNames"]:
        FileDir, FileName, FileExt = MiscUtil.ParseFileName(Infile)
        InfileRoot = FileName
        ChainsAndLigandInfo = PyMOLUtil.GetChainsAndLigandsInfo(Infile, InfileRoot)
        
        InfilesInfo["InfilesNames"].append(Infile)
        InfilesInfo["InfilesRoots"].append(InfileRoot)
        InfilesInfo["ChainsAndLigandsInfo"].append(ChainsAndLigandInfo)
    
    OptionsInfo["InfilesInfo"] = InfilesInfo

def ProcessChainIDs():
    """Process specified chain IDs for infiles."""
    
    OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"] = []
    OptionsInfo["InfilesInfo"]["ChainSelections"] = []
    
    for FileIndex in range(0, len(OptionsInfo["InfilesInfo"]["InfilesNames"])):
        MiscUtil.PrintInfo("\nProcessing specified chain IDs for input file %s..." % OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex])
        
        ChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["ChainsAndLigandsInfo"][FileIndex]
        SpecifiedChainsAndLigandsInfo = PyMOLUtil.ProcessChainsAndLigandsOptionsInfo(ChainsAndLigandsInfo, "-c, --chainIDs", OptionsInfo["ChainIDs"], None, None)

        # Setup chain selections...
        ChainSelections = None
        if not OptionsInfo["AllChains"]:
            Chains = []
            for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
                Chains.append("chain %s" % ChainID)
            ChainSelections = " or ".join(Chains)
            ChainSelections = "(%s)" % ChainSelections
        
        OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"].append(SpecifiedChainsAndLigandsInfo)
        OptionsInfo["InfilesInfo"]["ChainSelections"].append(ChainSelections)
        
        MiscUtil.PrintInfo("Specified chain IDs: %s" % (", ".join(SpecifiedChainsAndLigandsInfo["ChainIDs"])))
        
def ProcessKeepSelectionOptions():
    """Process keep selection options."""

    KeepSelections = []
    if not OptionsInfo["KeepSolvents"]:
        KeepSelections.append("(not solvent)")
    if  not OptionsInfo["KeepInorganics"]:
        KeepSelections.append("(not inorganic)")
    if not OptionsInfo["KeepLigands"]:
        KeepSelections.append("(not organic)")

    NonChainSelections = None
    if len(KeepSelections):
        NonChainSelections  = " and ".join(KeepSelections)
        NonChainSelections  = "(%s)" % NonChainSelections
    
    OptionsInfo["NonChainSelections"] = NonChainSelections
    
def ProcessOptions():
    """Process and validate command line arguments and options"""

    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()

    OptionsInfo["Addhydrogens"] = True if re.match("^Yes$", Options["--addHydrogens"], re.I) else False
    
    OptionsInfo["Infiles"] = Options["--infiles"]
    OptionsInfo["InfilesNames"] =  Options["--infilesNames"]
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    
    OptionsInfo["OutDelim"] = " "
    if MiscUtil.CheckFileExt(OptionsInfo["Outfile"], "csv"):
        OptionsInfo["OutDelim"] = ","
    elif MiscUtil.CheckFileExt(OptionsInfo["Outfile"], "tsv txt"):
        OptionsInfo["OutDelim"] = "\t"
    else:
        MiscUtil.PrintError("The file name specified , %s, for option \"--outfile\" is not valid. Supported file formats: csv tsv txt\n" % (OptionsInfo["Outfile"]))
    
    OptionsInfo["KeepInorganics"] = True if re.match("^Yes$", Options["--keepInorganics"], re.I) else False
    OptionsInfo["KeepLigands"] = True if re.match("^Yes$", Options["--keepLigands"], re.I) else False
    OptionsInfo["KeepSolvents"] = True if re.match("^Yes$", Options["--keepSolvents"], re.I) else False
    ProcessKeepSelectionOptions()

    OptionsInfo["Overwrite"] = Options["--overwrite"]
    OptionsInfo["Quiet"] = 1 if re.match("^Yes$", Options["--quiet"], re.I) else 0
    
    OptionsInfo["Precision"] = int(Options["--precision"])
    
    OptionsInfo["Mode"] = Options["--mode"]
    ProcessSpecifiedPropertyNames()
    
    RetrieveInfilesInfo()
    OptionsInfo["ChainIDs"] = Options["--chainIDs"]
    OptionsInfo["AllChains"] = True if re.match("^All$", Options["--chainIDs"], re.I) else False
    ProcessChainIDs()

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
    
    # Handle listing of property names...
    if  Options["--list"]:
        ProcessListPropertyNames()
        sys.exit(0)

def ValidateOptions():
    """Validate option values"""

    MiscUtil.ValidateOptionTextValue("-a, -addHydrogens", Options["--addHydrogens"], "yes no")
    MiscUtil.ValidateOptionTextValue("--keepInorganics", Options["--keepInorganics"], "yes no")
    MiscUtil.ValidateOptionTextValue("--keepLigands", Options["--keepLigands"], "yes no")
    MiscUtil.ValidateOptionTextValue("--keepSolvents", Options["--keepSolvents"], "yes no")
    
    # Expand infile names..
    InfilesNames = MiscUtil.ExpandFileNames(Options["--infiles"], ",")

    # Validate file extensions...
    for Infile in InfilesNames:
        MiscUtil.ValidateOptionFilePath("-i, --infiles", Infile)
        MiscUtil.ValidateOptionFileExt("-i, --infiles", Infile, "pdb cif")
    Options["--infilesNames"] = InfilesNames

    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "csv tsv txt")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    
    MiscUtil.ValidateOptionIntegerValue("-p, --precision", Options["--precision"], {">": 0})
    MiscUtil.ValidateOptionTextValue("--quiet", Options["--quiet"], "yes no")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLCalculateProperties.py - Calculate physicochemical properties

Usage:
    PyMOLCalculateProperties.py  [--addHydrogens <yes or no>]
                                 [--chainIDs <First, All or ID1,ID2...>] [--list] [--keepInorganics <yes or no>]
                                 [--keepLigands <yes or no>] [--keepSolvents <yes or no>]
                                 [--mode <All or Name1,Name2,Name3,...>]
                                 [--overwrite]  [--precision <number>] [--quiet <yes or no>]
                                 [-w <dir>] -i <infile1,infile2,infile3...> -o <outfile>
    PyMOLCalculateProperties.py -l | --list
    PyMOLCalculateProperties.py -h | --help | -e | --examples

Description:
    Calculate physicochemical properties for macromolecules. The properties may
    be calculated for the complete complex or a specified list of chain IDs. Ligands,
    inorganics, and solvents may be optionally excluded during the calculation
    of properties.

    The supported input  file format are: PDB (.pdb), mmCIF (.cif)

    The supported output file formats are:  CSV (.csv), TSV (.tsv, .txt)

Options:
    -a, --addHydrogens <yes or no>  [default: yes]
        Add hydrogens before calculating physiochemical properties.
    -c, --chainIDs <First, All or ID1,ID2...>  [default: All]
        List of chain IDs to use for calculating physicochemical properties. Possible
        values: First, All, or a comma delimited list of chain IDs. The default is to use
        all chain IDs in input file.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infiles <infile1,infile2,infile3...>
        A comma delimited list of input files. The wildcards are also allowed
        in file names.
    --keepInorganics <yes or no>  [default: yes]
        Keep inorganic molecules during calculation of physiochemical properties.
        The inorganic molecules are identified using inorganic selection operator
        available in PyMOL.
    --keepLigands <yes or no>  [default: yes]
        Keep ligand molecules during calculation of physiochemical properties.
        The ligand molecules are identified using organic selection operator
        available in PyMOL.
    --keepSolvents <yes or no>  [default: yes]
        Keep solvent molecules during calculation of physiochemical properties.
        The solvent molecules are identified using solvent selection operator
        available in PyMOL.
    -l, --list
        List available property names without performing any calculations.
    -m, --mode <All or Name1,Name2,Name3,...>  [default: All]
        Comma delimited lists of physicochemical properties to calculate. Default:
         'All'. The following properties may be calculated for macromolecules:
         
            CenterOfMass,MolecularWeight,MolecularSurfaceArea
            SumOfFormalCharges,SumOfPartialCharges,SolventAccessibleSurfaceArea
        
    -o, --outfile <outfile>
        Output file name for writing out calculated values. Supported text file extensions:
        csv, tsv or txt.
    --overwrite
        Overwrite existing files.
    -p, --precision <number>  [default: 3]
        Floating point precision for writing the calculated property values.
    -q, --quiet <yes or no>  [default: yes]
        Do not print information during the calculation of properties.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To calculate all available properties for all chains in input file along with all
    ligands, inorganics and solvents after adding hydrogens and write out a CSV
    file containing calculated values and PDB IDs, type:

        % PyMOLCalculateProperties.py  -i Sample3.pdb -o Sample3Out.csv 

    To calculate specified properties for all chains in input file along with all
    ligands, inorganics and solvents after adding hydrogens and write out a CSV
    file containing calculated values and PDB IDs, type:

        % PyMOLCalculateProperties.py  -m "MolecularWeight,CenterOfMass"
          -i Sample3.pdb -o Sample3Out.csv 

    To calculate all available properties for chain E in input file without including
    ligands, inorganics and solvents, and addition of hydrogens, and write out a
    TSV file containing calculated values and PDB IDs, type:

        % PyMOLCalculateProperties.py  --addHydrogens no -c E --keepLigands
          no --keepInorganics  no --keepSolvents no -i Sample3.pdb -o
          Sample3Out.tsv

    To calculate all available properties for all chains in multiple files along with all
    ligands, inorganics and solvents after adding hydrogens and write out a CSV
    file containing calculated values and PDB IDs, type:

        % PyMOLCalculateProperties.py  -i "Sample3.pdb,Sample4.pdb,Sample5.pdb"
          -o SampleOut.csv 

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
