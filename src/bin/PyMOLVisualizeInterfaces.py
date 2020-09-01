#!/usr/bin/env python
#
# File: PyMOLVisualizeInterfaces.py
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
    GenerateMacromolecularInterfacesVisualization()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def GenerateMacromolecularInterfacesVisualization():
    """Generate visualization for macromolecular interfaces."""

    Outfile = OptionsInfo["PMLOutfile"]
    OutFH = open(Outfile, "w")
    if OutFH is None:
        MiscUtil.PrintError("Failed to open output fie %s " % Outfile)
    
    MiscUtil.PrintInfo("\nGenerating file %s..." % Outfile)

    # Setup header...
    WritePMLHeader(OutFH, ScriptName)
    WritePyMOLParameters(OutFH)
    
    WriteComplexesAndChainsViews(OutFH)
    WriteInterfaceViews(OutFH)

    OutFH.close()

    # Generate PSE file as needed...
    if OptionsInfo["PSEOut"]:
        GeneratePyMOLSessionFile()

def WriteComplexesAndChainsViews(OutFH):
    """Write out PML for viewing complexes and chains in input files. """
    
    # Setup views for input file(s)...
    for FileIndex in range(0, len(OptionsInfo["InfilesInfo"]["InfilesNames"])):
        # Setup PyMOL object names...
        PyMOLObjectNamesInfo = OptionsInfo["InfilesInfo"]["PyMOLObjectNamesInfo"][FileIndex]

        # Setup complex view...
        WriteComplexView(OutFH, FileIndex, PyMOLObjectNamesInfo)

        # Setup chain views...
        SpecifiedChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]

        for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
            WriteChainView(OutFH, FileIndex, PyMOLObjectNamesInfo, ChainID)
            
            # Setup ligand views...
            for LigandID in SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID]:
                WriteChainLigandView(OutFH, FileIndex, PyMOLObjectNamesInfo, ChainID, LigandID)
                
                # Set up ligand level group...
                Enable, Action = [True, "close"]
                GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNamesInfo["Ligands"][ChainID][LigandID]["ChainLigandGroup"], PyMOLObjectNamesInfo["Ligands"][ChainID][LigandID]["ChainLigandGroupMembers"], Enable, Action)
            
            # Setup Chain level group...
            Enable, Action = [True, "open"]
            GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNamesInfo["Chains"][ChainID]["ChainGroup"], PyMOLObjectNamesInfo["Chains"][ChainID]["ChainGroupMembers"], Enable, Action)
    
        # Set up complex level group...
        Enable, Action = [True, "close"]
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNamesInfo["PDBGroup"], PyMOLObjectNamesInfo["PDBGroupMembers"], Enable, Action)
        
        # Delete empty PyMOL objects...
        DeleteEmptyPyMOLObjects(OutFH, FileIndex, PyMOLObjectNamesInfo)
    
def WriteInterfaceViews(OutFH):
    """Write out PML for viewing macromolecular interfaces among specified chains."""
    
    InterfaceChainsAndResiduesInfo = OptionsInfo["InfilesInfo"]["InterfaceChainsAndResiduesInfo"]
    if InterfaceChainsAndResiduesInfo is None:
        return

    PyMOLInterfaceObjectNamesInfo = OptionsInfo["InfilesInfo"]["PyMOLInterfaceObjectNamesInfo"]
    InterfaceChainPairsAndResiduesInfo = OptionsInfo["InfilesInfo"]["InterfaceChainPairsAndResiduesInfo"]
    
    FirstInterfaceID = True
    for InterfaceID in InterfaceChainsAndResiduesInfo["InterfaceIDs"]:
        WriteInterfacePolarContactsView(OutFH, InterfaceID, InterfaceChainPairsAndResiduesInfo, PyMOLInterfaceObjectNamesInfo)
        WriteInterfaceHydrophobicContactsView(OutFH, InterfaceID, InterfaceChainPairsAndResiduesInfo, PyMOLInterfaceObjectNamesInfo)
        
        for InterfaceChainID in InterfaceChainsAndResiduesInfo["InterfaceChainIDs"][InterfaceID]:
            ChainID = InterfaceChainsAndResiduesInfo["ChainIDs"][InterfaceID][InterfaceChainID]
            ResNums = InterfaceChainsAndResiduesInfo["ChainIDsResNums"][InterfaceID][InterfaceChainID]
            ComplexName = InterfaceChainsAndResiduesInfo["ChainIDsComplexNames"][InterfaceID][InterfaceChainID]
            FileIndex = InterfaceChainsAndResiduesInfo["ChainIDsInfileIndices"][InterfaceID][InterfaceChainID]

            WriteInterfaceChainView(OutFH, FileIndex, InterfaceID, InterfaceChainID, ChainID, ResNums, ComplexName, PyMOLInterfaceObjectNamesInfo)

            # Setup interface Chain level group...
            Enable, Action = [True, "open"]
            GenerateAndWritePMLForGroup(OutFH, PyMOLInterfaceObjectNamesInfo["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainGroup"], PyMOLInterfaceObjectNamesInfo["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainGroupMembers"], Enable, Action)

        if FirstInterfaceID:
            FirstInterfaceID = False
            Enable, Action = [True, "open"]
        else:
            Enable, Action = [True, "close"]
            
        GenerateAndWritePMLForGroup(OutFH, PyMOLInterfaceObjectNamesInfo["InterfaceIDs"][InterfaceID]["InterfaceIDGroup"], PyMOLInterfaceObjectNamesInfo["InterfaceIDs"][InterfaceID]["InterfaceIDGroupMembers"], Enable, Action)

    # Setup interfaces level group...
    Enable, Action = [True, "open"]
    GenerateAndWritePMLForGroup(OutFH, PyMOLInterfaceObjectNamesInfo["InterfacesGroup"], PyMOLInterfaceObjectNamesInfo["InterfacesGroupMembers"], Enable, Action)

    DeleteEmptyPyMOLInterfaceObjects(OutFH)
    
    # Setup orientation...
    OutFH.write("""\ncmd.orient("visible", animate = -1)\n""")

def WritePMLHeader(OutFH, ScriptName):
    """Write out PML setting up complex view"""

    HeaderInfo = PyMOLUtil.SetupPMLHeaderInfo(ScriptName)
    OutFH.write("%s\n" % HeaderInfo)

def WritePyMOLParameters(OutFH):
    """Write out PyMOL global parameters. """

    PMLCmds = []
    PMLCmds.append("""cmd.set("transparency", %.2f, "", 0)""" % (OptionsInfo["SurfaceTransparency"]))
    PMLCmds.append("""cmd.set("label_font_id", %s)""" % (OptionsInfo["LabelFontID"]))
    PML = "\n".join(PMLCmds)
    
    OutFH.write("""\n""\n"Setting up PyMOL gobal parameters..."\n""\n""")
    OutFH.write("%s\n" % PML)
    
def WriteComplexView(OutFH, FileIndex, PyMOLObjectNames):
    """Write out PML for viewing polymer complex."""

    # Setup complex...
    Infile = OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex]
    PML = PyMOLUtil.SetupPMLForPolymerComplexView(PyMOLObjectNames["Complex"], Infile, True)
    OutFH.write("""\n""\n"Loading %s and setting up view for complex..."\n""\n""" % Infile)
    OutFH.write("%s\n" % PML)

    if OptionsInfo["SurfaceComplex"]:
        # Setup hydrophobic surface...
        PML = PyMOLUtil.SetupPMLForHydrophobicSurfaceView(PyMOLObjectNames["ComplexHydrophobicSurface"], PyMOLObjectNames["Complex"], ColorPalette = OptionsInfo["SurfaceColorPalette"], Enable = False)
        OutFH.write("\n%s\n" % PML)
    
    # Setup complex group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["ComplexGroup"], PyMOLObjectNames["ComplexGroupMembers"], False, "close")

def WriteChainView(OutFH, FileIndex, PyMOLObjectNames, ChainID):
    """Write out PML for viewing chain."""
    
    OutFH.write("""\n""\n"Setting up views for chain %s..."\n""\n""" % ChainID)
    
    ChainComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
    
    # Setup chain complex group view...
    WriteChainComplexViews(OutFH, FileIndex, PyMOLObjectNames, ChainID)

    # Setup chain view...
    WriteChainAloneViews(OutFH, FileIndex, PyMOLObjectNames, ChainID)
    
    # Setup chain solvent view...
    PML = PyMOLUtil.SetupPMLForSolventView(PyMOLObjectNames["Chains"][ChainID]["Solvent"], ChainComplexName, False)
    OutFH.write("\n%s\n" % PML)

    # Setup chain inorganic view...
    PML = PyMOLUtil.SetupPMLForInorganicView(PyMOLObjectNames["Chains"][ChainID]["Inorganic"], ChainComplexName, False)
    OutFH.write("\n%s\n" % PML)

def WriteChainComplexViews(OutFH, FileIndex, PyMOLObjectNames, ChainID):
    """Write chain complex views. """

    # Setup chain complex...
    ChainComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
    PML = PyMOLUtil.SetupPMLForPolymerChainComplexView(ChainComplexName, PyMOLObjectNames["Complex"], ChainID, True)
    OutFH.write("%s\n" % PML)

    if OptionsInfo["SurfaceChainComplex"]:
        # Setup hydrophobic surface...
        PML = PyMOLUtil.SetupPMLForHydrophobicSurfaceView(PyMOLObjectNames["Chains"][ChainID]["ChainComplexHydrophobicSurface"], ChainComplexName, ColorPalette = OptionsInfo["SurfaceColorPalette"], Enable = False)
        OutFH.write("\n%s\n" % PML)
    
    # Setup chain complex group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroupMembers"], False, "close")
    
def WriteChainAloneViews(OutFH, FileIndex, PyMOLObjectNames, ChainID):
    """Write individual chain views. """

    ChainComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]

    # Setup chain view...
    ChainName = PyMOLObjectNames["Chains"][ChainID]["ChainAlone"]
    PML = PyMOLUtil.SetupPMLForPolymerChainView(ChainName, ChainComplexName, True)
    OutFH.write("\n%s\n" % PML)

    # Setup a non-interface chain view...
    NonInterfaceChainName = PyMOLObjectNames["Chains"][ChainID]["ChainAloneNonInterface"]
    
    InterfaceResNums = OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["InterfaceResNums"][ChainID]
    InterfaceResNumsSelection = "+".join(InterfaceResNums)
    
    Selection = "%s and chain %s and polymer and (not (resi %s))" % (ChainComplexName, ChainID, InterfaceResNumsSelection)
    PML = PyMOLUtil.SetupPMLForSelectionDisplayView(NonInterfaceChainName, Selection, "cartoon", Enable = False)
    OutFH.write("\n%s\n" % PML)
    
    if GetChainAloneContainsSurfacesStatus(FileIndex, ChainID):
        # Setup a generic colored surface...
        PML = PyMOLUtil.SetupPMLForSurfaceView(PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurface"], NonInterfaceChainName, Enable = False, Color = OptionsInfo["SurfaceNonInterfaceColor"])
        OutFH.write("\n%s\n" % PML)
        
        if GetChainAloneSurfaceChainStatus(FileIndex, ChainID):
            # Setup surface colored by hydrophobicity...
            PML = PyMOLUtil.SetupPMLForHydrophobicSurfaceView(PyMOLObjectNames["Chains"][ChainID]["ChainAloneHydrophobicSurface"], NonInterfaceChainName, ColorPalette = OptionsInfo["SurfaceColorPalette"], Enable = False)
            OutFH.write("\n%s\n" % PML)
            
            # Setup surface colored by hyrdophobicity and charge...
            PML = PyMOLUtil.SetupPMLForHydrophobicAndChargeSurfaceView(PyMOLObjectNames["Chains"][ChainID]["ChainAloneHydrophobicChargeSurface"], NonInterfaceChainName, OptionsInfo["AtomTypesColorNames"]["HydrophobicAtomsColor"], OptionsInfo["AtomTypesColorNames"]["NegativelyChargedAtomsColor"], OptionsInfo["AtomTypesColorNames"]["PositivelyChargedAtomsColor"], OptionsInfo["AtomTypesColorNames"]["OtherAtomsColor"], Enable = False, DisplayAs = None)
            OutFH.write("\n%s\n" % PML)
        
        if GetChainAloneSurfaceChainElectrostaticsStatus(FileIndex, ChainID):
            # Setup electrostatics surface...
            SelectionObjectName = NonInterfaceChainName
            ElectrostaticsGroupName = PyMOLObjectNames["Chains"][ChainID]["ChainAloneElectrostaticsGroup"]
            ElectrostaticsGroupMembers = PyMOLObjectNames["Chains"][ChainID]["ChainAloneElectrostaticsGroupMembers"]
            WriteSurfaceElectrostaticsView("Chain", OutFH, SelectionObjectName, ElectrostaticsGroupName, ElectrostaticsGroupMembers)

        # Setup surface group...
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroupMembers"], True, "open")

    # Setup a non-interface group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainAloneNonInterfaceGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainAloneNonInterfaceGroupMembers"], True, "open")
    
    # Setup chain group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"], True, "close")
    
def WriteChainLigandView(OutFH, FileIndex, PyMOLObjectNames, ChainID, LigandID):
    """Write out PML for viewing ligand in a chain."""

    GroupID = "Ligand"
    ComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
    LigandName = PyMOLObjectNames["Ligands"][ChainID][LigandID]["Ligand"]
    
    # Setup main object...
    GroupTypeObjectID = "%s" % (GroupID)
    GroupTypeObjectName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupTypeObjectID]
    
    OutFH.write("""\n""\n"Setting up views for ligand %s in chain %s..."\n""\n""" % (LigandID, ChainID))
    PML = PyMOLUtil.SetupPMLForLigandView(GroupTypeObjectName, ComplexName, LigandID, True)
    OutFH.write("%s\n" % PML)
    
    # Setup ball and stick view...
    BallAndStickNameID = "%sBallAndStick" % (GroupID)
    BallAndStickName = PyMOLObjectNames["Ligands"][ChainID][LigandID][BallAndStickNameID]
    PML = PyMOLUtil.SetupPMLForBallAndStickView(BallAndStickName, GroupTypeObjectName, Enable = False)
    OutFH.write("\n%s\n" % PML)
            
    # Setup group....
    GroupNameID = "%sGroup" % (GroupID)
    GroupMembersID = "%sGroupMembers" % (GroupID)
    GroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupNameID]
    GroupMembers = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID]

    Action = "open"
    Enable = True
    GenerateAndWritePMLForGroup(OutFH, GroupName, GroupMembers, Enable, Action)

def WriteInterfacePolarContactsView(OutFH, InterfaceID, InterfaceChainPairsAndResiduesInfo, PyMOLObjectNames):
    """Write out PML for viewing polar contacts between interface residues."""
    
    if not OptionsInfo["InterfacePolarContacts"]:
        return
    
    WriteInterfaceContactsView("InterfacePolarConatcts", OutFH, InterfaceID, InterfaceChainPairsAndResiduesInfo, PyMOLObjectNames)
    
def WriteInterfaceHydrophobicContactsView(OutFH, InterfaceID, InterfaceChainPairsAndResiduesInfo, PyMOLObjectNames):
    """Write out PML for viewing hydrophobic contacts between interface residues."""
    
    if not OptionsInfo["InterfaceHydrophobicContacts"]:
        return

    WriteInterfaceContactsView("InterfaceHydrophobicContacts", OutFH, InterfaceID, InterfaceChainPairsAndResiduesInfo, PyMOLObjectNames)
    
def WriteInterfaceContactsView(Mode, OutFH, InterfaceID, InterfaceChainPairsAndResiduesInfo, PyMOLObjectNames):
    """Write out PML for viewing polar or hydrophobic contacts between interface residues."""
    
    InterfacePolarContacts = True if re.match("^InterfacePolarConatcts$", Mode, re.I) else False
    
    ChainIDs1, ChainIDs2 =  InterfaceChainPairsAndResiduesInfo["ChainIDsPairs"][InterfaceID]
    ResNums1, ResNums2 = InterfaceChainPairsAndResiduesInfo["ChainIDsResNumsPairs"][InterfaceID]
    
    FileIndex1, FileIndex2 = InterfaceChainPairsAndResiduesInfo["InfileIndicesPairs"][InterfaceID]
    ComplexName1 = OptionsInfo["InfilesInfo"]["PyMOLObjectNamesInfo"][FileIndex1]["Complex"]
    ComplexName2 = OptionsInfo["InfilesInfo"]["PyMOLObjectNamesInfo"][FileIndex2]["Complex"]
    
    Selection1 = SetupSelectionForInterfaceContactsView(ComplexName1, ChainIDs1, ResNums1)
    Selection2 = SetupSelectionForInterfaceContactsView(ComplexName2, ChainIDs2, ResNums2)

    if InterfacePolarContacts:
        ContactsName = PyMOLObjectNames["InterfaceIDs"][InterfaceID]["PolarContacts"]
        ContactsColor = OptionsInfo["InterfacePolarContactsColor"]
        ContactsCutoff = OptionsInfo["InterfaceContactsCutoff"]
        
        PML = PyMOLUtil.SetupPMLForPolarContactsView(ContactsName, Selection1, Selection2, Enable = False, Color = ContactsColor, Cutoff = ContactsCutoff)
    else:
        ContactsName = PyMOLObjectNames["InterfaceIDs"][InterfaceID]["HydrophobicContacts"]
        ContactsColor = OptionsInfo["InterfaceHydrophobicContactsColor"]
        ContactsCutoff = OptionsInfo["InterfaceContactsCutoff"]
        
        PML = PyMOLUtil.SetupPMLForHydrophobicContactsView(ContactsName, Selection1, Selection2, Enable = False, Color = ContactsColor, Cutoff = ContactsCutoff)
    
    OutFH.write("\n%s\n" % PML)
    
    OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (ContactsColor, ContactsName))

def SetupSelectionForInterfaceContactsView(ComplexName, ChainIDs, ResNums):
    """Setup a selection for generating polar or hyrophobic contacts for an interface."""

    ChainSelections = []
    
    for ChainID in ChainIDs:
        ChainResNumsSelection = "+".join(ResNums["ResNums"][ChainID])
        Selection = "(%s and chain %s and polymer and (resi %s))" % (ComplexName, ChainID, ChainResNumsSelection)
        ChainSelections.append(Selection)

    Selection = " or ".join(ChainSelections)
    
    return Selection

def WriteInterfaceChainView(OutFH, FileIndex, InterfaceID, InterfaceChainID, ChainID, ResiduesNums, ComplexName, PyMOLObjectNames):
    """Write out PML for viewing interface residues in a chain."""

    OutFH.write("""\n""\n"Setting up interface views for interface in %s..."\n""\n""" % InterfaceChainID)
    
    ResiduesNumsSelection = "+".join(ResiduesNums)
    
    # Setup a chain for interface residues...
    ChainName = PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["Chain"]
    Selection = "%s and chain %s and polymer and (resi %s)" % (ComplexName, ChainID, ResiduesNumsSelection)

    PML = PyMOLUtil.SetupPMLForSelectionDisplayView(ChainName, Selection, "lines", Enable = True)
    OutFH.write("\n%s\n" % PML)

    OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (OptionsInfo["InterfaceLabelColor"], ChainName))

    WriteInterfaceChainResidueTypesView(OutFH, FileIndex, InterfaceID, InterfaceChainID, ChainID, PyMOLObjectNames)
    
    if GetInterfaceContainsSurfacesStatus(FileIndex, ChainID):
        # Setup generic colored surface...
        PML = PyMOLUtil.SetupPMLForSurfaceView(PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainSurface"], ChainName, Color = OptionsInfo["SurfaceInterfaceColor"], Enable = True, DisplayAs = "lines")
        OutFH.write("\n%s\n" % PML)
        
        if GetInterfaceSurfaceChainStatus(FileIndex, ChainID):
            # Setup surface colored by hydrophobicity...
            PML = PyMOLUtil.SetupPMLForHydrophobicSurfaceView(PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainHydrophobicSurface"], ChainName, ColorPalette = OptionsInfo["SurfaceColorPalette"], Enable = False, DisplayAs = "lines")
            OutFH.write("\n%s\n" % PML)
            
            # Setup surface colored by hyrdophobicity and charge...
            PML = PyMOLUtil.SetupPMLForHydrophobicAndChargeSurfaceView(PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainHydrophobicChargeSurface"], ChainName, OptionsInfo["AtomTypesColorNames"]["HydrophobicAtomsColor"], OptionsInfo["AtomTypesColorNames"]["NegativelyChargedAtomsColor"], OptionsInfo["AtomTypesColorNames"]["PositivelyChargedAtomsColor"], OptionsInfo["AtomTypesColorNames"]["OtherAtomsColor"], Enable = False, DisplayAs = "lines")
            OutFH.write("\n%s\n" % PML)
        
        if GetInterfaceSurfaceChainElectrostaticsStatus(FileIndex, ChainID):
            # Setup electrostatics surface...
            SelectionObjectName = ChainName
            ElectrostaticsGroupName = PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainElectrostaticsGroup"]
            ElectrostaticsGroupMembers = PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainElectrostaticsGroupMembers"]
            WriteSurfaceElectrostaticsView("InterfaceResidues", OutFH, SelectionObjectName, ElectrostaticsGroupName, ElectrostaticsGroupMembers)
        
        # Setup surface group...
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainSurfaceGroup"], PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainSurfaceGroupMembers"], True, "open")
        
    # Setup chain group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainGroup"], PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainGroupMembers"], True, "close")
    
def WriteInterfaceChainResidueTypesView(OutFH, FileIndex, InterfaceID, InterfaceChainID, ChainID, PyMOLObjectNames):
    """Write out PML for viewing interface residue types for a chain. """

    if not GetInterfaceResidueTypesStatus(FileIndex, ChainID):
        return
    
    ChainName = PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["Chain"]
    
    # Setup residue types objects...
    ResiduesGroupIDPrefix = "ChainResidues"
    for SubGroupType in ["Aromatic", "Hydrophobic", "Polar", "Positively_Charged", "Negatively_Charged", "Other"]:
        SubGroupID = re.sub("_", "", SubGroupType)

        ResiduesObjectID = "%s%sResidues" % (ResiduesGroupIDPrefix, SubGroupID)
        ResiduesObjectName = PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesObjectID]

        ResiduesSurfaceObjectID = "%s%sSurface" % (ResiduesGroupIDPrefix, SubGroupID)
        ResiduesSurfaceObjectName = PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesSurfaceObjectID]

        ResiduesColor = OptionsInfo["ResidueTypesParams"][SubGroupType]["Color"] 
        ResiduesNames = OptionsInfo["ResidueTypesParams"][SubGroupType]["Residues"]

        NegateResidueNames = True if re.match("^Other$", SubGroupType, re.I) else False
        WriteResidueTypesResiduesAndSurfaceView(OutFH, ChainName, ResiduesObjectName, ResiduesSurfaceObjectName, ResiduesColor, ResiduesNames, NegateResidueNames)

        # Setup residue type sub groups...
        ResiduesSubGroupID = "%s%sGroup" % (ResiduesGroupIDPrefix, SubGroupID)
        ResiduesSubGroupMembersID = "%s%sGroupMembers" % (ResiduesGroupIDPrefix, SubGroupID)

        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesSubGroupID], PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesSubGroupMembersID], True, "close")
        
    # Setup residue types group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainResiduesGroup"], PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainResiduesGroupMembers"], False, "close")

def WriteResidueTypesResiduesAndSurfaceView(OutFH, SelectionObjectName, Name, SurfaceName, ResiduesColor, ResiduesNames, NegateResidueNames):
    """Write residue types residues and surface view. """

    ResidueNamesSelection = "+".join(ResiduesNames)
    if NegateResidueNames:
        Selection = "%s and (not resn %s)" % (SelectionObjectName, ResidueNamesSelection)
    else:
        Selection = "%s and (resn %s)" % (SelectionObjectName, ResidueNamesSelection)

    # Setup residues...
    PML = PyMOLUtil.SetupPMLForSelectionDisplayView(Name, Selection, "lines", ResiduesColor, True)
    OutFH.write("\n%s\n" % PML)

    # Setup surface...
    PML = PyMOLUtil.SetupPMLForSelectionDisplayView(SurfaceName, Selection, "surface", ResiduesColor, True)
    OutFH.write("\n%s\n" % PML)
    
def WriteSurfaceElectrostaticsView(Mode, OutFH, SelectionObjectName, ElectrostaticsGroupName, ElectrostaticsGroupMembers):
    """Write out PML for viewing surface electrostatics. """

    if len(ElectrostaticsGroupMembers) == 5:
        Name, ContactPotentialName, MapName, LegendName, VolumeName = ElectrostaticsGroupMembers
    else:
        Name, ContactPotentialName, MapName, LegendName = ElectrostaticsGroupMembers
        VolumeName = None

    PMLCmds = []

    # Setup chain...
    PMLCmds.append("""cmd.create("%s", "(%s)")""" % (Name, SelectionObjectName))

    # Setup vacuum electrostatics surface along with associated objects...
    PMLCmds.append("""util.protein_vacuum_esp("%s", mode=2, quiet=0, _self=cmd)""" % (Name))
    
    PMLCmds.append("""cmd.set_name("%s_e_chg", "%s")""" % (Name, ContactPotentialName))
    if re.match("^Chain$", Mode, re.I):
        DisplayStyle = "cartoon"
    else:
        DisplayStyle = "lines"
    PMLCmds.append("""cmd.show("%s", "(%s)")""" % (DisplayStyle, ContactPotentialName))
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(ContactPotentialName, Enable = True))
    
    PMLCmds.append("""cmd.set_name("%s_e_map", "%s")""" % (Name, MapName))
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(MapName, Enable = False))
    
    PMLCmds.append("""cmd.set_name("%s_e_pot", "%s")""" % (Name, LegendName))
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(LegendName, Enable = False))

    if VolumeName is not None:
        PMLCmds.append("""cmd.volume("%s", "%s", "%s", "(%s)")""" % (VolumeName, MapName, "esp", Name))
        PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(VolumeName, Enable = False))
    
    # Delete name and take it out from the group membership. It is
    # is already part of ContactPotential object.
    PMLCmds.append("""cmd.delete("%s")""" % (Name))
    ElectrostaticsGroupMembers.pop(0)
    
    PML = "\n".join(PMLCmds)
    
    OutFH.write("\n%s\n" % PML)
    
    # Setup group...
    GenerateAndWritePMLForGroup(OutFH, ElectrostaticsGroupName, ElectrostaticsGroupMembers, False, "close")

def GenerateAndWritePMLForGroup(OutFH, GroupName, GroupMembers, Enable = False, Action = "close"):
    """Generate and write PML for group. """
    
    PML = PyMOLUtil.SetupPMLForGroup(GroupName, GroupMembers, Enable, Action)
    OutFH.write("""\n""\n"Setting up group %s..."\n""\n""" % GroupName)
    OutFH.write("%s\n" % PML)

def GeneratePyMOLSessionFile():
    """Generate PME file from PML file. """

    PSEOutfile = OptionsInfo["PSEOutfile"]
    PMLOutfile = OptionsInfo["PMLOutfile"]
    
    MiscUtil.PrintInfo("\nGenerating file %s..." % PSEOutfile)
    
    PyMOLUtil.ConvertPMLFileToPSEFile(PMLOutfile, PSEOutfile)
    
    if not os.path.exists(PSEOutfile):
        MiscUtil.PrintWarning("Failed to generate PSE file, %s..." % (PSEOutfile))
    
    if not OptionsInfo["PMLOut"]:
        MiscUtil.PrintInfo("Deleting file %s..." % PMLOutfile)
        os.remove(PMLOutfile)

def DeleteEmptyPyMOLObjects(OutFH, FileIndex, PyMOLObjectNames):
    """Delete empty PyMOL objects. """
    
    if OptionsInfo["AllowEmptyObjects"]:
        return
    
    SpecifiedChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
        OutFH.write("""\n""\n"Checking and deleting empty objects for chain %s..."\n""\n""" % (ChainID))
        
        # Delete any chain level objects...
        WritePMLToCheckAndDeleteEmptyObjects(OutFH, PyMOLObjectNames["Chains"][ChainID]["Solvent"])
        WritePMLToCheckAndDeleteEmptyObjects(OutFH, PyMOLObjectNames["Chains"][ChainID]["Inorganic"])

def DeleteEmptyPyMOLInterfaceObjects(OutFH):
    """Delete empty PyMOL interface objects. """
    
    if OptionsInfo["AllowEmptyObjects"]:
        return
    
    InterfaceChainsAndResiduesInfo = OptionsInfo["InfilesInfo"]["InterfaceChainsAndResiduesInfo"]
    PyMOLInterfaceObjectNamesInfo = OptionsInfo["InfilesInfo"]["PyMOLInterfaceObjectNamesInfo"]
    
    if InterfaceChainsAndResiduesInfo is None:
        return
    
    for InterfaceID in InterfaceChainsAndResiduesInfo["InterfaceIDs"]:
        for InterfaceChainID in InterfaceChainsAndResiduesInfo["InterfaceChainIDs"][InterfaceID]:
            ChainID = InterfaceChainsAndResiduesInfo["ChainIDs"][InterfaceID][InterfaceChainID]
            FileIndex = InterfaceChainsAndResiduesInfo["ChainIDsInfileIndices"][InterfaceID][InterfaceChainID]
            
            # Delete interface residue type objects...
            DeleteEmptyInterfaceChainResidueTypesObjects(OutFH, FileIndex, InterfaceID, InterfaceChainID, ChainID, PyMOLInterfaceObjectNamesInfo)
            
def DeleteEmptyInterfaceChainResidueTypesObjects(OutFH, FileIndex, InterfaceID, InterfaceChainID, ChainID, PyMOLObjectNames):
    """Delete empty interface chain residue objects. """
    
    if not GetInterfaceResidueTypesStatus(FileIndex, ChainID):
        return
    
    ResiduesGroupIDPrefix = "ChainResidues"
    for GroupType in ["Aromatic", "Hydrophobic", "Polar", "Positively_Charged", "Negatively_Charged", "Other"]:
        GroupID = re.sub("_", "", GroupType)
        
        ResiduesGroupID = "%s%sGroup" % (ResiduesGroupIDPrefix, GroupID)
        GroupName = PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesGroupID]
        
        GroupObjectNamesList = []
        
        ResiduesObjectID = "%s%sResidues" % (ResiduesGroupIDPrefix, GroupID)
        ResiduesObjectName = PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesObjectID]
        GroupObjectNamesList.append(ResiduesObjectName)
        
        ResiduesSurfaceObjectID = "%s%sSurface" % (ResiduesGroupIDPrefix, GroupID)
        ResiduesSurfaceObjectName = PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesSurfaceObjectID]
        GroupObjectNamesList.append(ResiduesSurfaceObjectName)
        
        GroupObjectNames = ",".join(GroupObjectNamesList)
        WritePMLToCheckAndDeleteEmptyObjects(OutFH, GroupObjectNames, GroupName)
    
def WritePMLToCheckAndDeleteEmptyObjects(OutFH, ObjectName, ParentObjectName = None):
    """Write PML to check and delete empty PyMOL objects. """
    
    if ParentObjectName is None:
        PML = """CheckAndDeleteEmptyObjects("%s")""" % (ObjectName)
    else:
        PML = """CheckAndDeleteEmptyObjects("%s", "%s")""" % (ObjectName, ParentObjectName)
    
    OutFH.write("%s\n" % PML)
    
def RetrieveInfilesInfo():
    """Retrieve information for input files."""

    InfilesInfo = {}
    
    InfilesInfo["InfilesNames"] = []
    InfilesInfo["InfilesRoots"] = []
    InfilesInfo["ChainsAndLigandsInfo"] = []
    InfilesInfo["SpecifiedChainsAndLigandsInfo"] = []
    
    InfilesInfo["PyMOLObjectNamesInfo"] = []
    
    InfilesInfo["SingleInfileMode"] = False
    InfilesInfo["InterfaceChainsAndResiduesInfo"] = None
    InfilesInfo["InterfaceChainPairsAndResiduesInfo"] = None
    
    InfilesInfo["PyMOLInterfaceObjectNamesInfo"] = None

    InfilesCount = 0
    for Infile in OptionsInfo["InfilesNames"]:
        InfilesCount += 1
        FileDir, FileName, FileExt = MiscUtil.ParseFileName(Infile)
        InfileRoot = FileName
        
        ChainsAndLigandInfo = PyMOLUtil.GetChainsAndLigandsInfo(Infile, InfileRoot)
        
        InfilesInfo["InfilesNames"].append(Infile)
        InfilesInfo["InfilesRoots"].append(InfileRoot)
        InfilesInfo["ChainsAndLigandsInfo"].append(ChainsAndLigandInfo)

    if InfilesCount > 2:
        MiscUtil.PrintError("Number of input files, %s, specified using \"-i, --infiles\" option is not valid. Number of allowed files: 1 or 2" % (InfilesCount))
    
    InfilesInfo["SingleInfileMode"] = True if InfilesCount == 1 else False
    
    OptionsInfo["InfilesInfo"] = InfilesInfo

def ProcessInterfaceChainIDs():
    """Process specified interface chain IDs for input files."""

    ValidateInterfaceChainIDs()

    SetupChainsAndLigandsInfo()
    SetupPyMOLObjectNamesInfo()
    
    SetupInterfaceChainPairsAndResiduesInfo()
    ProcessInterfaceChainPairsAndResiduesInfo()

    SetupPyMOLInterfaceObjectNamesInfo()

def ValidateInterfaceChainIDs():
    """Check for the presence of interface IDs in input file(s). """

    if  re.match("^auto$", OptionsInfo["InterfaceChainIDs"], re.I):
        AutoAssignInterfaceChainIDs()
        SetupInfilesInterfaceChainIDsLists()
        return

    MiscUtil.PrintInfo("\nValidating interface chain IDs...")
    
    # Check for the presences of interface chain IDs across input files..
    SingleInfileMode = OptionsInfo["InfilesInfo"]["SingleInfileMode"]
    Infile1ChainIDs = OptionsInfo["InfilesInfo"]["ChainsAndLigandsInfo"][0]["ChainIDs"]
    Infile2ChainIDs = []
    if not SingleInfileMode:
        Infile2ChainIDs = OptionsInfo["InfilesInfo"]["ChainsAndLigandsInfo"][1]["ChainIDs"]

    InterfaceChainIDsList = OptionsInfo["InterfaceChainIDsList"]

    for Index in range(0, len(InterfaceChainIDsList), 2):
        ChainIDs1 = InterfaceChainIDsList[Index]
        ChainIDs2 = InterfaceChainIDsList[Index + 1]

        for ChainID in ChainIDs1:
            if not ChainID in Infile1ChainIDs:
                MiscUtil.PrintError("The chain ID, %s, specified using \"-c, --chainIDs\" for a chain IDs pairs, \"%s,%s\", must be present in first input file." % (ChainID, "+".join(ChainIDs1), "+".join(ChainIDs2)))
        
        for ChainID in ChainIDs2:
            if  SingleInfileMode:
                if not ChainID in Infile1ChainIDs:
                    MiscUtil.PrintError("The chain ID, %s, specified using \"-c, --chainIDs\" for a chain IDs pairs, \"%s,%s\", must be present in first input file." % (ChainID, "+".join(ChainIDs1), "+".join(ChainIDs2)))
            else:
                if not ChainID in Infile2ChainIDs:
                    MiscUtil.PrintError("The chain ID, %s, specified using \"-c, --chainIDs\" for a chain IDs pairs, \"%s,%s\", must be present in second input file." % (ChainID, "+".join(ChainIDs1), "+".join(ChainIDs2)))

    # Check for any duplicate interface chain IDs specifications...
    CanonicalInterfaceIDsMap = {}
    
    for Index in range(0, len(InterfaceChainIDsList), 2):
        ChainIDs1 = InterfaceChainIDsList[Index]
        ChainIDs2 = InterfaceChainIDsList[Index + 1]
        InterfaceID = "%s,%s" % ("+".join(ChainIDs1), "+".join(ChainIDs2))
        
        SortedChainIDs1 = sorted(ChainIDs1)
        SortedChainIDs2 = sorted(ChainIDs2)
        CanonicalInterfaceID = "%s,%s" % ("+".join(SortedChainIDs1), "+".join(SortedChainIDs2))

        if CanonicalInterfaceID in CanonicalInterfaceIDsMap:
            MiscUtil.PrintError("The chain ID pair, \"%s\", using \"-c, --chainIDs\", option has been specified multiple times, \"%s\"." % (CanonicalInterfaceIDsMap[CanonicalInterfaceID], OptionsInfo["InterfaceChainIDs"]))
        else:
            CanonicalInterfaceIDsMap[CanonicalInterfaceID] =  InterfaceID

    SetupInfilesInterfaceChainIDsLists()
    
def SetupInfilesInterfaceChainIDsLists():
    """Setup interface chain IDs list for infiles. """
    
    Infie1InterfaceChainIDsList = []
    Infie2InterfaceChainIDsList = []
    
    SingleInfileMode = OptionsInfo["InfilesInfo"]["SingleInfileMode"]
    InterfaceChainIDsList = OptionsInfo["InterfaceChainIDsList"]
    
    for Index in range(0, len(InterfaceChainIDsList), 2):
        ChainIDs1 = InterfaceChainIDsList[Index]
        ChainIDs2 = InterfaceChainIDsList[Index + 1]

        Infie1InterfaceChainIDsList.extend(ChainIDs1)
        if SingleInfileMode:
            Infie1InterfaceChainIDsList.extend(ChainIDs2)
        else:
            Infie2InterfaceChainIDsList.extend(ChainIDs2)

    Infie1InterfaceChainIDsList = sorted(Infie1InterfaceChainIDsList)
    Infie2InterfaceChainIDsList = sorted(Infie2InterfaceChainIDsList)
    
    OptionsInfo["InfilesInterfaceChainIDsList"] = [Infie1InterfaceChainIDsList, Infie2InterfaceChainIDsList]

def AutoAssignInterfaceChainIDs():
    """Handle automatic assignment of interface chain IDs. """
    
    Infile1ChainIDs = OptionsInfo["InfilesInfo"]["ChainsAndLigandsInfo"][0]["ChainIDs"]
    Infile2ChainIDs = []
    if not OptionsInfo["InfilesInfo"]["SingleInfileMode"]:
        Infile2ChainIDs = OptionsInfo["InfilesInfo"]["ChainsAndLigandsInfo"][1]["ChainIDs"]
    
    InterfaceChainIDsList = []
    if OptionsInfo["InfilesInfo"]["SingleInfileMode"]:
        # Take first two chains from first input file...
        if len(Infile1ChainIDs) < 2:
            MiscUtil.PrintError("Failed to automatically set interface chain IDs. Number of chains, %s, in input file specified using \"-i, --infiles\" option must be >= 2. " % (len(Infile1ChainIDs)))
        InterfaceChainIDsList.append([Infile1ChainIDs[0]])
        InterfaceChainIDsList.append([Infile1ChainIDs[1]])
    else:
        # Take first chain from each input file...
        if len(Infile1ChainIDs) < 1:
            MiscUtil.PrintError("Failed to automatically set interface chain IDs. Number of chains, %s, in first input file specified using \"-i, --infiles\" option must be >= 1. " % (len(Infile1ChainIDs)))
        if len(Infile2ChainIDs) < 1:
            MiscUtil.PrintError("Failed to automatically set interface chain IDs. Number of chains, %s, in second input file specified using \"-i, --infiles\" option must be >= 1. " % (len(Infile1ChainIDs)))
        InterfaceChainIDsList.append([Infile1ChainIDs[0]])
        InterfaceChainIDsList.append([Infile2ChainIDs[1]])
    
    OptionsInfo["InterfaceChainIDsList"] = []
    OptionsInfo["InterfaceChainIDsList"].extend(InterfaceChainIDsList)
    
    return

def SetupChainsAndLigandsInfo():
    """Setup chains and ligands info for input files to visualize macromolecules. """

    OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"] = []
    for FileIndex in range(0, len(OptionsInfo["InfilesInfo"]["InfilesNames"])):
        Infile = OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex]
        MiscUtil.PrintInfo("\nSetting up chain and ligand information for input file %s..." % Infile)
        
        ChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["ChainsAndLigandsInfo"][FileIndex]
        InfilesInterfaceChainIDsList = OptionsInfo["InfilesInterfaceChainIDsList"][FileIndex]
        
        SpecifiedChainsAndLigandsInfo = {}
        SpecifiedChainsAndLigandsInfo["ChainIDs"] = []
        SpecifiedChainsAndLigandsInfo["InterfaceResNums"] = {}
        SpecifiedChainsAndLigandsInfo["LigandIDs"] = {}

        if len(InfilesInterfaceChainIDsList):
            # Add unique interface IDs to the chain IDs list for visualization. Interface chain IDs
            # may contain duplicate chain IDs due to the presence # of same chain in multiple
            # interfaces.
            for ChainID in InfilesInterfaceChainIDsList:
                if not ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
                    SpecifiedChainsAndLigandsInfo["ChainIDs"].append(ChainID)

            for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
                # Initialize interface residue numbers to be assigned later...
                SpecifiedChainsAndLigandsInfo["InterfaceResNums"][ChainID] = []
                
                # Setup ligand IDs...
                SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID] = SetupSpecifiedLigandIDs(FileIndex, ChainID, ChainsAndLigandsInfo)
        
        ProcessResidueTypesAndSurfaceOptions(FileIndex, SpecifiedChainsAndLigandsInfo)
        
        OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"].append(SpecifiedChainsAndLigandsInfo)

def SetupSpecifiedLigandIDs(FileIndex, ChainID, ChainsAndLigandsInfo):
    """Setup specified ligand IDs for input file. """

    LigandIDs = []
    
    if re.match("^All$", OptionsInfo["LigandIDs"], re.I):
        LigandIDs = ChainsAndLigandsInfo["LigandIDs"][ChainID]
        return LigandIDs
    elif re.match("^Largest|Auto$", OptionsInfo["LigandIDs"], re.I):
        LargestLigandID = ChainsAndLigandsInfo["LigandIDs"][ChainID][0] if (len(ChainsAndLigandsInfo["LigandIDs"][ChainID])) else None
        if LargestLigandID is not None:
            LigandIDs.append(LargestLigandID)
        return LigandIDs
    elif re.match("^None$", OptionsInfo["LigandIDs"], re.I):
        return LigandIDs

    Infile = OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex]
    ValidLigandIDs = ChainsAndLigandsInfo["LigandIDs"][ChainID]
    
    SpecifiedLigandIDs = re.sub(" ", "", OptionsInfo["LigandIDs"])
    if not SpecifiedLigandIDs:
        MiscUtil.PrintError("No valid value specified using \-l, --ligandIDs\" option.")

    LigandIDsWords = SpecifiedLigandIDs.split(",")
    for LigandID in LigandIDsWords:
        if not LigandID in ValidLigandIDs:
            LigandIDsListNames = ",".join(ValidLigandIDs) if len(ValidLigandIDs) else "None"
            MiscUtil.PrintWarning("The ligand ID, %s, specified using \"-l, --ligandiIDs\" option is not valid for chain, %s, in input file, %s. It'll be ignored. Valid ligand IDs: %s" % (LigandID, ChainID, Infile, LigandIDsListNames))
            continue
        
        if LigandID in LigandIDs:
            MiscUtil.PrintWarning("The ligand ID, %s, has already been specified using \"-l, --ligandIDs\" option for chain, %s, in input file, %s. It'll be ignored." % (LigandID, ChainID, Infile))
            continue
    
        LigandIDs.append(LigandID)
    
    if not len(LigandIDs):
            MiscUtil.PrintWarning("No valid ligand IDs \"%s\" specified using \"-l, --ligandIDs\" option for chain ID, %s, in input file, %s." % (OptionsInfo["LigandIDs"], ChainID, Infile))
    
    return LigandIDs
        
def SetupInterfaceChainPairsAndResiduesInfo():
    """Setup chain and residue pairs corresponding to interfaces."""

    MiscUtil.PrintInfo("\nIdentifying interface residues...")
    
    SingleInfileMode = OptionsInfo["InfilesInfo"]["SingleInfileMode"]
    InterfaceChainIDsList = OptionsInfo["InterfaceChainIDsList"]

    if not len(InterfaceChainIDsList):
        MiscUtil.PrintError("Failed to identify interface residues: No valid chain ID pairs available for interfaces")

    # Load infiles to identify interface residues...
    Infile1 = OptionsInfo["InfilesInfo"]["InfilesNames"][0]
    MolName1 = OptionsInfo["InfilesInfo"]["InfilesRoots"][0]
    pymol.cmd.load(Infile1, MolName1)
    if SingleInfileMode:
        MolName2 = MolName1
    else:
        Infile2 = OptionsInfo["InfilesInfo"]["InfilesNames"][1]
        MolName2 = OptionsInfo["InfilesInfo"]["InfilesRoots"][1]
        pymol.cmd.load(Infile2, MolName2)

    # Initialize data...
    InterfaceChainPairsAndResiduesInfo = {}
    InterfaceChainPairsAndResiduesInfo["InterfaceIDs"] = []
    InterfaceChainPairsAndResiduesInfo["ChainIDsPairs"] = {}
    InterfaceChainPairsAndResiduesInfo["ChainIDsResNumsPairs"] = {}
    InterfaceChainPairsAndResiduesInfo["InfileIndicesPairs"] = {}
    
    Method = OptionsInfo["Method"]
    MethodCutoff = OptionsInfo["MethodCutoff"]
    MiscUtil.PrintInfo("Methodology: %s; Cutoff: %.2f" % (Method, MethodCutoff))

    for Index in range(0, len(InterfaceChainIDsList), 2):
        ChainIDs1 = sorted(InterfaceChainIDsList[Index])
        ChainIDs2 = sorted(InterfaceChainIDsList[Index + 1])

        ChainIDs1Prefix = "Chains" if len(ChainIDs1) > 1 else "Chain"
        ChainIDs2Prefix = "Chains" if len(ChainIDs2) > 1 else "Chain"
        InterfaceID = "%s%s_%s%s" % (ChainIDs1Prefix, "+".join(ChainIDs1), ChainIDs2Prefix, "+".join(ChainIDs2))
        
        ListInterfaceID = "%s_%s" % ("+".join(ChainIDs1), "+".join(ChainIDs2))

        FileIndex1 = 0
        FileIndex2 = 0 if SingleInfileMode else 1

        if InterfaceID in InterfaceChainPairsAndResiduesInfo["InterfaceIDs"]:
            MiscUtil.PrintInfo("Ignoring interface ID %s: It has already been defined" % InterfaceID)
            continue

        # Identify and list interface chains and residues...
        ChainsAndResiduesInfo1, ChainsAndResiduesInfo2 = GetInterfaceChainsAndResiduesInfo(MolName1, ChainIDs1, MolName2, ChainIDs2, Method, MethodCutoff)
        ListInterfaceChainsAndResidues(ListInterfaceID, MolName1, ChainIDs1, ChainsAndResiduesInfo1, MolName2, ChainIDs2, ChainsAndResiduesInfo2)

        # Check presence of interface residues...
        if (not (len(ChainsAndResiduesInfo1["ChainIDs"]) and len(ChainsAndResiduesInfo2["ChainIDs"]))):
            MiscUtil.PrintWarning("Ignoring interface ID %s: Failed to identify any interface residues. PyMOL groups and objects won't be created." % InterfaceID)
            continue
        
        # Collect residue numbers for set of interface residues in each chain...
        InterfaceChainIDs1, InterfaceChainResidueNums1 = GetInterfaceChainsAndResidueNumbers(ChainIDs1, ChainsAndResiduesInfo1)
        InterfaceChainIDs2, InterfaceChainResidueNums2 = GetInterfaceChainsAndResidueNumbers(ChainIDs2, ChainsAndResiduesInfo2)

        InterfaceChainPairsAndResiduesInfo["InterfaceIDs"].append(InterfaceID)
        InterfaceChainPairsAndResiduesInfo["ChainIDsPairs"][InterfaceID] = [InterfaceChainIDs1, InterfaceChainIDs2]
        InterfaceChainPairsAndResiduesInfo["ChainIDsResNumsPairs"][InterfaceID] = [InterfaceChainResidueNums1, InterfaceChainResidueNums2]
        InterfaceChainPairsAndResiduesInfo["InfileIndicesPairs"][InterfaceID] = [FileIndex1, FileIndex2]

    InterfaceChainPairsAndResiduesInfo["InterfaceIDs"] = sorted(InterfaceChainPairsAndResiduesInfo["InterfaceIDs"])
    OptionsInfo["InfilesInfo"]["InterfaceChainPairsAndResiduesInfo"] = InterfaceChainPairsAndResiduesInfo

    # Delete loaded objects...
    pymol.cmd.delete(MolName1)
    if not SingleInfileMode:
        pymol.cmd.delete(MolName2)
    
def ProcessInterfaceChainPairsAndResiduesInfo():
    """Process chain and residue pairs for visualizing interfaces."""

    InterfaceChainsAndResiduesInfo = {}
    InterfaceChainsAndResiduesInfo["InterfaceIDs"] = []
    InterfaceChainsAndResiduesInfo["InterfaceChainIDs"] = {}
    InterfaceChainsAndResiduesInfo["ChainIDs"] = {}
    InterfaceChainsAndResiduesInfo["ChainIDsResNums"] = {}
    InterfaceChainsAndResiduesInfo["ChainIDsComplexNames"] = {}
    InterfaceChainsAndResiduesInfo["ChainIDsInfileIndices"] = {}
    
    InterfaceChainPairsAndResiduesInfo = OptionsInfo["InfilesInfo"]["InterfaceChainPairsAndResiduesInfo"]
    
    for InterfaceID in InterfaceChainPairsAndResiduesInfo["InterfaceIDs"]:
        # Flatten pairwise interface chain and residues info...
        FlattenedInterfaceChainIDs, FlattenedChainIDs, FlattenedResNums, FlattenedComplexNames, FlattenedInfileIndices = FlattenIntferfaceChainPairsAndResiduesInfo(InterfaceChainPairsAndResiduesInfo, InterfaceID)
        
        if not InterfaceID in InterfaceChainsAndResiduesInfo["InterfaceIDs"]:
            InterfaceChainsAndResiduesInfo["InterfaceIDs"].append(InterfaceID)
            
            InterfaceChainsAndResiduesInfo["InterfaceChainIDs"][InterfaceID] = []
            InterfaceChainsAndResiduesInfo["ChainIDs"][InterfaceID] = {}
            InterfaceChainsAndResiduesInfo["ChainIDsResNums"][InterfaceID] = {}
            InterfaceChainsAndResiduesInfo["ChainIDsComplexNames"][InterfaceID] = {}
            InterfaceChainsAndResiduesInfo["ChainIDsInfileIndices"][InterfaceID] = {}

        # Track interface information for visualization...
        for Index in range(0, len(FlattenedInterfaceChainIDs)):
            InterfaceChainID = FlattenedInterfaceChainIDs[Index]
            ChainID = FlattenedChainIDs[Index]
            ResNums = FlattenedResNums[Index]
            ComplexName = FlattenedComplexNames[Index]
            FileIndex = FlattenedInfileIndices[Index]

            if not len(ResNums):
                continue

            InterfaceChainsAndResiduesInfo["InterfaceChainIDs"][InterfaceID].append(InterfaceChainID)
            InterfaceChainsAndResiduesInfo["ChainIDs"][InterfaceID][InterfaceChainID] = ChainID
            InterfaceChainsAndResiduesInfo["ChainIDsResNums"][InterfaceID][InterfaceChainID] = ResNums
            InterfaceChainsAndResiduesInfo["ChainIDsComplexNames"][InterfaceID][InterfaceChainID] = ComplexName
            InterfaceChainsAndResiduesInfo["ChainIDsInfileIndices"][InterfaceID][InterfaceChainID] = FileIndex

            OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["InterfaceResNums"][ChainID].extend(ResNums)

    OptionsInfo["InfilesInfo"]["InterfaceChainsAndResiduesInfo"] = InterfaceChainsAndResiduesInfo

def FlattenIntferfaceChainPairsAndResiduesInfo(InterfaceChainPairsAndResiduesInfo, InterfaceID):
    """Flatten interface chain and residues info. """

    SingleInfileMode = OptionsInfo["InfilesInfo"]["SingleInfileMode"]
    
    ChainIDs1, ChainIDs2 =  InterfaceChainPairsAndResiduesInfo["ChainIDsPairs"][InterfaceID]
    ResNums1, ResNums2 = InterfaceChainPairsAndResiduesInfo["ChainIDsResNumsPairs"][InterfaceID]
    FileIndex1, FileIndex2 = InterfaceChainPairsAndResiduesInfo["InfileIndicesPairs"][InterfaceID]

    # Set complex names..
    PDBGroup1 = OptionsInfo["InfilesInfo"]["PyMOLObjectNamesInfo"][FileIndex1]["PDBGroup"]
    ComplexName1 = OptionsInfo["InfilesInfo"]["PyMOLObjectNamesInfo"][FileIndex1]["Complex"]
    
    PDBGroup2 = OptionsInfo["InfilesInfo"]["PyMOLObjectNamesInfo"][FileIndex2]["PDBGroup"]
    ComplexName2 = OptionsInfo["InfilesInfo"]["PyMOLObjectNamesInfo"][FileIndex2]["Complex"]
    
    ChainIDs = []
    InterfaceChainIDs = []
    ResNums = []
    ComplexNames = []
    FileIndices = []
    
    for ChainID in ChainIDs1:
        ChainIDs.append(ChainID)
        InterfaceChainID = ("Chain%s" % ChainID) if SingleInfileMode else ("Chain%s_%s" % (ChainID, PDBGroup1))
        InterfaceChainIDs.append(InterfaceChainID)
        ResNums.append(ResNums1["ResNums"][ChainID])
        ComplexNames.append(ComplexName1)
        FileIndices.append(FileIndex1)
    
    for ChainID in ChainIDs2:
        ChainIDs.append(ChainID)
        InterfaceChainID = ("Chain%s" % ChainID) if SingleInfileMode else ("Chain%s_%s" % (ChainID, PDBGroup2))
        InterfaceChainIDs.append(InterfaceChainID)
        ResNums.append(ResNums2["ResNums"][ChainID])
        ComplexNames.append(ComplexName2)
        FileIndices.append(FileIndex2)

    return InterfaceChainIDs, ChainIDs, ResNums, ComplexNames, FileIndices
    
def GetInterfaceChainsAndResiduesInfo(MolName1, ChainIDs1, MolName2, ChainIDs2, Method, Cutoff):
    """Get interface chains and residues info for chains using a specified methodology."""

    InterfaceChainsResiduesInfo1 = None
    InterfaceChainsResiduesInfo2 = None

    ChainNames1 = ",".join(ChainIDs1)
    ChainNames2 = ",".join(ChainIDs2)
    
    if re.match("^BySASAChange$", Method, re.I):
        InterfaceChainsResiduesInfo1, InterfaceChainsResiduesInfo2 = PyMOLUtil.GetInterfaceChainsResiduesBySASAChange(MolName1, ChainNames1, MolName2, ChainNames2, Cutoff)
    elif re.match("^ByHeavyAtomsDistance$", Method, re.I):
        InterfaceChainsResiduesInfo1, InterfaceChainsResiduesInfo2 = PyMOLUtil.GetnterfaceChainsResiduesByHeavyAtomsDistance(MolName1, ChainNames1, MolName2, ChainNames2, Cutoff)
    elif re.match("^ByCAlphaAtomsDistance$", Method, re.I):
        InterfaceChainsResiduesInfo1, InterfaceChainsResiduesInfo2 = PyMOLUtil.GetInterfaceChainsResiduesByCAlphaAtomsDistance(MolName1, ChainNames1, MolName2, ChainNames2, Cutoff)
    else:
        MiscUtil.PrintError("Failed to retrieve interface residues information: Method %s is not valid..." % Method)

    return InterfaceChainsResiduesInfo1, InterfaceChainsResiduesInfo2

def ListInterfaceChainsAndResidues(InterfaceID, MolName1, ChainIDs1, ChainsAndResiduesInfo1, MolName2, ChainIDs2, ChainsAndResiduesInfo2):
    """List interface chains and residues for an interface."""
    
    ChainNames1, ResiduesInfo1 = PrepareInterfaceChainsAndResiduesInfo(ChainsAndResiduesInfo1)
    ChainNames2, ResiduesInfo2 = PrepareInterfaceChainsAndResiduesInfo(ChainsAndResiduesInfo2)
    
    if len(ChainNames1) and len(ChainNames2):
        MiscUtil.PrintInfo("\nListing interface residues for interface chain IDs: %s" % (InterfaceID))
        
        ListInterfaceChainsAndResiduesInfo(InterfaceID,  MolName1, ChainIDs1, ChainNames1, ResiduesInfo1)
        ListInterfaceChainsAndResiduesInfo(InterfaceID,  MolName2, ChainIDs2, ChainNames2, ResiduesInfo2)
    else:
        MiscUtil.PrintInfo("\nListing interface residues for interface chain IDs: %s" % (InterfaceID))
        MiscUtil.PrintInfo("Interface chain IDs: None; ChainID: None; Number of interface residues: 0")

def PrepareInterfaceChainsAndResiduesInfo(ChainsAndResiduesInfo):
    """Prepare interface chains and residues info for listing."""

    ChainNames = []
    ResiduesInfo = []

    for ChainID in ChainsAndResiduesInfo["ChainIDs"]:
        ChainNames.append(ChainID)
        
        # Setup distribution of residues...
        LineWords = []
        ResiduesCount = 0
        SortedResNames = sorted(ChainsAndResiduesInfo["ResNames"][ChainID], key = lambda ResName: ChainsAndResiduesInfo["ResCount"][ChainID][ResName], reverse = True)
        for ResName in SortedResNames:
            ResCount = ChainsAndResiduesInfo["ResCount"][ChainID][ResName]
            LineWords.append("%s - %s" % (ResName, ResCount))
            ResiduesCount += ResCount
        
        ResiduesDistribution = "; ".join(LineWords) if len(LineWords) else None
        
        # Setup residue IDs sorted by residue numbers...
        ResNumMap = {}
        for ResName in ChainsAndResiduesInfo["ResNames"][ChainID]:
            for ResNum in ChainsAndResiduesInfo["ResNum"][ChainID][ResName]:
                ResNumMap[ResNum] = ResName
        
        LineWords = []
        for ResNum in sorted(ResNumMap, key = int):
            ResName = ResNumMap[ResNum]
            ResID = "%s_%s" % (ResName, ResNum)
            LineWords.append(ResID)
        ResiduesIDs = ", ".join(LineWords) if len(LineWords) else None

        ResiduesInfo.append([ResiduesCount, ResiduesDistribution, ResiduesIDs])
                
    return ChainNames, ResiduesInfo

def ListInterfaceChainsAndResiduesInfo(InterfaceID, MolName, InterfaceChainIDs, ChainNames, ResiduesInfo):
    """List interface chains and residues. """

    for ChainID in InterfaceChainIDs:
        if not ChainID in ChainNames:
            MiscUtil.PrintWarning("Interface chain IDs: %s; MoleculeID: %s; ChainID: %s; Number of interface residues: 0. PyMOL groups and objects related to chain won't be created." % (InterfaceID, MolName, ChainID))
            continue
        
        for Index in range(0, len(ChainNames)):
            ChainName = ChainNames[Index]
            ChainResiduesInfo = ResiduesInfo[Index]
            if re.match(ChainID, ChainName, re.I):
                MiscUtil.PrintInfo("\nInterface chain IDs: %s; MoleculeID: %s; ChainID: %s; Number of interface residues: %d" % (InterfaceID, MolName, ChainName, ChainResiduesInfo[0]))
                MiscUtil.PrintInfo("Residue distribution: %s\nResidue IDs: %s" % (ChainResiduesInfo[1], ChainResiduesInfo[2]))
                continue

def GetInterfaceChainsAndResidueNumbers(ChainIDs, ChainsAndResiduesInfo):
    """Collect interface residue numbers for chains. """

    ChainResidueNums = {}
    ChainResidueNums["ChainIDs"] = []
    ChainResidueNums["ResNums"] = {}
    
    for ChainID in ChainIDs:
        if not ChainID in ChainsAndResiduesInfo["ChainIDs"]:
            continue
        
        ChainResidueNums["ChainIDs"].append(ChainID)
        ChainResidueNums["ResNums"][ChainID] = []

        ResNums = []
        for ResName in ChainsAndResiduesInfo["ResNames"][ChainID]:
            ResNums.extend(ChainsAndResiduesInfo["ResNum"][ChainID][ResName])
        
        ChainResidueNums["ResNums"][ChainID] = sorted(ResNums, key = int)

    InterfaceChainIDs = ChainResidueNums["ChainIDs"]
    
    return InterfaceChainIDs, ChainResidueNums

def SetupPyMOLObjectNamesInfo():
    """Setup PyMOL object names for displaying macromolecules. """

    OptionsInfo["InfilesInfo"]["PyMOLObjectNamesInfo"] = []
    
    for FileIndex in range(0, len(OptionsInfo["InfilesInfo"]["InfilesNames"])):
        PyMOLObjectNamesInfo = SetupPyMOLObjectNames(FileIndex)
        OptionsInfo["InfilesInfo"]["PyMOLObjectNamesInfo"].append(PyMOLObjectNamesInfo)
    
def SetupPyMOLObjectNames(FileIndex):
    """Setup hierarchy of PyMOL groups and objects for ligand centric views of
    chains and ligands present in input file.
    """

    PyMOLObjectNames = {}
    PyMOLObjectNames["Chains"] = {}
    PyMOLObjectNames["Ligands"] = {}

    # Setup groups and objects for complex...
    SetupPyMOLObjectNamesForComplex(FileIndex, PyMOLObjectNames)
    
    # Setup groups and objects for chain...
    SpecifiedChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
        SetupPyMOLObjectNamesForChain(FileIndex, PyMOLObjectNames, ChainID)
        
        # Setup groups and objects for ligand...
        for LigandID in SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID]:
            SetupPyMOLObjectNamesForLigand(FileIndex, PyMOLObjectNames, ChainID, LigandID)

    return PyMOLObjectNames

def SetupPyMOLObjectNamesForComplex(FileIndex, PyMOLObjectNames):
    """Setup groups and objects for complex. """
    
    PDBFileRoot = OptionsInfo["InfilesInfo"]["InfilesRoots"][FileIndex]
    
    PDBGroupName = "%s" % PDBFileRoot
    PyMOLObjectNames["PDBGroup"] = PDBGroupName
    PyMOLObjectNames["PDBGroupMembers"] = []

    ComplexGroupName = "%s.Complex" % PyMOLObjectNames["PDBGroup"]
    PyMOLObjectNames["ComplexGroup"] = ComplexGroupName
    PyMOLObjectNames["PDBGroupMembers"].append(ComplexGroupName)
    
    PyMOLObjectNames["Complex"] = "%s.Complex" % ComplexGroupName
    if OptionsInfo["SurfaceComplex"]:
        PyMOLObjectNames["ComplexHydrophobicSurface"] = "%s.Surface" % ComplexGroupName

    PyMOLObjectNames["ComplexGroupMembers"] = []
    PyMOLObjectNames["ComplexGroupMembers"].append(PyMOLObjectNames["Complex"])
    if OptionsInfo["SurfaceComplex"]:
        PyMOLObjectNames["ComplexGroupMembers"].append(PyMOLObjectNames["ComplexHydrophobicSurface"])

def SetupPyMOLObjectNamesForChain(FileIndex, PyMOLObjectNames, ChainID):
    """Setup groups and objects for chain."""

    PDBGroupName = PyMOLObjectNames["PDBGroup"]
    
    PyMOLObjectNames["Chains"][ChainID] = {}
    PyMOLObjectNames["Ligands"][ChainID] = {}
    
    # Set up chain group and chain objects...
    ChainGroupName = "%s.Chain%s" % (PDBGroupName, ChainID)
    PyMOLObjectNames["Chains"][ChainID]["ChainGroup"] = ChainGroupName
    PyMOLObjectNames["PDBGroupMembers"].append(ChainGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"] = []
    
    # Setup chain complex group and objects...
    ChainComplexGroupName = "%s.Complex" % (ChainGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroup"] = ChainComplexGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"].append(ChainComplexGroupName)

    PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroupMembers"] = []
    
    Name = "%s.Complex" % (ChainComplexGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainComplex"] = Name
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroupMembers"].append(Name)
    
    if OptionsInfo["SurfaceChainComplex"]:
        Name = "%s.Surface" % (ChainComplexGroupName)
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexHydrophobicSurface"] = Name
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroupMembers"].append(Name)

    # Setup up a group for individual chains...
    ChainAloneGroupName = "%s.Chain" % (ChainGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroup"] = ChainAloneGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"].append(ChainAloneGroupName)
        
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"] = []
        
    Name = "%s.Chain" % (ChainAloneGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainAlone"] = Name
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"].append(Name)

    # Setup a group for non-interface residues...
    NonInterfaceGroupName = "%s.NonInterface" % (ChainAloneGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneNonInterfaceGroup"] = NonInterfaceGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"].append(NonInterfaceGroupName)
    
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneNonInterfaceGroupMembers"] = []

    # Setup a chain for non-interface residues...
    Name = "%s.Chain" % (NonInterfaceGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneNonInterface"] = Name
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneNonInterfaceGroupMembers"].append(Name)

    if GetChainAloneContainsSurfacesStatus(FileIndex, ChainID):
        # Setup a surface group and add it to chain alone non-interface group...
        SurfaceGroupName = "%s.Surface" % (NonInterfaceGroupName)
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroup"] = SurfaceGroupName
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneNonInterfaceGroupMembers"].append(SurfaceGroupName)
        
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroupMembers"] = []

        # Setup a generic colored surface...
        Name = "%s.Surface" % (SurfaceGroupName)
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurface"] = Name
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroupMembers"].append(Name)
        
        if GetChainAloneSurfaceChainStatus(FileIndex, ChainID):
            # Setup hydrophobicity surface...
            Name = "%s.Hydrophobicity" % (SurfaceGroupName)
            PyMOLObjectNames["Chains"][ChainID]["ChainAloneHydrophobicSurface"] = Name
            PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroupMembers"].append(Name)
            
            # Setup hydrophobicity and charge surface...
            Name = "%s.Hydrophobicity_Charge" % (SurfaceGroupName)
            PyMOLObjectNames["Chains"][ChainID]["ChainAloneHydrophobicChargeSurface"] = Name
            PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroupMembers"].append(Name)
    
        if GetChainAloneSurfaceChainElectrostaticsStatus(FileIndex, ChainID):
            # Setup electrostatics group...
            GroupName = "%s.Vacuum_Electrostatics" % (SurfaceGroupName)
            PyMOLObjectNames["Chains"][ChainID]["ChainAloneElectrostaticsGroup"] = GroupName
            PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroupMembers"].append(GroupName)
            
            # Setup electrostatics group members...
            PyMOLObjectNames["Chains"][ChainID]["ChainAloneElectrostaticsGroupMembers"] = []
            
            for MemberType in ["Chain", "Contact_Potential", "Map", "Legend", "Volume"]:
                MemberID = re.sub("_", "", MemberType)
                
                Name = "%s.%s" % (GroupName, MemberType)
                NameID = "ChainAloneElectrostaticsSurface%s" % MemberID
                
                PyMOLObjectNames["Chains"][ChainID][NameID] = Name
                PyMOLObjectNames["Chains"][ChainID]["ChainAloneElectrostaticsGroupMembers"].append(Name)
            
    # Setup solvent and inorganic objects for chain...
    for NameID in ["Solvent", "Inorganic"]:
        Name = "%s.%s" % (ChainGroupName, NameID)
        PyMOLObjectNames["Chains"][ChainID][NameID] = Name
        PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"].append(Name)

def SetupPyMOLObjectNamesForLigand(FileIndex, PyMOLObjectNames, ChainID, LigandID):
    """Stetup groups and objects for ligand."""

    PyMOLObjectNames["Ligands"][ChainID][LigandID] = {}
    
    ChainGroupName = PyMOLObjectNames["Chains"][ChainID]["ChainGroup"]
    
    # Setup a chain level ligand group...
    ChainLigandGroupName = "%s.Ligand%s" % (ChainGroupName, LigandID)
    PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroup"] = ChainLigandGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"].append(ChainLigandGroupName)
    
    PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroupMembers"] = []

    # Set up ligand group and its members...
    GroupName = "%s.Ligand" % (ChainLigandGroupName)
    GroupNameID = "LigandGroup"
    GroupMembersID = "LigandGroupMembers"
    
    PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupNameID] = GroupName
    PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroupMembers"].append(GroupName)
    
    LigandName = "%s.Ligand" % (GroupName)
    LigandNameID = "Ligand"
    PyMOLObjectNames["Ligands"][ChainID][LigandID][LigandNameID] = LigandName
    
    PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID] = []
    PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(LigandName)
    
    # Add ball and stick...
    BallAndStickName = "%s.BallAndStick" % (GroupName)
    BallAndStickID = "LigandBallAndStick"
    PyMOLObjectNames["Ligands"][ChainID][LigandID][BallAndStickID] = BallAndStickName
    PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(BallAndStickName)

def SetupPyMOLInterfaceObjectNamesInfo():
    """Setup PyMOL object names for displaying macromolecular interfaces. """

    InterfaceChainsAndResiduesInfo = OptionsInfo["InfilesInfo"]["InterfaceChainsAndResiduesInfo"]
    
    PyMOLObjectNames = {}
    PyMOLObjectNames["InterfaceIDs"] = {}
    PyMOLObjectNames["InterfaceChainIDs"] = {}

    # Setup top level interfaces group...
    InterfacesGroupName = "Interfaces"
    PyMOLObjectNames["InterfacesGroup"] = InterfacesGroupName
    PyMOLObjectNames["InterfacesGroupMembers"] = []

    for InterfaceID in InterfaceChainsAndResiduesInfo["InterfaceIDs"]:
        PyMOLObjectNames["InterfaceIDs"][InterfaceID] = {}
        PyMOLObjectNames["InterfaceChainIDs"][InterfaceID] = {}
        
        # Setup an interface group...
        InterfaceIDGroupName = "%s.%s" % (InterfacesGroupName, InterfaceID)
        PyMOLObjectNames["InterfaceIDs"][InterfaceID]["InterfaceIDGroup"] = InterfaceIDGroupName
        PyMOLObjectNames["InterfacesGroupMembers"].append(InterfaceIDGroupName)
        PyMOLObjectNames["InterfaceIDs"][InterfaceID]["InterfaceIDGroupMembers"] = []

        # Setup a polar contact group...
        if OptionsInfo["InterfacePolarContacts"]:
            PolarContactName = "%s.Polar_Contacts" % (InterfaceIDGroupName)
            PyMOLObjectNames["InterfaceIDs"][InterfaceID]["PolarContacts"] = PolarContactName
            PyMOLObjectNames["InterfaceIDs"][InterfaceID]["InterfaceIDGroupMembers"].append(PolarContactName)
        
        # Setup a hydrophobic contact group...
        if OptionsInfo["InterfaceHydrophobicContacts"]:
            HydrophobicContactsName = "%s.Hydrophobic_Contacts" % (InterfaceIDGroupName)
            PyMOLObjectNames["InterfaceIDs"][InterfaceID]["HydrophobicContacts"] = HydrophobicContactsName
            PyMOLObjectNames["InterfaceIDs"][InterfaceID]["InterfaceIDGroupMembers"].append(HydrophobicContactsName)
        
        for InterfaceChainID in InterfaceChainsAndResiduesInfo["InterfaceChainIDs"][InterfaceID]:
            SetupPyMOLInterfaceObjectNamesForChain(PyMOLObjectNames, InterfaceID, InterfaceChainID)
            
    OptionsInfo["InfilesInfo"]["PyMOLInterfaceObjectNamesInfo"] = PyMOLObjectNames

def SetupPyMOLInterfaceObjectNamesForChain(PyMOLObjectNames, InterfaceID, InterfaceChainID):
    """Setup PyMOL interface object names for a chain. """
    
    InterfaceChainsAndResiduesInfo = OptionsInfo["InfilesInfo"]["InterfaceChainsAndResiduesInfo"]
    FileIndex = InterfaceChainsAndResiduesInfo["ChainIDsInfileIndices"][InterfaceID][InterfaceChainID]
    ChainID = InterfaceChainsAndResiduesInfo["ChainIDs"][InterfaceID][InterfaceChainID]
    
    InterfaceIDGroupName = PyMOLObjectNames["InterfaceIDs"][InterfaceID]["InterfaceIDGroup"]
    
    PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID] = {}
    
    # Setup a chain group...
    ChainGroupName = "%s.%s" % (InterfaceIDGroupName, InterfaceChainID)
    PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainGroup"] = ChainGroupName
    PyMOLObjectNames["InterfaceIDs"][InterfaceID]["InterfaceIDGroupMembers"].append(ChainGroupName)
    PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainGroupMembers"] = []
    
    # Setup chain...
    Name = "%s.Chain" % (ChainGroupName)
    PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["Chain"] = Name
    PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainGroupMembers"].append(Name)
    
    if GetInterfaceResidueTypesStatus(FileIndex, ChainID):
        # Setup residue type group and its subgroups...
        ResiduesGroupName = "%s.Residues" % (ChainGroupName)
        
        ResiduesGroupIDPrefix = "ChainResidues"
        ResiduesGroupID = "%sGroup" % ResiduesGroupIDPrefix
        
        # Add residue group to chain group...
        PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesGroupID] = ResiduesGroupName
        PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainGroupMembers"].append(ResiduesGroupName)
        
        # Initialize residue group members...
        ResiduesGroupMembersID = "%sGroupMembers" % ResiduesGroupIDPrefix
        PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesGroupMembersID] = []
        
        # Setup residues sub groups and its members...
        for SubGroupType in ["Aromatic", "Hydrophobic", "Polar", "Positively_Charged", "Negatively_Charged", "Other"]:
            SubGroupID = re.sub("_", "", SubGroupType)

            ResiduesSubGroupName = "%s.%s" % (ResiduesGroupName, SubGroupType)
            ResiduesSubGroupID = "%s%sGroup" % (ResiduesGroupIDPrefix, SubGroupID)

            # Add sub group to residues group...
            PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesSubGroupID] = ResiduesSubGroupName
            PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesGroupMembersID].append(ResiduesSubGroupName)
            
            # Initialize sub group members...
            ResiduesSubGroupMembersID = "%s%sGroupMembers" % (ResiduesGroupIDPrefix, SubGroupID)
            PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesSubGroupMembersID] = []
            
            # Add sub group members to subgroup...
            for MemberType in ["Residues", "Surface"]:
                MemberID = re.sub("_", "", MemberType)

                SubGroupMemberName = "%s.%s" % (ResiduesSubGroupName, MemberType)
                SubGroupMemberID = "%s%s%s" % (ResiduesGroupIDPrefix, SubGroupID, MemberID)
                
                PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][SubGroupMemberID] = SubGroupMemberName
                PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][ResiduesSubGroupMembersID].append(SubGroupMemberName)
                
    if GetInterfaceContainsSurfacesStatus(FileIndex, ChainID):
        # Setup a surface group and add it to chain group...
        SurfaceGroupName = "%s.Surface" % (ChainGroupName)
        PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainSurfaceGroup"] = SurfaceGroupName
        PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainGroupMembers"].append(SurfaceGroupName)
        
        PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainSurfaceGroupMembers"] = []

        # Setup a generic colored surface...
        Name = "%s.Surface" % (SurfaceGroupName)
        PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainSurface"] = Name
        PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainSurfaceGroupMembers"].append(Name)
        
        if GetInterfaceSurfaceChainStatus(FileIndex, ChainID):
            # Setup hydrophobicity surface...
            Name = "%s.Hydrophobicity" % (SurfaceGroupName)
            PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainHydrophobicSurface"] = Name
            PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainSurfaceGroupMembers"].append(Name)
            
            # Setup hydrophobicity and charge surface...
            Name = "%s.Hydrophobicity_Charge" % (SurfaceGroupName)
            PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainHydrophobicChargeSurface"] = Name
            PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainSurfaceGroupMembers"].append(Name)
            
        if GetInterfaceSurfaceChainElectrostaticsStatus(FileIndex, ChainID):
            # Setup electrostatics group...
            GroupName = "%s.Vacuum_Electrostatics" % (SurfaceGroupName)
            PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainElectrostaticsGroup"] = GroupName
            PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainSurfaceGroupMembers"].append(GroupName)
            
            # Setup electrostatics group members...
            PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainElectrostaticsGroupMembers"] = []
            
            for MemberType in ["Chain", "Contact_Potential", "Map", "Legend", "Volume"]:
                MemberID = re.sub("_", "", MemberType)
                
                Name = "%s.%s" % (GroupName, MemberType)
                NameID = "ChainElectrostaticsSurface%s" % MemberID
                
                PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID][NameID] = Name
                PyMOLObjectNames["InterfaceChainIDs"][InterfaceID][InterfaceChainID]["ChainElectrostaticsGroupMembers"].append(Name)

def ProcessResidueTypesAndSurfaceOptions(FileIndex, SpecifiedChainsAndLigandsInfo):
    """Process residue types and surface options for chains and interfaces."""

    SpecifiedChainsAndLigandsInfo["ChainSurfaces"] = {}
    SpecifiedChainsAndLigandsInfo["SurfaceChain"] = {}
    SpecifiedChainsAndLigandsInfo["SurfaceChainElectrostatics"] = {}

    SpecifiedChainsAndLigandsInfo["InterfaceSurfaces"] = {}
    SpecifiedChainsAndLigandsInfo["SurfaceInterface"] = {}
    SpecifiedChainsAndLigandsInfo["SurfaceInterfaceElectrostatics"] = {}
    
    SpecifiedChainsAndLigandsInfo["ResidueTypesInterface"] = {}

    # Load infile...
    Infile = OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex]
    MolName = OptionsInfo["InfilesInfo"]["InfilesRoots"][FileIndex]
    pymol.cmd.load(Infile, MolName)
    
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
        AminoAcidsPresent = PyMOLUtil.AreAminoAcidResiduesPresent(MolName, ChainID)

        # Process surfaces for chains...
        if re.match("^auto$", OptionsInfo["SurfaceChain"], re.I):
            SurfaceChain = True if AminoAcidsPresent else False
        else:
            SurfaceChain = True if re.match("^yes$", OptionsInfo["SurfaceChain"], re.I) else False
        SpecifiedChainsAndLigandsInfo["SurfaceChain"][ChainID] = SurfaceChain

        if re.match("^auto$", OptionsInfo["SurfaceChainElectrostatics"], re.I):
            SurfaceChainElectrostatics = True if AminoAcidsPresent else False
        else:
            SurfaceChainElectrostatics = True if re.match("^yes$", OptionsInfo["SurfaceChainElectrostatics"], re.I) else False
        SpecifiedChainsAndLigandsInfo["SurfaceChainElectrostatics"][ChainID] = SurfaceChainElectrostatics

        # A colored surface is always created by default...
        ChainSurfaces = True
        SpecifiedChainsAndLigandsInfo["ChainSurfaces"][ChainID] = ChainSurfaces
        
        # Process residue types and surfaces for interfaces...
        if re.match("^auto$", OptionsInfo["InterfaceSurface"], re.I):
            InterfaceSurface = True if AminoAcidsPresent else False
        else:
            InterfaceSurface = True if re.match("^yes$", OptionsInfo["InterfaceSurface"], re.I) else False
        SpecifiedChainsAndLigandsInfo["SurfaceInterface"][ChainID] = InterfaceSurface
        
        if re.match("^auto$", OptionsInfo["InterfaceSurfaceElectrostatics"], re.I):
            InterfaceSurfaceElectrostatics = True if AminoAcidsPresent else False
        else:
            InterfaceSurfaceElectrostatics = True if re.match("^yes$", OptionsInfo["InterfaceSurfaceElectrostatics"], re.I) else False
        SpecifiedChainsAndLigandsInfo["SurfaceInterfaceElectrostatics"][ChainID] = InterfaceSurfaceElectrostatics
        
        # A colored surface is always created by default...
        InterfaceSurfaces = True
        SpecifiedChainsAndLigandsInfo["InterfaceSurfaces"][ChainID] = InterfaceSurfaces
        
        if re.match("^auto$", OptionsInfo["InterfaceResidueTypes"], re.I):
            InterfaceResidueTypes = True if AminoAcidsPresent else False
        else:
            InterfaceResidueTypes = True if re.match("^yes$", OptionsInfo["InterfaceResidueTypes"], re.I) else False
        SpecifiedChainsAndLigandsInfo["ResidueTypesInterface"][ChainID] = InterfaceResidueTypes
    
    # Delete loaded object...
    pymol.cmd.delete(MolName)

def GetInterfaceResidueTypesStatus(FileIndex, ChainID):
    """Get status of residue types for an interface."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["ResidueTypesInterface"][ChainID]

def GetChainAloneContainsSurfacesStatus(FileIndex, ChainID):
    """Get status of surfaces present in chain alone object."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["ChainSurfaces"][ChainID]

def GetInterfaceContainsSurfacesStatus(FileIndex, ChainID):
    """Get status of surfaces present in an interface."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["InterfaceSurfaces"][ChainID]

def GetChainAloneSurfaceChainStatus(FileIndex, ChainID):
    """Get status of hydrophobic surfaces for chain alone object."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["SurfaceChain"][ChainID]

def GetChainAloneSurfaceChainElectrostaticsStatus(FileIndex, ChainID):
    """Get status of electrostatics surfaces for chain alone object."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["SurfaceChainElectrostatics"][ChainID]

def GetInterfaceSurfaceChainStatus(FileIndex, ChainID):
    """Get status of hydrophobic surfaces for an interface."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["SurfaceInterface"][ChainID]

def GetInterfaceSurfaceChainElectrostaticsStatus(FileIndex, ChainID):
    """Get status of hydrophobic surfaces for an interface."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["SurfaceInterfaceElectrostatics"][ChainID]

def CheckPresenceOfValidLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo):
    """Check presence of valid ligand IDs."""

    MiscUtil.PrintInfo("\nSpecified chain IDs: %s" % (", ".join(SpecifiedChainsAndLigandsInfo["ChainIDs"])))
    
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
        if len (SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID]):
            MiscUtil.PrintInfo("Chain ID: %s; Specified LigandIDs: %s" % (ChainID, ", ".join(SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID])))
        else:
            MiscUtil.PrintInfo("Chain IDs: %s; Specified LigandIDs: None" % (ChainID))
            MiscUtil.PrintWarning("No valid ligand IDs found for chain ID, %s. PyMOL groups and objects related to ligand and binding pockect won't be created." % (ChainID))

def RetrieveFirstChainID(FileIndex):
    """Get first chain ID."""
    
    ChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["ChainsAndLigandsInfo"][FileIndex]
    
    FirstChainID = None
    if len(ChainsAndLigandsInfo["ChainIDs"]):
        FirstChainID = ChainsAndLigandsInfo["ChainIDs"][0]
    
    return FirstChainID

def ProcessResidueTypes():
    """Process residue types. """

    ResidueTypesNamesInfo, ResidueTypesParamsInfo = PyMOLUtil.ProcessResidueTypesOptionsInfo("-r, --residueTypes", OptionsInfo["ResidueTypes"])
    OptionsInfo["ResidueTypesNames"] = ResidueTypesNamesInfo
    OptionsInfo["ResidueTypesParams"] = ResidueTypesParamsInfo
    
def ProcessSurfaceAtomTypesColors():
    """Process surface atom types colors. """

    AtomTypesColorNamesInfo = PyMOLUtil.ProcessSurfaceAtomTypesColorsOptionsInfo("--surfaceAtomTypesColors", OptionsInfo["SurfaceAtomTypesColors"])
    OptionsInfo["AtomTypesColorNames"] = AtomTypesColorNamesInfo

def ProcessOptions():
    """Process and validate command line arguments and options"""

    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["AllowEmptyObjects"] = True if re.match("^Yes$", Options["--allowEmptyObjects"], re.I) else False

    OptionsInfo["Infiles"] = Options["--infiles"]
    OptionsInfo["InfilesNames"] =  Options["--infileNames"]

    OptionsInfo["Overwrite"] = Options["--overwrite"]
    OptionsInfo["PMLOut"] = True if re.match("^Yes$", Options["--PMLOut"], re.I) else False
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Outfile"])
    OptionsInfo["PSEOut"] = False 
    if re.match("^pml$", FileExt, re.I):
        OptionsInfo["PMLOutfile"] = OptionsInfo["Outfile"] 
        OptionsInfo["PMEOutfile"] = re.sub(".pml$", ".pme", OptionsInfo["Outfile"]) 
    elif re.match("^pse$", FileExt, re.I):
        OptionsInfo["PSEOut"] = True 
        OptionsInfo["PSEOutfile"] = OptionsInfo["Outfile"] 
        OptionsInfo["PMLOutfile"] = re.sub(".pse$", ".pml", OptionsInfo["Outfile"]) 
        if os.path.exists(OptionsInfo["PMLOutfile"]) and (not OptionsInfo["Overwrite"]):
            MiscUtil.PrintError("The intermediate output file to be generated, %s, already exist. Use option \"--ov\" or \"--overwrite\" and try again." % OptionsInfo["PMLOutfile"] )

    OptionsInfo["InterfaceLabelColor"] = Options["--interfaceLabelColor"]
    
    OptionsInfo["InterfaceContactsCutoff"] = float(Options["--interfaceContactsCutoff"])
    
    OptionsInfo["InterfaceHydrophobicContactsColor"] = Options["--interfaceHydrophobicContactsColor"]
    OptionsInfo["InterfaceHydrophobicContacts"] = True if re.match("^Yes$", Options["--interfaceHydrophobicContacts"], re.I) else False
    
    OptionsInfo["InterfacePolarContactsColor"] = Options["--interfacePolarContactsColor"]
    OptionsInfo["InterfacePolarContacts"] = True if re.match("^Yes$", Options["--interfacePolarContacts"], re.I) else False
    
    OptionsInfo["InterfaceResidueTypes"] = Options["--interfaceResidueTypes"]
    OptionsInfo["InterfaceSurface"] = Options["--interfaceSurface"]
    OptionsInfo["InterfaceSurfaceElectrostatics"] = Options["--interfaceSurfaceElectrostatics"]
    
    OptionsInfo["LabelFontID"] = int(Options["--labelFontID"])
    
    Method = Options["--method"]
    MethodCutoff = Options["--methodCutoff"]
    if re.match("^auto$", MethodCutoff, re.I):
        if re.match("BySASAChange", Method, re.I):
            MethodCutoff = 1.0
        elif re.match("ByHeavyAtomsDistance", Method, re.I):
            MethodCutoff = 5.0
        elif re.match("ByCAlphaAtomsDistance", Method, re.I):
            MethodCutoff = 8.0
        else:
            MiscUtil.PrintError("The method, %s, specified using \"-m, --method\" option  is  not a valid method." % (Method))
    else:
        MethodCutoff = float(Options["--methodCutoff"])
    OptionsInfo["Method"] = Method
    OptionsInfo["MethodCutoff"] = MethodCutoff
    
    OptionsInfo["ResidueTypes"] = Options["--residueTypes"]
    ProcessResidueTypes()
    
    OptionsInfo["SurfaceChain"] = Options["--surfaceChain"]
    OptionsInfo["SurfaceChainElectrostatics"] = Options["--surfaceChainElectrostatics"]
    
    OptionsInfo["SurfaceChainComplex"] = True if re.match("^Yes$", Options["--surfaceChainComplex"], re.I) else False
    OptionsInfo["SurfaceComplex"] = True if re.match("^Yes$", Options["--surfaceComplex"], re.I) else False

    # Retrieve surface colors...
    SurfaceColors = re.sub(" ", "", Options["--surfaceColors"])
    SurfaceColorsWords = SurfaceColors.split(",")
    if len(SurfaceColorsWords) != 2:
        MiscUtil.PrintError("The number of comma delinited color names, %d, specified using \"--surfaceColors\" option, \"%s\",  must be a 2." % (len(SurfaceColorsWords), Options["--surfaceColors"]))
    OptionsInfo["SurfaceColors"] = SurfaceColors
    OptionsInfo["SurfaceInterfaceColor"] = SurfaceColorsWords[0]
    OptionsInfo["SurfaceNonInterfaceColor"] = SurfaceColorsWords[1]
    
    OptionsInfo["SurfaceColorPalette"] = Options["--surfaceColorPalette"]
    OptionsInfo["SurfaceAtomTypesColors"] = Options["--surfaceAtomTypesColors"]
    ProcessSurfaceAtomTypesColors()
    
    OptionsInfo["SurfaceTransparency"] = float(Options["--surfaceTransparency"])
    
    RetrieveInfilesInfo()

    # Retrieve interface chain IDs...
    InterfaceChainIDs = re.sub(" ", "", Options["--chainIDs"])
    InterfaceChainIDsList = []
    CanonicalInterfaceChainIDsMap = {}
    if not re.match("^auto$", InterfaceChainIDs, re.I):
        InterfaceChainIDsWords = InterfaceChainIDs.split(",")
        if len(InterfaceChainIDsWords) % 2:
            MiscUtil.PrintError("The number of comma delimited chain IDs, %d, specified using \"-c, --chainIDs\" option, \"%s\",  must be a multple of 2." % (len(InterfaceChainIDsWords), Options["--chainIDs"]))
        for ChainID in InterfaceChainIDsWords:
            if not len(ChainID):
                MiscUtil.PrintError("A chain ID specified, \"%s\" using \"-c, --chainIDs\" option is empty." % (Options["--chainIDs"]))
            ChainIDWords = ChainID.split("+")
            InterfaceChainIDsList.append(ChainIDWords)
    OptionsInfo["InterfaceChainIDs"]  = InterfaceChainIDs
    OptionsInfo["InterfaceChainIDsList"]  = InterfaceChainIDsList
    
    # Process interface chain IDs...
    OptionsInfo["LigandIDs"] = Options["--ligandIDs"]
    ProcessInterfaceChainIDs()
    
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
    
    MiscUtil.ValidateOptionTextValue("--allowEmptyObjects", Options["--allowEmptyObjects"], "yes no")

    # Expand infiles to handle presence of multiple input files...
    InfileNames = MiscUtil.ExpandFileNames(Options["--infiles"], ",")
    InfilesCount = len(InfileNames)
    if not InfilesCount:
        MiscUtil.PrintError("No input files specified for \"-i, --infiles\" option")

    # Validate file extensions...
    for Infile in InfileNames:
        MiscUtil.ValidateOptionFilePath("-i, --infiles", Infile)
        MiscUtil.ValidateOptionFileExt("-i, --infiles", Infile, "pdb cif")
        MiscUtil.ValidateOptionsDistinctFileNames("-i, --infiles", Infile, "-o, --outfile", Options["--outfile"])
    
    # Validate file count...
    if InfilesCount > 2:
        MiscUtil.PrintError("Number of input files, %s, specified using \"-i, --infiles\" option is not valid. Number of allowed files: 1 or 2" % (InfilesCount))

    # Validate distinct file names...
    if InfilesCount == 2:
        Infile1Name, Infile2Name = InfileNames
        Infile1Pattern = r"^" + re.escape(Infile1Name) + r"$";
        if re.match(Infile1Pattern, Infile2Name, re.I):
            MiscUtil.PrintError("The file names specified, \"%s\", for option \"-i, --infiles\" must be different.\n" % (Options["--infiles"]))
        
    Options["--infileNames"] = InfileNames
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "pml pse")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])

    MiscUtil.ValidateOptionFloatValue("--interfaceContactsCutoff", Options["--interfaceContactsCutoff"], {">": 0.0})
    MiscUtil.ValidateOptionTextValue("--interfaceHydrophobicContacts", Options["--interfaceHydrophobicContacts"], "yes no")
    MiscUtil.ValidateOptionTextValue("--interfacePolarContacts", Options["--interfacePolarContacts"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--interfaceResidueTypes", Options["--interfaceResidueTypes"], "yes no auto")
    MiscUtil.ValidateOptionTextValue("--interfaceSurface", Options["--interfaceSurface"], "yes no auto")
    MiscUtil.ValidateOptionTextValue("--interfaceSurfaceElectrostatics", Options["--interfaceSurfaceElectrostatics"], "yes no auto")
                       
    MiscUtil.ValidateOptionIntegerValue("--labelFontID", Options["--labelFontID"], {})
    
    MiscUtil.ValidateOptionTextValue("--method", Options["--method"], "BySASAChange ByHeavyAtomsDistance ByCAlphaAtomsDistance")
    if not re.match("^auto$", Options["--methodCutoff"], re.I):
        MiscUtil.ValidateOptionFloatValue("--methodCutoff", Options["--methodCutoff"], {">": 0.0})
    
    MiscUtil.ValidateOptionTextValue("--PMLOut", Options["--PMLOut"], "yes no")

    MiscUtil.ValidateOptionTextValue("--surfaceChain", Options["--surfaceChain"], "yes no auto")
    MiscUtil.ValidateOptionTextValue("--surfaceChainElectrostatics", Options["--surfaceChainElectrostatics"], "yes no auto")
    MiscUtil.ValidateOptionTextValue("--surfaceComplex", Options["--surfaceComplex"], "yes no")
    MiscUtil.ValidateOptionTextValue("--surfaceChainComplex", Options["--surfaceChainComplex"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--surfaceColorPalette", Options["--surfaceColorPalette"], "RedToWhite WhiteToGreen")
    MiscUtil.ValidateOptionFloatValue("--surfaceTransparency", Options["--surfaceTransparency"], {">=": 0.0, "<=": 1.0})
    
# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLVisualizeInterfaces.py - Visualize macromolecular interfaces

Usage:
    PyMOLVisualizeInterfaces.py [--allowEmptyObjects <yes or no>] [--chainIDs <ChainID1 or ChainID1,ChainID2>]
                                              [--interfaceLabelColor <text>] [ --interfaceContactsCutoff <number>]
                                              [--interfaceHydrophobicContacts <yes or no>] [--interfaceHydrophobicContactsColor <text>]
                                              [--interfacePolarContacts <yes or no>] [--interfacePolarContactsColor <text>]
                                              [--interfaceResidueTypes <yes or no>] [--interfaceSurface <yes or no>]
                                              [--interfaceSurfaceElectrostatics <yes or no>] [--labelFontID <number>]
                                              [--ligandIDs <Largest, All, None or ID1,ID2...>] [--method <text>] [--methodCutoff <number>]
                                              [--PMLOut <yes or no>] [--residueTypes <Type,Color,ResNames,...>] [--surfaceChain <yes or no>]
                                              [--surfaceChainElectrostatics <yes or no>] [--surfaceChainComplex <yes or no>]
                                              [--surfaceComplex <yes or no>] [--surfaceColors <ColorName1,ColorName2>]
                                              [--surfaceColorPalette <RedToWhite or WhiteToGreen>]
                                              [--surfaceAtomTypesColors <ColorType,ColorSpec,...>] [--surfaceTransparency <number>]
                                              [--overwrite] [-w <dir>] -i <infile1,...> -o <outfile>
    PyMOLVisualizeInterfaces.py -h | --help | -e | --examples

Description:
    Generate PyMOL visualization files for viewing interfaces between macromolecules
    including proteins and nucleic acids. The interfaces may be generated between
    pairs of chains in a single file or across two different files.

    The supported input file format are: PDB (.pdb), CIF (.cif)

    The supported output file formats are: PyMOL script file (.pml), PyMOL session
    file (.pse)

    A variety of PyMOL groups and objects may be  created for visualization of
    macromolecular interfaces. These groups and objects correspond to complexes,
    surfaces, chains, ligands, and interfaces. A complete hierarchy of all possible
    PyMOL groups and objects is shown below:

        <PDBFileRoot>
            .Complex
                .Complex
                .Surface
            .Chain<ID>
                .Complex
                    .Complex
                    .Surface
                .Chain
                    .Chain
                    .NonInterface
                        .Chain
                        .Surface
                            .Surface
                            .Hydrophobicity
                            .Hydrophobicity_Charge
                            .Vacuum_Electrostatics
                                .Contact_Potentials
                                .Map
                                .Legend
                                .Volume
                .Solvent
                .Inorganic
                .Ligand<ID>
                    .Ligand
                        .Ligand
                        .BallAndStick
                .Ligand<ID>
                    .Ligand
                        ... ... ...
            .Chain<ID>
                ... ... ...
                .Ligand<ID>
                    ... ... ...
                .Ligand<ID>
                    ... ... ...
            .Chain<ID>
                ... ... ...
        <PDBFileRoot>
            .Complex
                ... ... ...
            .Chain<ID>
                ... ... ...
                .Ligand<ID>
                    ... ... ...
                .Ligand<ID>
                    ... ... ...
            .Chain<ID>
                ... ... ...
        <Interfaces>
            .Chain<IDs1>_Chain<IDs2>
                .Polar_Contacts
                .Hydrophobic_Contacts
                .Chain<ID> or Chain<ID>_<PDBFileRoot>
                    .Chain
                    .Residues
                        .Aromatic
                            .Residues
                            .Surface
                        .Hydrophobic
                            .Residues
                            .Surface
                        .Polar
                            .Residues
                            .Surface
                        .Positively_Charged
                            .Residues
                            .Surface
                        .Negatively_Charged
                            .Residues
                            .Surface
                        .Other
                            .Residues
                            .Surface
                    .Surface
                        .Surface
                        .Hydrophobicity
                        .Hydrophobicity_Charge
                        .Vacuum_Electrostatics
                            .Contact_Potentials
                            .Map
                            .Legend
                            .Volume
                .Chain<ID> or <PDBFileRoot>_Chain<ID>
                    .Chain
                    .Residues
                        ... ... ...
                    .Surface
                        ... ... ...
            .Chain<IDs>_Chain<IDs>
                .Polar_Contacts
                .Hydrophobic_Contacts
                .Chain<ID> or Chain<ID>_<PDBFileRoot>
                    .Chain
                    .Residues
                        ... ... ...
                    .Surface
                        ... ... ...
                .Chain<ID> or Chain<ID>_<PDBFileRoot>
                    .Chain
                    .Residues
                        ... ... ...
                    .Surface
                        ... ... ...
    
    The hydrophobic and electrostatic surfaces are not created for complete complex
    and chain complex in input file(s) by default. A word to the wise: The creation of
    surface objects may slow down loading of PML file and generation of PSE file, based
    on the size of input complexes. The generation of PSE file may also fail.

Options:
    --allowEmptyObjects <yes or no>  [default: no]
        Allow creation of empty PyMOL objects corresponding to interface,
        solvent, and inorganic atom selections across chains and ligands in
        input file(s). By default, the empty objects are marked for deletion.
    -c, --chainIDs <ChainID1,ChainD2,...>  [default: Auto]
        Pairwise comma delimited list of chain IDs for the identification of
        macromolecular interfaces. All chain IDs must be present in the
        same file for a single input file. Otherwise, the first and second
        chain ID(s) in a pair belong to the first and second input file.
        
        The default values for interface chain IDs depend on the number
        of input files as shown below:
        
        One input file: First two chains
        Two input files: First chain in each input file
        
        Each chain may contain multiple chain IDs delimited by a plus sign. For
        example, A+B,C+D chain pair specifies interface between chain complexes
        A+B and C+D in first input file or across two input files.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infiles <infile or infile1,infile2>
        Name of an input file or a comma delmited list of names for two input
        files.
    --interfaceLabelColor <text>  [default: magenta]
        Color for drawing residue or atom level labels for residues in an interface.
        The specified value must be valid color. No validation is performed.
    --interfaceContactsCutoff <number>  [default: 4.0]
        Distance in Angstroms for identifying polar and hyrdophobic contacts
        between atoms in interface reisudes.
    --interfaceHydrophobicContacts <yes or no>  [default: yes]
        Hydrophobic contacts between residues in an interface. The hydrophobic
        contacts are shown between pairs of carbon atoms not connected to
        hydrogen bond donor or acceptors atoms as identified by PyMOL.
    --interfaceHydrophobicContactsColor <text>  [default: purpleblue]
        Color for drawing hydrophobic contacts between residues in an interface.
        The specified value must be valid color. No validation is performed.
    --interfacePolarContacts <yes or no>  [default: yes]
        Polar contacts between residues in an interface.
    --interfacePolarContactsColor <text>  [default: orange]
        Color for drawing polar contacts between residues in an interface.
        The specified value must be valid color. No validation is performed.
    --interfaceResidueTypes <yes or no>  [default: auto]
        Interface residue types. The residue groups are generated using residue types,
        colors, and names specified by '--residueTypes' option. It is only valid for
        amino acids.  By default, the residue type groups are automatically created
        for interfaces containing amino acids and skipped for chains only containing
        nucleic acids.
    --interfaceSurface <yes or no>  [default: auto]
        Surfaces around interface residues colored by hydrophobicity alone and
        both hydrophobicity and charge. The hydrophobicity surface is colored
        at residue level using Eisenberg hydrophobicity scale for residues and color
        gradient specified by '--surfaceColorPalette' option. The  hydrophobicity and
        charge surface is colored [ Ref 140 ] at atom level using colors specified for
        groups of atoms by '--surfaceAtomTypesColors' option. This scheme allows
        simultaneous mapping of hyrophobicity and charge values on the surfaces.
        
        This option is only valid for amino acids. By default, both surfaces are
        automatically created for pockets containing amino acids and skipped for
        pockets containing only nucleic acids.
        
        In addition, generic surfaces colored by '--surfaceColors' are always created
        for interface residues containing amino acids and nucleic acids.
    --interfaceSurfaceElectrostatics <yes or no>  [default: no]
        Vacuum electrostatics contact potential surface around interface residues.
        A word to the wise from PyMOL documentation: The computed protein
        contact potentials are only qualitatively useful, due to short cutoffs,
        truncation, and lack of solvent "screening".
        
        This option is only valid for amino acids. By default, the electrostatics surface
        is automatically created for chains containing amino acids and skipped for chains
        containing only nucleic acids.
    --labelFontID <number>  [default: 7]
        Font ID for drawing labels. Default: 7 (Sans Bold). Valid values: 5 to 16.
        The specified value must be a valid PyMOL font ID. No validation is
        performed. The complete lists of valid font IDs is available at:
        pymolwiki.org/index.php/Label_font_id. Examples: 5 - Sans;
        7 - Sans Bold; 9 - Serif; 10 - Serif Bold.
    -l, --ligandIDs <Largest, All, None or ID1,ID2...>  [default: All]
        List of ligand IDs to show in chains during visualization of interfaces. Possible
        values: Largest, All, None, or a comma delimited list of ligand IDs. The
        default is to show all ligands present in chains involved in interfaces.
        
        Ligands are identified using organic selection operator available in PyMOL.
        It'll also  identify buffer molecules as ligands. The largest ligand contains
        the highest number of heavy atoms.
    -m, --method <text>  [default: BySASAChange]
        Methodology for the identification of interface residues between a pair
        of chains in an input file. The interface residues may be identified by
        change in solvent accessible surface area (SASA) for a residue between
        a chain and chains complex, distance between heavy atoms
        in two chains, or distance between CAlpha atoms. Possible values:
        BySASAChange, ByHeavyAtomsDistance, or ByCAlphaAtomsDistance. 
    --methodCutoff <number>  [default: auto]
        Cutoff value used by different methodologies during the identification of
        interface residues between a pair of chains. The default values are
        shown below:
        
            BySASAChange: 1.0; Units: Angstrom**2 [ Ref 141 ]
            ByHeavyAtomsDistance: 5.0; Units: Angstrom [ Ref 142 ]
            ByCAlphaAtomsDistance: 8.0; Units: Angstrom [ Ref 143 ]
        
    -o, --outfile <outfile>
        Output file name.
    -p, --PMLOut <yes or no>  [default: yes]
        Save PML file during generation of PSE file.
    -r, --residueTypes <Type,Color,ResNames,...>  [default: auto]
        Residue types, colors, and names to generate for residue groups during
        and '--residueTypesChain' option. It is only valid for amino acids.
        
        It is a triplet of comma delimited list of amino acid residues type, residues
        color, and a space delimited list three letter residue names. 
        
        The default values for residue type, color, and name triplets  are shown
        below:
            
            Aromatic,brightorange,HIS PHE TRP TYR,
            Hydrophobic,orange,ALA GLY VAL LEU ILE PRO MET,
            Polar,palegreen,ASN GLN SER THR CYS,
            Positively_Charged,marine,ARG LYS,
            Negatively_Charged,red,ASP GLU
            
        The color name must be a valid PyMOL name. No validation is performed.
        An amino acid name may appear across multiple residue types. All other
        residues are grouped under 'Other'.
    --surfaceChain <yes or no>  [default: auto]
        Surfaces around non-interface residues in individual  chain colored by
        hydrophobicity alone and both hydrophobicity and charge. The hydrophobicity
        surface is colored at residue level using Eisenberg hydrophobicity scale for residues
        and color gradient specified by '--surfaceColorPalette' option. The  hydrophobicity
        and charge surface is colored [ Ref 140 ] at atom level using colors specified for
        groups of atoms by '--surfaceAtomTypesColors' option. This scheme allows
        simultaneous mapping of hyrophobicity and charge values on the surfaces.
        
        This option is only valid for amino acids. By default, both surfaces are
        automatically created for chains containing amino acids and skipped for
        chains containing only nucleic acids.
        
        In addition, generic surfaces colored by '--surfaceColors' are always created
        for non-interface residues containing amino acids and nucleic acids.
    --surfaceChainElectrostatics <yes or no>  [default: no]
        Vacuum electrostatics contact potential surface and volume around non-interface
        residues in individual chain. A word to the wise from PyMOL documentation: The
        computed protein contact potentials are only qualitatively useful, due to short cutoffs,
        truncation, and lack of solvent "screening".
        
        This option is only valid for amino acids. By default, the electrostatics surface
        and volume are automatically created for chains containing amino acids and
        skipped for chains containing only nucleic acids.
    --surfaceChainComplex <yes or no>  [default: no]
        Hydrophobic surface around chain complex. The  surface is colored by
        hydrophobicity. It is only valid for amino acids.
    --surfaceComplex <yes or no>  [default: no]
        Hydrophobic surface around complete complex. The  surface is colored by
        hydrophobicity. It is only valid for amino acids.
    --surfaceColors <ColorName1,ColorName2>  [default: salmon,lightblue]
        Color names for surfaces around interface residues and non-interface
        residues in chains. These colors are not used for surfaces colored by
        hydrophobicity and charge. The color names must be valid PyMOL names.
    --surfaceColorPalette <RedToWhite or WhiteToGreen>  [default: RedToWhite]
        Color palette for hydrophobic surfaces around chains and interface residues
        in proteins. Possible values: RedToWhite or WhiteToGreen from most
        hydrophobic amino acid to least hydrophobic. The colors values for amino
        acids are taken from color_h script available as part of the Script Library at
        PyMOL Wiki.
    --surfaceAtomTypesColors <ColorType,ColorSpec,...>  [default: auto]
        Atom colors for generating surfaces colored by hyrophobicity and charge
        around chains and interface residues in proteins. It's a pairwise comma
        delimited list of atom color type and color specification for goups of atoms.
        
        The default values for color types [ Ref 140 ] along wth color specifications
        are shown below: 
            
            HydrophobicAtomsColor, yellow,
            NegativelyChargedAtomsColor, red,
            PositivelyChargedAtomsColor, blue,
            OtherAtomsColor, gray90
            
        The color names must be valid PyMOL names.
        
        The color values may also be specified as space delimited RGB triplets:
             
            HydrophobicAtomsColor, 0.95 0.78 0.0,
            NegativelyChargedAtomsColor, 1.0 0.4 0.4,
            PositivelyChargedAtomsColor, 0.2 0.5 0.8,
            OtherAtomsColor, 0.95 0.95 0.95
            
    --surfaceTransparency <number>  [default: 0.25]
        Surface transparency for molecular surfaces.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To visualize interface residues between the first two chains in a PDB file,
    using default methodology to identify interfaces, and and generate a PML
    file, type: 

        % PyMOLVisualizeInterfaces.py -i Sample8.pdb -o Sample8.pml

    To visualize interface residues between a pair of specific chains in a PDB
    file using a specific methodology and cutoff value to identify interfaces, and
    generate a PML file, type: 

        % PyMOLVisualizeInterfaces.py -m BySASAChange --methodCutoff 1.0
        -c "A,B" -i Sample8.pdb -o Sample8.pml

    To visualize interface residues between multiple pairs of specified chains in
    a PDB file using a specific methodology and cutoff value to identify interfaces,
    and generate a PML file, type: 

        % PyMOLVisualizeInterfaces.py -m ByHeavyAtomsDistance
        --methodCutoff 5.0 -c "A,B,B,D" -i Sample8.pdb -o Sample8.pml

    To visualize interface residues between a pair of specified chains, each member
    containing multiple chains, a PDB file using a specific methodology and cutoff
    value to identify interfaces, and generate a PML file, type: 

        % PyMOLVisualizeInterfaces.py -m ByCAlphaAtomsDistance
        --methodCutoff 8.0 -c "A+C,B+D" -i Sample8.pdb -o Sample8.pml

    To visualize interface residues between a pair of specific chains across two PDB
    files using a specific methodology and cutoff value to identify interfaces, and
    generate a PML file, type: 

        % PyMOLVisualizeInterfaces.py -m BySASAChange --methodCutoff 1.0 
        -c "A,B" -i Sample8Part1.pdb,Sample8Part2.pdb
        -o Sample8.pml

    To visualize interface residues between multiple pairs of specified chains across
    two PDB files using a specific methodology and cutoff value to identify interfaces,
    and generate a PML file, type: 

        % PyMOLVisualizeInterfaces.py -m ByHeavyAtomsDistance
        --methodCutoff 5.0  -c "A,B,C,B" -i Sample8Part1.pdb,Sample8Part2.pdb
        -o Sample8.pml

    To visualize interface residues between a pair of specified chains, each member
    containing multiple chains, across two PDB files using a specific methodology
    and cutoff value to identify interfaces, and generate a PML file, type:

        % PyMOLVisualizeInterfaces.py -m ByCAlphaAtomsDistance
        --methodCutoff 8.0  -c "A+C,B+D" -i "Sample8Part1.pdb,Sample8Part2.pdb"
        -o Sample8.pml

Author:
    Manish Sud(msud@san.rr.com)

See also:
    DownloadPDBFiles.pl,  PyMOLVisualizeCryoEMDensity.py,
    PyMOLVisualizeElectronDensity.py, PyMOLVisualizeMacromolecules.py,
    PyMOLVisualizeSurfaceAndBuriedResidues.py

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
