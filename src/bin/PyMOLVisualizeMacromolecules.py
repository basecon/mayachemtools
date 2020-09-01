#!/usr/bin/env python
#
# File: PyMOLVisualizeMacromolecules.py
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
    GenerateMacromolecularVisualization()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def GenerateMacromolecularVisualization():
    """Generate macromolecular visualization."""

    Outfile = OptionsInfo["PMLOutfile"]
    OutFH = open(Outfile, "w")
    if OutFH is None:
        MiscUtil.PrintError("Failed to open output fie %s " % Outfile)
    
    MiscUtil.PrintInfo("\nGenerating file %s..." % Outfile)

    # Setup header...
    WritePMLHeader(OutFH, ScriptName)
    WritePyMOLParameters(OutFH)
    
    # Load reffile for alignment..
    if OptionsInfo["Align"]:
        WriteAlignReference(OutFH)

    # Setup view for each input file...
    FirstComplex = True
    FirstComplexFirstChainName = None
    for FileIndex in range(0, len(OptionsInfo["InfilesInfo"]["InfilesNames"])):
        # Setup PyMOL object names...
        PyMOLObjectNames = SetupPyMOLObjectNames(FileIndex)

        # Setup complex view...
        WriteComplexView(OutFH, FileIndex, PyMOLObjectNames, FirstComplex)
        
        SpecifiedChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]
        FirstChain = True
        for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
            if FirstComplex and FirstChain:
                FirstComplexFirstChainName = PyMOLObjectNames["Chains"][ChainID]["ChainAlone"]
                
            WriteChainView(OutFH, FileIndex, PyMOLObjectNames, ChainID)
            
            # Setup ligand views...
            FirstLigand = True
            for LigandID in SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID]:
                WriteChainLigandView(OutFH, FileIndex, PyMOLObjectNames, ChainID, LigandID)
                
                # Set up ligand level group...
                Enable, Action = [False, "close"]
                if FirstLigand:
                    FirstLigand = False
                    Enable, Action = [True, "open"]
                GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroup"], PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroupMembers"], Enable, Action)
            
            # Setup Chain level group...
            Enable, Action = [False, "close"]
            if FirstChain:
                FirstChain = False
                Enable, Action = [True, "open"]
            GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"], Enable, Action)
    
        # Set up complex level group...
        Enable, Action = [False, "close"]
        if FirstComplex:
            FirstComplex = False
            Enable, Action = [True, "open"]
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["PDBGroup"], PyMOLObjectNames["PDBGroupMembers"], Enable, Action)
        
        # Delete empty PyMOL objects...
        DeleteEmptyPyMOLObjects(OutFH, FileIndex, PyMOLObjectNames)
        
    if OptionsInfo["Align"]:
        DeleteAlignReference(OutFH)

    if FirstComplexFirstChainName is not None:
        OutFH.write("""\ncmd.orient("%s", animate = -1)\n""" % FirstComplexFirstChainName)
    else:
        OutFH.write("""\ncmd.orient("visible", animate = -1)\n""")
    
    OutFH.close()

    # Generate PSE file as needed...
    if OptionsInfo["PSEOut"]:
        GeneratePyMOLSessionFile()

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
    
def WriteAlignReference(OutFH):
    """Setup object for alignment reference """

    RefFileInfo = OptionsInfo["RefFileInfo"]
    RefFile = RefFileInfo["RefFileName"]
    RefName = RefFileInfo["PyMOLObjectName"]
    
    PMLCmds = []
    PMLCmds.append("""cmd.load("%s", "%s")""" % (RefFile, RefName))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (RefName))
    PMLCmds.append("""cmd.disable("%s")""" % (RefName))
    PML = "\n".join(PMLCmds)
    
    OutFH.write("""\n""\n"Loading %s and setting up view for align reference..."\n""\n""" % RefFile)
    OutFH.write("%s\n" % PML)
    
def WriteAlignComplex(OutFH, FileIndex, PyMOLObjectNames):
    """Setup alignment of complex to reference"""

    RefFileInfo = OptionsInfo["RefFileInfo"]
    RefName = RefFileInfo["PyMOLObjectName"]
    
    ComplexName = PyMOLObjectNames["Complex"]
    
    if re.match("^FirstChain$", OptionsInfo["AlignMode"], re.I):
        RefFirstChainID = RefFileInfo["ChainsAndLigandsInfo"]["ChainIDs"][0]
        RefAlignSelection = "%s and chain %s" % (RefName, RefFirstChainID)
        
        ComplexFirstChainID = RetrieveFirstChainID(FileIndex)
        ComplexAlignSelection = "%s and chain %s" % (ComplexName, ComplexFirstChainID)
    else:
        RefAlignSelection = RefName
        ComplexAlignSelection = ComplexName

    PML = PyMOLUtil.SetupPMLForAlignment(OptionsInfo["AlignMethod"], RefAlignSelection, ComplexAlignSelection)
    OutFH.write("""\n""\n"Aligning %s against reference %s ..."\n""\n""" % (ComplexAlignSelection, RefAlignSelection))
    OutFH.write("%s\n" % PML)
    
def DeleteAlignReference(OutFH):
    """Delete alignment reference object."""
    
    RefName = OptionsInfo["RefFileInfo"]["PyMOLObjectName"]
    OutFH.write("""\n""\n"Deleting alignment reference object %s..."\n""\n""" % RefName)
    OutFH.write("""cmd.delete("%s")\n""" % RefName)

def WriteComplexView(OutFH, FileIndex, PyMOLObjectNames, FirstComplex):
    """Write out PML for viewing polymer complex."""

    # Setup complex...
    Infile = OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex]
    PML = PyMOLUtil.SetupPMLForPolymerComplexView(PyMOLObjectNames["Complex"], Infile, True)
    OutFH.write("""\n""\n"Loading %s and setting up view for complex..."\n""\n""" % Infile)
    OutFH.write("%s\n" % PML)

    if OptionsInfo["Align"]:
        # No need to align complex on to itself...
        if not (re.match("^FirstInputFile$", OptionsInfo["AlignRefFile"], re.I) and FirstComplex):
            WriteAlignComplex(OutFH, FileIndex, PyMOLObjectNames)
    
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

    WriteChainAloneChainSelectionsView(OutFH,  FileIndex, PyMOLObjectNames, ChainID)
    WriteChainAloneResidueTypesView(OutFH,  FileIndex, PyMOLObjectNames, ChainID)
        
    if GetChainAloneContainsSurfacesStatus(FileIndex, ChainID):
        # Setup a generic color surface...
        PML = PyMOLUtil.SetupPMLForSurfaceView(PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurface"], ChainName,  Enable = False, Color = OptionsInfo["SurfaceColor"])
        OutFH.write("\n%s\n" % PML)
        
        if GetChainAloneSurfaceChainStatus(FileIndex, ChainID):
            # Setup surface colored by hydrophobicity...
            PML = PyMOLUtil.SetupPMLForHydrophobicSurfaceView(PyMOLObjectNames["Chains"][ChainID]["ChainAloneHydrophobicSurface"], ChainName, ColorPalette = OptionsInfo["SurfaceColorPalette"], Enable = False)
            OutFH.write("\n%s\n" % PML)
            
            # Setup surface colored by hyrdophobicity and charge...
            PML = PyMOLUtil.SetupPMLForHydrophobicAndChargeSurfaceView(PyMOLObjectNames["Chains"][ChainID]["ChainAloneHydrophobicChargeSurface"], ChainName, OptionsInfo["AtomTypesColorNames"]["HydrophobicAtomsColor"], OptionsInfo["AtomTypesColorNames"]["NegativelyChargedAtomsColor"], OptionsInfo["AtomTypesColorNames"]["PositivelyChargedAtomsColor"], OptionsInfo["AtomTypesColorNames"]["OtherAtomsColor"], Enable = False, DisplayAs = None)
            OutFH.write("\n%s\n" % PML)
        
        if GetChainAloneSurfaceChainElectrostaticsStatus(FileIndex, ChainID):
            # Setup electrostatics surface...
            SelectionObjectName = ChainName
            ElectrostaticsGroupName = PyMOLObjectNames["Chains"][ChainID]["ChainAloneElectrostaticsGroup"]
            ElectrostaticsGroupMembers = PyMOLObjectNames["Chains"][ChainID]["ChainAloneElectrostaticsGroupMembers"]
            WriteSurfaceElectrostaticsView("Surface", OutFH, SelectionObjectName, ElectrostaticsGroupName, ElectrostaticsGroupMembers, DisplayAs = "cartoon")

        # Setup surface group...
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroupMembers"], True, "open")

    # Setup chain group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"], True, "close")
    
def WriteChainLigandView(OutFH, FileIndex, PyMOLObjectNames, ChainID, LigandID):
    """Write out PML for viewing ligand in a chain."""
    
    for GroupID in ["Ligand", "Pocket", "PocketSolvent", "PocketInorganic"]:
        ComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
        LigandName = PyMOLObjectNames["Ligands"][ChainID][LigandID]["Ligand"]
        
        # Setup main object...
        GroupTypeObjectID = "%s" % (GroupID)
        GroupTypeObjectName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupTypeObjectID]
        
        if re.match("^Ligand$", GroupID, re.I):
            OutFH.write("""\n""\n"Setting up views for ligand %s in chain %s..."\n""\n""" % (LigandID, ChainID))
            PML = PyMOLUtil.SetupPMLForLigandView(GroupTypeObjectName, ComplexName, LigandID, True)
            OutFH.write("%s\n" % PML)
        elif re.match("^Pocket$", GroupID, re.I):
            OutFH.write("""\n""\n"Setting up views for pocket around ligand %s in chain %s..."\n""\n""" % (LigandID, ChainID))
            PML = PyMOLUtil.SetupPMLForLigandPocketView(GroupTypeObjectName, ComplexName, LigandName, OptionsInfo["PocketDistanceCutoff"], True)
            OutFH.write("%s\n" % PML)
            OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (OptionsInfo["PocketLabelColor"], GroupTypeObjectName))
        elif re.match("^PocketSolvent$", GroupID, re.I):
            OutFH.write("""\n""\n"Setting up views for solvent in pockect around ligand %s in chain %s..."\n""\n""" % (LigandID, ChainID))
            PML = PyMOLUtil.SetupPMLForLigandPocketSolventView(GroupTypeObjectName, ComplexName, LigandName, OptionsInfo["PocketDistanceCutoff"], Enable = True)
            OutFH.write("%s\n" % PML)
        elif re.match("^PocketInorganic$", GroupID, re.I):
            OutFH.write("""\n""\n"Setting up views for inorganic in pockect around ligand %s in chain %s..."\n""\n""" % (LigandID, ChainID))
            PML = PyMOLUtil.SetupPMLForLigandPocketInorganicView(GroupTypeObjectName, ComplexName, LigandName, OptionsInfo["PocketDistanceCutoff"], Enable = True)
            OutFH.write("%s\n" % PML)
        
        # Set up polar contacts...
        if re.match("^(Pocket|PocketSolvent|PocketInorganic)$", GroupID, re.I):
            PolarContactsID = "%sPolarContacts" % (GroupID)
            PolarContactsName = PyMOLObjectNames["Ligands"][ChainID][LigandID][PolarContactsID]
            
            PolarContactsColor = OptionsInfo["PocketContactsLigandColor"]
            if re.match("^PocketSolvent$", GroupID, re.I):
                PolarContactsColor = OptionsInfo["PocketContactsSolventColor"]
            elif re.match("^PocketInorganic$", GroupID, re.I):
                PolarContactsColor = OptionsInfo["PocketContactsInorganicColor"]

            PML = PyMOLUtil.SetupPMLForPolarContactsView(PolarContactsName, LigandName, GroupTypeObjectName, Enable = False, Color = PolarContactsColor, Cutoff = OptionsInfo["PocketContactsCutoff"])
            OutFH.write("\n%s\n" % PML)
            
            OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (PolarContactsColor, PolarContactsName))
            
        if re.match("^Pocket$", GroupID, re.I):
            # Setup hydrophobic contacts...
            HydrophobicContactsID = "%sHydrophobicContacts" % (GroupID)
            HydrophobicContactsName = PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicContactsID]
            HydrophobicContactsColor = OptionsInfo["PocketContactsLigandHydrophobicColor"]
            
            PML = PyMOLUtil.SetupPMLForHydrophobicContactsView(HydrophobicContactsName, LigandName, GroupTypeObjectName, Enable = False, Color = HydrophobicContactsColor, Cutoff = OptionsInfo["PocketContactsCutoff"])
            OutFH.write("\n%s\n" % PML)
            OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (HydrophobicContactsColor, HydrophobicContactsName))

            # Setup pocket residues...
            WritePocketResidueTypesView(OutFH,  FileIndex, PyMOLObjectNames, ChainID, LigandID, GroupTypeObjectID)
            WritePocketSurfacesTypesView(OutFH,  FileIndex, PyMOLObjectNames, ChainID, LigandID, GroupTypeObjectID)
        
        if re.match("^Ligand$", GroupID, re.I):
            # Setup ball and stick view...
            BallAndStickNameID = "%sBallAndStick" % (GroupID)
            BallAndStickName = PyMOLObjectNames["Ligands"][ChainID][LigandID][BallAndStickNameID]
            PML = PyMOLUtil.SetupPMLForBallAndStickView(BallAndStickName, GroupTypeObjectName, Enable = False)
            OutFH.write("\n%s\n" % PML)
            
        # Setup group...
        GroupNameID = "%sGroup" % (GroupID)
        GroupMembersID = "%sGroupMembers" % (GroupID)
        GroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupNameID]
        GroupMembers = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID]

        Action = "close"
        Enable = False
        if  re.match("^(Ligand|Pocket)$", GroupID, re.I):
            Action = "open"
            Enable = True
        GenerateAndWritePMLForGroup(OutFH, GroupName, GroupMembers, Enable, Action)

def WritePocketSurfacesTypesView(OutFH,  FileIndex, PyMOLObjectNames, ChainID, LigandID, PocketObjectID):
    """Write out PML for viewing surfaces for a ligand pocket."""

    if not GetPocketContainsSurfaceStatus(FileIndex, ChainID, LigandID):
        return
    
    PocketObjectName = PyMOLObjectNames["Ligands"][ChainID][LigandID][PocketObjectID]
    
    SurfacesGroupID = "%sSurfacesGroup" % (PocketObjectID)
    SurfacesGroupMembersID = "%sSurfacesGroupMembers" % (PocketObjectID)

    # Cavity modes: 1 or 2.  1: Cavity surfaces; 2: Culled cavity surfaces...
    CavityMode = 2
    
    # Setup surfaces subgroup and its members...
    for SubGroupType in ["Surface", "Cavity"]:
        ProcessingCavity = True if re.match("^Cavity$", SubGroupType, re.I) else False
        
        SubGroupID = re.sub("_", "", SubGroupType)
        SurfacesSubGroupID = "%s%sGroup" % (SurfacesGroupID, SubGroupID)
        SurfacesSubGroupMembersID = "%sMembers" % (SurfacesSubGroupID)

        # Turn off lines display for cavity surfaces...
        DisplayStyle = None if ProcessingCavity else "lines"

        # Setup a generic color surface...
        SurfaceID = "%sSurface" % (SurfacesSubGroupID)
        SurfaceName = PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfaceID]
        PML = PyMOLUtil.SetupPMLForSurfaceView(SurfaceName, PocketObjectName, Enable = False, Color = OptionsInfo["SurfaceColor"], DisplayAs = DisplayStyle)
        OutFH.write("\n%s\n" % PML)
        
        if ProcessingCavity:
            OutFH.write("""cmd.set("surface_cavity_mode", %d, "%s")\n""" % (CavityMode, SurfaceName))
        
        OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (OptionsInfo["PocketLabelColor"], SurfaceName))
        
        if GetPocketSurfaceChainStatus(FileIndex, ChainID, LigandID):
            # Setup a surface colored by hydrphobicity...
            HydrophobicSurfaceID = "%sHydrophobicSurface" % (SurfacesSubGroupID)
            HydrophobicSurfaceName = PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicSurfaceID]
            PML = PyMOLUtil.SetupPMLForHydrophobicSurfaceView(HydrophobicSurfaceName, PocketObjectName, ColorPalette = OptionsInfo["SurfaceColorPalette"], Enable = False, DisplayAs = DisplayStyle)
            OutFH.write("\n%s\n" % PML)
            
            if ProcessingCavity:
                OutFH.write("""cmd.set("surface_cavity_mode", %d, "%s")\n""" % (CavityMode, HydrophobicSurfaceName))
            
            OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (OptionsInfo["PocketLabelColor"], HydrophobicSurfaceName))
            
            # Setup a surface colored by hydrphobicity and charge...
            HydrophobicChargeSurfaceID = "%sHydrophobicChargeSurface" % (SurfacesSubGroupID)
            HydrophobicChargeSurfaceName = PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicChargeSurfaceID]
            PML = PyMOLUtil.SetupPMLForHydrophobicAndChargeSurfaceView(HydrophobicChargeSurfaceName, PocketObjectName,  OptionsInfo["AtomTypesColorNames"]["HydrophobicAtomsColor"], OptionsInfo["AtomTypesColorNames"]["NegativelyChargedAtomsColor"], OptionsInfo["AtomTypesColorNames"]["PositivelyChargedAtomsColor"], OptionsInfo["AtomTypesColorNames"]["OtherAtomsColor"], Enable = False, DisplayAs = DisplayStyle)
            OutFH.write("\n%s\n" % PML)
            
            if ProcessingCavity:
                OutFH.write("""cmd.set("surface_cavity_mode", %d, "%s")\n""" % (CavityMode, HydrophobicChargeSurfaceName))
            
            OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (OptionsInfo["PocketLabelColor"], HydrophobicChargeSurfaceName))
            
            if GetPocketSurfaceChainElectrostaticsStatus(FileIndex, ChainID, LigandID):
                # Set up a electrostatics surface...
                ElectrostaticsGroupID = "%sElectrostaticsGroup" % (SurfacesSubGroupID)
                ElectrostaticsGroupMembersID = "%sElectrostaticsGroupMembers" % (SurfacesSubGroupID)
                ElectrostaticsGroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][ElectrostaticsGroupID]
                ElectrostaticsGroupMembers = PyMOLObjectNames["Ligands"][ChainID][LigandID][ElectrostaticsGroupMembersID]
                WriteSurfaceElectrostaticsView(SubGroupType, OutFH, PocketObjectName, ElectrostaticsGroupName, ElectrostaticsGroupMembers, DisplayAs = DisplayStyle, SurfaceCavityMode = CavityMode)
                
            # Setup surfaces sub group...
            GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesSubGroupID], PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesSubGroupMembersID], True, "open")
        
    # Setup surface group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesGroupID], PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesGroupMembersID], True, "open")
    
    
def WritePocketResidueTypesView(OutFH,  FileIndex, PyMOLObjectNames, ChainID, LigandID, PocketObjectID):
    """Write out PML for viewing residue types for a ligand pocket."""
    
    if not GetPocketResidueTypesStatus(FileIndex, ChainID, LigandID):
        return

    PocketObjectName = PyMOLObjectNames["Ligands"][ChainID][LigandID][PocketObjectID]
    
    ResiduesGroupID = "%sResiduesGroup" % (PocketObjectID)
    ResiduesGroupMembersID = "%sMembers" % (ResiduesGroupID)
    
    # Setup residue types objects...
    for SubGroupType in ["Aromatic", "Hydrophobic", "Polar", "Positively_Charged", "Negatively_Charged", "Other"]:
        SubGroupID = re.sub("_", "", SubGroupType)
        
        ResiduesSubGroupID = "%s%sGroup" % (ResiduesGroupID, SubGroupID)
        ResiduesSubMembersGroupID = "%sMembers" % (ResiduesSubGroupID)
    
        SubGroupMemberID = "%sResidues" % (ResiduesSubGroupID)
        ResiduesObjectName = PyMOLObjectNames["Ligands"][ChainID][LigandID][SubGroupMemberID]

        SubGroupMemberID = "%sSurface" % (ResiduesSubGroupID)
        ResiduesSurfaceObjectName = PyMOLObjectNames["Ligands"][ChainID][LigandID][SubGroupMemberID]

        ResiduesColor = OptionsInfo["ResidueTypesParams"][SubGroupType]["Color"] 
        ResiduesNames = OptionsInfo["ResidueTypesParams"][SubGroupType]["Residues"]

        NegateResidueNames = True if re.match("^Other$", SubGroupType, re.I) else False
        WriteResidueTypesResiduesAndSurfaceView(OutFH, PocketObjectName, ResiduesObjectName, ResiduesSurfaceObjectName, ResiduesColor, ResiduesNames, NegateResidueNames)
        
        # Setup residue type sub groups...
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesSubGroupID], PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesSubMembersGroupID], True, "close")
    
    # Setup residue types group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesGroupID], PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesGroupMembersID], False, "close")

def WriteChainAloneChainSelectionsView(OutFH,  FileIndex, PyMOLObjectNames, ChainID):
    """Write out PML for viewing selections for a chain. """

    if not GetChainAloneContainsSelectionsStatus(FileIndex, ChainID):
        return
    
    ChainName = PyMOLObjectNames["Chains"][ChainID]["ChainAlone"]
    SelectionsGroupIDPrefix = "ChainAloneSelections"
    
    for Index in range(0, len(OptionsInfo["ChainSelectionsInfo"]["Names"])):
        SelectionName = OptionsInfo["ChainSelectionsInfo"]["Names"][Index]
        Selection = OptionsInfo["ChainSelectionsInfo"]["Selections"][Index]
        
        SelectionNameGroupID = SelectionName
        
        # Setup a selection object...
        SelectionObjectID = "%s%sSelection" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        SelectionObjectName = PyMOLObjectNames["Chains"][ChainID][SelectionObjectID]
        SelectionCmd = "(%s and (%s))" % (ChainName, Selection)
        PML = PyMOLUtil.SetupPMLForSelectionDisplayView(SelectionObjectName, SelectionCmd, OptionsInfo["SelectionsChainStyle"], Enable = True)
        OutFH.write("\n%s\n" % PML)
        
        if GetChainAloneContainsChainSelectionSurfacesStatus(FileIndex, ChainID):
            # Display style for selection objects in surfaces...
            DisplayStyle = "lines"
            
            # Setup a generic color surface...
            SurfaceID = "%s%s%sSurface" % (SelectionsGroupIDPrefix, SelectionNameGroupID, "Surface")
            SurfaceName = PyMOLObjectNames["Chains"][ChainID][SurfaceID]
            PML = PyMOLUtil.SetupPMLForSurfaceView(SurfaceName, SelectionObjectName, Enable = False, Color = OptionsInfo["SurfaceColor"], DisplayAs = DisplayStyle)
            OutFH.write("\n%s\n" % PML)
            
            if GetChainAloneSurfaceChainSelectionStatus(FileIndex, ChainID):
                # Setup a surface colored by hydrphobicity...
                HydrophobicSurfaceID = "%s%sSurfaceHydrophobicity" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
                HydrophobicSurfaceName = PyMOLObjectNames["Chains"][ChainID][HydrophobicSurfaceID]
                PML = PyMOLUtil.SetupPMLForHydrophobicSurfaceView(HydrophobicSurfaceName, SelectionObjectName, ColorPalette = OptionsInfo["SurfaceColorPalette"], Enable = False, DisplayAs = DisplayStyle)
                OutFH.write("\n%s\n" % PML)
                
                # Setup a surface colored by hydrphobicity and charge...
                HydrophobicChargeSurfaceID = "%s%sSurfaceHydrophobicityCharge" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
                HydrophobicChargeSurfaceName = PyMOLObjectNames["Chains"][ChainID][HydrophobicChargeSurfaceID]
                PML = PyMOLUtil.SetupPMLForHydrophobicAndChargeSurfaceView(HydrophobicChargeSurfaceName, SelectionObjectName,  OptionsInfo["AtomTypesColorNames"]["HydrophobicAtomsColor"], OptionsInfo["AtomTypesColorNames"]["NegativelyChargedAtomsColor"], OptionsInfo["AtomTypesColorNames"]["PositivelyChargedAtomsColor"], OptionsInfo["AtomTypesColorNames"]["OtherAtomsColor"], Enable = False, DisplayAs = DisplayStyle)
                OutFH.write("\n%s\n" % PML)

            # Setup group for surfaces...
            SurfaceGroupID = "%s%sSurfaceGroup" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
            SurfaceGroupMembersID = "%s%sSurfaceGroupMembers" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
            GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID][SurfaceGroupID], PyMOLObjectNames["Chains"][ChainID][SurfaceGroupMembersID], True, "close")

        # Setup groups for named selections...
        SelectionsNameGroupID = "%s%sGroup" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        SelectionsNameGroupMembersID = "%s%sGroupMembers" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID][SelectionsNameGroupID], PyMOLObjectNames["Chains"][ChainID][SelectionsNameGroupMembersID], True, "open")

    # Setup a group for selections...
    SelectionsGroupID = "%sGroup" % (SelectionsGroupIDPrefix)
    SelectionsGroupMembersID = "%sGroupMembers" % (SelectionsGroupIDPrefix)
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID][SelectionsGroupID], PyMOLObjectNames["Chains"][ChainID][SelectionsGroupMembersID], False, "close")
    
def WriteChainAloneResidueTypesView(OutFH,  FileIndex, PyMOLObjectNames, ChainID):
    """Write out PML for viewing residue types for a chain. """

    if not GetChainAloneResidueTypesStatus(FileIndex, ChainID):
        return
    
    ChainName = PyMOLObjectNames["Chains"][ChainID]["ChainAlone"]
    
    # Setup residue types objects...
    ResiduesGroupIDPrefix = "ChainAloneResidues"
    for SubGroupType in ["Aromatic", "Hydrophobic", "Polar", "Positively_Charged", "Negatively_Charged", "Other"]:
        SubGroupID = re.sub("_", "", SubGroupType)

        ResiduesObjectID = "%s%sResidues" % (ResiduesGroupIDPrefix, SubGroupID)
        ResiduesObjectName = PyMOLObjectNames["Chains"][ChainID][ResiduesObjectID]

        ResiduesSurfaceObjectID = "%s%sSurface" % (ResiduesGroupIDPrefix, SubGroupID)
        ResiduesSurfaceObjectName = PyMOLObjectNames["Chains"][ChainID][ResiduesSurfaceObjectID]

        ResiduesColor = OptionsInfo["ResidueTypesParams"][SubGroupType]["Color"] 
        ResiduesNames = OptionsInfo["ResidueTypesParams"][SubGroupType]["Residues"]

        NegateResidueNames = True if re.match("^Other$", SubGroupType, re.I) else False
        WriteResidueTypesResiduesAndSurfaceView(OutFH, ChainName, ResiduesObjectName, ResiduesSurfaceObjectName, ResiduesColor, ResiduesNames, NegateResidueNames)

        # Setup sub groups for residue types..
        ResiduesSubGroupID = "%s%sGroup" % (ResiduesGroupIDPrefix, SubGroupID)
        ResiduesSubGroupMembersID = "%s%sGroupMembers" % (ResiduesGroupIDPrefix, SubGroupID)

        # Setup residue type sub groups...
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID][ResiduesSubGroupID], PyMOLObjectNames["Chains"][ChainID][ResiduesSubGroupMembersID], True, "close")
        
    # Setup residue types group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainAloneResiduesGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainAloneResiduesGroupMembers"], False, "close")

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
    
def WriteSurfaceElectrostaticsView(Mode, OutFH, SelectionObjectName, ElectrostaticsGroupName, ElectrostaticsGroupMembers, DisplayAs = None, SurfaceCavityMode = 2):
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
    
    if DisplayAs is not None:
        PMLCmds.append("""cmd.show("%s", "(%s)")""" % (DisplayAs, ContactPotentialName))

    if re.match("^Cavity$", Mode, re.I):
        if SurfaceCavityMode is not None:
            PMLCmds.append("""cmd.set("surface_cavity_mode", %d, "%s")\n""" % (SurfaceCavityMode, ContactPotentialName))
        
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

        # Delete chain selection objects...
        DeleteEmptyChainSelectionsObjects(OutFH, FileIndex, PyMOLObjectNames, ChainID)
        
        # Delete residue type objects...
        DeleteEmptyChainResidueTypesObjects(OutFH, FileIndex, PyMOLObjectNames, ChainID)
        
        for LigandID in SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID]:
            # Delete ligand level objects...
            for GroupID in ["Pocket", "PocketSolvent", "PocketInorganic"]:
                GroupNameID = "%sGroup" % (GroupID)
                GroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupNameID]

                GroupTypeObjectID = "%s" % (GroupID)
                GroupTypeObjectName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupTypeObjectID]
                
                WritePMLToCheckAndDeleteEmptyObjects(OutFH, GroupTypeObjectName, GroupName)
                
                if re.match("^Pocket$", GroupID, re.I):
                    DeleteEmptyPocketResidueTypesObjects(OutFH, FileIndex, PyMOLObjectNames, ChainID, LigandID, GroupTypeObjectID)

def DeleteEmptyChainSelectionsObjects(OutFH, FileIndex, PyMOLObjectNames, ChainID):
    """Delete empty chain selection objects. """

    if not GetChainAloneContainsSelectionsStatus(FileIndex, ChainID):
        return
    
    SelectionsGroupIDPrefix = "ChainAloneSelections"

    for SelectionName in OptionsInfo["ChainSelectionsInfo"]["Names"]:
        SelectionNameGroupID = SelectionName
        
        # Delete surface objects and surface group...
        SurfaceGroupID = "%s%sSurfaceGroup" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        SurfaceGroupMembersID = "%s%sSurfaceGroupMembers" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        WritePMLToCheckAndDeleteEmptyObjects(OutFH, ",".join(PyMOLObjectNames["Chains"][ChainID][SurfaceGroupMembersID]), PyMOLObjectNames["Chains"][ChainID][SurfaceGroupID])
        
        # Delete Selection object and selection name group...
        SelectionObjectID = "%s%sSelection" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        SelectionsGroupID = "%s%sGroup" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        WritePMLToCheckAndDeleteEmptyObjects(OutFH, PyMOLObjectNames["Chains"][ChainID][SelectionObjectID], PyMOLObjectNames["Chains"][ChainID][SelectionsGroupID])

def DeleteEmptyChainResidueTypesObjects(OutFH, FileIndex, PyMOLObjectNames, ChainID):
    """Delete empty chain residue objects. """
    
    if not GetChainAloneResidueTypesStatus(FileIndex, ChainID):
        return
    
    ResiduesGroupIDPrefix = "ChainAloneResidues"
    for GroupType in ["Aromatic", "Hydrophobic", "Polar", "Positively_Charged", "Negatively_Charged", "Other"]:
        GroupID = re.sub("_", "", GroupType)
        
        ResiduesGroupID = "%s%sGroup" % (ResiduesGroupIDPrefix, GroupID)
        GroupName = PyMOLObjectNames["Chains"][ChainID][ResiduesGroupID]
        
        GroupObjectNamesList = []
        
        ResiduesObjectID = "%s%sResidues" % (ResiduesGroupIDPrefix, GroupID)
        ResiduesObjectName = PyMOLObjectNames["Chains"][ChainID][ResiduesObjectID]
        GroupObjectNamesList.append(ResiduesObjectName)
        
        ResiduesSurfaceObjectID = "%s%sSurface" % (ResiduesGroupIDPrefix, GroupID)
        ResiduesSurfaceObjectName = PyMOLObjectNames["Chains"][ChainID][ResiduesSurfaceObjectID]
        GroupObjectNamesList.append(ResiduesSurfaceObjectName)
        
        GroupObjectNames = ",".join(GroupObjectNamesList)
        WritePMLToCheckAndDeleteEmptyObjects(OutFH, GroupObjectNames, GroupName)

def DeleteEmptyPocketResidueTypesObjects(OutFH, FileIndex, PyMOLObjectNames, ChainID, LigandID, PocketObjectID):
    """Delete empty chain residue objects. """
    
    if not GetPocketResidueTypesStatus(FileIndex, ChainID, LigandID):
        return
    
    ResiduesGroupID = "%sResiduesGroup" % (PocketObjectID)
    ResiduesGroupMembersID = "%sMembers" % (ResiduesGroupID)
    
    for SubGroupType in ["Aromatic", "Hydrophobic", "Polar", "Positively_Charged", "Negatively_Charged", "Other"]:
        SubGroupID = re.sub("_", "", SubGroupType)
        
        ResiduesSubGroupID = "%s%sGroup" % (ResiduesGroupID, SubGroupID)
        ResiduesSubMembersGroupID = "%sMembers" % (ResiduesSubGroupID)

        SubGroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesSubGroupID]
        SubGroupObjectNames = ",".join(PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesSubMembersGroupID])
        
        WritePMLToCheckAndDeleteEmptyObjects(OutFH, SubGroupObjectNames, SubGroupName)

def WritePMLToCheckAndDeleteEmptyObjects(OutFH, ObjectName, ParentObjectName = None):
    """Write PML to check and delete empty PyMOL objects. """
    
    if ParentObjectName is None:
        PML = """CheckAndDeleteEmptyObjects("%s")""" % (ObjectName)
    else:
        PML = """CheckAndDeleteEmptyObjects("%s", "%s")""" % (ObjectName, ParentObjectName)
    
    OutFH.write("%s\n" % PML)
    
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

    if GetChainAloneContainsSelectionsStatus(FileIndex, ChainID):
        # Setup selections group and its subgroups..
        SelectionsGroupName = "%s.Selections" % (ChainAloneGroupName)
        
        SelectionsGroupIDPrefix = "ChainAloneSelections"
        SelectionsGroupID = "%sGroup" % SelectionsGroupIDPrefix
        
        # Add selections group to chain alone group...
        PyMOLObjectNames["Chains"][ChainID][SelectionsGroupID] = SelectionsGroupName
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"].append(SelectionsGroupName)
        
        # Initialize selections group members...
        SelectionsGroupMembersID = "%sGroupMembers" % SelectionsGroupIDPrefix
        PyMOLObjectNames["Chains"][ChainID][SelectionsGroupMembersID] = []

        # Setup selections name sub group and its members...
        for SelectionName in OptionsInfo["ChainSelectionsInfo"]["Names"]:
            SelectionNameGroupID = SelectionName
            
            SelectionsNameGroupName = "%s.%s" % (SelectionsGroupName, SelectionName)
            SelectionsNameGroupID = "%s%sGroup" % (SelectionsGroupIDPrefix, SelectionNameGroupID)

            # Add selections name sub group to selections group...
            PyMOLObjectNames["Chains"][ChainID][SelectionsNameGroupID] = SelectionsNameGroupName
            PyMOLObjectNames["Chains"][ChainID][SelectionsGroupMembersID].append(SelectionsNameGroupName)
            
            # Initialize selections names sub group members...
            SelectionsNameGroupMembersID = "%s%sGroupMembers" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
            PyMOLObjectNames["Chains"][ChainID][SelectionsNameGroupMembersID] = []

            # Add selection member to selections name group...
            SubGroupMemberName = "%s.Selection" % (SelectionsNameGroupName)
            SubGroupMemberID = "%s%sSelection" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
            
            PyMOLObjectNames["Chains"][ChainID][SubGroupMemberID] = SubGroupMemberName
            PyMOLObjectNames["Chains"][ChainID][SelectionsNameGroupMembersID].append(SubGroupMemberName)
            
            if GetChainAloneContainsChainSelectionSurfacesStatus(FileIndex, ChainID):
                # Setup a surface sub group and add it to selections name group...
                
                SelectionsNameSurfaceGroupName = "%s.Surface" % (SelectionsNameGroupName)
                SelectionsNameSurfaceGroupID = "%s%sSurfaceGroup" % (SelectionsGroupIDPrefix, SelectionNameGroupID)

                # Add selection surface group to selectio name group...
                PyMOLObjectNames["Chains"][ChainID][SelectionsNameSurfaceGroupID] = SelectionsNameSurfaceGroupName
                PyMOLObjectNames["Chains"][ChainID][SelectionsNameGroupMembersID].append(SelectionsNameSurfaceGroupName)

                # Initialize surface names sub group members...
                SelectionsNameSurfaceGroupMembersID = "%s%sSurfaceGroupMembers" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
                PyMOLObjectNames["Chains"][ChainID][SelectionsNameSurfaceGroupMembersID] = []
                
                # Setup a generic color surface...
                SubGroupMemberName = "%s.Surface" % (SelectionsNameSurfaceGroupName)
                SubGroupMemberID = "%s%s%sSurface" % (SelectionsGroupIDPrefix, SelectionNameGroupID, "Surface")
                PyMOLObjectNames["Chains"][ChainID][SubGroupMemberID] = SubGroupMemberName
                PyMOLObjectNames["Chains"][ChainID][SelectionsNameSurfaceGroupMembersID].append(SubGroupMemberName)
                
                if GetChainAloneSurfaceChainSelectionStatus(FileIndex, ChainID):
                    # Setup surfaces...
                    for MemberType in ["Hydrophobicity", "Hydrophobicity_Charge"]:
                        MemberID = re.sub("_", "", MemberType)
                        
                        SubGroupMemberName = "%s.%s" % (SelectionsNameSurfaceGroupName, MemberType)
                        SubGroupMemberID = "%s%s%s%s" % (SelectionsGroupIDPrefix, SelectionNameGroupID, "Surface", MemberID)
                        
                        PyMOLObjectNames["Chains"][ChainID][SubGroupMemberID] = SubGroupMemberName
                        PyMOLObjectNames["Chains"][ChainID][SelectionsNameSurfaceGroupMembersID].append(SubGroupMemberName)

    if GetChainAloneResidueTypesStatus(FileIndex, ChainID):
        # Setup residue type group and its subgroups...
        ResiduesGroupName = "%s.Residues" % (ChainAloneGroupName)
        
        ResiduesGroupIDPrefix = "ChainAloneResidues"
        ResiduesGroupID = "%sGroup" % ResiduesGroupIDPrefix

        # Add residue group to chain alone group...
        PyMOLObjectNames["Chains"][ChainID][ResiduesGroupID] = ResiduesGroupName
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"].append(ResiduesGroupName)
        
        # Initialize residue group members...
        ResiduesGroupMembersID = "%sGroupMembers" % ResiduesGroupIDPrefix
        PyMOLObjectNames["Chains"][ChainID][ResiduesGroupMembersID] = []

        # Setup residues sub groups and its members...
        for SubGroupType in ["Aromatic", "Hydrophobic", "Polar", "Positively_Charged", "Negatively_Charged", "Other"]:
            SubGroupID = re.sub("_", "", SubGroupType)

            ResiduesSubGroupName = "%s.%s" % (ResiduesGroupName, SubGroupType)
            ResiduesSubGroupID = "%s%sGroup" % (ResiduesGroupIDPrefix, SubGroupID)

            # Add sub group to residues group...
            PyMOLObjectNames["Chains"][ChainID][ResiduesSubGroupID] = ResiduesSubGroupName
            PyMOLObjectNames["Chains"][ChainID][ResiduesGroupMembersID].append(ResiduesSubGroupName)

            # Initialize sub group members...
            ResiduesSubGroupMembersID = "%s%sGroupMembers" % (ResiduesGroupIDPrefix, SubGroupID)
            PyMOLObjectNames["Chains"][ChainID][ResiduesSubGroupMembersID] = []
            
            # Add sub group members to subgroup...
            for MemberType in ["Residues", "Surface"]:
                MemberID = re.sub("_", "", MemberType)

                SubGroupMemberName = "%s.%s" % (ResiduesSubGroupName, MemberType)
                SubGroupMemberID = "%s%s%s" % (ResiduesGroupIDPrefix, SubGroupID, MemberID)
                
                PyMOLObjectNames["Chains"][ChainID][SubGroupMemberID] = SubGroupMemberName
                PyMOLObjectNames["Chains"][ChainID][ResiduesSubGroupMembersID].append(SubGroupMemberName)

    if GetChainAloneContainsSurfacesStatus(FileIndex, ChainID):
        # Setup a surface group and add it to chain alone group...
        SurfaceGroupName = "%s.Surface" % (ChainAloneGroupName)
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroup"] = SurfaceGroupName
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"].append(SurfaceGroupName)
        
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneSurfaceGroupMembers"] = []

        # Setup a generic color surface...
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

    # Set up groups and objects for a specific ligand group...
    for GroupType in ["Ligand", "Pocket", "Pocket_Solvent", "Pocket_Inorganic"]:
        GroupID = re.sub("_", "", GroupType)
        GroupName = "%s.%s" % (ChainLigandGroupName, GroupType)
                
        GroupNameID = "%sGroup" % (GroupID)
        GroupMembersID = "%sGroupMembers" % (GroupID)
        
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupNameID] = GroupName
        PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroupMembers"].append(GroupName)
        
        GroupTypeObjectName = "%s.%s" % (GroupName, GroupType)
        GroupTypeObjectID = "%s" % (GroupID)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupTypeObjectID] = GroupTypeObjectName
        
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID] = []
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(GroupTypeObjectName)
        
        if re.match("^Ligand$", GroupType, re.I):
            # Only need to add ball and stick...
            BallAndStickName = "%s.BallAndStick" % (GroupName)
            BallAndStickID = "%sBallAndStick" % (GroupID)
            PyMOLObjectNames["Ligands"][ChainID][LigandID][BallAndStickID] = BallAndStickName
            PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(BallAndStickName)
        
        if re.match("^(Pocket|Pocket_Solvent|Pocket_Inorganic)$", GroupType, re.I):
            PolarContactsName = "%s.Polar_Contacts" % (GroupName)
            PolarContactsID = "%sPolarContacts" % (GroupID)
            PyMOLObjectNames["Ligands"][ChainID][LigandID][PolarContactsID] = PolarContactsName
            PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(PolarContactsName)
                
        if re.match("^Pocket$", GroupType, re.I):
            HydrophobicContactsName = "%s.Hydrophobic_Contacts" % (GroupName)
            HydrophobicContactsID = "%sHydrophobicContacts" % (GroupID)
            PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicContactsID] = HydrophobicContactsName
            PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(HydrophobicContactsName)
            
            if GetPocketResidueTypesStatus(FileIndex, ChainID, LigandID):
                # Setup residue type group and add to pocket group...
                ResiduesGroupName = "%s.Residues" % (GroupName)
                ResiduesGroupID = "%sResiduesGroup" % (GroupID)
                ResiduesGroupMembersID = "%sMembers" % (ResiduesGroupID)

                PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesGroupID] = ResiduesGroupName
                PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(ResiduesGroupName)
                PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesGroupMembersID] = []
                
                # Setup residue group subgroup and its members...
                for SubGroupType in ["Aromatic", "Hydrophobic", "Polar", "Positively_Charged", "Negatively_Charged", "Other"]:
                    SubGroupID = re.sub("_", "", SubGroupType)
                    ResiduesSubGroupName = "%s.%s" % (ResiduesGroupName, SubGroupType)
                    ResiduesSubGroupID = "%s%sGroup" % (ResiduesGroupID, SubGroupID)
                    ResiduesSubMembersGroupID = "%sMembers" % (ResiduesSubGroupID)
                    
                    # Add sub group to residues group...
                    PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesSubGroupID] = ResiduesSubGroupName
                    PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesGroupMembersID].append(ResiduesSubGroupName)
                    PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesSubMembersGroupID] = []
                    
                    # Add sub group members to subgroup...
                    for MemberType in ["Residues", "Surface"]:
                        MemberID = re.sub("_", "", MemberType)
                        SubGroupMemberName = "%s.%s" % (ResiduesSubGroupName, MemberType)
                        SubGroupMemberID = "%s%s" % (ResiduesSubGroupID, MemberID)
                        
                        PyMOLObjectNames["Ligands"][ChainID][LigandID][SubGroupMemberID] = SubGroupMemberName
                        PyMOLObjectNames["Ligands"][ChainID][LigandID][ResiduesSubMembersGroupID].append(SubGroupMemberName)

            if GetPocketContainsSurfaceStatus(FileIndex, ChainID, LigandID):
                # Setup a surfaces group and add it to pocket group...
                SurfacesGroupName = "%s.Surfaces" % (GroupName)
                SurfacesGroupID = "%sSurfacesGroup" % (GroupID)
                SurfacesGroupMembersID = "%sMembers" % (SurfacesGroupID)

                PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesGroupID] = SurfacesGroupName
                PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(SurfacesGroupName)
                PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesGroupMembersID] = []

                # Setup surfaces subgroup and its members...
                for SubGroupType in ["Surface", "Cavity"]:
                    SubGroupID = re.sub("_", "", SubGroupType)
                    SurfacesSubGroupName = "%s.%s" % (SurfacesGroupName, SubGroupType)
                    SurfacesSubGroupID = "%s%sGroup" % (SurfacesGroupID, SubGroupID)
                    SurfacesSubGroupMembersID = "%sMembers" % (SurfacesSubGroupID)
                    
                    # Add sub group to surfaces group...
                    PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesSubGroupID] = SurfacesSubGroupName
                    PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesGroupMembersID].append(SurfacesSubGroupName)
                    PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesSubGroupMembersID] = []
                    
                    # Setup a generic color surface...
                    SurfaceName = "%s.Surface" % (SurfacesSubGroupName)
                    SurfaceID = "%sSurface" % (SurfacesSubGroupID)
                    PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfaceID] = SurfaceName
                    PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesSubGroupMembersID].append(SurfaceName)
                    
                    if GetPocketSurfaceChainStatus(FileIndex, ChainID, LigandID):
                        # Surface colored by hydrophobicity...
                        HydrophobicSurfaceName = "%s.Hydrophobicity" % (SurfacesSubGroupName)
                        HydrophobicSurfaceID = "%sHydrophobicSurface" % (SurfacesSubGroupID)
                        PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicSurfaceID] = HydrophobicSurfaceName
                        PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesSubGroupMembersID].append(HydrophobicSurfaceName)
                        
                        # Surface colored by hydrophobicity and charge...
                        HydrophobicChargeSurfaceName = "%s.Hydrophobicity_Charge" % (SurfacesSubGroupName)
                        HydrophobicChargeSurfaceID = "%sHydrophobicChargeSurface" % (SurfacesSubGroupID)
                        PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicChargeSurfaceID] = HydrophobicChargeSurfaceName
                        PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesSubGroupMembersID].append(HydrophobicChargeSurfaceName)
                
                    if GetPocketSurfaceChainElectrostaticsStatus(FileIndex, ChainID, LigandID):
                        ElectrostaticsGroupName = "%s.Vacuum_Electrostatics" % (SurfacesSubGroupName)
                        ElectrostaticsGroupID = "%sElectrostaticsGroup" % (SurfacesSubGroupID)
                        ElectrostaticsGroupMembersID = "%sElectrostaticsGroupMembers" % (SurfacesSubGroupID)
                        
                        PyMOLObjectNames["Ligands"][ChainID][LigandID][ElectrostaticsGroupID] = ElectrostaticsGroupName
                        PyMOLObjectNames["Ligands"][ChainID][LigandID][SurfacesSubGroupMembersID].append(ElectrostaticsGroupName)
                        
                        # Setup electrostatics group members without the volume object...
                        PyMOLObjectNames["Ligands"][ChainID][LigandID][ElectrostaticsGroupMembersID] = []
                        
                        for MemberType in ["Pocket", "Contact_Potential", "Map", "Legend"]:
                            MemberID = re.sub("_", "", MemberType)
                            
                            Name = "%s.%s" % (ElectrostaticsGroupName, MemberType)
                            NameID = "%s%s" % (ElectrostaticsGroupID, MemberID)
                            
                            PyMOLObjectNames["Ligands"][ChainID][LigandID][NameID] = Name
                            PyMOLObjectNames["Ligands"][ChainID][LigandID][ElectrostaticsGroupMembersID].append(Name)
    
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

def RetrieveRefFileInfo():
    """Retrieve information for ref file."""

    RefFileInfo = {}
    if not OptionsInfo["Align"]:
        OptionsInfo["RefFileInfo"] = RefFileInfo
        return

    RefFile = OptionsInfo["RefFileName"]
    
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(RefFile)
    RefFileRoot = FileName
    
    if re.match("^FirstInputFile$", OptionsInfo["AlignRefFile"], re.I):
        ChainsAndLigandInfo = OptionsInfo["InfilesInfo"]["ChainsAndLigandsInfo"][0]
    else:
        MiscUtil.PrintInfo("\nRetrieving chain and ligand information for alignment reference file %s..." % RefFile)
        ChainsAndLigandInfo = PyMOLUtil.GetChainsAndLigandsInfo(RefFile, RefFileRoot)

    RefFileInfo["RefFileName"] = RefFile
    RefFileInfo["RefFileRoot"] = RefFileRoot
    RefFileInfo["PyMOLObjectName"] = "AlignRef_%s" % RefFileRoot
    RefFileInfo["ChainsAndLigandsInfo"] = ChainsAndLigandInfo
    
    OptionsInfo["RefFileInfo"] = RefFileInfo

def ProcessChainAndLigandIDs():
    """Process specified chain and ligand IDs for infiles."""
    
    OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"] = []
    
    for FileIndex in range(0, len(OptionsInfo["InfilesInfo"]["InfilesNames"])):
        MiscUtil.PrintInfo("\nProcessing specified chain and ligand IDs for input file %s..." % OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex])
        
        ChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["ChainsAndLigandsInfo"][FileIndex]
        SpecifiedChainsAndLigandsInfo = PyMOLUtil.ProcessChainsAndLigandsOptionsInfo(ChainsAndLigandsInfo, "-c, --chainIDs", OptionsInfo["ChainIDs"], "-l, --ligandIDs", OptionsInfo["LigandIDs"])
        ProcessResidueTypesAndSurfaceAndChainSelectionsOptions(FileIndex, SpecifiedChainsAndLigandsInfo)
        
        OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"].append(SpecifiedChainsAndLigandsInfo)
        
        CheckPresenceOfValidLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo)
        
def ProcessResidueTypesAndSurfaceAndChainSelectionsOptions(FileIndex, SpecifiedChainsAndLigandsInfo):
    """Process residue types, surface, and chian selections options for chains and pockets."""

    SpecifiedChainsAndLigandsInfo["ChainSurfaces"] = {}
    SpecifiedChainsAndLigandsInfo["SurfaceChain"] = {}
    SpecifiedChainsAndLigandsInfo["SurfaceChainElectrostatics"] = {}
    
    SpecifiedChainsAndLigandsInfo["PocketSurfaces"] = {}
    SpecifiedChainsAndLigandsInfo["SurfacePocket"] = {}
    SpecifiedChainsAndLigandsInfo["SurfacePocketElectrostatics"] = {}
    
    SpecifiedChainsAndLigandsInfo["ResidueTypesChain"] = {}
    SpecifiedChainsAndLigandsInfo["ResidueTypesPocket"] = {}

    SpecifiedChainsAndLigandsInfo["ChainSelections"] = {}
    SpecifiedChainsAndLigandsInfo["ChainSelectionsSurfaces"] = {}
    SpecifiedChainsAndLigandsInfo["SurfaceChainSelections"] = {}
    
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

        # A generic color surafce is always created...
        ChainSurfaces = True
        SpecifiedChainsAndLigandsInfo["ChainSurfaces"][ChainID] = ChainSurfaces
        
        # Process residue types for chains...
        if re.match("^auto$", OptionsInfo["ResidueTypesChain"], re.I):
            ResidueTypesChain = True if AminoAcidsPresent else False
        else:
            ResidueTypesChain = True if re.match("^yes$", OptionsInfo["ResidueTypesChain"], re.I) else False
        SpecifiedChainsAndLigandsInfo["ResidueTypesChain"][ChainID] = ResidueTypesChain

        # Process chain selections...
        ChainSelections = True if len(OptionsInfo["ChainSelectionsInfo"]["Names"]) else False
        SpecifiedChainsAndLigandsInfo["ChainSelections"][ChainID] = ChainSelections
        
        # Process surfaces for chain selections...
        if re.match("^auto$", OptionsInfo["SelectionsChainSurface"], re.I):
            SurfaceChainSelections = True if AminoAcidsPresent else False
        else:
            SurfaceChainSelections = True if re.match("^yes$", OptionsInfo["SelectionsChainSurface"], re.I) else False
        SpecifiedChainsAndLigandsInfo["SurfaceChainSelections"][ChainID] = SurfaceChainSelections
        
        # A generic color surafce is always created...
        ChainSelectionsSurfaces = True
        SpecifiedChainsAndLigandsInfo["ChainSelectionsSurfaces"][ChainID] = ChainSelectionsSurfaces
        
        # Process residue types and surfaces for pockets...
        SpecifiedChainsAndLigandsInfo["PocketSurfaces"][ChainID] = {}
        SpecifiedChainsAndLigandsInfo["SurfacePocket"][ChainID] = {}
        SpecifiedChainsAndLigandsInfo["SurfacePocketElectrostatics"][ChainID] = {}
        
        SpecifiedChainsAndLigandsInfo["ResidueTypesPocket"][ChainID] = {}
        
        for LigandID in SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID]:
            if re.match("^auto$", OptionsInfo["PocketSurface"], re.I):
                PocketSurface = True if AminoAcidsPresent else False
            else:
                PocketSurface = True if re.match("^yes$", OptionsInfo["PocketSurface"], re.I) else False
            SpecifiedChainsAndLigandsInfo["SurfacePocket"][ChainID][LigandID] = PocketSurface
            
            if re.match("^auto$", OptionsInfo["PocketSurfaceElectrostatics"], re.I):
                PocketSurfaceElectrostatics = True if AminoAcidsPresent else False
            else:
                PocketSurfaceElectrostatics = True if re.match("^yes$", OptionsInfo["PocketSurfaceElectrostatics"], re.I) else False
            SpecifiedChainsAndLigandsInfo["SurfacePocketElectrostatics"][ChainID][LigandID] = PocketSurfaceElectrostatics

            # A generic color surafce is always displayed...
            PocketSurfaces = True
            SpecifiedChainsAndLigandsInfo["PocketSurfaces"][ChainID][LigandID] = PocketSurfaces
        
            if re.match("^auto$", OptionsInfo["PocketResidueTypes"], re.I):
                PocketResidueTypes = True if AminoAcidsPresent else False
            else:
                PocketResidueTypes = True if re.match("^yes$", OptionsInfo["PocketResidueTypes"], re.I) else False
            SpecifiedChainsAndLigandsInfo["ResidueTypesPocket"][ChainID][LigandID] = PocketResidueTypes
    
    # Delete loaded object...
    pymol.cmd.delete(MolName)

def GetChainAloneResidueTypesStatus(FileIndex, ChainID):
    """Get status of residue types for chain alone object."""

    Status = OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["ResidueTypesChain"][ChainID]
    
    return Status

def GetPocketResidueTypesStatus(FileIndex, ChainID, LigandID):
    """Get status of residue types for a pocket."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["ResidueTypesPocket"][ChainID][LigandID]

def GetChainAloneContainsSurfacesStatus(FileIndex, ChainID):
    """Get status of surfaces present in chain alone object."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["ChainSurfaces"][ChainID]

def GetPocketContainsSurfaceStatus(FileIndex, ChainID, LigandID):
    """Get status of surfaces present in a pocket."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["PocketSurfaces"][ChainID][LigandID]

def GetChainAloneSurfaceChainStatus(FileIndex, ChainID):
    """Get status of hydrophobic surfaces for chain alone object."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["SurfaceChain"][ChainID]

def GetChainAloneSurfaceChainElectrostaticsStatus(FileIndex, ChainID):
    """Get status of electrostatics surfaces for chain alone object."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["SurfaceChainElectrostatics"][ChainID]

def GetPocketSurfaceChainStatus(FileIndex, ChainID, LigandID):
    """Get status of hydrophobic surfaces for a pocket."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["SurfacePocket"][ChainID][LigandID]

def GetPocketSurfaceChainElectrostaticsStatus(FileIndex, ChainID, LigandID):
    """Get status of hydrophobic surfaces for a pocket."""

    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["SurfacePocketElectrostatics"][ChainID][LigandID]

def GetChainAloneContainsSelectionsStatus(FileIndex, ChainID):
    """Get status of selections present in chain alone object."""
    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["ChainSelections"][ChainID]

def GetChainAloneContainsChainSelectionSurfacesStatus(FileIndex, ChainID):
    """Get status of chain selections surfaces present in chain alone object."""
    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["ChainSelectionsSurfaces"][ChainID]

def GetChainAloneSurfaceChainSelectionStatus(FileIndex, ChainID):
    """Get status of hydrophobic surfaces for chain alone object."""
    return OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]["SurfaceChainSelections"][ChainID]

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

def ProcessChainSelections():
    """Process custom selections for chains. """

    ChainSelectionsInfo = PyMOLUtil.ProcessChainSelectionsOptionsInfo("--selectionsChain", OptionsInfo["SelectionsChain"])
    OptionsInfo["ChainSelectionsInfo"] = ChainSelectionsInfo

def ProcessOptions():
    """Process and validate command line arguments and options"""

    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Align"] = True if re.match("^Yes$", Options["--align"], re.I) else False
    OptionsInfo["AlignMethod"] = Options["--alignMethod"].lower()
    OptionsInfo["AlignMode"] = Options["--alignMode"]
    
    OptionsInfo["AllowEmptyObjects"] = True if re.match("^Yes$", Options["--allowEmptyObjects"], re.I) else False

    OptionsInfo["Infiles"] = Options["--infiles"]
    OptionsInfo["InfilesNames"] =  Options["--infileNames"]

    OptionsInfo["AlignRefFile"] = Options["--alignRefFile"]
    if re.match("^FirstInputFile$", Options["--alignRefFile"], re.I):
        OptionsInfo["RefFileName"] = OptionsInfo["InfilesNames"][0]
    else:
        OptionsInfo["RefFileName"] = Options["--alignRefFile"]
    
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

    OptionsInfo["LabelFontID"] = int(Options["--labelFontID"])
    
    OptionsInfo["PocketContactsLigandColor"] = Options["--pocketContactsLigandColor"]
    OptionsInfo["PocketContactsLigandHydrophobicColor"] = Options["--pocketContactsLigandHydrophobicColor"]
    OptionsInfo["PocketContactsSolventColor"] = Options["--pocketContactsSolventColor"]
    OptionsInfo["PocketContactsInorganicColor"] = Options["--pocketContactsInorganicColor"]

    OptionsInfo["PocketContactsCutoff"] = float(Options["--pocketContactsCutoff"])
    OptionsInfo["PocketDistanceCutoff"] = float(Options["--pocketDistanceCutoff"])
    
    OptionsInfo["PocketLabelColor"] = Options["--pocketLabelColor"]
    
    OptionsInfo["PocketResidueTypes"] = Options["--pocketResidueTypes"]
    OptionsInfo["PocketSurface"] = Options["--pocketSurface"]
    OptionsInfo["PocketSurfaceElectrostatics"] = Options["--pocketSurfaceElectrostatics"]
    
    OptionsInfo["ResidueTypesChain"] = Options["--residueTypesChain"]
    OptionsInfo["ResidueTypes"] = Options["--residueTypes"]
    ProcessResidueTypes()

    OptionsInfo["SelectionsChain"] = Options["--selectionsChain"]
    OptionsInfo["SelectionsChainSurface"] = Options["--selectionsChainSurface"]
    OptionsInfo["SelectionsChainStyle"] = Options["--selectionsChainStyle"]
    ProcessChainSelections()
    
    OptionsInfo["SurfaceChain"] = Options["--surfaceChain"]
    OptionsInfo["SurfaceChainElectrostatics"] = Options["--surfaceChainElectrostatics"]
    
    OptionsInfo["SurfaceChainComplex"] = True if re.match("^Yes$", Options["--surfaceChainComplex"], re.I) else False
    OptionsInfo["SurfaceComplex"] = True if re.match("^Yes$", Options["--surfaceComplex"], re.I) else False
    
    OptionsInfo["SurfaceColor"] = Options["--surfaceColor"]
    OptionsInfo["SurfaceColorPalette"] = Options["--surfaceColorPalette"]
    OptionsInfo["SurfaceAtomTypesColors"] = Options["--surfaceAtomTypesColors"]
    ProcessSurfaceAtomTypesColors()
    
    OptionsInfo["SurfaceTransparency"] = float(Options["--surfaceTransparency"])
    
    RetrieveInfilesInfo()
    RetrieveRefFileInfo()
    
    OptionsInfo["ChainIDs"] = Options["--chainIDs"]
    OptionsInfo["LigandIDs"] = Options["--ligandIDs"]
    
    ProcessChainAndLigandIDs()

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
    
    MiscUtil.ValidateOptionTextValue("--align", Options["--align"], "yes no")
    MiscUtil.ValidateOptionTextValue("--alignMethod", Options["--alignMethod"], "align cealign super")
    MiscUtil.ValidateOptionTextValue("--alignMode", Options["--alignMode"], "FirstChain Complex")
    
    MiscUtil.ValidateOptionTextValue("--allowEmptyObjects", Options["--allowEmptyObjects"], "yes no")

    # Expand infiles to handle presence of multiple input files...
    InfileNames = MiscUtil.ExpandFileNames(Options["--infiles"], ",")
    if not len(InfileNames):
        MiscUtil.PrintError("No input files specified for \"-i, --infiles\" option")

    # Validate file extensions...
    for Infile in InfileNames:
        MiscUtil.ValidateOptionFilePath("-i, --infiles", Infile)
        MiscUtil.ValidateOptionFileExt("-i, --infiles", Infile, "pdb cif")
        MiscUtil.ValidateOptionsDistinctFileNames("-i, --infiles", Infile, "-o, --outfile", Options["--outfile"])
    Options["--infileNames"] = InfileNames
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "pml pse")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])

    if re.match("^yes$", Options["--align"], re.I):
        if not re.match("^FirstInputFile$", Options["--alignRefFile"], re.I):
            AlignRefFile = Options["--alignRefFile"]
            MiscUtil.ValidateOptionFilePath("--alignRefFile", AlignRefFile)
            MiscUtil.ValidateOptionFileExt("--alignRefFile", AlignRefFile, "pdb cif")
            MiscUtil.ValidateOptionsDistinctFileNames("--AlignRefFile", AlignRefFile, "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionTextValue("--PMLOut", Options["--PMLOut"], "yes no")
    MiscUtil.ValidateOptionIntegerValue("--labelFontID", Options["--labelFontID"], {})

    MiscUtil.ValidateOptionFloatValue("--pocketContactsCutoff", Options["--pocketContactsCutoff"], {">": 0.0})
    MiscUtil.ValidateOptionFloatValue("--pocketDistanceCutoff", Options["--pocketDistanceCutoff"], {">": 0.0})
    if (float(Options["--pocketContactsCutoff"]) > float(Options["--pocketDistanceCutoff"])):
        MiscUtil.PrintError("The value, %s, specified using option \"--pocketContactsCutoff\" must be less than value, %s, specified using \"-pocketDistanceCutoff\" option." % (Options["--pocketContactsCutoff"], Options["--pocketDistanceCutoff"]))
        
    MiscUtil.ValidateOptionTextValue("--pocketResidueTypes", Options["--pocketResidueTypes"], "yes no auto")
    MiscUtil.ValidateOptionTextValue("--pocketSurface", Options["--pocketSurface"], "yes no auto")
    MiscUtil.ValidateOptionTextValue("--pocketSurfaceElectrostatics", Options["--pocketSurfaceElectrostatics"], "yes no auto")
    
    MiscUtil.ValidateOptionTextValue("--residueTypesChain", Options["--residueTypesChain"], "yes no auto")
    
    MiscUtil.ValidateOptionTextValue("--surfaceComplex", Options["--surfaceComplex"], "yes no")
    MiscUtil.ValidateOptionTextValue("--surfaceChainComplex", Options["--surfaceChainComplex"], "yes no")
    MiscUtil.ValidateOptionTextValue("--surfaceChain", Options["--surfaceChain"], "yes no auto")
    MiscUtil.ValidateOptionTextValue("--surfaceChainElectrostatics", Options["--surfaceChainElectrostatics"], "yes no auto")
    
    MiscUtil.ValidateOptionTextValue("--surfaceColorPalette", Options["--surfaceColorPalette"], "RedToWhite WhiteToGreen")
    MiscUtil.ValidateOptionFloatValue("--surfaceTransparency", Options["--surfaceTransparency"], {">=": 0.0, "<=": 1.0})
    
# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLVisualizeMacromolecules.py - Visualize macromolecules

Usage:
    PyMOLVisualizeMacromolecules.py [--align <yes or no>] [--alignMethod <align, cealign, super>]
                                    [--alignMode <FirstChain or Complex>] [--alignRefFile <filename>]
                                    [--allowEmptyObjects <yes or no>] [--chainIDs <First, All or ID1,ID2...>]
                                    [--ligandIDs <Largest, All or ID1,ID2...>] [--labelFontID <number>]
                                    [--PMLOut <yes or no>] [--pocketContactsInorganicColor <text>]
                                    [--pocketContactsLigandColor <text>] [--pocketContactsLigandHydrophobicColor <text>]
                                    [--pocketContactsSolventColor <text>] [--pocketContactsCutoff <number>]
                                    [--pocketDistanceCutoff <number>] [--pocketLabelColor <text>] [--pocketResidueTypes <yes or no>]
                                    [--pocketSurface <yes or no>] [--pocketSurfaceElectrostatics <yes or no>]
                                    [--residueTypes <Type,Color,ResNames,...>] [--residueTypesChain <yes or no>]
                                    [--selectionsChain <ObjectName,SelectionSpec,...>] [--selectionsChainSurface <yes or no>]
                                    [--selectionsChainStyle <DisplayStyle>] [--surfaceChain <yes or no>] [--surfaceChainElectrostatics <yes or no>]
                                    [--surfaceChainComplex <yes or no>] [--surfaceComplex <yes or no>]
                                    [--surfaceColor <ColorName>] [--surfaceColorPalette <RedToWhite or WhiteToGreen>]
                                    [--surfaceAtomTypesColors <ColorType,ColorSpec,...>]
                                    [--surfaceTransparency <number>] [--overwrite] [-w <dir>] -i <infile1,infile2,infile3...> -o <outfile>
    PyMOLVisualizeMacromolecules.py -h | --help | -e | --examples

Description:
    Generate PyMOL visualization files for viewing surfaces, chains, ligands, ligand
    binding pockets, and interactions between ligands and binding pockets in
    macromolecules including proteins and nucleic acids.

    The supported input file format are: PDB (.pdb), CIF (.cif)

    The supported output file formats are: PyMOL script file (.pml), PyMOL session
    file (.pse)

    A variety of PyMOL groups and objects may be  created for visualization of
    macromolecules. These groups and objects correspond to complexes, surfaces,
    chains, ligands, inorganics, ligand binding pockets, pocket, polar interactions,
    and pocket hydrophobic surfaces. A complete hierarchy of all possible PyMOL
    groups and objects is shown below:
    
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
                    .Selections
                        .<Name1>
                            .Selection
                            .Surface
                                .Surface
                                .Hydrophobicity
                                .Hydrophobicity_Charge
                        .<Name2>
                            ... ... ..
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
                .Solvent
                .Inorganic
                .Ligand<ID>
                    .Ligand
                        .Ligand
                        .BallAndStick
                    .Pocket
                        .Pocket
                        .Polar_Contacts
                        .Hydrophobic_Contacts
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
                        .Surfaces
                            .Surface
                                .Surface
                                .Hydrophobicity
                                .Hydrophobicity_Charge
                                .Vacuum_Electrostatics
                                    .Contact_Potentials
                                    .Map
                                    .Legend
                            .Cavity
                                .Surface
                                .Hydrophobicity
                                .Hydrophobicity_Charge
                                .Vacuum_Electrostatics
                                    .Contact_Potentials
                                    .Map
                                    .Legend
                    .Pocket_Solvent
                        .Pocket_Solvent
                        .Polar_Contacts
                    .Pocket_Inorganic
                        .Pocket_Inorganic
                        .Polar_Contacts
                .Ligand<ID>
                    .Ligand
                        ... ... ...
                    .Pocket
                        ... ... ...
                    .Pocket_Solvent
                        ... ... ...
                    .Pocket_Inorganic
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
    
    The hydrophobic and electrostatic surfaces are not created for complete complex
    and chain complex in input file(s) by default. A word to the wise: The creation of
    surface objects may slow down loading of PML file and generation of PSE file, based
    on the size of input complexes. The generation of PSE file may also fail.

Options:
    -a, --align <yes or no>  [default: no]
        Align input files to a reference file before visualization.
    --alignMethod <align, cealign, super>  [default: super]
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
    --allowEmptyObjects <yes or no>  [default: no]
        Allow creation of empty PyMOL objects corresponding to solvent and
        inorganic atom selections across chains and ligands in input file(s). By
        default, the empty objects are marked for deletion.
    -c, --chainIDs <First, All or ID1,ID2...>  [default: First]
        List of chain IDs to use for visualizing macromolecules. Possible values:
        First, All, or a comma delimited list of chain IDs. The default is to use the
        chain ID for the first chain in each input file.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infiles <infile1,infile2,infile3...>
        Input file names.
    -l, --ligandIDs <Largest, All or ID1,ID2...>  [default: Largest]
        List of ligand IDs present in chains for visualizing macromolecules to
        highlight ligand interactions. Possible values: Largest, All, or a comma
        delimited list of ligand IDs. The default is to use the largest ligand present
        in all or specified chains in each input file.
        
        Ligands are identified using organic selection operator available in PyMOL.
        It'll also  identify buffer molecules as ligands. The largest ligand contains
        the highest number of heavy atoms.
    --labelFontID <number>  [default: 7]
        Font ID for drawing labels. Default: 7 (Sans Bold). Valid values: 5 to 16.
        The specified value must be a valid PyMOL font ID. No validation is
        performed. The complete lists of valid font IDs is available at:
        pymolwiki.org/index.php/Label_font_id. Examples: 5 - Sans;
        7 - Sans Bold; 9 - Serif; 10 - Serif Bold.
    -o, --outfile <outfile>
        Output file name.
    -p, --PMLOut <yes or no>  [default: yes]
        Save PML file during generation of PSE file.
    --pocketContactsInorganicColor <text>  [default: deepsalmon]
        Color for drawing polar contacts between inorganic and pocket residues.
        The specified value must be valid color. No validation is performed.
    --pocketContactsLigandColor <text>  [default: orange]
        Color for drawing polar contacts between ligand and pocket residues.
        The specified value must be valid color. No validation is performed.
    --pocketContactsLigandHydrophobicColor <text>  [default: purpleblue]
        Color for drawing hydrophobic contacts between ligand and pocket residues.
        The specified value must be valid color. No validation is performed. The
        hydrophobic contacts are shown between pairs of carbon atoms not
        connected to hydrogen bond donor or acceptors atoms as identified
        by PyMOL.
    --pocketContactsSolventColor <text>  [default: marine]
        Color for drawing polar contacts between solvent and pocket residues..
        The specified value must be valid color. No validation is performed.
    --pocketContactsCutoff <number>  [default: 4.0]
        Distance in Angstroms for identifying polar and hyrdophobic contacts
        between atoms in pocket residues and ligands.
    --pocketDistanceCutoff <number>  [default: 5.0]
        Distance in Angstroms for identifying pocket residues around ligands.
    --pocketLabelColor <text>  [default: magenta]
        Color for drawing residue or atom level labels for a pocket. The specified
        value must be valid color. No validation is performed.
    --pocketResidueTypes <yes or no>  [default: auto]
        Pocket residue types. The residue groups are generated using residue types,
        colors, and names specified by '--residueTypes' option. It is only valid for
        amino acids.  By default, the residue type groups are automatically created
        for pockets containing amino acids and skipped for chains only containing
        nucleic acids.
    --pocketSurface <yes or no>  [default: auto]
        Surfaces around pocket residues colored by hydrophobicity alone and
        both hydrophobicity and charge. The hydrophobicity surface is colored
        at residue level using Eisenberg hydrophobicity scale for residues and color
        gradient specified by '--surfaceColorPalette' option. The  hydrophobicity and
        charge surface is colored [ Ref 140 ] at atom level using colors specified for
        groups of atoms by '--surfaceAtomTypesColors' option. This scheme allows
        simultaneous mapping of hyrophobicity and charge values on the surfaces.
        
        The cavity surfaces around ligands are also generated. These surfaces are
        colored by hydrophobicity along and both hydrophobicity and charge.
        
        This option is only valid for amino acids. By default, both surfaces are
        automatically created for pockets containing amino acids and skipped for
        pockets containing only nucleic acids.
        
        In addition, generic pockect and cavity surfaces colored by '--surfaceColor'
        are always created for pocket residues containing amino acids and nucleic
        acids.
    --pocketSurfaceElectrostatics <yes or no>  [default: no]
        Vacuum electrostatics contact potential surface around pocket residues.
        A word to the wise from PyMOL documentation: The computed protein
        contact potentials are only qualitatively useful, due to short cutoffs,
        truncation, and lack of solvent "screening".
        
        The cavity surface around ligands is also generated. This surface is
        colored by vacuum electrostatics contact potential.
        
        This option is only valid for amino acids. By default, the electrostatics surface
        is automatically created for chains containing amino acids and skipped for chains
        containing only nucleic acids.
    -r, --residueTypes <Type,Color,ResNames,...>  [default: auto]
        Residue types, colors, and names to generate for residue groups during
        '--pocketResidueTypes' and '--residueTypesChain' option. It is only
        valid for amino acids.
        
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
    --residueTypesChain <yes or no>  [default: auto]
        Chain residue types. The residue groups are generated using residue types,
        colors, and names specified by '--residueTypes' option. It is only valid for
        amino acids.  By default, the residue type groups are automatically created
        for chains containing amino acids and skipped for chains only containing
        nucleic acids.
    --selectionsChain <ObjectName,SelectionSpec,...>  [default: None]
        Custom selections for chains. It is a pairwise list of comma delimited values
        corresponding to PyMOL object names and selection specifications.  The
        selection specification must be a valid PyMOL specification. No validation is
        performed.
        
        The PyMOL objects are created for each chain corresponding to the
        specified selections. The display style for PyMOL objects is set using
        value of '--selectionsChainStyle' option.
        
        The specified selection specification is automatically appended to appropriate
        chain specification before creating PyMOL objects.
        
        For example, the following specification for '--selectionsChain' option will
        generate PyMOL objects for chains containing Cysteines and Serines:
            
            Cysteines,resn CYS,Serines,resn SER
            
    --selectionsChainSurface <yes or no>  [default: auto]
        Surfaces around individual chain selections colored by hydrophobicity alone
        and both hydrophobicity and charge. This option is similar to '--surafceChain'
        options for creating surfaces for chain. Additional details are available in the
        documentation section for '--surfaceChain' options.
    --selectionsChainStyle <DisplayStyle>  [default: sticks]
        Display style for PyMOL objects created for '--selectionsChain' option. It
        must be a valid PyMOL display style. No validation is performed.
    --surfaceChain <yes or no>  [default: auto]
        Surfaces around individual chain colored by hydrophobicity alone and
        both hydrophobicity and charge. The hydrophobicity surface is colored
        at residue level using Eisenberg hydrophobicity scale for residues and color
        gradient specified by '--surfaceColorPalette' option. The  hydrophobicity and
        charge surface is colored [ Ref 140 ] at atom level using colors specified for
        groups of atoms by '--surfaceAtomTypesColors' option. This scheme allows
        simultaneous mapping of hyrophobicity and charge values on the surfaces.
        
        This option is only valid for amino acids. By default, both surfaces are
        automatically created for chains containing amino acids and skipped for
        chains containing only nucleic acids.
        
        In addition, generic surfaces colored by '--surfaceColor' are always created
        for chain residues containing amino acids and nucleic acids.
    --surfaceChainElectrostatics <yes or no>  [default: no]
        Vacuum electrostatics contact potential surface and volume around individual
        chain. A word to the wise from PyMOL documentation: The computed protein
        contact potentials are only qualitatively useful, due to short cutoffs,
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
    --surfaceAtomTypesColors <ColorType,ColorSpec,...>  [default: auto]
        Atom colors for generating surfaces colored by hyrophobicity and charge
        around chains and pockets in proteins. It's a pairwise comma delimited list
        of atom color type and color specification for goups of atoms.
        
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
            
    --surfaceColor <ColorName>  [default: lightblue]
        Color name for surfaces around chains and pockets. This color is not used
        for surfaces colored by hydrophobicity and charge. The color name must be
        a valid PyMOL name.
    --surfaceColorPalette <RedToWhite or WhiteToGreen>  [default: RedToWhite]
        Color palette for hydrophobic surfaces around chains and pockets in proteins.
        Possible values: RedToWhite or WhiteToGreen from most hydrophobic amino
        acid to least hydrophobic. The colors values for amino acids are taken from
        color_h script available as part of the Script Library at PyMOL Wiki.
    --surfaceTransparency <number>  [default: 0.25]
        Surface transparency for molecular surfaces.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To visualize the first chain, the largest ligand in the first chain, and ligand
    binding pockets to highlight ligand interaction with pocket resiudes, solvents
    and inorganics, in a PDB file, and generate a PML file, type:

        % PyMOLVisualizeMacromolecules.py -i Sample4.pdb -o Sample4.pml

    To visualize the first chain along with all cysteines and serines, the largest
    ligand in the first chain, and ligand binding pockets to highlight ligand
    interaction with pocket resiudes, solvents and inorganics, in a PDB file,
    and generate a PML file, type:

        % PyMOLVisualizeMacromolecules.py -i Sample4.pdb -o Sample4.pml
          --selectionsChain "Cysteines,resn cys,Serines,resn ser"

    To visualize all chains, all ligands in all chains, and all ligand binding pockets to
    highlight ligand interaction with pocket resiudes, solvents and inorganics, in a
    PDB file, and generate a PML file, type:

        % PyMOLVisualizeMacromolecules.py -c All -l All -i Sample4.pdb -o
          Sample4.pml

    To visualize all chains, ligands, and ligand binding pockets along with displaying
    all hydrophibic surfaces and chain electrostatic surface, in a PDB file, and
    generate a PML file, type:

        % PyMOLVisualizeMacromolecules.py -c All -l All
          --surfaceChainElectrostatics yes --surfaceChainComplex yes
          --surfaceComplex yes -i Sample4.pdb -o Sample4.pml

    To visualize chain E, ligand ADP in chain E, and ligand binding pockets to
    highlight ligand interaction with pocket resiudes, solvents and inorganics,
    in a PDB file, and generate a PML file, type:

        % PyMOLVisualizeMacromolecules.py -c E -l ADP -i Sample3.pdb
          -o Sample3.pml

    To visualize chain E, ligand ADP in chain E, and ligand binding pockets to
    highlight ligand interaction with pocket resiudes, solvents and inorganics,
    in a PDB file, and generate a PSE file, type:

        % PyMOLVisualizeMacromolecules.py -c E -l ADP -i Sample3.pdb
          -o Sample3.pse

    To visualize the first chain, the largest ligand in the first chain, and ligand
    binding pockets to highlight ligand interaction with pocket resiudes, solvents
    and inorganics, in PDB files, along with aligning first chain in each input file to
    the first chain in first input file, and generate a PML file, type:

        % PyMOLVisualizeMacromolecules.py --align yes -i
          "Sample5.pdb,Sample6.pdb,Sample7.pdb" -o SampleOut.pml

    To visualize all chains, all ligands in all chains, and all ligand binding pockets to
    highlight ligand interaction with pocket resiudes, solvents and inorganics, in
    PDB files, along with aligning first chain in each input file to the first chain in
    first input file, and generate a PML file, type:

        % PyMOLVisualizeMacromolecules.py --align yes  -c All -l All -i
          "Sample5.pdb,Sample6.pdb,Sample7.pdb" -o SampleOut.pml

    To visualize all chains, all ligands in all chains, and all ligand binding pockets to
    highlight ligand interaction with pocket resiudes, solvents and inorganics, in
    PDB files, along with aligning first chain in each input file to the first chain in a
    specified PDB file using a specified alignment method, and generate a PML
    file, type:

        % PyMOLVisualizeMacromolecules.py --align yes  --alignMode FirstChain
          --alignRefFile Sample5.pdb --alignMethod super   -c All  -l All -i
          "Sample5.pdb,Sample6.pdb,Sample7.pdb" -o SampleOut.pml

Author:
    Manish Sud(msud@san.rr.com)

See also:
    DownloadPDBFiles.pl, PyMOLVisualizeCavities.py,
    PyMOLVisualizeCryoEMDensity.py, PyMOLVisualizeElectronDensity.py,
    PyMOLVisualizeInterfaces.py, PyMOLVisualizeSurfaceAndBuriedResidues.py

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
