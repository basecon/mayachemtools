#!/bin/env python
#
# File: PyMOLVisualizeCryoEMDensity.py
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
import xml.etree.ElementTree as ElementTree

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
    GenerateCryoEMDensityVisualization()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def GenerateCryoEMDensityVisualization():
    """Generate cryo-EM density visualization."""
    
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
    PMLCmds.append("""cmd.set("mesh_width", %.2f)""" % (OptionsInfo["MeshWidth"]))
    PMLCmds.append("""cmd.set("transparency", %.2f, "", 0)""" % (OptionsInfo["SurfaceTransparency"]))
    PMLCmds.append("""cmd.set("label_font_id", %s)""" % (OptionsInfo["LabelFontID"]))

    if OptionsInfo["VolumeColorRampCreate"]:
        ColorRampName = OptionsInfo["VolumeColorRampName"]
        ContourLevel = OptionsInfo["VolumeColorRampContourLevel"]
        LowerContourLevel = ContourLevel - 0.3
        UpperContourLevel = ContourLevel + 0.3
        PMLCmds.append("""cmd.volume_ramp_new("%s", "%.2f blue 0.00 %.2f cyan 0.20 %.2f blue 0.00")""" % (ColorRampName, LowerContourLevel, ContourLevel, UpperContourLevel))
        
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
    """Write out PML for viewing polymer complex along with cryo-EM density."""

    # Setup complex...
    Infile = OptionsInfo["InfilesInfo"]["InfilesNames"][FileIndex]
    PML = PyMOLUtil.SetupPMLForPolymerComplexView(PyMOLObjectNames["Complex"], Infile, True)
    OutFH.write("""\n""\n"Loading %s and setting up view for complex..."\n""\n""" % Infile)
    OutFH.write("%s\n" % PML)

    if OptionsInfo["Align"]:
        # No need to align complex on to itself...
        if not (re.match("^FirstInputFile$", OptionsInfo["AlignRefFile"], re.I) and FirstComplex):
            WriteAlignComplex(OutFH, FileIndex, PyMOLObjectNames)

    # Setup cryo-EM density maps and meshes...
    DensityMapFile = OptionsInfo["DensityMapFilesNames"][FileIndex]
    ContourLevel = OptionsInfo["MeshLevels"][FileIndex]
    WriteComplexCryoEMDensityMapView(OutFH, PyMOLObjectNames, DensityMapFile, ContourLevel)

    # Setup complex group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["ComplexGroup"], PyMOLObjectNames["ComplexGroupMembers"], False, "close")

def WriteComplexCryoEMDensityMapView(OutFH, PyMOLObjectNames, MapFileName, ContourLevel):
    """Write out PML for viewing cryoEM density map."""

    # Load cryo-EM density map and setup mesh views...
    Info = """\
""
"Loading cryo-EM density map %s and setting up mesh view for complex..."
"" """ % MapFileName
    OutFH.write("\n%s\n" % Info)

    MapName = PyMOLObjectNames["ComplexCryoEMMap"]
    ComplexName = PyMOLObjectNames["Complex"]
    
    Color = OptionsInfo["MeshColor"]
    VolumeColorRamp = OptionsInfo["VolumeColorRampName"]
    
    VolumeName = PyMOLObjectNames["ComplexCryoEMVolume"]
    MeshName = PyMOLObjectNames["ComplexCryoEMMesh"]
    SurfaceName = PyMOLObjectNames["ComplexCryoEMSurface"]
    
    AlignMapToObjectName = ComplexName if OptionsInfo["Align"] else None
    EnableMap = True
    PML = SetupPMLForCryoEMDensityMap(MapFileName, MapName, AlignMapToObjectName, EnableMap)
    OutFH.write("%s\n" % PML)

    EnableMesh = OptionsInfo["MeshComplex"]
    
    EnableVolume = OptionsInfo["VolumeComplex"]
    if EnableVolume and EnableMesh:
        EnableVolume = False
        
    EnableSurface = OptionsInfo["SurfaceComplex"]
    if EnableSurface and (EnableVolume or EnableMesh):
        EnableSurface = False
    
    if OptionsInfo["VolumeComplex"]:
        PML = SetupPMLForCryoEMDensityVolume(MapName, VolumeName, VolumeColorRamp, Enable = EnableVolume, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)
    
    if OptionsInfo["MeshComplex"]:
        PML = SetupPMLForCryoEMDensityMesh(MapName, MeshName, ContourLevel, Color, Enable = EnableMesh, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)

    if OptionsInfo["SurfaceComplex"]:
        PML = SetupPMLForCryoEMDensitySurface(MapName, SurfaceName, ContourLevel, Color, Enable = EnableSurface, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)

    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["ComplexCryoEMGroup"], PyMOLObjectNames["ComplexCryoEMGroupMembers"], True, "close")
    
def WriteChainView(OutFH, FileIndex, PyMOLObjectNames, ChainID):
    """Write out PML for viewing chain."""
    
    OutFH.write("""\n""\n"Setting up views for chain %s..."\n""\n""" % ChainID)
    
    ChainComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
    
    # Setup chain complex group view...
    WriteChainComplexAndMeshViews(OutFH, FileIndex, PyMOLObjectNames, ChainID)

    # Setup chain view...
    WriteChainAloneViews(OutFH, FileIndex, PyMOLObjectNames, ChainID)
    
    # Setup chain solvent view...
    PML = PyMOLUtil.SetupPMLForSolventView(PyMOLObjectNames["Chains"][ChainID]["Solvent"], ChainComplexName, False)
    OutFH.write("\n%s\n" % PML)

    # Setup chain inorganic view...
    PML = PyMOLUtil.SetupPMLForInorganicView(PyMOLObjectNames["Chains"][ChainID]["Inorganic"], ChainComplexName, False)
    OutFH.write("\n%s\n" % PML)

def WriteChainComplexAndMeshViews(OutFH, FileIndex, PyMOLObjectNames, ChainID):
    """Write chain complex and mesh views. """
    
    # Setup chain complex...
    ChainComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
    PML = PyMOLUtil.SetupPMLForPolymerChainComplexView(ChainComplexName, PyMOLObjectNames["Complex"], ChainID, True)
    OutFH.write("%s\n" % PML)

    SpecifiedChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]
    
    MeshChainComplex = SpecifiedChainsAndLigandsInfo["MeshChainComplex"][ChainID]
    VolumeChainComplex = SpecifiedChainsAndLigandsInfo["VolumeChainComplex"][ChainID]
    SurfaceChainComplex = SpecifiedChainsAndLigandsInfo["SurfaceChainComplex"][ChainID]
    
    EnableVolumeChainComplex = SpecifiedChainsAndLigandsInfo["EnableVolumeChainComplex"][ChainID]
    EnableMeshChainComplex = SpecifiedChainsAndLigandsInfo["EnableMeshChainComplex"][ChainID]
    EnableSurfaceChainComplex = SpecifiedChainsAndLigandsInfo["EnableSurfaceChainComplex"][ChainID]
    
    if MeshChainComplex or VolumeChainComplex or SurfaceChainComplex:
        # Set up cryoEM mesh and group...
        MapName = PyMOLObjectNames["ComplexCryoEMMap"]
        ContourLevel = OptionsInfo["MeshLevels"][FileIndex]
        Color = OptionsInfo["MeshColor"]
        
        VolumeColorRamp = OptionsInfo["VolumeColorRampName"]
        
        VolumeName = PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMVolume"]
        MeshName = PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMMesh"]
        SurfaceName = PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMSurface"]
        
        if VolumeChainComplex:
            PML = SetupPMLForCryoEMDensityVolume(MapName, VolumeName, VolumeColorRamp, Enable = EnableVolumeChainComplex, Selection = ChainComplexName)
            OutFH.write("\n%s\n" % PML)
            
        if MeshChainComplex:
            PML = SetupPMLForCryoEMDensityMesh(MapName, MeshName, ContourLevel, Color, Enable = EnableMeshChainComplex, Selection = ChainComplexName)
            OutFH.write("\n%s\n" % PML)
        
        if SurfaceChainComplex:
            PML = SetupPMLForCryoEMDensitySurface(MapName, SurfaceName, ContourLevel, Color, Enable = EnableSurfaceChainComplex, Selection = ChainComplexName)
            OutFH.write("\n%s\n" % PML)
        
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroupMembers"], True, "close")
        
    # Setup chain complex group...
    EnableChainComplexGroup = SpecifiedChainsAndLigandsInfo["EnableChainComplexGroup"][ChainID]
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroupMembers"], EnableChainComplexGroup, "close")
    
def WriteChainAloneViews(OutFH, FileIndex, PyMOLObjectNames, ChainID):
    """Write individual chain views. """

    ChainComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
    
    # Setup chain view...
    ChainName = PyMOLObjectNames["Chains"][ChainID]["ChainAlone"]
    PML = PyMOLUtil.SetupPMLForPolymerChainView(ChainName, ChainComplexName, Enable = True)
    OutFH.write("\n%s\n" % PML)

    # Setup chain putty by B-factor view...
    if OptionsInfo["BFactorChainCartoonPutty"]:
        BFactorPuttyName = PyMOLObjectNames["Chains"][ChainID]["ChainAloneBFactorPutty"]
        PML = PyMOLUtil.SetupPMLForBFactorPuttyView(BFactorPuttyName, ChainName, ColorPalette = OptionsInfo["BFactorColorPalette"], Enable = False)
        OutFH.write("\n%s\n" % PML)
        
    # Setup chain selections view...
    SetupChainSelectionsView(OutFH, FileIndex, PyMOLObjectNames, ChainID)
    
    # Setup chain group...
    SpecifiedChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]
    EnableChainAloneGroup = SpecifiedChainsAndLigandsInfo["EnableChainAloneGroup"][ChainID]
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"], EnableChainAloneGroup, "close")
    
def SetupChainSelectionsView(OutFH, FileIndex, PyMOLObjectNames, ChainID):
    """Setup chain selectons view. """

    if not OptionsInfo["ChainSelections"]:
        return
    
    ChainName = PyMOLObjectNames["Chains"][ChainID]["ChainAlone"]
    SelectionsGroupIDPrefix = "ChainAloneSelections"
    
    for Index in range(0, len(OptionsInfo["ChainSelectionsInfo"]["Names"])):
        SelectionName = OptionsInfo["ChainSelectionsInfo"]["Names"][Index]
        SpecifiedSelection = OptionsInfo["ChainSelectionsInfo"]["Selections"][Index]
        
        SelectionNameGroupID = SelectionName
        
        # Setup selection object...
        SelectionObjectID = "%s%sSelection" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        SelectionObjectName = PyMOLObjectNames["Chains"][ChainID][SelectionObjectID]
        SelectionCmd = "(%s and (%s))" % (ChainName, SpecifiedSelection)
        PML = PyMOLUtil.SetupPMLForSelectionDisplayView(SelectionObjectName, SelectionCmd, OptionsInfo["SelectionsChainStyle"], Enable = True)
        OutFH.write("\n%s\n" % PML)
        
        # Set up cryo-EM mesh and group...
        CryoEMVolumeID = "%s%sCryoEMVolume" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        CryoEMMeshID = "%s%sCryoEMMesh" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        CryoEMSurfaceID = "%s%sCryoEMSurface" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        CryoEMMeshGroupID = "%s%sCryoEMGroup" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        CryoEMMeshGroupMembersID = "%s%sCryoEMGroupMembers" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        
        CryoEMVolumeName = PyMOLObjectNames["Chains"][ChainID][CryoEMVolumeID]
        CryoEMMeshName = PyMOLObjectNames["Chains"][ChainID][CryoEMMeshID]
        CryoEMSurfaceName = PyMOLObjectNames["Chains"][ChainID][CryoEMSurfaceID]
        CryoEMMeshGroupName = PyMOLObjectNames["Chains"][ChainID][CryoEMMeshGroupID]
        CryoEMMeshGroupMembers = PyMOLObjectNames["Chains"][ChainID][CryoEMMeshGroupMembersID]

        MapName = PyMOLObjectNames["ComplexCryoEMMap"]
        ContourLevel = OptionsInfo["MeshLevels"][FileIndex]
        Color = OptionsInfo["MeshColor"]
        
        PML = SetupPMLForCryoEMDensityVolume(MapName, CryoEMVolumeName, OptionsInfo["VolumeColorRampName"], Enable = False, Selection = SelectionObjectName)
        OutFH.write("\n%s\n" % PML)
        
        PML = SetupPMLForCryoEMDensityMesh(MapName, CryoEMMeshName, ContourLevel, Color, Enable = True, Selection = SelectionObjectName)
        OutFH.write("\n%s\n" % PML)
        
        PML = SetupPMLForCryoEMDensitySurface(MapName, CryoEMSurfaceName, ContourLevel, Color, Enable = False, Selection = SelectionObjectName)
        OutFH.write("\n%s\n" % PML)
        
        GenerateAndWritePMLForGroup(OutFH, CryoEMMeshGroupName, CryoEMMeshGroupMembers, True, "close")
        
        # Setup groups for named selections...
        SelectionsNameGroupID = "%s%sGroup" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        SelectionsNameGroupMembersID = "%s%sGroupMembers" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID][SelectionsNameGroupID], PyMOLObjectNames["Chains"][ChainID][SelectionsNameGroupMembersID], True, "open")
    
    # Setup a group for selections...
    SelectionsGroupID = "%sGroup" % (SelectionsGroupIDPrefix)
    SelectionsGroupMembersID = "%sGroupMembers" % (SelectionsGroupIDPrefix)
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID][SelectionsGroupID], PyMOLObjectNames["Chains"][ChainID][SelectionsGroupMembersID], False, "close")
        
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
        
        # Set up cryoEM mesh and group...
        CryoEMMeshGroupID = "%sCryoEMMeshGroup" % (GroupID)
        CryoEMMeshGroupMembersID = "%sCryoEMMeshGroupMembers" % (GroupID)
        CryoEMVolumeID = "%sCryoEMVolume" % (GroupID)
        CryoEMMeshID = "%sCryoEMMesh" % (GroupID)
        CryoEMSurfaceID = "%sCryoEMSurface" % (GroupID)

        CryoEMMapName = PyMOLObjectNames["ComplexCryoEMMap"]
        CryoEMVolumeName = PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMVolumeID]
        CryoEMMeshName = PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshID]
        CryoEMSurfaceName = PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMSurfaceID]
        CryoEMMeshGroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupID]
        CryoEMMeshGroupMembers = PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupMembersID]
        
        PML = SetupPMLForCryoEMDensityVolume(CryoEMMapName, CryoEMVolumeName, OptionsInfo["VolumeColorRampName"], Enable = False, Selection = GroupTypeObjectName)
        OutFH.write("\n%s\n" % PML)

        ContourLevel = OptionsInfo["MeshLevels"][FileIndex]
        PML = SetupPMLForCryoEMDensityMesh(CryoEMMapName, CryoEMMeshName, ContourLevel, OptionsInfo["MeshColor"], Enable = True, Selection = GroupTypeObjectName)
        OutFH.write("\n%s\n" % PML)
        
        PML = SetupPMLForCryoEMDensitySurface(CryoEMMapName, CryoEMSurfaceName, ContourLevel, OptionsInfo["MeshColor"], Enable = False, Selection = GroupTypeObjectName)
        OutFH.write("\n%s\n" % PML)
        
        GenerateAndWritePMLForGroup(OutFH, CryoEMMeshGroupName, CryoEMMeshGroupMembers, True, "close")
        
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
            
        # Set up hydrophobic contacts...
        if re.match("^Pocket$", GroupID, re.I):
            HydrophobicContactsID = "%sHydrophobicContacts" % (GroupID)
            HydrophobicContactsName = PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicContactsID]
            HydrophobicContactsColor = OptionsInfo["PocketContactsLigandHydrophobicColor"]
            
            PML = PyMOLUtil.SetupPMLForHydrophobicContactsView(HydrophobicContactsName, LigandName, GroupTypeObjectName, Enable = False, Color = HydrophobicContactsColor, Cutoff = OptionsInfo["PocketContactsCutoff"])
            OutFH.write("\n%s\n" % PML)
            OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (HydrophobicContactsColor, HydrophobicContactsName))
            
        # Set up hydrophobic surface...
        if re.match("^Pocket$", GroupID, re.I) and OptionsInfo["PocketSurface"]:
            HydrophobicSurfaceID = "%sHydrophobicSurface" % (GroupID)
            HydrophobicSurfaceName = PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicSurfaceID]
            PML = PyMOLUtil.SetupPMLForHydrophobicSurfaceView(HydrophobicSurfaceName, GroupTypeObjectName, ColorPalette = "RedToWhite", Enable = False)
            OutFH.write("\n%s\n" % PML)
            
            OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (OptionsInfo["PocketLabelColor"], HydrophobicSurfaceName))
        
        # Setup group....
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

def GenerateAndWritePMLForGroup(OutFH, GroupName, GroupMembers, Enable = False, Action = "close"):
    """Generate and write PML for group. """
    
    PML = PyMOLUtil.SetupPMLForGroup(GroupName, GroupMembers, Enable, Action)
    OutFH.write("""\n""\n"Setting up group %s..."\n""\n""" % GroupName)
    OutFH.write("%s\n" % PML)

def SetupPMLForCryoEMDensityMap(MapFileName, MapName, AlignMapToObjectName = None, Enable = True):
    """Setup PML for loading and viewing cryo-EM density map. """

    PMLCmds = []
    PMLCmds.append("""cmd.load("%s", "%s")""" % (MapFileName, MapName))
    if AlignMapToObjectName is not None:
        PMLCmds.append("""cmd.matrix_copy("%s", "%s")""" % (AlignMapToObjectName, MapName))
        
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(MapName, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML
    
def SetupPMLForCryoEMDensityMesh(MapName, MeshName, SigmaLevel, Color, Enable = True, Selection = None):
    """Setup PML for cryo-EM density mesh. """

    Carve = OptionsInfo["MeshCarveRadius"]
    
    PMLCmds = []
    if Selection is None:
        PMLCmds.append("""cmd.isomesh("%s", "%s", %.1f)""" % (MeshName, MapName, SigmaLevel))
    else:
        PMLCmds.append("""cmd.isomesh("%s", "%s", %.1f, "(%s)", carve = %.1f)""" % (MeshName, MapName, SigmaLevel, Selection, Carve))
    PMLCmds.append("""util.color_deep("%s", "%s")""" % (Color, MeshName))
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(MeshName, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForCryoEMDensityVolume(MapName, VolumeName, VolumeColorRamp, Enable = True, Selection = None):
    """Setup PML for cryo-EM density volume. """

    Carve = OptionsInfo["VolumeCarveRadius"]
    
    PMLCmds = []
    if Selection is None:
        PMLCmds.append("""cmd.volume("%s", "%s", "%s")""" % (VolumeName, MapName, VolumeColorRamp))
    else:
        PMLCmds.append("""cmd.volume("%s", "%s", "%s", "(%s)", carve = %.1f)""" % (VolumeName, MapName, VolumeColorRamp, Selection, Carve))
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(VolumeName, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForCryoEMDensitySurface(MapName, SurfaceName, SigmaLevel, Color, Enable = True, Selection = None):
    """Setup PML for cryo-EM density surface. """

    Carve = OptionsInfo["MeshCarveRadius"]
    
    PMLCmds = []
    if Selection is None:
        PMLCmds.append("""cmd.isosurface("%s", "%s", %.1f)""" % (SurfaceName, MapName, SigmaLevel))
    else:
        PMLCmds.append("""cmd.isosurface("%s", "%s", %.1f, "(%s)", carve = %.1f)""" % (SurfaceName, MapName, SigmaLevel, Selection, Carve))
    PMLCmds.append("""util.color_deep("%s", "%s")""" % (Color, SurfaceName))
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(SurfaceName, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

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
        
        for LigandID in SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID]:
            # Delete ligand level objects...
            for GroupID in ["Pocket", "PocketSolvent", "PocketInorganic"]:
                GroupNameID = "%sGroup" % (GroupID)
                GroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupNameID]

                GroupTypeObjectID = "%s" % (GroupID)
                GroupTypeObjectName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupTypeObjectID]
                
                WritePMLToCheckAndDeleteEmptyObjects(OutFH, GroupTypeObjectName, GroupName)

def WritePMLToCheckAndDeleteEmptyObjects(OutFH, ObjectName, ParentObjectName = None):
    """Write PML to check and delete empty PyMOL objects. """
    
    if ParentObjectName is None:
        PML = """CheckAndDeleteEmptyObjects("%s")""" % (ObjectName)
    else:
        PML = """CheckAndDeleteEmptyObjects("%s", "%s")""" % (ObjectName, ParentObjectName)
    
    OutFH.write("%s\n" % PML)

def SetupPyMOLObjectNames(FileIndex):
    """Setup hierarchy of PyMOL groups and objects for ligand centric views of
    cryo-EM density for chains and ligands present in input file.
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
    """Stetup groups and objects for complex. """
    
    PDBFileRoot = OptionsInfo["InfilesInfo"]["InfilesRoots"][FileIndex]
    
    PDBGroupName = "%s" % PDBFileRoot
    PyMOLObjectNames["PDBGroup"] = PDBGroupName
    PyMOLObjectNames["PDBGroupMembers"] = []

    ComplexGroupName = "%s.Complex" % PyMOLObjectNames["PDBGroup"]
    PyMOLObjectNames["ComplexGroup"] = ComplexGroupName
    PyMOLObjectNames["PDBGroupMembers"].append(ComplexGroupName)
    
    PyMOLObjectNames["Complex"] = "%s.Complex" % ComplexGroupName

    CryoEMMeshGroupName = "%s.CryoEM" % (ComplexGroupName)
    CryoEMMapName = "%s.Map" % (CryoEMMeshGroupName)
    CryoEMVolumeName = "%s.Volume" % (CryoEMMeshGroupName)
    CryoEMMeshName = "%s.Mesh" % (CryoEMMeshGroupName)
    CryoEMSurfaceName = "%s.Surface" % (CryoEMMeshGroupName)
    
    PyMOLObjectNames["ComplexCryoEMGroup"] = CryoEMMeshGroupName
    PyMOLObjectNames["ComplexCryoEMMap"] = CryoEMMapName
    PyMOLObjectNames["ComplexCryoEMVolume"] = CryoEMVolumeName
    PyMOLObjectNames["ComplexCryoEMMesh"] = CryoEMMeshName
    PyMOLObjectNames["ComplexCryoEMSurface"] = CryoEMSurfaceName

    PyMOLObjectNames["ComplexCryoEMGroupMembers"] = []
    PyMOLObjectNames["ComplexCryoEMGroupMembers"].append(CryoEMMapName)
    if OptionsInfo["VolumeComplex"]:
        PyMOLObjectNames["ComplexCryoEMGroupMembers"].append(CryoEMVolumeName)
    if OptionsInfo["MeshComplex"]:
        PyMOLObjectNames["ComplexCryoEMGroupMembers"].append(CryoEMMeshName)
    if OptionsInfo["SurfaceComplex"]:
        PyMOLObjectNames["ComplexCryoEMGroupMembers"].append(CryoEMSurfaceName)
    
    PyMOLObjectNames["ComplexGroupMembers"] = []
    PyMOLObjectNames["ComplexGroupMembers"].append(PyMOLObjectNames["Complex"])
    PyMOLObjectNames["ComplexGroupMembers"].append(PyMOLObjectNames["ComplexCryoEMGroup"])
    
def SetupPyMOLObjectNamesForChain(FileIndex, PyMOLObjectNames, ChainID):
    """Setup groups and objects for chain."""
    
    PDBGroupName = PyMOLObjectNames["PDBGroup"]
    
    SpecifiedChainsAndLigandsInfo = OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"][FileIndex]
    MeshChainComplex = SpecifiedChainsAndLigandsInfo["MeshChainComplex"][ChainID]
    VolumeChainComplex = SpecifiedChainsAndLigandsInfo["VolumeChainComplex"][ChainID]
    SurfaceChainComplex = SpecifiedChainsAndLigandsInfo["SurfaceChainComplex"][ChainID]
    
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
    
    PyMOLObjectNames["Chains"][ChainID]["ChainComplex"] = "%s.Complex" % (ChainComplexGroupName)
    
    CryoEMMeshGroupName = "%s.CryoEM" % (ChainComplexGroupName)
    CryoEMVolumeName = "%s.Volume" % (CryoEMMeshGroupName)
    CryoEMMeshName = "%s.Mesh" % (CryoEMMeshGroupName)
    CryoEMSurfaceName = "%s.Surface" % (CryoEMMeshGroupName)
    
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroup"] = CryoEMMeshGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMVolume"] = CryoEMVolumeName
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMMesh"] = CryoEMMeshName
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMSurface"] = CryoEMSurfaceName
    
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroupMembers"] = []
    if VolumeChainComplex:
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroupMembers"].append(CryoEMVolumeName)
    if MeshChainComplex:
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroupMembers"].append(CryoEMMeshName)
    if SurfaceChainComplex:
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroupMembers"].append(CryoEMSurfaceName)
    
    NameIDs = ["ChainComplex"]
    if MeshChainComplex or VolumeChainComplex or SurfaceChainComplex :
        NameIDs.append("ChainComplexCryoEMGroup")
    
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroupMembers"] = []
    for NameID in NameIDs:
        Name = PyMOLObjectNames["Chains"][ChainID][NameID]
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroupMembers"].append(Name)

    # Setup up a group for individual chains...
    ChainAloneGroupName = "%s.Chain" % (ChainGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroup"] = ChainAloneGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"].append(ChainAloneGroupName)
        
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"] = []
        
    Name = "%s.Chain" % (ChainAloneGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainAlone"] = Name
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"].append(Name)
        
    if OptionsInfo["BFactorChainCartoonPutty"]:
        Name = "%s.BFactor" % (ChainAloneGroupName)
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneBFactorPutty"] = Name
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"].append(Name)
    
    if OptionsInfo["ChainSelections"]:
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

            # Add selection object to selections name group...
            SelectionObjectName = "%s.Selection" % (SelectionsNameGroupName)
            SelectionObjectID = "%s%sSelection" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
            
            PyMOLObjectNames["Chains"][ChainID][SelectionObjectID] = SelectionObjectName
            PyMOLObjectNames["Chains"][ChainID][SelectionsNameGroupMembersID].append(SelectionObjectName)
            
            # Setup cryo-EM mesh group and add it to selections name group...
            CryoEMMeshGroupName = "%s.CryoEM" % (SelectionsNameGroupName)
            CryoEMMeshGroupID = "%s%sCryoEMGroup" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
            
            PyMOLObjectNames["Chains"][ChainID][CryoEMMeshGroupID] = CryoEMMeshGroupName
            PyMOLObjectNames["Chains"][ChainID][SelectionsNameGroupMembersID].append(CryoEMMeshGroupName)
            
            # Initialize cryo-EM mesh group members...
            CryoEMMeshGroupMembersID = "%s%sCryoEMGroupMembers" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
            PyMOLObjectNames["Chains"][ChainID][CryoEMMeshGroupMembersID] = []

            # Setup members of cryo-EM mesh group...
            CryoEMVolumeName = "%s.Volume" % (CryoEMMeshGroupName)
            CryoEMVolumeID = "%s%sCryoEMVolume" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
            CryoEMMeshName = "%s.Mesh" % (CryoEMMeshGroupName)
            CryoEMMeshID = "%s%sCryoEMMesh" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
            CryoEMSurfaceName = "%s.Surface" % (CryoEMMeshGroupName)
            CryoEMSurfaceID = "%s%sCryoEMSurface" % (SelectionsGroupIDPrefix, SelectionNameGroupID)
    
            PyMOLObjectNames["Chains"][ChainID][CryoEMVolumeID] = CryoEMVolumeName
            PyMOLObjectNames["Chains"][ChainID][CryoEMMeshID] = CryoEMMeshName
            PyMOLObjectNames["Chains"][ChainID][CryoEMSurfaceID] = CryoEMSurfaceName

            PyMOLObjectNames["Chains"][ChainID][CryoEMMeshGroupMembersID].append(CryoEMVolumeName)
            PyMOLObjectNames["Chains"][ChainID][CryoEMMeshGroupMembersID].append(CryoEMMeshName)
            PyMOLObjectNames["Chains"][ChainID][CryoEMMeshGroupMembersID].append(CryoEMSurfaceName)
    
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
        
        CryoEMMeshGroupName = "%s.CryoEM" % (GroupName)
        CryoEMVolumeName = "%s.Volume" % (CryoEMMeshGroupName)
        CryoEMMeshName = "%s.Mesh" % (CryoEMMeshGroupName)
        CryoEMSurfaceName = "%s.Surface" % (CryoEMMeshGroupName)
                
        CryoEMMeshGroupID = "%sCryoEMMeshGroup" % (GroupID)
        CryoEMMeshGroupMembersID = "%sCryoEMMeshGroupMembers" % (GroupID)
        CryoEMVolumeID = "%sCryoEMVolume" % (GroupID)
        CryoEMMeshID = "%sCryoEMMesh" % (GroupID)
        CryoEMSurfaceID = "%sCryoEMSurface" % (GroupID)
        
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupID] = CryoEMMeshGroupName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMVolumeID] = CryoEMVolumeName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshID] = CryoEMMeshName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMSurfaceID] = CryoEMSurfaceName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupMembersID] = []
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupMembersID].append(CryoEMVolumeName)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupMembersID].append(CryoEMMeshName)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupMembersID].append(CryoEMSurfaceName)
                
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID] = []
        NameIDs = [GroupTypeObjectID, CryoEMMeshGroupID]
        
        for NameID in NameIDs:
            Name = PyMOLObjectNames["Ligands"][ChainID][LigandID][NameID]
            PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(Name)
        
        if re.match("^Ligand$", GroupType, re.I):
            # No other object needed for Ligand group...
            continue
        
        PolarContactsName = "%s.Polar_Contacts" % (GroupName)
        PolarContactsID = "%sPolarContacts" % (GroupID)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][PolarContactsID] = PolarContactsName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(PolarContactsName)
                
        if not re.match("^Pocket$", GroupType, re.I):
            # No other object needed for any other group besides Pocket...
            continue
        
        if not OptionsInfo["PocketSurface"]:
            continue

        HydrophobicContactsName = "%s.Hydrophobic_Contacts" % (GroupName)
        HydrophobicContactsID = "%sHydrophobicContacts" % (GroupID)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicContactsID] = HydrophobicContactsName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(HydrophobicContactsName)
        
        HydrophobicSurfaceName = "%s.Surface" % (GroupName)
        HydrophobicSurfaceID = "%sHydrophobicSurface" % (GroupID)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicSurfaceID] = HydrophobicSurfaceName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(HydrophobicSurfaceName)

def ProcessDensityMapFiles():
    """Process density map files."""
    
    DensityMapFiles = OptionsInfo["DensityMapFiles"]
    if re.match("^auto$", DensityMapFiles, re.I):
        ProcessAutoDensityMapFiles()
    else:
        ProcessSpecifiedDensityMapFiles()
        
def ProcessSpecifiedDensityMapFiles():
    """Process specified density map files."""
    
    OptionsInfo["DensityMapFilesNames"] = []
    
    MiscUtil.PrintInfo("\nProcessing cryo-EM density file names...")
    
    DensityMapFiles = re.sub(" ", "", OptionsInfo["DensityMapFiles"])
    if not DensityMapFiles:
        MiscUtil.PrintError("No valid value specified using \"--densityMapFiles\" option.")
    
    DensityMapFilesWords = DensityMapFiles.split(",")
    DensityMapFilesWordsCount = len(DensityMapFilesWords)
    
    InfilesNamesCount = len(OptionsInfo["InfilesNames"])
    
    if DensityMapFilesWordsCount != InfilesNamesCount:
        MiscUtil.PrintError("The number of comma delimited cryo-EM density map files, %d, specified using \"--DensityMapFiles\" must be equal to the number of input files, %s, specified using \"-i, --infiles\" option." % (DensityMapFilesWordsCount, InfilesNamesCount))
    
    for Index in range(0, InfilesNamesCount):
        Infile = OptionsInfo["InfilesNames"][Index]
        DensityMapFile = DensityMapFilesWords[Index]
        
        if not os.path.exists(DensityMapFile):
            MiscUtil.PrintError("The cryo-EM density ED map file, %s, specified using option \"--DensityMapFiles\", corresponding to input file, %s,  doesn't exist.\n" % (DensityMapFile, Infile))
        
        OptionsInfo["DensityMapFilesNames"].append(DensityMapFile)
    
def ProcessAutoDensityMapFiles():
    """Set up and process name of density map files."""
    
    OptionsInfo["DensityMapFilesNames"] = []
    InfilesNamesCount = len(OptionsInfo["InfilesNames"])

    MiscUtil.PrintInfo("\nSetting cryo-EM density file names...")
    
    for Index in range(0, InfilesNamesCount):
        Infile = OptionsInfo["InfilesNames"][Index]
        
        EMDBID = RetrieveEMDBID(Infile)
        if EMDBID is None:
            MiscUtil.PrintError("Failed to retrieve EMDB ID from input file %s to automatically set density map file name. Use option \"-d, --densityMapFiles \" to specify density map file name and try again." % Infile)
        
        DensityMapFile = None
        MapFileRoot = "emd_%s" % EMDBID
        MapFile1 = "%s.map.gz" % MapFileRoot
        MapFile2 = "%s.map" % MapFileRoot
        if os.path.exists(MapFile1):
            DensityMapFile = MapFile1
        elif os.path.exists(MapFile2):
            DensityMapFile = MapFile2
        else:
            MiscUtil.PrintError("Density map files %s or %s don't exist. Use option \"-d, --densityMapFiles \" to specify density map file name and try again" % (MapFile1, MapFile2))
    
        MiscUtil.PrintInfo("Setting density map file name as %s for input file %s..." % (DensityMapFile, Infile))
        OptionsInfo["DensityMapFilesNames"].append(DensityMapFile)
        
def RetrieveRecommededContourLevel(Infile):
    """Retrieve recommened contour level."""

    if Infile in OptionsInfo["InfilesRecommededContourLevels"]:
        RecommendedContourLevel = OptionsInfo["InfilesRecommededContourLevels"][Infile]
        return RecommendedContourLevel
    
    RecommendedContourLevel = None
    EMDBID = RetrieveEMDBID(Infile)
    if EMDBID is None:
        MiscUtil.PrintWarning("Failed to retrieve EMDB ID from input file %s to detect local header file already downloaded from EMDB server..." % Infile)
        OptionsInfo["InfilesRecommededContourLevels"][Infile] = RecommendedContourLevel
        return RecommendedContourLevel

    MetadataHeaderFile = "emd-%s.xml" % (EMDBID)
    if not os.path.exists(MetadataHeaderFile):
        MiscUtil.PrintWarning("Failed to find a local header file, %s, for EMDB ID %s..." % (MetadataHeaderFile, EMDBID))
        OptionsInfo["InfilesRecommededContourLevels"][Infile] = RecommendedContourLevel
        return RecommendedContourLevel

    MiscUtil.PrintInfo("\nRetrieving recommeded contour level from header file, %s, for input file, %s..." % (MetadataHeaderFile, Infile))

    ContourLevel = None
    Source = None
    XMLTree = ElementTree.parse(MetadataHeaderFile)
    XMLRoot = XMLTree.getroot()

    MapElement = XMLTree.find("map")
    if MapElement is not None:
        ContourLevelElement = MapElement.find("contourLevel")
        if ContourLevelElement is not None:
            ContourLevel = ContourLevelElement.text
            Source = ContourLevelElement.get("source")

    if ContourLevel is not None:
        if Source is None:
            Source = "NA"
        MiscUtil.PrintInfo("Retrieved recommended (Source: %s) contour level %s..." % (Source, ContourLevel))
        RecommendedContourLevel = ContourLevel
    
    OptionsInfo["InfilesRecommededContourLevels"][Infile] = RecommendedContourLevel
    
    return RecommendedContourLevel

def RetrieveEMDBID(Infile):
    """Retrieve EMDB ID from input file. """

    if Infile in OptionsInfo["InfilesEMDBIDs"]:
        EMDBID = OptionsInfo["InfilesEMDBIDs"][Infile]
        return EMDBID
    
    EMDBID = None
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Infile)

    if re.match("^pdb$", FileExt, re.I):
        EMDBID = RetriveEMDBIDFromPDBFile(Infile)
    elif re.match("^cif$", FileExt, re.I):
        EMDBID = RetriveEMDBIDFromCIFFile(Infile)
    else:
        EMDBID = None

    OptionsInfo["InfilesEMDBIDs"][Infile] = EMDBID
    
    return EMDBID

def RetriveEMDBIDFromPDBFile(Infile):
    """Retrieve EMDB ID from PDB file. """

    EMDBID = None
    InfileFH = open(Infile, "r")
    if InfileFH is None:
        MiscUtil.PrintError("Couldn't open input file: %s.\n" % (Infile))

    MiscUtil.PrintInfo("\nRetrieving EMDB ID from input file %s..." % Infile)
    
    EMDBID = None
    for Line in InfileFH:
        Line = Line.rstrip()
        if re.match("^REMARK", Line, re.I):
            if re.search("DB: EMDB", Line, re.I):
                for Word in Line.split(" "):
                    # Retrieve string with EMD-
                    if re.search("EMD-", Word, re.I):
                        Word = Word.strip()
                        EMDBID = re.sub("EMD-", "", Word)
                        break
                break
    InfileFH.close()
    
    return EMDBID

def RetriveEMDBIDFromCIFFile(Infile):
    """Retrieve EMDB ID from CIF file. """

    InfileFH = open(Infile, "r")
    if InfileFH is None:
        MiscUtil.PrintError("Couldn't open input file: %s.\n" % (Infile))

    MiscUtil.PrintInfo("\nRetrieving EMDB ID from input file %s..." % Infile)
    
    EMDBID = None
    for Line in InfileFH:
        Line = Line.rstrip()
        if re.match("^EMDB  EMD", Line, re.I):
            for Word in Line.split(" "):
                # Retrieve string with EMD-
                if re.search("EMD-", Word, re.I):
                    Word = Word.strip()
                    EMDBID = re.sub("EMD-", "", Word)
                    break
            break
    InfileFH.close()
    
    return EMDBID

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
        ProcessChainMeshesVolumesAndSurfacesOptions(SpecifiedChainsAndLigandsInfo)
        OptionsInfo["InfilesInfo"]["SpecifiedChainsAndLigandsInfo"].append(SpecifiedChainsAndLigandsInfo)
        
        CheckPresenceOfValidLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo)
        
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

def ProcessChainMeshesVolumesAndSurfacesOptions(SpecifiedChainsAndLigandsInfo):
    """Process options to create meshes and surfaces for chains."""

    SpecifiedChainsAndLigandsInfo["VolumeChainComplex"] = {}
    SpecifiedChainsAndLigandsInfo["MeshChainComplex"] = {}
    SpecifiedChainsAndLigandsInfo["SurfaceChainComplex"] = {}

    SpecifiedChainsAndLigandsInfo["EnableVolumeChainComplex"] = {}
    SpecifiedChainsAndLigandsInfo["EnableMeshChainComplex"] = {}
    SpecifiedChainsAndLigandsInfo["EnableSurfaceChainComplex"] = {}
    
    SpecifiedChainsAndLigandsInfo["EnableChainComplexGroup"] = {}
    SpecifiedChainsAndLigandsInfo["EnableChainAloneGroup"] = {}
    
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
        LigandsPresent = True if len(SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID]) else False

        # Create and enable mesh  or volume in auto mode...
        if re.match("^auto$", OptionsInfo["MeshChainComplex"], re.I):
            MeshChainComplex = False if LigandsPresent else True
            EnableMeshChainComplex = False if LigandsPresent else True
        else:
            MeshChainComplex = True if re.match("^Yes$", OptionsInfo["MeshChainComplex"], re.I) else False
            EnableMeshChainComplex = True if re.match("^Yes$", OptionsInfo["MeshChainComplex"], re.I) else False
        
        if re.match("^auto$", OptionsInfo["VolumeChainComplex"], re.I):
            VolumeChainComplex = False if LigandsPresent else True
            EnableVolumeChainComplex = False if LigandsPresent else True
        else:
            VolumeChainComplex = True if re.match("^Yes$", OptionsInfo["VolumeChainComplex"], re.I) else False
            EnableVolumeChainComplex = True if re.match("^Yes$", OptionsInfo["VolumeChainComplex"], re.I) else False
        
        if MeshChainComplex and EnableMeshChainComplex:
            EnableVolumeChainComplex = False
        
        # Create and enable surface in auto mode based on the status of mesh and volume...
        if re.match("^auto$", OptionsInfo["SurfaceChainComplex"], re.I):
            SurfaceChainComplex = False if LigandsPresent else True
            EnableSurfaceChainComplex = False if LigandsPresent else True
            
            if MeshChainComplex or VolumeChainComplex:
                SurfaceChainComplex = False
                EnableSurfaceChainComplex = False
        else:
            SurfaceChainComplex = True if re.match("^Yes$", OptionsInfo["SurfaceChainComplex"], re.I) else False
            EnableSurfaceChainComplex = True if re.match("^Yes$", OptionsInfo["SurfaceChainComplex"], re.I) else False
        
        if (MeshChainComplex and EnableMeshChainComplex) or (VolumeChainComplex or EnableVolumeChainComplex):
            EnableSurfaceChainComplex = False
            
        if LigandsPresent:
            EnableChainComplexGroup = False
            EnableChainAloneGroup = True
        else:
            EnableChainComplexGroup = True
            EnableChainAloneGroup = False

        SpecifiedChainsAndLigandsInfo["VolumeChainComplex"][ChainID] = VolumeChainComplex
        SpecifiedChainsAndLigandsInfo["MeshChainComplex"][ChainID] = MeshChainComplex
        SpecifiedChainsAndLigandsInfo["SurfaceChainComplex"][ChainID] = SurfaceChainComplex
        
        SpecifiedChainsAndLigandsInfo["EnableVolumeChainComplex"][ChainID] = EnableVolumeChainComplex
        SpecifiedChainsAndLigandsInfo["EnableMeshChainComplex"][ChainID] = EnableMeshChainComplex
        SpecifiedChainsAndLigandsInfo["EnableSurfaceChainComplex"][ChainID] = EnableSurfaceChainComplex
        
        SpecifiedChainsAndLigandsInfo["EnableChainComplexGroup"][ChainID] = EnableChainComplexGroup
        SpecifiedChainsAndLigandsInfo["EnableChainAloneGroup"][ChainID] = EnableChainAloneGroup

def ProcessChainSelections():
    """Process custom selections for chains. """

    ChainSelectionsInfo = PyMOLUtil.ProcessChainSelectionsOptionsInfo("--selectionsChain", OptionsInfo["SelectionsChain"])
    OptionsInfo["ChainSelectionsInfo"] = ChainSelectionsInfo
    
    ChainSelections = True if len(OptionsInfo["ChainSelectionsInfo"]["Names"]) else False
    OptionsInfo["ChainSelections"] = ChainSelections
    
def ProcessMeshLevel():
    """Process mesh level."""
    
    MeshLevel = OptionsInfo["MeshLevel"]
    
    if re.match("^auto$", MeshLevel, re.I):
        ProcessAutoMeshLevel()
    else:
        ProcessSpecifiedMeshLevel()

def ProcessAutoMeshLevel():
    """Process auto mesh level."""
    
    OptionsInfo["MeshLevels"] = []
    InfilesNamesCount = len(OptionsInfo["InfilesNames"])
    
    MiscUtil.PrintInfo("\nSetting mesh levels...")
    
    for Index in range(0, InfilesNamesCount):
        Infile = OptionsInfo["InfilesNames"][Index]
        
        RecommededContourLevel = RetrieveRecommededContourLevel(Infile)
        if RecommededContourLevel is None:
            MiscUtil.PrintWarning("Failed to retrieve recommended mesh contour level from header. It's being set to 1.0. Use \"--meshLevel\" option to specify a different contour mesh level.")
            MeshLevel = 1.0
        else:
            MeshLevel = float(RecommededContourLevel)
        OptionsInfo["MeshLevels"].append(MeshLevel)
        
def ProcessSpecifiedMeshLevel():
    """Process specified mesh level."""

    MiscUtil.PrintInfo("\nProcessing mesh levels...")
    
    OptionsInfo["MeshLevels"] = []
    InfilesNamesCount = len(OptionsInfo["InfilesNames"])
    
    MeshLevels = re.sub(" ", "", OptionsInfo["MeshLevel"])
    if not MeshLevels:
        MiscUtil.PrintError("No valid value specified using \"--meshLevel\" option.")
    
    MeshLevelWords = MeshLevels.split(",")
    MeshLevelWordsCount = len(MeshLevelWords)
    
    if MeshLevelWordsCount != InfilesNamesCount:
        MiscUtil.PrintError("The number of comma delimited mesh levels, %d, specified using \"--meshLevel\" must be equal to the number of input files, %s, specified using \"-i, --infiles\" option." % (MeshLevelWordsCount, InfilesNamesCount))
    
    for Index in range(0, InfilesNamesCount):
        Infile = OptionsInfo["InfilesNames"][Index]
        MeshLevel = MeshLevelWords[Index]
        if not MiscUtil.IsFloat(MeshLevel):
            MiscUtil.PrintError("The mesh level, %s, specified using \"--meshLevel\" for input file, %s, must be a float." % (MeshLevel, Infile))
        
        MeshLevel = float(MeshLevel)
        OptionsInfo["MeshLevels"].append(MeshLevel)
    
def ProcessVolumeColorRamp():
    """Process volume color ramp."""
    
    MiscUtil.PrintInfo("\nProcessing volume color ramp...")
    
    ColorRamp = Options["--volumeColorRamp"]
    ColorRampName = None
    CreateColorRamp = False
    ColorRampContourLevel = None
    if re.match("^auto$", ColorRamp, re.I):
        FirstInfile = OptionsInfo["InfilesNames"][0]
        RecommededContourLevel = RetrieveRecommededContourLevel(FirstInfile)
        if RecommededContourLevel is None:
            ColorRampName = "default"
            MiscUtil.PrintWarning("Failed to retrieve recommended contour level from header file corresponding to input file, %s, to create a new ramp. Using PyMOL default volume color ramp." % (FirstInfile))
        else:
            ColorRampName = "CryoEMAuto"
            CreateColorRamp = True
            ColorRampContourLevel = float(RecommededContourLevel)
            MiscUtil.PrintInfo("\nSetting up a  new volume color ramp, %s, at recommended contour level, %.2f, from input file, %s..." % (ColorRampName, ColorRampContourLevel, FirstInfile))
    else:
        ColorRampName = ColorRamp
        
    OptionsInfo["VolumeColorRamp"] = ColorRamp
    OptionsInfo["VolumeColorRampName"] = ColorRampName
    OptionsInfo["VolumeColorRampCreate"] = CreateColorRamp
    OptionsInfo["VolumeColorRampContourLevel"] = ColorRampContourLevel

def ProcessOptions():
    """Process and validate command line arguments and options."""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Align"] = True if re.match("^Yes$", Options["--align"], re.I) else False
    OptionsInfo["AlignMethod"] = Options["--alignMethod"].lower()
    OptionsInfo["AlignMode"] = Options["--alignMode"]
    
    OptionsInfo["AllowEmptyObjects"] = True if re.match("^Yes$", Options["--allowEmptyObjects"], re.I) else False
    
    OptionsInfo["BFactorChainCartoonPutty"] = True if re.match("^Yes$", Options["--BFactorChainCartoonPutty"], re.I) else False
    OptionsInfo["BFactorColorPalette"] = Options["--BFactorColorPalette"]
    
    OptionsInfo["Infiles"] = Options["--infiles"]
    OptionsInfo["InfilesNames"] =  Options["--infileNames"]
    
    OptionsInfo["InfilesEMDBIDs"] =  {}
    OptionsInfo["InfilesRecommededContourLevels"] =  {}

    OptionsInfo["AlignRefFile"] = Options["--alignRefFile"]
    if re.match("^FirstInputFile$", Options["--alignRefFile"], re.I):
        OptionsInfo["RefFileName"] = OptionsInfo["InfilesNames"][0]
    else:
        OptionsInfo["RefFileName"] = Options["--alignRefFile"]
    
    OptionsInfo["DensityMapFiles"] = Options["--densityMapFiles"]
    ProcessDensityMapFiles()

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
    
    # Process mesh parameters...
    OptionsInfo["MeshCarveRadius"] = float(Options["--meshCarveRadius"])
    OptionsInfo["MeshComplex"] = True if re.match("^Yes$", Options["--meshComplex"], re.I) else False
    OptionsInfo["MeshChainComplex"] = Options["--meshChainComplex"]
    
    OptionsInfo["MeshWidth"] = float(Options["--meshWidth"])
    OptionsInfo["MeshColor"] = Options["--meshColor"]
    
    OptionsInfo["MeshLevel"] = Options["--meshLevel"]
    ProcessMeshLevel()
    
    OptionsInfo["SurfaceComplex"] = True if re.match("^Yes$", Options["--surfaceComplex"], re.I) else False
    OptionsInfo["SurfaceChainComplex"] = Options["--surfaceChainComplex"]
    OptionsInfo["SurfaceTransparency"] = float(Options["--surfaceTransparency"])
    
    OptionsInfo["PocketContactsLigandColor"] = Options["--pocketContactsLigandColor"]
    OptionsInfo["PocketContactsLigandHydrophobicColor"] = Options["--pocketContactsLigandHydrophobicColor"]
    OptionsInfo["PocketContactsSolventColor"] = Options["--pocketContactsSolventColor"]
    OptionsInfo["PocketContactsInorganicColor"] = Options["--pocketContactsInorganicColor"]
    
    OptionsInfo["PocketContactsCutoff"] = float(Options["--pocketContactsCutoff"])
    OptionsInfo["PocketDistanceCutoff"] = float(Options["--pocketDistanceCutoff"])
    
    OptionsInfo["PocketLabelColor"] = Options["--pocketLabelColor"]
    OptionsInfo["PocketSurface"] = True if re.match("^Yes$", Options["--pocketSurface"], re.I) else False

    OptionsInfo["SelectionsChain"] = Options["--selectionsChain"]
    OptionsInfo["SelectionsChainStyle"] = Options["--selectionsChainStyle"]
    ProcessChainSelections()
    
    OptionsInfo["VolumeCarveRadius"] = float(Options["--volumeCarveRadius"])
    OptionsInfo["VolumeComplex"] = True if re.match("^Yes$", Options["--volumeComplex"], re.I) else False
    OptionsInfo["VolumeChainComplex"] = Options["--volumeChainComplex"]

    ProcessVolumeColorRamp()

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
    
    MiscUtil.ValidateOptionTextValue("--BFactorChainCartoonPutty", Options["--BFactorChainCartoonPutty"], "yes no")

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

    MiscUtil.ValidateOptionFloatValue("--meshCarveRadius", Options["--meshCarveRadius"], {">": 0.0})
    MiscUtil.ValidateOptionTextValue("--meshComplex", Options["--meshComplex"], "yes no")
    MiscUtil.ValidateOptionTextValue("--meshChainComplex", Options["--meshChainComplex"], "yes no auto")
    MiscUtil.ValidateOptionFloatValue("--meshWidth", Options["--meshWidth"], {">": 0.0})
    
    MiscUtil.ValidateOptionTextValue("--PMLOut", Options["--PMLOut"], "yes no")
    
    MiscUtil.ValidateOptionFloatValue("--pocketContactsCutoff", Options["--pocketContactsCutoff"], {">": 0.0})
    MiscUtil.ValidateOptionFloatValue("--pocketDistanceCutoff", Options["--pocketDistanceCutoff"], {">": 0.0})
    if (float(Options["--pocketContactsCutoff"]) > float(Options["--pocketDistanceCutoff"])):
        MiscUtil.PrintError("The value, %s, specified using option \"--pocketContactsCutoff\" must be less than value, %s, specified using \"-pocketDistanceCutoff\" option." % (Options["--pocketContactsCutoff"], Options["--pocketDistanceCutoff"]))
        
    MiscUtil.ValidateOptionTextValue("--pocketSurface", Options["--pocketSurface"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--surfaceComplex", Options["--surfaceComplex"], "yes no")
    MiscUtil.ValidateOptionTextValue("--surfaceChainComplex", Options["--surfaceChainComplex"], "yes no auto")
    MiscUtil.ValidateOptionFloatValue("--surfaceTransparency", Options["--surfaceTransparency"], {">=": 0.0, "<=": 1.0})
    
    MiscUtil.ValidateOptionFloatValue("--volumeCarveRadius", Options["--volumeCarveRadius"], {">": 0.0})
    MiscUtil.ValidateOptionTextValue("--volumeComplex", Options["--volumeComplex"], "yes no")
    MiscUtil.ValidateOptionTextValue("--volumeChainComplex", Options["--volumeChainComplex"], "yes no auto")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLVisualizeCryoEMDensity.py - Visualize cryo-EM density

Usage:
    PyMOLVisualizeCryoEMDensity.py  [--align <yes or no>] [--alignMethod <align, cealign, super>]
                                   [--alignMode <FirstChain or Complex>] [--alignRefFile <filename>]
                                   [--allowEmptyObjects <yes or no>] [--BFactorChainCartoonPutty <yes or no>]
                                   [--BFactorColorPalette <text> ] [--chainIDs <First, All or ID1,ID2...>]
                                   [--densityMapFiles <file1,file2,file3,...>]
                                   [--ligandIDs <Largest, All or ID1,ID2...>] [--labelFontID <number>]
                                   [--meshCarveRadius <number>] [--meshComplex <yes or no>]
                                   [--meshChainComplex <yes, no, or auto>] [--meshColor <text>]
                                   [--meshLevel <number>] [--meshWidth <number>] [--PMLOut <yes or no>]
                                   [--pocketContactsLigandColor <text>] [--pocketContactsLigandHydrophobicColor <text>]
                                   [--pocketContactsSolventColor <text>]  [--pocketContactsCutoff <number>]
                                   [--pocketContactsInorganicColor <text>] [--pocketDistanceCutoff <number>]
                                   [--pocketLabelColor <text>] [--pocketSurface <yes or no>]
                                   [--selectionsChain <ObjectName,SelectionSpec,...>] [--selectionsChainStyle <DisplayStyle>]
                                   [--surfaceComplex <yes or no>] [--surfaceChainComplex <yes, no or auto>]
                                   [--surfaceTransparency <number>] [--volumeCarveRadius <number>]
                                   [--volumeComplex <yes or no>] [--volumeChainComplex <yes, no, or auto>]
                                   [--volumeColorRamp <text>]   [--overwrite] [-w <dir>] -i <infile1,infile2,...> -o <outfile>
    PyMOLVisualizeCryoEMDensity.py -h | --help | -e | --examples

Description:
    Generate PyMOL visualization files for viewing electron microscopy (EM) or
    cryo-EM density around chains, ligands, and ligand binding pockets in
    macromolecules including proteins and nucleic acids.

    The supported input file formats are: Macromolecule - PDB (.pdb) or CIF(.cif),
    Cryo-EM Density - Collaborative Computational Project Number 4 (CCP4) ( .map)

    The supported output file formats are: PyMOL script file (.pml), PyMOL session
    file (.pse)

    The cryo-EM density and header files along with PDB files may be downloaded
    from appropriate servers using DownloadPDBFiles.pl script.

    A variety of PyMOL groups and objects may be  created for visualization of
    cryo-EM density present in map files. These groups and objects correspond to
    maps, volumes, meshes, surfaces,chains, ligands, inorganics, ligand binding
    pockets, polar interactions, and pocket hydrophobic surfaces. A complete
    hierarchy of all possible PyMOL groups and objects is shown below:
    
        <PDBFileRoot>
            .Complex
                .Complex
                .CryoEM
                    .Map
                    .Volume
                    .Mesh
                    .Surface
            .Chain<ID>
                .Complex
                    .Complex
                    .CryoEM
                        .Volume
                        .Mesh
                        .Surface
                .Chain
                    .Chain
                    .BFactor
                    .Selections
                        .<Name1>
                            .Selection
                            .CryoEM
                                .Volume
                                .Mesh
                                .Surface
                        .<Name2>
                            ... ... ..
                .Solvent
                .Inorganic
                .Ligand<ID>
                    .Ligand
                        .Ligand
                        .CryoEM
                            .Volume
                            .Mesh
                            .Surface
                    .Pocket
                        .Pocket
                        .CryoEM
                            .Volume
                            .Mesh
                            .Surface
                        .Polar_Contacts
                        .Hydrophobic_Contacts
                        .Surface
                    .Pocket_Solvent
                        .Pocket_Solvent
                        .CryoEM
                            .Volume
                            .Mesh
                            .Surface
                        .Polar_Contacts
                    .Pocket_Inorganic
                        .Pocket_Inorganic
                        .CryoEM
                            .Volume
                            .Mesh
                            .Surface
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
    
    The meshes, volumes, and surfaces  are not created for complete complex in input
    files by default. A word to the wise: The creation of these mesh, volume, and surface
    objects may slow down loading of PML file and generation of PSE file, based on the
    size of input complex and map files. The generation of PSE file may also fail. In 
    addition, you may want to interactively manipulate the contour level for meshes,
    volumes, and surfaces. The recommended value for contour level is automatically
    retrieved from header files available from EM density server. The recommended
    value may not always work.

Options:
    -a, --align <yes or no>  [default: no]
        Align input files to a reference file before visualization along with
        available cryo-EM density map files.
    --alignMethod <align, cealign, super>  [default: super]
        Alignment methodology to use for aligning input files to a reference
        file.
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
        inorganic atom selections across chains, ligands, and ligand binding pockets
        in input file(s).
    -b, --BFactorChainCartoonPutty <yes or no>  [default: yes]
        A cartoon putty around individual chains colored by B factors. The minimum
        and maximum values for B factors are automatically detected. These values
        indicate spread of cryo-EM density around atoms. The 'blue_white_red' color
        palette is deployed for coloring the cartoon putty.
    --BFactorColorPalette <text>  [default: blue_white_red]
        Color palette for coloring cartoon putty around chains generated using B
        factors. Any valid PyMOL color palette name is allowed. No validation is
        performed. The complete list of valid color palette names is a available
        at: pymolwiki.org/index.php/Spectrum. Examples: blue_white_red,
        blue_white_magenta, blue_red, green_white_red, green_red.
    -c, --chainIDs <First, All or ID1,ID2...>  [default: First]
        List of chain IDs to use for visualizing cryo-EM density. Possible values:
        First, All, or a comma delimited list of chain IDs. The default is to use the
        chain ID for the first chain in each input file.
    -d, --densityMapFiles <file1,file2,file3,...>  [default: auto]
        CryoEM density map file names. The EMDB ID is retrieved from PDB and CIF
        file to set the cryo-EM density file name during automatic detection of
        density files. The format of the file name is as follows:
        
            emd_<EMDBID>.map.gz or emd_<EMDBID>.map
         
        The density files must be present in the working directory.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infiles <infile1,infile2,infile3...>
        Input file names.
    -l, --ligandIDs <Largest, All or ID1,ID2...>  [default: Largest]
        List of ligand IDs present in chains for visualizing cryo-EM density across
        ligands and ligand binding pockets. Possible values: Largest, All, or a comma
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
    --meshCarveRadius <number>  [default: 1.6]
        Radius in Angstroms around atoms for including cryo-EM density.
    --meshComplex <yes or no>  [default: no]
        Create meshes for complete complex in each input file using corresponding
        density map file.
    --meshChainComplex <yes, no, or auto>  [default: auto]
        Create meshes for individual chain complex in each input file using
        corresponding density map file. By default, the meshes are automatically
        created for chain complexes without any ligands. 
    --meshColor <text>  [default: blue]
        Line color for meshes corresponding to density maps.. The specified value
        must be valid color. No validation is performed.
    --meshLevel <number1,number2,...>  [default: auto]
        Comma delimited list of contour levels in sigma units for generating meshes
        for each input file using corresponding density map file. The default is to
        automatically retrieve the recommended contour levels for each input 
        file. The header file emd-<EMDBID>.xml corresponding to an input file
        must be present in the working directory  to automatically retrieve
        recommended value for mesh contour level. Otherwise, the default contour
        level is set to 1.
        
        You may want to interactively manipulate the contour level for meshes and
        surfaces. The default recommended value may not always work.
    --meshWidth <number>  [default: 0.5]
        Line width for mesh lines corresponding to density maps.
    -o, --outfile <outfile>
        Output file name.
    -p, --PMLOut <yes or no>  [default: yes]
        Save PML file during generation of PSE file.
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
        Color for drawing polar contacts between solvent and pocket residues.
        The specified value must be valid color. No validation is performed.
    --pocketContactsInorganicColor <text>  [default: deepsalmon]
        Color for drawing polar contacts between inorganic and pocket residues.
        The specified value must be valid color. No validation is performed.
    --pocketContactsCutoff <number>  [default: 4.0]
        Distance in Angstroms for identifying polar and hyrdophobic contacts
        between atoms in pocket residues and ligands.
    --pocketDistanceCutoff <number>  [default: 5.0]
        Distance in Angstroms for identifying pocket residues around ligands.
    --pocketLabelColor <text>  [default: magenta]
        Color for drawing residue or atom level labels for a pocket. The specified
        value must be valid color. No validation is performed.
    --pocketSurface <yes or no>  [default: yes]
        Hydrophobic surface around pocket. The pocket surface is colored by
        hydrophobicity. It is only valid for proteins. The color of amino acids is
        set using the Eisenberg hydrophobicity scale. The color varies from red
        to white, red being the most hydrophobic amino acid.
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
            
    --selectionsChainStyle <DisplayStyle>  [default: sticks]
        Display style for PyMOL objects created for '--selectionsChain' option. It
        must be a valid PyMOL display style. No validation is performed.
    --surfaceComplex <yes or no>  [default: no]
        Create surfaces for complete complex in input file(s) corresponding to density
        map.
    --surfaceChainComplex <yes, no or auto>  [default: auto]
        Create surfaces for individual chain complexes in each input file using corresponding
        density map file. By default, the surfaces are automatically created for chain complexes
        without any ligands.
    --surfaceTransparency <number>  [default: 0.25]
        Surface transparency for molecular and cryo-EM density surfaces.
    --overwrite
        Overwrite existing files.
    --volumeCarveRadius <number>  [default: 1.6]
        Radius in Angstroms around atoms for including cryo-EM density.
    --volumeComplex <yes or no>  [default: no]
        Create volumes for complete complex in each input file using corresponding density
        map file.
    --volumeChainComplex <yes, no, or auto>  [default: auto]
        Create volumes for individual chain complex in each input file using corresponding
        density map file. By default, the volumes are automatically created for chain
        complexes without any ligands.
    --volumeColorRamp <text>  [default: auto]
        Name of a volume color ramp for density map files. The specified value must
        be a valid name. No validation is performed. The following volume color ramps
        are currently available in PyMOL: default, 2fofc, fofc, rainbow, and rainbow2.
        
        The default is to automatically create a new volume color ramp for the first
        input file using recommended contour level with an offset of 0.3 around this value.
        The header file emd-<EMDBID>.xml must be present in the working directory  to
        automatically retrieve recommended contour level and generate a  volume color ramp.
        Otherwise, PyMOL default volume color ramp is employed to color volumes.
        
        The volume color ramp automatically created for the first input file is used for all
        other input files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To download structure and cryo-EM data for 5K12, 5UMD, 5W81, and 5UAK
    before running the following examples, type:

        % DownloadPDBFiles.pl --DensityMap yes 5K12,5UMD,5W81,5UAK

    To visualize cryo-EM density at recommended contour level for the first
    chain complex in a PDB file using corresponding density map and header
    file, and generate a PML file type:

        % PyMOLVisualizeCryoEMDensity.py -i 5K12.pdb -o 5K12.pml

    To visualize cryo-EM density at recommended contour level for the first
    chain complex in a PDB file and highlighting densities for all cysteines and
    serines  using corresponding density map and header file, and generate
    a PML file type:

        % PyMOLVisualizeCryoEMDensity.py -i 5K12.pdb -o 5K12.pml
          --selectionsChain "Csysteines,resn cys,Serines,resn ser"

    To visualize electron density for the largest ligand in  chain K, and ligand
    binding pocket to highlight ligand interactions with pockect residues,
    solvents and inorganics, in a PDB and using corresponding map files, and
    generate a PML file, type:

        % PyMOLVisualizeCryoEMDensity.py -c K -i 5UMD.cif -o 5UMD.pml

    To visualize cryo-EM density for all  chains along with any solvents in a
    PDB file and using corresponding map files, and generate a PML file, type:

        % PyMOLVisualizeCryoEMDensity.py -c all -i 5K12.pdb -o 5K12.pml

    To visualize cryo-EM density at a specific contour level for the first chain
    complex along with volume and surface in a PDB file using corresponding
    to a specific density map file, and generate a PML file, type:

        % PyMOLVisualizeCryoEMDensity.py -d emd_8194.map.gz --meshLevel 1.0
          --surfaceChainComplex yes --volumeChainComplex yes -i 5K12.pdb
          -o 5K12.pml

    To align and visualize cryo-EM density at recommended contour levels for the
    largest ligand in the first chain along with pockets or the first chain complex
    in input files using corresponding maps and header files, type:

        % PyMOLVisualizeCryoEMDensity.py -a yes -i "5W81.pdb,5UAK.pdb"
          -o SampleOut.pml

    To align and visualize cryo-EM density at recommended contour levels for all
    chains and ligands in input files using specified density files, type:
    in input files using corresponding maps and header files, type:

        % PyMOLVisualizeCryoEMDensity.py -a yes -i "5W81.pdb,5UAK.pdb"
          -o SampleOut.pml -c all -l all -d "emd_8782.map.gz,emd_8516.map.gz"

Author:
    Manish Sud(msud@san.rr.com)

See also:
    DownloadPDBFiles.pl, PyMOLVisualizeCavities.py,
    PyMOLVisualizeElectronDensity.py, PyMOLVisualizeInterfaces.py,
    PyMOLVisualizeMacromolecules.py, PyMOLVisualizeSurfaceAndBuriedResidues.py

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
