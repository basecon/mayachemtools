#!/usr/bin/env python
#
# File: PyMOLInfoMacromolecules.py
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
    ListInfo()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def ListInfo():
    """List information for macromolecules."""

    for Infile in OptionsInfo["InfilesNames"]:
        MiscUtil.PrintInfo("\nProcessing file %s..." % Infile)
        ListFileInfo(Infile)

def ListFileInfo(Infile):
    """List information for macromolecules in a file."""
    
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Infile)
    MolName = FileName

    # Load infile...
    pymol.cmd.load(Infile, MolName)
    
    ChainIDs = PyMOLUtil.GetChains(MolName)
    ListHeaderInfo(Infile)
    
    ListChainsInfo(MolName, ChainIDs)
    ListChainsResiduesInfo(MolName, ChainIDs)
    
    ListLigandsInfo(MolName, ChainIDs)
    ListSolventsInfo(MolName, ChainIDs)
    ListInorganicsInfo(MolName, ChainIDs)

    ListPocketsInfo(MolName, ChainIDs)
    
    ListInterfaceResiduesInfo(MolName, ChainIDs)
    ListSurfaceResiduesInfo(MolName, ChainIDs)
    ListPhiPsiAnglesInfo(MolName, ChainIDs)
    
    ListBoundingBoxInfo(MolName)
    
    # Delete infile object...
    pymol.cmd.delete(MolName)
    
    ListFileSizeAndModificationInfo(Infile)

def ListHeaderInfo(Infile):
    """List header information."""

    if not (OptionsInfo["All"] or OptionsInfo["Header"]):
        return

    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Infile)

    Classification, DepositionDate, IDCode, ExperimentalTechnique, Resolution = ["Not Available"] * 5
    if re.match("^pdb$", FileExt, re.I):
         Classification, DepositionDate, IDCode, ExperimentalTechnique, Resolution = RetriveHeadAndExperimentalInfoFromPDBFile(Infile)
    elif re.match("^cif$", FileExt, re.I):
         Classification, DepositionDate, IDCode, ExperimentalTechnique, Resolution = RetriveHeadAndExperimentalInfoFromCIFFile(Infile)

    MiscUtil.PrintInfo("\nID: %s\nClassification: %s\nDeposition date: %s" % (IDCode, Classification, DepositionDate))
    MiscUtil.PrintInfo("\nExperimental technique: %s\nResolution: %s" % (ExperimentalTechnique, Resolution))

def RetriveHeadAndExperimentalInfoFromPDBFile(Infile):
    """Retrieve header and experimental information from PDB file. """

    Classification, DepositionDate, IDCode, ExperimentalTechnique, Resolution = ["Not Available"] * 5

    Lines = MiscUtil.GetTextLines(Infile)

    # Retrieve header info...
    for Line in Lines:
        if re.match("^HEADER", Line, re.I):
            # Format: 10x40s9s3x4s
            FormatSize = 66
            Line = PrepareLineForFormatSize(Line, FormatSize)
            
            Classification = Line[10:50]
            DepositionDate = Line[50:59]
            IDCode = Line[62:66]
            
            Classification = Classification.strip() if len(Classification.strip()) else "Not Available"
            DepositionDate = DepositionDate.strip() if len(DepositionDate.strip()) else "Not Available"
            IDCode = IDCode.strip() if len(IDCode.strip()) else "Not Available"
            
            break

    # Retrieve experimental info...
    for Line in Lines:
        if re.match("^EXPDTA", Line, re.I):
            ExperimentalTechnique = re.sub("^EXPDTA", "", Line, re.I)
            ExperimentalTechnique = ExperimentalTechnique.strip()
        elif re.match("^REMARK   2 RESOLUTION.", Line, re.I):
            if re.search("NOT APPLICABLE", Line, re.I):
                Resolution = "NOT APPLICABLE"
            else:
                FormatSize = 70
                Line = PrepareLineForFormatSize(Line, FormatSize)
                Resolution = Line[22:70]
                Resolution = Resolution.strip() if len(Resolution.strip()) else "Not Available"
        elif re.match("^(ATOM|HETATM)", Line, re.I):
            break

    return Classification, DepositionDate, IDCode, ExperimentalTechnique, Resolution

def RetriveHeadAndExperimentalInfoFromCIFFile(Infile):
    """Retrieve header information from CIF file. """

    Classification, DepositionDate, IDCode, ExperimentalTechnique, Resolution = ["Not Available"] * 5
    
    Lines = MiscUtil.GetTextLines(Infile)
    
    # IDCode...
    for Line in Lines:
        if re.match("^_struct_keywords.entry_id", Line, re.I):
            IDCode = re.sub("^_struct_keywords.entry_id", "", Line, re.I)
            IDCode = IDCode.strip() if len(IDCode.strip()) else "Not Available"
            break

    # Classification...
    for Line in Lines:
        if re.match("^_struct_keywords.pdbx_keywords", Line, re.I):
            Classification = re.sub("^_struct_keywords.pdbx_keywords", "", Line, re.I)
            Classification = Classification.strip() if len(Classification.strip()) else "Not Available"
            break
    
    # Deposition date...
    for Line in Lines:
        if re.match("^_pdbx_database_status.recvd_initial_deposition_date", Line, re.I):
            DepositionDate = re.sub("^_pdbx_database_status.recvd_initial_deposition_date", "", Line, re.I)
            DepositionDate = DepositionDate.strip() if len(DepositionDate.strip()) else "Not Available"
            break
    
    # Experimental technique...
    for Line in Lines:
        if re.match("^_exptl.method", Line, re.I):
            ExperimentalTechnique = re.sub("(_exptl.method|')", "", Line, flags = re.I)
            ExperimentalTechnique = ExperimentalTechnique.strip() if len(ExperimentalTechnique.strip()) else "Not Available"
            break
        
    # Resolution...
    for Line in Lines:
        if re.match("^_reflns.d_resolution_high", Line, re.I):
            Resolution = re.sub("^_reflns.d_resolution_high", "", Line, re.I)
            Resolution = Resolution.strip() if len(Resolution.strip()) else "Not Available"
            break
        
    return Classification, DepositionDate, IDCode, ExperimentalTechnique, Resolution

def PrepareLineForFormatSize(Text, FormatSize):
    """Prepare text string for format size padding or truncation to
    alter its size.
    """

    TextLen = len(Text)
    if TextLen < FormatSize:
        PaddingLen = FormatSize - TextLen
        TextPadding = " " * PaddingLen
        Text = Text + TextPadding
    elif TextLen > FormatSize:
        Text = Text[:FormatSize]
    
    return Text

def ListChainsInfo(MolName, ChainIDs):
    """List chains information across all chains."""

    if not (OptionsInfo["All"] or OptionsInfo["Chains"]):
        return

    ChainsInfo = ", ".join(ChainIDs) if len(ChainIDs) else "None"

    MiscUtil.PrintInfo("\nNumber of chains: %s" % len(ChainIDs))
    MiscUtil.PrintInfo("ChainIDs: %s" % ChainsInfo)

def ListChainsResiduesInfo(MolName, ChainIDs):
    """List polymer chain residue information across all chains."""
    
    if not (OptionsInfo["All"] or OptionsInfo["CountResidues"]):
        return
    
    ListSelectionResiduesInfo(MolName, ChainIDs, "Chains")

    # List information for non-standard amino acids...
    ListSelectionResiduesInfo(MolName, ChainIDs, "NonStandardAminoAcids")

def ListLigandsInfo(MolName, ChainIDs):
    """List ligand information across all chains."""
    
    if not (OptionsInfo["All"] or OptionsInfo["Ligands"]):
        return
    
    ListSelectionResiduesInfo(MolName, ChainIDs, "Ligands")

def ListSolventsInfo(MolName, ChainIDs):
    """List solvents information across all chains."""
    
    if not (OptionsInfo["All"] or OptionsInfo["Solvents"]):
        return

    ListSelectionResiduesInfo(MolName, ChainIDs, "Solvents")
    
def ListInorganicsInfo(MolName, ChainIDs):
    """List inorganics information across all chains."""
    
    if not (OptionsInfo["All"] or OptionsInfo["Inorganics"]):
        return
    
    ListSelectionResiduesInfo(MolName, ChainIDs, "Inorganics")

def ListSelectionResiduesInfo(MolName, ChainIDs, SelectionType):
    """List residues information for a specified selection type. """

    Lines = []
    TotalResCount = 0
    SelectionLabel = None
    
    if not len(ChainIDs):
        SelectionInfo, SelectionLabel = GetSelectionResiduesInfo(MolName, None, SelectionType)
        MiscUtil.PrintInfo("\nNumber of %s residues: %s" % (SelectionLabel, TotalResCount))
        return
        
    for ChainID in ChainIDs:
        SelectionInfo, SelectionLabel = GetSelectionResiduesInfo(MolName, ChainID, SelectionType)

        ChainResCount = 0
        LineWords = []
        
        SortedResNames = sorted(SelectionInfo["ResNames"], key = lambda ResName: SelectionInfo["ResCount"][ResName], reverse = True)
        for ResName in SortedResNames:
            ResCount = SelectionInfo["ResCount"][ResName]
            LineWords.append("%s - %s" % (ResName, ResCount))

            ChainResCount += ResCount
            TotalResCount += ResCount
        
        Line = "; ".join(LineWords) if len(LineWords) else None
        Lines.append("Chain ID: %s; Count: %s; Names:  %s" % (ChainID, ChainResCount, Line))
    
    MiscUtil.PrintInfo("\nNumber of %s residues: %s" % (SelectionLabel, TotalResCount))
    for Line in Lines:
        MiscUtil.PrintInfo("%s" % Line)

def GetSelectionResiduesInfo(MolName, ChainID, SelectionType):
    """Get residues info for a specified selection type. """

    SelectionInfo = None
    SelectionLabel = None
    
    if re.match("^Ligands$", SelectionType, re.I):
        SelectionLabel = "ligand"
        SelectionInfo = PyMOLUtil.GetLigandResiduesInfo(MolName, ChainID) if ChainID is not None else None
    elif re.match("^Solvents$", SelectionType, re.I):
        SelectionLabel = "solvent"
        SelectionInfo = PyMOLUtil.GetSolventResiduesInfo(MolName, ChainID) if ChainID is not None else None
    elif re.match("^Inorganics$", SelectionType, re.I):
        SelectionLabel = "inorganic"
        SelectionInfo = PyMOLUtil.GetInorganicResiduesInfo(MolName, ChainID) if ChainID is not None else None
    elif re.match("^Chains$", SelectionType, re.I):
        SelectionLabel = "polymer chain"
        SelectionInfo = PyMOLUtil.GetPolymerResiduesInfo(MolName, ChainID) if ChainID is not None else None
    elif re.match("^NonStandardAminoAcids$", SelectionType, re.I):
        SelectionLabel = "non-standard amino acids"
        SelectionInfo = PyMOLUtil.GetAminoAcidResiduesInfo(MolName, ChainID, "NonStandard") if ChainID is not None else None
    else:
        MiscUtil.PrintError("Failed to retrieve residues information: Selection type %s is not valid..." % SelectionType)

    return SelectionInfo, SelectionLabel

def ListPocketsInfo(MolName, ChainIDs):
    """List pockect residues information across all chains."""
    
    if not (OptionsInfo["All"] or OptionsInfo["PocketLigands"] or OptionsInfo["PocketSolvents"] or OptionsInfo["PocketInorganics"]) :
        return
    
    for ChainID in ChainIDs:
        MiscUtil.PrintInfo("\nListing ligand pockets information for chain %s..." % (ChainID))
        
        LigandsInfo = PyMOLUtil.GetLigandResiduesInfo(MolName, ChainID)
        if not len(LigandsInfo["ResNames"]):
            MiscUtil.PrintInfo("\nNumber of residues in ligand pocket: None")
            MiscUtil.PrintInfo("Chain ID: %s; Ligands: None" % (ChainID))
            continue
        
        for LigandResName in sorted(LigandsInfo["ResNames"]):
            for LigandResNum in LigandsInfo["ResNum"][LigandResName]:
                ListPocketsPolymerInfo(MolName, ChainID, LigandResName, LigandResNum)
                ListPocketsSolventsInfo(MolName, ChainID, LigandResName, LigandResNum)
                ListPocketsInorganicsInfo(MolName, ChainID, LigandResName, LigandResNum)

def ListPocketsPolymerInfo(MolName, ChainID, LigandResName, LigandResNum):
    """List pockect residues information across all chains."""

    if not (OptionsInfo["All"] or OptionsInfo["PocketLigands"]):
        return
    
    ListPocketSelectionResiduesInfo(MolName, ChainID, LigandResName, LigandResNum, "Pockets")

def ListPocketsSolventsInfo(MolName, ChainID, LigandResName, LigandResNum):
    """List pockect solvent residues information across all chains."""
    
    if not (OptionsInfo["All"] or OptionsInfo["PocketSolvents"]):
        return
    
    ListPocketSelectionResiduesInfo(MolName, ChainID, LigandResName, LigandResNum, "PocketSolvents")

def ListPocketsInorganicsInfo(MolName, ChainID, LigandResName, LigandResNum):
    """List pockect inorganic residues information across all chains."""
    
    if not (OptionsInfo["All"] or OptionsInfo["PocketInorganics"]):
        return
    
    ListPocketSelectionResiduesInfo(MolName, ChainID, LigandResName, LigandResNum, "PocketInorganics")

def ListPocketSelectionResiduesInfo(MolName, ChainID, LigandResName, LigandResNum, SelectionType):
    """List residues information for a specified pocket selection type. """

    PocketDistanceCutoff  = OptionsInfo["PocketDistanceCutoff"]
    SelectionLabel = GetPocketSelectionResiduesInfoLabel(SelectionType)
    
    SelectionInfo = GetPocketSelectionResiduesInfo(MolName, ChainID, LigandResName, LigandResNum, PocketDistanceCutoff, SelectionType)

    ResiduesCount, ResiduesDistribution, ResiduesIDs = SetupSelectionResiduesInfo(SelectionInfo)
    MiscUtil.PrintInfo("\nNumber of %s residues in ligand pocket: %s" % (SelectionLabel, ResiduesCount))
    MiscUtil.PrintInfo("Chain ID: %s; Ligand ID: %s_%s\nResidue distribution: %s\nResidue IDs: %s" % (ChainID, LigandResName, LigandResNum, ResiduesDistribution, ResiduesIDs))

def GetPocketSelectionResiduesInfo(MolName, ChainID, LigandResName, LigandResNum, PocketDistanceCutoff, SelectionType):
    """Get pocket residues info for a specified selection type. """

    SelectionInfo = None
    
    if re.match("^Pockets$", SelectionType, re.I):
        SelectionInfo = PyMOLUtil.GetPocketPolymerResiduesInfo(MolName, ChainID, LigandResName, LigandResNum, PocketDistanceCutoff)
    elif re.match("^PocketSolvents$", SelectionType, re.I):
        SelectionInfo = PyMOLUtil.GetPocketSolventResiduesInfo(MolName, ChainID, LigandResName, LigandResNum, PocketDistanceCutoff)
    elif re.match("^PocketInorganics$", SelectionType, re.I):
        SelectionInfo = PyMOLUtil.GetPocketInorganicResiduesInfo(MolName, ChainID, LigandResName, LigandResNum, PocketDistanceCutoff)
    else:
        MiscUtil.PrintError("Failed to retrieve pocket residues information: Selection type %s is not valid..." % SelectionType)

    return SelectionInfo

def ListSurfaceResiduesInfo(MolName, ChainIDs):
    """List surface and buried residues for polymer chains with in a molecule."""
    
    if not (OptionsInfo["All"] or OptionsInfo["SurfaceResidues"]):
        return

    MiscUtil.PrintInfo("\nListing surface and buried residues information...")
    if not len(ChainIDs):
        MiscUtil.PrintInfo("\nNumber of surface residues: None\nNumber of buried residues: None")
        return

    TotalSurfaceResidues, TotalBuriedResidues = [0] * 2
    
    for ChainID in ChainIDs:
        SurfaceResiduesSelectionInfo, BuriedResiduesSelectionInfo = PyMOLUtil.GetSurfaceAndBuriedResiduesInfo(MolName, ChainID, OptionsInfo["SurfaceResiduesCutoff"])
        
        SurfaceResiduesCount, SurfaceResiduesDistribution, SurfaceResiduesIDs = SetupSelectionResiduesInfo(SurfaceResiduesSelectionInfo)
        BuriedResiduesCount, BuriedResiduesDistribution, BuriedResiduesIDs = SetupSelectionResiduesInfo(BuriedResiduesSelectionInfo)

        TotalSurfaceResidues += SurfaceResiduesCount
        MiscUtil.PrintInfo("\nChainID: %s; Number of surface residues: %d" % (ChainID, SurfaceResiduesCount))
        MiscUtil.PrintInfo("Residue distribution: %s" % (SurfaceResiduesDistribution))
        if OptionsInfo["SurfaceResiduesIDs"]:
            MiscUtil.PrintInfo("Residue IDs: %s" % (SurfaceResiduesIDs))
        
        TotalBuriedResidues += BuriedResiduesCount
        MiscUtil.PrintInfo("\nChainID: %s; Number of buried residues: %d" % (ChainID, BuriedResiduesCount))
        MiscUtil.PrintInfo("Residue distribution: %s" % (BuriedResiduesDistribution))
        if OptionsInfo["SurfaceResiduesIDs"]:
            MiscUtil.PrintInfo("Residue IDs: %s" % (BuriedResiduesIDs))
    
    MiscUtil.PrintInfo("\nTotal number of surface residues: %d\nTotal number of buried residues: %s" % (TotalSurfaceResidues, TotalBuriedResidues))
        
def ListPhiPsiAnglesInfo(MolName, ChainIDs):
    """List phi and psi torsion angles for polymer chains with in a molecule."""
    
    if not (OptionsInfo["All"] or OptionsInfo["PhiPsi"]):
        return

    MiscUtil.PrintInfo("\nListing phi and psi angles information...")
    
    if not len(ChainIDs):
        MiscUtil.PrintInfo("\nNumber of phi and psi angles: None\n")
        return

    for ChainID in ChainIDs:
        if re.match("^Categories$", OptionsInfo["PhiPsiMode"], re.I):
            # Retrieve phi and psi angles by categories used for Ramachandran plots
            GeneralPhiPsiInfo, GlyPhiPsiInfo, ProPhiPsiInfo, PreProPhiPsiInfo = PyMOLUtil.GetPhiPsiCategoriesResiduesInfo(MolName, ChainID)
            
            SetupAndListPhiPsiResiduesInfo(GeneralPhiPsiInfo, "General (All residues except glycine, proline, or pre-proline)", ChainID)
            SetupAndListPhiPsiResiduesInfo(GlyPhiPsiInfo, "Glycine (Only glycine residues)", ChainID)
            SetupAndListPhiPsiResiduesInfo(ProPhiPsiInfo, "Proline (Only proline residues)", ChainID)
            SetupAndListPhiPsiResiduesInfo(PreProPhiPsiInfo, "Pre-Proline (Only residues before proline not including glycine or proline)", ChainID)
        else:
            PhiPsiResiduesInfo = PyMOLUtil.GetPhiPsiResiduesInfo(MolName, ChainID, Categorize = True)
            SetupAndListPhiPsiResiduesInfo(PhiPsiResiduesInfo, "All", ChainID)

def SetupAndListPhiPsiResiduesInfo(PhiPsiResiduesInfo, AnglesType, ChainID):
    """Setup and  list information for phi and psi torsion angles."""

    MiscUtil.PrintInfo("\nChainID: %s; Phi and Psi angles: %s" %(ChainID, AnglesType))

    PhiPsiCount = len(PhiPsiResiduesInfo["ResNums"])
    if not PhiPsiCount:
        MiscUtil.PrintInfo("Number of angles: None")
        return

    MiscUtil.PrintInfo("Number of phi and psi angles: %d" % PhiPsiCount)

    CategoriesMode = True if re.match("^Categories$", OptionsInfo["PhiPsiMode"], re.I) else False

    if CategoriesMode:
        MiscUtil.PrintInfo("\nResNum, ResName, Phi, Psi")
    else:
        MiscUtil.PrintInfo("\nResNum, ResName, Phi, Psi, Category")
    
    Precision = OptionsInfo["PhiPsiPrecision"]
    for ResNum in PhiPsiResiduesInfo["ResNums"]:
        ResName = PhiPsiResiduesInfo["ResName"][ResNum]
        Phi = PhiPsiResiduesInfo["Phi"][ResNum]
        Psi = PhiPsiResiduesInfo["Psi"][ResNum]
        if CategoriesMode:
            MiscUtil.PrintInfo("%s, %s, %.*f, %.*f" % (ResNum, ResName, Precision, Phi, Precision, Psi))
        else:
            Category = PhiPsiResiduesInfo["Category"][ResNum]
            MiscUtil.PrintInfo("%s, %s, %.*f, %.*f, %s " % (ResNum, ResName, Precision, Phi, Precision, Psi, Category))

def ListInterfaceResiduesInfo(MolName, ChainIDs):
    """List interface residues between all unique pairs of polymer chains
    with in a molecule."""
    
    if not (OptionsInfo["All"] or OptionsInfo["InterfaceResidues"]):
        return

    MiscUtil.PrintInfo("\nListing interface residues information...")

    # Get chain ID pairs for identifying interface residues...
    ChainIDsPairs = GetChainIDsPairsForInterfaceResidues(ChainIDs)
    if not len(ChainIDsPairs):
        MiscUtil.PrintInfo("Valid interface chain ID pairs: None")
        return

    InterfaceResiduesMethod = OptionsInfo["InterfaceResiduesMethod"]
    InterfaceResiduesCutoff = OptionsInfo["InterfaceResiduesCutoff"]
    MiscUtil.PrintInfo("Methodology: %s; Cutoff: %.2f" % (InterfaceResiduesMethod, InterfaceResiduesCutoff))
    
    for ChainIDsPair in ChainIDsPairs:
        InterfaceChainIDs1, InterfaceChainIDs2 = ChainIDsPair
        InterfaceID = "%s-%s" % ("+".join(InterfaceChainIDs1), "+".join(InterfaceChainIDs2))

        InterfaceResiduesInfo1, InterfaceResiduesInfo2 = GetInterfaceChainsAndResiduesInfo(MolName, InterfaceChainIDs1, MolName, InterfaceChainIDs2, InterfaceResiduesMethod, InterfaceResiduesCutoff)

        ChainNames1, ResiduesInfo1 = SetupSelectionChainsResiduesInfo(InterfaceResiduesInfo1)
        ChainNames2, ResiduesInfo2 = SetupSelectionChainsResiduesInfo(InterfaceResiduesInfo2)
        
        if len(ChainNames1) and len(ChainNames2):
            MiscUtil.PrintInfo("\nListing interface residues for interface chain IDs: %s" % (InterfaceID))
            
            ListInterfaceChainsAndResiduesInfo(InterfaceID,  InterfaceChainIDs1, ChainNames1, ResiduesInfo1)
            ListInterfaceChainsAndResiduesInfo(InterfaceID,  InterfaceChainIDs2, ChainNames2, ResiduesInfo2)
        else:
            MiscUtil.PrintInfo("\nListing interface residues for interface chain IDs: %s" % (InterfaceID))
            MiscUtil.PrintInfo("Interface chain IDs: None; ChainID: None; Number of interface residues: 0")

def ListInterfaceChainsAndResiduesInfo(InterfaceID, InterfaceChainIDs, ChainNames, ResiduesInfo):
    """List interface chains and residues. """

    for ChainID in InterfaceChainIDs:
        if not ChainID in ChainNames:
            MiscUtil.PrintInfo("\nInterface chain IDs: %s; ChainID: %s; Number of interface residues: 0" % (InterfaceID, ChainID))
            continue
        
        for Index in range(0, len(ChainNames)):
            ChainName = ChainNames[Index]
            ChainResiduesInfo = ResiduesInfo[Index]
            if re.match(ChainID, ChainName, re.I):
                MiscUtil.PrintInfo("\nInterface chain IDs: %s; ChainID: %s; Number of interface residues: %d" % (InterfaceID, ChainName, ChainResiduesInfo[0]))
                MiscUtil.PrintInfo("Residue distribution: %s\nResidue IDs: %s" % (ChainResiduesInfo[1], ChainResiduesInfo[2]))
                continue

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

def GetChainIDsPairsForInterfaceResidues(ChainIDs):
    """Get chain IDs pairs for identifying interface residues. """

    ChainIDsPairsList = []
    
    InterfaceResiduesChainsList = OptionsInfo["InterfaceResiduesChainsList"]
    if not len(InterfaceResiduesChainsList):
        # Use first two chain IDs...
        if len(ChainIDs) >= 2:
            ChainIDsPair = [ChainIDs[0], ChainIDs[1]]
            ChainIDsPairsList.append(ChainIDsPair)
        return ChainIDsPairsList

    # Validate specified pairwise chain IDs...
    for Index in range(0, len(InterfaceResiduesChainsList), 2):
        ChainIDs1 = InterfaceResiduesChainsList[Index]
        ChainIDs2 = InterfaceResiduesChainsList[Index + 1]
        
        ValidChainIDs = True
        SpecifiedChainIDs = []
        SpecifiedChainIDs.extend(ChainIDs1)
        SpecifiedChainIDs.extend(ChainIDs2)
        
        for ChainID in (SpecifiedChainIDs):
            if not ChainID in ChainIDs:
                ValidChainIDs = False
                MiscUtil.PrintWarning("The chain ID, %s, specified using \"--interfaceResiduesChains\" for a chain IDs pairs is not a valid chain ID." % (ChainID))
        
        if not ValidChainIDs:
            MiscUtil.PrintWarning("Ignoring chain IDs pair: %s, %s" % ("+".join(ChainIDs1), "+".join(ChainIDs2)))
            continue
        
        ChainIDsPair = [ChainIDs1, ChainIDs2]
        ChainIDsPairsList.append(ChainIDsPair)

    return ChainIDsPairsList
    
def SetupSelectionResiduesInfo(SelectionInfo):
    """Setup residues info."""
    
    # Setup distribution of residues...
    LineWords = []
    ResiduesCount = 0
    SortedResNames = sorted(SelectionInfo["ResNames"], key = lambda ResName: SelectionInfo["ResCount"][ResName], reverse = True)
    for ResName in SortedResNames:
        ResCount = SelectionInfo["ResCount"][ResName]
        LineWords.append("%s - %s" % (ResName, ResCount))
        ResiduesCount += ResCount
    
    ResiduesDistribution = "; ".join(LineWords) if len(LineWords) else None
    
    # Setup residue IDs sorted by residue numbers...
    ResNumMap = {}
    for ResName in SelectionInfo["ResNames"]:
        for ResNum in SelectionInfo["ResNum"][ResName]:
            ResNumMap[ResNum] = ResName
    
    LineWords = []
    for ResNum in sorted(ResNumMap, key = int):
        ResName = ResNumMap[ResNum]
        ResID = "%s_%s" % (ResName, ResNum)
        LineWords.append(ResID)
    ResiduesIDs = ", ".join(LineWords) if len(LineWords) else None
                
    return ResiduesCount, ResiduesDistribution, ResiduesIDs

def SetupSelectionChainsResiduesInfo(SelectionInfo):
    """Setup chains and residues info."""

    ChainNames = []
    ResiduesInfo = []

    for ChainID in SelectionInfo["ChainIDs"]:
        ChainNames.append(ChainID)
        
        # Setup distribution of residues...
        LineWords = []
        ResiduesCount = 0
        SortedResNames = sorted(SelectionInfo["ResNames"][ChainID], key = lambda ResName: SelectionInfo["ResCount"][ChainID][ResName], reverse = True)
        for ResName in SortedResNames:
            ResCount = SelectionInfo["ResCount"][ChainID][ResName]
            LineWords.append("%s - %s" % (ResName, ResCount))
            ResiduesCount += ResCount
        
        ResiduesDistribution = "; ".join(LineWords) if len(LineWords) else None
        
        # Setup residue IDs sorted by residue numbers...
        ResNumMap = {}
        for ResName in SelectionInfo["ResNames"][ChainID]:
            for ResNum in SelectionInfo["ResNum"][ChainID][ResName]:
                ResNumMap[ResNum] = ResName
        
        LineWords = []
        for ResNum in sorted(ResNumMap, key = int):
            ResName = ResNumMap[ResNum]
            ResID = "%s_%s" % (ResName, ResNum)
            LineWords.append(ResID)
        ResiduesIDs = ", ".join(LineWords) if len(LineWords) else None

        ResiduesInfo.append([ResiduesCount, ResiduesDistribution, ResiduesIDs])
                
    return ChainNames, ResiduesInfo

def ListBoundingBoxInfo(MolName):
    """List bounding box information. """
    
    if not (OptionsInfo["All"] or OptionsInfo["BoundingBox"]):
        return

    MolSelection = "(%s)" % MolName
    MolExtents = pymol.cmd.get_extent(MolSelection)
    
    XMin, YMin, ZMin = MolExtents[0]
    XMax, YMax, ZMax = MolExtents[1]

    XSize = abs(XMax - XMin)
    YSize = abs(YMax - YMin)
    ZSize = abs(ZMax - ZMin)

    MiscUtil.PrintInfo("\nBounding box coordinates: <XMin, XMax> - <%.3f, %.3f>; <YMin, YMax> - <%.3f, %.3f>; <ZMin, ZMax> - <%.3f, %.3f>" % (XMin, XMax, YMin, YMax, ZMin, ZMax))
    MiscUtil.PrintInfo("Bounding box size in Angstroms: XSize - %.3f; YSize - %.3f; ZSize - %.3f" % (XSize, YSize, ZSize))
    
def ListFileSizeAndModificationInfo(Infile):
    """List file size and modification time info."""

    MiscUtil.PrintInfo("\nFile size: %s" % MiscUtil.GetFormattedFileSize(Infile))
    MiscUtil.PrintInfo("Last modified: %s" % time.ctime(os.path.getmtime(Infile)))
    MiscUtil.PrintInfo("Created: %s" % time.ctime(os.path.getctime(Infile)))
    
def GetPocketSelectionResiduesInfoLabel(SelectionType):
    """Get pocket residues info label for a specified selection type. """

    SelectionLabel = None
    
    if re.match("^Pockets$", SelectionType, re.I):
        SelectionLabel = "polymer"
    elif re.match("^PocketSolvents$", SelectionType, re.I):
        SelectionLabel = "solvent"
    elif re.match("^PocketInorganics$", SelectionType, re.I):
        SelectionLabel = "inorganic"
    else:
        MiscUtil.PrintError("Failed to retrieve pocket residues label information: Selection type %s is not valid..." % SelectionType)

    return SelectionLabel

def ProcessOptions():
    """Process and validate command line arguments and options"""

    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()

    OptionsInfo["All"] = Options["--all"]
    OptionsInfo["BoundingBox"] = Options["--boundingBox"]
    
    OptionsInfo["Chains"] = Options["--chains"]
    
    OptionsInfo["CountResidues"] = Options["--countResidues"]
    OptionsInfo["Header"] = Options["--header"]
    
    OptionsInfo["Infiles"] = Options["--infiles"]
    OptionsInfo["InfilesNames"] =  Options["--infilesNames"]
    
    OptionsInfo["Inorganics"] = Options["--inorganics"]

    OptionsInfo["InterfaceResidues"] = Options["--interfaceResidues"]
    OptionsInfo["InterfaceResiduesMethod"] = Options["--interfaceResiduesMethod"]
    
    InterfaceResiduesChains = Options["--interfaceResiduesChains"]
    InterfaceResiduesChainsList = []
    if not re.match("^auto$", InterfaceResiduesChains, re.I):
        InterfaceResiduesChains = re.sub(" ", "", Options["--interfaceResiduesChains"])
        InterfaceResiduesChainsWords = InterfaceResiduesChains.split(",")
        if len(InterfaceResiduesChainsWords) % 2:
            MiscUtil.PrintError("The number of comma delimited chain IDs, %d, specified using \"--interfaceResiduesChains\" option, \"%s\",  must be a multple of 2." % (len(InterfaceResiduesChainsWords), Options["--interfaceResiduesChains"]))
        for ChainID in InterfaceResiduesChainsWords:
            ChainIDWords = ChainID.split("+")
            InterfaceResiduesChainsList.append(ChainIDWords)
    OptionsInfo["InterfaceResiduesChains"] = InterfaceResiduesChains
    OptionsInfo["InterfaceResiduesChainsList"] = InterfaceResiduesChainsList
    
    InterfaceResiduesCutoff = Options["--interfaceResiduesCutoff"]
    InterfaceResiduesMethod = OptionsInfo["InterfaceResiduesMethod"]
    if re.match("^auto$", InterfaceResiduesCutoff, re.I):
        if re.match("^BySASAChange$", InterfaceResiduesMethod, re.I):
            InterfaceResiduesCutoff = 1.0
        elif re.match("^ByHeavyAtomsDistance$", InterfaceResiduesMethod, re.I):
            InterfaceResiduesCutoff = 5.0
        elif re.match("^ByCAlphaAtomsDistance$", InterfaceResiduesMethod, re.I):
            InterfaceResiduesCutoff = 8.0
        else:
            MiscUtil.PrintError("The specified value, %s, for option \"--interfaceResiduesMethod\" is not supported." % (InterfaceResiduesMethod))
    else:
        InterfaceResiduesCutoff = float(InterfaceResiduesCutoff)
    OptionsInfo["InterfaceResiduesCutoff"] = InterfaceResiduesCutoff
    
    OptionsInfo["Ligands"] = Options["--ligands"]
    
    OptionsInfo["PocketLigands"] = Options["--pocketLigands"]
    OptionsInfo["PocketDistanceCutoff"] = float(Options["--pocketDistanceCutoff"])
    OptionsInfo["PocketSolvents"] = Options["--pocketSolvents"]
    OptionsInfo["PocketInorganics"] = Options["--pocketInorganics"]
    
    OptionsInfo["PhiPsi"] = Options["--phiPsi"]
    OptionsInfo["PhiPsiMode"] = Options["--phiPsiMode"]
    OptionsInfo["PhiPsiPrecision"] = int(Options["--phiPsiPrecision"])
    
    OptionsInfo["Solvents"] = Options["--solvents"]
    
    OptionsInfo["SurfaceResidues"] = Options["--surfaceResidues"]
    OptionsInfo["SurfaceResiduesCutoff"] = float(Options["--surfaceResiduesCutoff"])
    OptionsInfo["SurfaceResiduesIDs"] = True if re.match("^Yes$", Options["--surfaceResiduesIDs"], re.I) else False

    # Always list header, chains, and ligands...
    OptionNames = ["Chains", "Header", "Ligands"]
    for Name in OptionNames:
        if Name in OptionsInfo:
            OptionsInfo[Name] = True
        else:
            MiscUtil.PrintError("Option name %s is not a valid name..." % s)
    
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

    # Expand infile names..
    InfilesNames = MiscUtil.ExpandFileNames(Options["--infiles"], ",")

    # Validate file extensions...
    for Infile in InfilesNames:
        MiscUtil.ValidateOptionFilePath("-i, --infiles", Infile)
        MiscUtil.ValidateOptionFileExt("-i, --infiles", Infile, "pdb cif")
    Options["--infilesNames"] = InfilesNames
    
    if not re.match("^auto$", Options["--interfaceResiduesCutoff"], re.I):
        MiscUtil.ValidateOptionFloatValue("--interfaceResiduesCutoff", Options["--interfaceResiduesCutoff"], {">": 0.0})
    MiscUtil.ValidateOptionTextValue("--interfaceResiduesMethod", Options["--interfaceResiduesMethod"], "BySASAChange ByHeavyAtomsDistance ByCAlphaAtomsDistance")
    
    MiscUtil.ValidateOptionFloatValue("--pocketDistanceCutoff", Options["--pocketDistanceCutoff"], {">": 0.0})
    MiscUtil.ValidateOptionTextValue("--phiPsiMode", Options["--phiPsiMode"], "All  Categories")
    MiscUtil.ValidateOptionIntegerValue("--phiPsiPrecision", Options["--phiPsiPrecision"], {">": 0})
    MiscUtil.ValidateOptionFloatValue("--surfaceResiduesCutoff", Options["--surfaceResiduesCutoff"], {">": 0.0})
    MiscUtil.ValidateOptionTextValue("--surfaceResiduesIDs", Options["--surfaceResiduesIDs"], "Yes No")

# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLInfoMacromolecules.py - List information about macromolecules

Usage:
    PyMOLInfoMacromolecules.py [--all] [--boundingBox] [--chains] [--countResidues] 
                               [--header] [--inorganics] [--interfaceResidues]
                               [--interfaceResiduesChains <ChainID1,ChainD2,...>] [--interfaceResiduesMethod <text>]
                               [--interfaceResiduesCutoff <number>] [--ligands] [--pocketLigands]
                               [--pocketDistanceCutoff  <number>] [--pocketSolvents] [--pocketInorganics]
                               [--phiPsi] [--phiPsiMode <All or Categories>] [--phiPsiPrecision <number>]
                               [--surfaceResidues] [--surfaceResiduesCutoff <number>] [--surfaceResiduesIDs <yes or no>]
                               [--solvents] [-w <dir>] -i <infile1,infile2,infile3...>
    PyMOLInfoMacromolecules.py -h | --help | -e | --examples

Description:
    List information regarding  ID, classification, experimental technique, chains,
    solvents, inorganics, ligands, and ligand binding pockets in macromolecules
    present including proteins and nucleic acids.

    The supported input  file format are: PDB (.pdb), mmCIF (.cif)

Options:
    -a, --all
        All available information.
    -b, --boundingBox
        Min and max coordinates for bounding box along with its size.
    -c, --chains
        Number of chains and their IDs. This is also default behavior.
     --countResidues
        Number of residues across chains. The chain residues are identified
        using polymer selection operator available in PyMOL. In addition,
        the non-standard amino acid residues are listed.
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
     --header
        Header information including experimental technique information
        along with any available resolution. This is also default behavior.
    -i, --infiles <infile1,infile2,infile3...>
        A comma delimited list of input files. The wildcards are also allowed
        in file names.
    --inorganics
        Inorganic residues across chains. The inorganic residues are identified
        using inorganic selection operator available in PyMOL.
    --interfaceResidues
        Interface residues between specified pairs of chains.
    --interfaceResiduesChains <ChainID1,Chain1D2,...>  [default: Auto]
        Pairwise comma delimited list of chain IDs for the identification of
        interface residues. Each chain ID may contain mutiple chain IDs
        delimited by a plus sign. For example: A+B,C+D chain pair specifies
        interface between chain complexes A+B and C+D.
        
        The interface residues are identified between first two chains in
        input files by default.
    --interfaceResiduesMethod <text>  [default: BySASAChange]
        Methodology for the identification of interface residues between a pair
        of chains in an input file. The interface residues may be identified by
        change in solvent accessible surface area (SASA) for a residue between
        a chain and chains complex, distance between heavy atoms
        in two chains, or distance between CAlpha atoms. Possible values:
        BySASAChange, ByHeavyAtomsDistance, or ByCAlphaAtomsDistance. 
    --interfaceResiduesCutoff <number>  [default: auto]
        Cutoff value used by different methodologies during identification of
        interface residues between a pair of chains. The default values are
        shown below:
        
            BySASAChange: 1.0; Units: Angstrom**2 [ Ref 141 ]
            ByHeavyAtomsDistance: 5.0; Units: Angstrom [ Ref 142 ]
            ByCAlphaAtomsDistance: 8.0; Units: Angstrom [ Ref 143 ]
        
    -l, --ligands
        Ligands across chains. This is also default behavior. The ligands
        residues are identified using organic selection operator available
        in PyMOL.
    -p, --pocketLigands
        Chain residues in ligand pockets.
    --pocketDistanceCutoff <number>  [default: 5.0]
        Distance in Angstroms for identifying pocket residues around ligands.
    --pocketSolvents
        Solvent residues in ligand pockets. The solvent residues are identified
        using solvent selection operator available in PyMOL.
    --pocketInorganics
        Inorganic residues in ligand pockets. The inorganic residues are identified
        using Inorganic selection operator available in PyMOL.
     --phiPsi
        Phi and psi torsion angles across chains in macromolecules containing
        amino acids.
     --phiPsiMode <All or Categories>  [default: Categories]
        List all phi and psi torsion angles for residues as a single group or split
        them into the following categories corresponding to four types of
        Ramachandran plots:
        
            General: All residues except glycine, proline, or pre-proline
            Glycine: Only glycine residues
            Proline: Only proline residues
            Pre-Proline: Only residues before proline not including glycine
                or proline
        
     --phiPsiPrecision <number>  [default: 2]
        Precision for listing phi and psi torsion angles.
    -s, --solvents
        Solvent residues across chains. The solvent residues are identified
        using solvent selection operator available in PyMOL.
    --surfaceResidues
        Surface and buried residues in chains.
    --surfaceResiduesCutoff <number>  [default: 2.5]
        Solvenet Accessible Surface Area (SASA) cutoff value in Angstroms**2
        for surface and buried resiudes in chains. The residues with SASA less than
        the cutoff value correspond to burried residues.
    --surfaceResiduesIDs <yes or no>  [default: No]
        List residue IDs for surface and buried residues during listing of the
        distribution of these residues for '--surfaceResidues' option.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To list header, chains, and ligand information for macromolecules in input
    file, type:

        % PyMOLInfoMacromolecules.py  -i Sample3.pdb

    To list all available information for macromolecules in input files, type:

        % PyMOLInfoMacromolecules.py  -a  -i "Sample3.pdb,Sample4.pdb"

    To list pockets residues information along with other default information
    for marcomolecules in input file, type:

        % PyMOLInfoMacromolecules.py  -p --pocketDistanceCutoff 4.5 
        --pocketSolvents  --pocketInorganics -i Sample3.pdb

    To list chain residues information along with other default information
    for marcomolecules in input file, type:

        % PyMOLInfoMacromolecules.py  -c --countResidues --solvents
        --inorganics -i "Sample3.pdb,Sample4.pdb" 

    To list interface residues between first two chains by SASA change for
    marcomolecules in input file, type:

        % PyMOLInfoMacromolecules.py  --interfaceResidues
        -i Sample3.pdb

    To list interface residues between chains E and I  by heay atoms
    distance for marcomolecules in input file, type:

        % PyMOLInfoMacromolecules.py  --interfaceResidues
        --interfaceResiduesChains E,I  --interfaceResiduesMethod
        ByHeavyAtomsDistance --interfaceResiduesCutoff  5 -i Sample3.pdb

    To list interface residues between two sets of chains by SASA change for
    marcomolecules in input file, type:

        % PyMOLInfoMacromolecules.py  --interfaceResidues
        --interfaceResiduesChains "A+B,C+D" -i Sample8.pdb 

Author:
    Manish Sud(msud@san.rr.com)

See also:
    DownloadPDBFiles.pl, PyMOLSplitChainsAndLigands.py,
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
