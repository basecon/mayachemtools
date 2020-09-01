#
# File: PyMOLUtil.py
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2020 Manish Sud. All rights reserved.
#
# The functionality available in this file is implemented using PyMOL, a
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

import os
import sys
import re

from pymol import cmd, stored, CmdException

import MiscUtil

__all__ = ["AreAminoAcidResiduesPresent", "AreNucleicAcidResiduesPresent", "CalculateCenterOfMass", "ConvertFileFormat", "ConvertPMLFileToPSEFile", "GetChains", "GetChainsAndLigandsInfo", "GetAminoAcidResiduesInfo", "GetInorganicResiduesInfo", "GetInterfaceChainsResiduesByCAlphaAtomsDistance", "GetInterfaceChainsResiduesByHeavyAtomsDistance", "GetInterfaceChainsResiduesBySASAChange", "GetLargestLigand", "GetLigandResiduesInfo", "GetLigands", "GetMolecules", "GetNucleicAcidResiduesInfo", "GetPocketInorganicResiduesInfo", "GetPocketPolymerResiduesInfo", "GetPocketSolventResiduesInfo",  "GetPolymerResiduesInfo", "GetPhiPsiResiduesInfo", "GetPhiPsiChainsAndResiduesInfo", "GetPhiPsiCategoriesResiduesInfo", "GetSelectionResiduesInfo", "GetSolventResiduesInfo", "GetSurfaceResiduesInfo", "ProcessChainsAndLigandsOptionsInfo", "ProcessResidueTypesOptionsInfo", "ProcessChainSelectionsOptionsInfo", "ProcessSurfaceAtomTypesColorsOptionsInfo", "SetupPMLForAlignment", "SetupPMLForBFactorPuttyView", "SetupPMLForBallAndStickView", "SetupPMLForEnableDisable", "SetupPMLForGroup", "SetupPMLForHydrophobicSurfaceView", "SetupPMLForHydrophobicAndChargeSurfaceView", "SetupPMLForInorganicView", "SetupPMLForLigandPocketInorganicView", "SetupPMLForLigandPocketSolventView", "SetupPMLForLigandPocketView", "SetupPMLForLigandView", "SetupPMLForPolarContactsView", "SetupPMLForHydrophobicContactsView", "SetupPMLForPolymerChainComplexView", "SetupPMLForPolymerChainView", "SetupPMLForPolymerComplexView", "SetupPMLForSolventView", "SetupPMLForSurfaceView", "SetupPMLForSelectionDisplayView", "SetupPMLHeaderInfo"]

def GetMolecules(Selection = "all"):
    """Get names of molecule objects in a selection or all molecule objects.

    Arguments:
        Selection: (str): A PyMOL selection.

    Returns:
        list: Names of molecule objects.

    """
    Names = cmd.get_object_list('(' + Selection + ')')

    return Names

def GetChains(MoleculeName, RemoveEmpty = True):
    """Get chain identifiers present in a molecule.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        RemoveEmpty (bool): Remove empty chain ID from the list of chain IDs
            returned by PyMOL.

    Returns:
        list: Names of chains present in a molecule, sorted alphabetically in a
            ascending order.

    """
    if not len(MoleculeName):
        return None

    ChainIDs = []
    try:
        ChainIDs = cmd.get_chains('model %s' % MoleculeName)
    except CmdException as ErrMsg:
        MiscUtil.PrintWarning("PyMOLUtil.GetChains: Invalid molecule name: %s" % MoleculeName)

    if not len(ChainIDs):
        return None

    # Remove empty Chain IDs from the list...
    if RemoveEmpty:
        NonEmptyChainIDs = []
        for ChainID in ChainIDs:
            if len(ChainID):
                NonEmptyChainIDs.append(ChainID)
        if len(NonEmptyChainIDs) != len(ChainIDs):
            MiscUtil.PrintInfo("PyMOLUtil.GetChains: Removing non-empty chain IDs from the list of chain IDs...")
            
        ChainIDs = NonEmptyChainIDs
        
    return ChainIDs

def GetChainsAndLigandsInfo(Infile, MolName, Quite = False, LigandSortBy = "Size", LigandSortOrder = "Auto", LigandIgnoreHydrogens = "Yes"):
    """Get chain identifiers present in a molecule along with names of the
    ligands present in chains. Ligands are identified using PyMOL 'organic'
    selection.

    Arguments:
        Infile (str) : Name of a file.
        MolName (str) : Name to use for PyMOL molecule object.
        Quite (bool) : Flag 
        LigandSortBy (str): Sort ligand names alphabetically or by size. Possible
            values: Alphabetical or Size
        LigandSortOrder (str): Sort order for sorting ligands. Possible values:
            Ascending, Descending, Auto. The 'Auto' value implies automatic
            determination of sort order based on the value of 'SortBy'.
            Automatic defaults: Descending for SortBy value of Size; Ascending
            for SortBy value of Alphabetical.
        LigandIgnoreHydrogens (str): Ignore hydrogens during determination of ligand
            size.

    Returns:
        dict: A dictionary containing list of chain identifiers and dictionaries
            of chains containing lists of ligand names for each chain. Names of
            ligands present in chain for a molecule sorted by size or
            alphabetically.

    Examples:

        ChainsAndLigandsInfo = GetChainsAndLigandsInfo(Infile, MolName)
        for ChainID in ChainsAndLigandsInfo["ChainIDs"]:
            for LigandID in ChainsAndLigandsInfo["LigandIDs"][ChainID]:
                MiscUtil.PrintInfo("ChainID: %s; LigandID: %s" % (ChainID,
                    LigandID))

    """
    if not Quite:
        MiscUtil.PrintInfo("\nRetrieving chain and ligand information from input file %s..." % Infile)
    
    # Collect chains and ligands information with ligands sorted by size to be used for
    # identification of largest ligand at the top...
    cmd.load(Infile, MolName)
    ChainsAndLigandsInfo = _GetChainsAndLigands(MolName, LigandSortBy = "size", LigandSortOrder = "descending")
    cmd.delete(MolName)

    # Print out chain and ligand IDs...
    if not Quite:
        ChainIDs = ", ".join(ChainsAndLigandsInfo["ChainIDs"]) if len(ChainsAndLigandsInfo["ChainIDs"]) else "None"
        MiscUtil.PrintInfo("Chain IDs: %s" % ChainIDs)
        
        for ChainID in ChainsAndLigandsInfo["ChainIDs"]:
            LigandIDs = ", ".join(ChainsAndLigandsInfo["LigandIDs"][ChainID]) if len(ChainsAndLigandsInfo["LigandIDs"][ChainID]) else "None"
            MiscUtil.PrintInfo("Chain ID: %s; LigandIDs: %s" % (ChainID, LigandIDs))
    
    return ChainsAndLigandsInfo

def GetLigands(MoleculeName, ChainName, SortBy = "Size", SortOrder = "Auto", IgnoreHydrogens = "Yes"):
    """Get names of ligands present in a chain of a  molecule. Ligands are
    identified using PyMOL 'organic' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        SortBy (str): Sort ligand names alphabetically or by size. Possible
            values: Alphabetical or Size
        SortOrder (str): Sort order for sorting ligands. Possible values:
            Ascending, Descending, Auto. The 'Auto' value implies automatic
            determination of sort order based on the value of 'SortBy'.
            Automatic defaults: Descending for SortBy value of Size; Ascending
            for SortBy value of Alphabetical.
        IgnoreHydrogens (str): Ignore hydrogens during determination of ligand
            size.

    Returns:
        list: Names of ligands present in chain for a molecule sorted by size
            or alphabetically.

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    LigandsInfoMap = _GetLigandsInfo(MoleculeName, ChainName, SortBy, SortOrder, IgnoreHydrogens)
    
    LigandIDs = LigandsInfoMap["LigandResNames"]
    if not len(LigandIDs):
        LigandIDs = None
    
    return LigandIDs

def GetLargestLigand(MoleculeName, ChainName, IgnoreHydrogens = 'Yes'):
    """Get name of the largest ligand for a chain present in a molecule. Ligands
    are identified using PyMOL 'organic' selection.

    Arguments:
        IgnoreHydrogens (str): Ignore hydrogens during determination of ligand
            size.

    Returns:
        str: Name of the largest ligand present in a chain. 

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SortBy = "Size"
    SortOrder = "Descending"
    LigandsInfoMap = _GetLigandsInfo(MoleculeName, ChainName, SortBy, SortOrder, IgnoreHydrogens)
    LigandIDs = LigandsInfoMap["LigandResNames"]
    
    if len(LigandIDs):
        LigandID = LigandIDs[0]
    else:
        LigandID = None

    return LigandID

def _GetChainsAndLigands(MoleculeName, LigandSortBy = "Size", LigandSortOrder = "Auto", LigandIgnoreHydrogens = "Yes"):
    """Get chain identifiers in a molecule along with names of the ligands
    present in chains. Ligands are identified using PyMOL 'organic' selection.
    """
    if not len(MoleculeName):
        return None

    ChainIDs = GetChains(MoleculeName)
    if ChainIDs is None:
        return None
        
    ChainsAndLigandsMap = {}
    ChainsAndLigandsMap["ChainIDs"] = []
    ChainsAndLigandsMap["LigandIDs"] = {}
    
    for ChainID in ChainIDs:
        ChainsAndLigandsMap["ChainIDs"].append(ChainID)
        ChainsAndLigandsMap["LigandIDs"][ChainID] = []
        
        LigandIDs = GetLigands(MoleculeName, ChainID, SortBy = LigandSortBy, SortOrder = LigandSortOrder, IgnoreHydrogens = LigandIgnoreHydrogens)
        if LigandIDs is not None:
            ChainsAndLigandsMap["LigandIDs"][ChainID] = LigandIDs
        
    return ChainsAndLigandsMap

def _GetLigandsInfo(MoleculeName, ChainName, SortBy = "Size", SortOrder = "Auto", IgnoreHydrogens = "Yes"):
    """Retrieve information about ligands present in a chain of a molecule."""

    if not MiscUtil.CheckTextValue(SortBy, "Size Alphabetical"):
        MiscUtil.PrintError("PyMOLUtil._GetLigandsInfo: The value specified, %s, for parameter SortBy is not valid. SupportedValues: Size Alphabetical" % SortBy)
        
    if not MiscUtil.CheckTextValue(SortOrder, "Ascending Descending Auto"):
        MiscUtil.PrintError("PyMOLUtil._GetLigandsInfo: The value specified, %s, for parameter SortOrder is not valid. SupportedValues: Ascending Descending Auto" % SortOrder)
        
    if not MiscUtil.CheckTextValue(IgnoreHydrogens, "Yes No"):
        MiscUtil.PrintError("PyMOLUtil._GetLigandsInfo: The value specified, %s, for parameter IgnoreHydrogens is not valid. SupportedValues: Yes No" % IgnoreHydrogens)
        
    SortBySize = True if re.match("^Size$", SortBy, re.I) else False
    if re.match("^Auto$", SortOrder, re.I):
        SortOrderDescending = True if re.match("^Size$", SortBy, re.I) else False
    else:
        SortOrderDescending = True if re.match("^Descending$", SortOrder, re.I) else False
    IgnoreHydrogenAtoms = True if re.match("^Yes$", IgnoreHydrogens, re.I) else False

    # Set up a command to retrieve all appropriate ligand atoms in organic ligands...
    SelectionCmd = "%s and chain %s and organic" % (MoleculeName, ChainName)
    if IgnoreHydrogenAtoms:
        SelectionCmd = "%s and not hydro" % (SelectionCmd)

    # Retrieve atoms...
    stored.LigandsInfo = []
    cmd.iterate(SelectionCmd, "stored.LigandsInfo.append([resi, resn])")

    # Retrieve ligands...
    LigandsInfoMap = {}
    LigandsInfoMap["LigandResNames"] = []
    LigandsInfoMap["LigandAtomCount"] = {}
    LigandsInfoMap["LigandResNumber"] = {}
    
    for LigandResNum, LigandResName in stored.LigandsInfo:
        if LigandResName in LigandsInfoMap["LigandResNames"]:
            LigandsInfoMap["LigandAtomCount"][LigandResName] += 1
        else:
            LigandsInfoMap["LigandResNames"].append(LigandResName)
            LigandsInfoMap["LigandAtomCount"][LigandResName] = 1
            LigandsInfoMap["LigandResNumber"][LigandResName] = LigandResNum

    if not len(LigandsInfoMap["LigandResNames"]):
        return LigandsInfoMap
        
    # Sort ligand names...
    ReverseOrder = True if SortOrderDescending else False
    if SortBySize:
        SortedLigandResNames = sorted(LigandsInfoMap["LigandResNames"], key = lambda LigandResName: LigandsInfoMap["LigandAtomCount"][LigandResName], reverse = ReverseOrder)
    else:
        # Sort alphabetically...
        SortedLigandResNames = sorted(LigandsInfoMap["LigandResNames"], reverse = ReverseOrder)

    LigandsInfoMap["LigandResNames"] = SortedLigandResNames
    
    return LigandsInfoMap

def AreAminoAcidResiduesPresent(MoleculeName, ChainName, Type = "Any"):
    """Check for the presence of amino acid residues in a chain of a
    molecule. Chains are identified using PyMOL 'polymer' selection.
    Nonstandard amino acid residues correspond to all residues other than
    the standard amino acids and nucleic acids. Any amino acid residues cover
    all residues other than the standard nucleic acids.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        Type (str): Types of amino acids: Standard, NonStandard, Any

    Returns:
        boolean: True or False.

    """
    
    if not (len(MoleculeName) and len(ChainName)):
        return False

    if not re.match("^(Standard|NonStandard|Any)$", Type, re.I):
        MiscUtil.PrintError("PyMOLUtil.AreAminoAcidResiduesPresent: Invalid amino acid type: %s" % Type)

    SelectionCmd = _SetupAminoAcidResiduesSelectionCmd(MoleculeName, ChainName, Type)
    
    Status = True if cmd.count_atoms(SelectionCmd) else False
    
    return Status
    
def AreNucleicAcidResiduesPresent(MoleculeName, ChainName):
    """Check for the presence of nucleic acid residues in a chain of a 
    molecule. Chains are identified using PyMOL 'polymer' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        boolean: True or False.

    """
    
    if not (len(MoleculeName) and len(ChainName)):
        return False
    
    ResidueNames = _GetNucleicAcidResidueNames()
    if not len(ResidueNames):
        return False
    
    ResidueNamesSelection = "+".join(ResidueNames)
    SelectionCmd = "(%s and chain %s and polymer and (resn %s))" % (MoleculeName, ChainName, ResidueNamesSelection)

    Status = True if cmd.count_atoms(SelectionCmd) else False
    
    return Status
    
def _SetupAminoAcidResiduesSelectionCmd(MoleculeName, ChainName, Type):
    """Set up amino acids selection command for PyMOL. """
    
    AminoAcidsSelection = "+".join(_GetAminoAcidResidueNames())
    NucleicAcidsSelection = "+".join(_GetNucleicAcidResidueNames())
    
    if re.match("^Standard$", Type, re.I):
        SelectionCmd = "(%s and chain %s and polymer and (resn %s))" % (MoleculeName, ChainName, AminoAcidsSelection)
    elif re.match("^NonStandard$", Type, re.I):
        SelectionCmd = "(%s and chain %s and polymer and (not ((resn %s) or (resn %s))))" % (MoleculeName, ChainName, AminoAcidsSelection, NucleicAcidsSelection)
    else:
        # Any amino acid resiudes...
        SelectionCmd = "(%s and chain %s and polymer and (not (resn %s)))" % (MoleculeName, ChainName, NucleicAcidsSelection)

    return SelectionCmd

def _GetAminoAcidResidueNames():
    """Get list of amino acid residue names.
    """
    
    ResidueNames = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    
    return ResidueNames

def _GetNucleicAcidResidueNames():
    """Get list of nucleic acid residue names.
    """
    
    ResidueNames = ["A", "G", "T", "U", "C", "DA", "DG", "DT", "DU", "DC"]
    
    return ResidueNames

def GetPolymerResiduesInfo(MoleculeName, ChainName):
    """Get information for residues present in a chain of a  molecule.
    Chains are identified using PyMOL 'polymer' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetPolymerResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SelectionCmd = "%s and chain %s and polymer" % (MoleculeName, ChainName)
    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetAminoAcidResiduesInfo(MoleculeName, ChainName, Type = "Standard"):
    """Get information for amino acid residues present in a chain of a
    molecule. Chains are identified using PyMOL 'polymer' selection.
    Nonstandard amino acid residues correspond to all residues other than
    the standard amino acids and nucleic acids. Any amino acid residues cover
    all residues other than the standard nucleic acids.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        Type (str): Types of amino acids: Standard, NonStandard, Any

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetPolymerResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    
    if not (len(MoleculeName) and len(ChainName)):
        return None

    if not re.match("^(Standard|NonStandard|Any)$", Type, re.I):
        MiscUtil.PrintError("PyMOLUtil.GetAminoAcidResiduesInfo: Invalid amino acid type: %s" % Type)
    
    SelectionCmd = _SetupAminoAcidResiduesSelectionCmd(MoleculeName, ChainName, Type)
    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetNucleicAcidResiduesInfo(MoleculeName, ChainName):
    """Get information for nucleic acid residues present in a chain of
    a  molecule. Chains are identified using PyMOL 'polymer' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetPolymerResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    
    if not (len(MoleculeName) and len(ChainName)):
        return None

    ResidueNames = _GetNucleicAcidResidueNames()
    
    ResidueNamesSelection = "+".join(ResidueNames)
    SelectionCmd = "(%s and chain %s and polymer and (resn %s))" % (MoleculeName, ChainName, ResidueNamesSelection)
    
    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetSelectionResiduesInfo(SelectionCmd):
    """Get information for residues in a chain specified by a selection command.

    Arguments:
        SelectionCmd (str): PyMOL selection command.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetSelectionResiduesInfo(SelectionCmd)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not len(SelectionCmd):
        return None

    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetSolventResiduesInfo(MoleculeName, ChainName):
    """Get information for solvent residues present in a chain of a  molecule.
    Solvents are identified using PyMOL 'solvent' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetSolventResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SelectionCmd = "%s and chain %s and solvent" % (MoleculeName, ChainName)
    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetInorganicResiduesInfo(MoleculeName, ChainName):
    """Get information for inorganic residues present in a chain of a  molecule.
    Inorganic residues are identified using PyMOL 'inorganic' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetInorganicResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SelectionCmd = "%s and chain %s and inorganic" % (MoleculeName, ChainName)
    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetLigandResiduesInfo(MoleculeName, ChainName):
    """Get information for ligand residues present in a chain of a  molecule.
    Ligands are identified using PyMOL 'organic' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetLigandResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SelectionCmd = "%s and chain %s and organic" % (MoleculeName, ChainName)
    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetPocketPolymerResiduesInfo(MoleculeName, ChainName, LigandResName, LigandResNum, PocketDistanceCutoff):
    """Get information for chain residues present in a pocket around a ligand
    in a molecule. Polymer residues are identified using negation of PyMOL
    selection operators 'organic', 'solvent', and 'inorganic'.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        LigandResName (str): Residue name of a ligand in a chain.
        LigandResNum (str): Residue number of a ligand in a chain.
        PocketDistanceCutoff (float): Distance around ligand to identify pocket
            residues.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetPocketPolymerResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName) and len(LigandResName) and len(LigandResNum)):
        return None

    LigandSelection = "%s and chain %s and organic and resn %s and resi %s" % (MoleculeName, ChainName, LigandResName, LigandResNum)
    MoleculeSelection = "%s and chain %s" % (MoleculeName, ChainName)
    SelectionCmd = "((byresidue (%s) within %.1f of (%s))  and (not solvent) and (not inorganic) and (not organic))" % (MoleculeSelection, PocketDistanceCutoff, LigandSelection)

    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetPocketSolventResiduesInfo(MoleculeName, ChainName, LigandResName, LigandResNum, PocketDistanceCutoff):
    """Get information for solvent residues present in a pocket around a ligand
    in a molecule. Solvent residues are identified using PyMOL 'solvent'
    selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        LigandResName (str): Residue name of a ligand in a chain.
        LigandResNum (str): Residue number of a ligand in a chain.
        PocketDistanceCutoff (float): Distance around ligand to identify pocket
            residues.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetPocketSolventResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName) and len(LigandResName) and len(LigandResNum)):
        return None

    LigandSelection = "%s and chain %s and organic and resn %s and resi %s" % (MoleculeName, ChainName, LigandResName, LigandResNum)
    MoleculeSelection = "%s and chain %s" % (MoleculeName, ChainName)
    SelectionCmd = "((byresidue (%s) within %.1f of (%s))  and solvent)" % (MoleculeSelection, PocketDistanceCutoff, LigandSelection)

    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetPocketInorganicResiduesInfo(MoleculeName, ChainName, LigandResName, LigandResNum, PocketDistanceCutoff):
    """Get information for inorganic residues present in a pocket around a
    ligand in a molecule. Inorganic residues are identified using PyMOL
    'inorganic' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        LigandResName (str): Residue name of a ligand in a chain.
        LigandResNum (str): Residue number of a ligand in a chain.
        PocketDistanceCutoff (float): Distance around a ligand to identify
            pocket residues.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetPocketInorganicResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName) and len(LigandResName) and len(LigandResNum)):
        return None

    LigandSelection = "%s and chain %s and organic and resn %s and resi %s" % (MoleculeName, ChainName, LigandResName, LigandResNum)
    MoleculeSelection = "%s and chain %s" % (MoleculeName, ChainName)
    SelectionCmd = "((byresidue (%s) within %.1f of (%s))  and inorganic)" % (MoleculeSelection, PocketDistanceCutoff, LigandSelection)

    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetPhiPsiResiduesInfo(MoleculeName, ChainName, Categorize = True):
    """Get phi and psi torsion angle information for residues in a chain of a
    molecule containing amino acids.
    
    The phi and psi angles are optionally categorized into the following groups
    corresponding to four types of  Ramachandran plots:
    
    General: All residues except glycine, proline, or pre-proline
    Glycine: Only glycine residues
    Proline: Only proline residues
    Pre-Proline: Only residues before proline not including glycine or proline
    
    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict: A dictionary containing sorted list of residue numbers and
            dictionaries of residue names, phi and psi angles for each residue
            number.

    Examples:

        PhiPsiInfoMap = GetPhiPsiResiduesInfo(MolName, ChainName, True)
        for ResNum in PhiPsiInfoMap["ResNums"]:
            ResName = PhiPsiInfoMap["ResName"][ResNum]
            Phi = PhiPsiInfoMap["Phi"][ResNum]
            Psi = PhiPsiInfoMap["Psi"][ResNum]
            Category = PhiPsiInfoMap["Category"][ResNum]
            MiscUtil.PrintInfo("ResNum: %s; ResName: %s; Phi: %8.2f;
                Psi: %8.2f; Category: %s" % (ResNum, ResName, Phi, Psi,
                Categorize))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None
    
    SelectionCmd = "%s and chain %s" % (MoleculeName, ChainName)
    PhiPsiResiduesInfoMap = _GetSelectionPhiPsiResiduesInfo(SelectionCmd, Categorize)

    return PhiPsiResiduesInfoMap

def GetPhiPsiChainsAndResiduesInfo(MoleculeName, Categorize = True):
    """Get phi and psi torsion angle information for residues across chains in
    a molecule containing amino acids.

    The phi and psi angles are optionally categorized into the following groups
    corresponding to four types of  Ramachandran plots:
    
    General: All residues except glycine, proline, or pre-proline
    Glycine: Only glycine residues
    Proline: Only proline residues
    Pre-Proline: Only residues before proline not including glycine or proline
    
    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.

    Returns:
        dict: A dictionary containing sorted list of residue numbers for each
            chain and dictionaries of residue names, phi and psi angles for each
            residue number.

    Examples:

        PhiPsiInfoMap = GetPhiPsiChainsAndResiduesInfo(MolName)
        for ChainID in PhiPsiInfoMap["ChainIDs"]:
            for ResNum in PhiPsiInfoMap["ResNums"][ChainID]:
                ResName = PhiPsiInfoMap["ResName"][ChainID][ResNum]
                Phi = PhiPsiInfoMap["Phi"][ChainID][ResNum]
                Psi = PhiPsiInfoMap["Psi"][ChainID][ResNum]
                Category = PhiPsiInfoMap["Category"][ChainID][ResNum]
                MiscUtil.PrintInfo("ChainID: %s; ResNum: %s; ResName: %s; Phi: %8.2f;
                    Psi: %8.2f; Category: %s" % (ChainID, ResNum, ResName, Phi,
                    Psi, Category))

    """
    if not len(MoleculeName):
        return None

    SelectionCmd = "%s" % (MoleculeName)
    PhiPsiResiduesInfoMap = _GetSelectionPhiPsiChainsAndResiduesInfo(SelectionCmd, Categorize)

    return PhiPsiResiduesInfoMap

def GetPhiPsiCategoriesResiduesInfo(MoleculeName, ChainName):
    """Get phi and psi torsion angle information for residues in a chain of a
    molecule containing amino acids.

    The phi and psi angles are optionally categorized into the following groups
    corresponding to four types of  Ramachandran plots:
    
    General: All residues except glycine, proline, or pre-proline
    Glycine: Only glycine residues
    Proline: Only proline residues
    Pre-Proline: Only residues before proline not including glycine or proline
    
    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict1: Phi and psi angle information for residues in General category.
            It's a dictionary containing sorted list of residue numbers and
            dictionaries of residue names, phi and psi angles for each residue
            number.
        dict2: Phi and psi angle information for residues in Gly category.
        dict3: Phi and psi angle information for residues in Pro category.
        dict2: Phi and psi angle information for residues in Pre-Pro category.

    Examples:

        GeneralPhiPsiInfo, GlyPhiPsiInfo, ProPhiPsiInfo, PreProPhiPsiInfo =
            GetPhiPsiCategoriesResiduesInfo(MolName, ChainID)
        for ResNum in GeneralPhiPsiInfo["ResNums"]:
            ResName = GeneralPhiPsiInfo["ResName"][ResNum]
            Phi = GeneralPhiPsiInfo["Phi"][ResNum]
            Psi = GeneralPhiPsiInfo["Psi"][ResNum]
            MiscUtil.PrintInfo("ResNum: %s; ResName: %s;
                Phi: %8.2f; Psi: %8.2f" % (ResNum, ResName, Phi, Psi))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SelectionCmd = "%s and chain %s" % (MoleculeName, ChainName)
    GeneralPhiPsiInfo, GlyPhiPsiInfo, ProPhiPsiInfo, PreProPhiPsiInfo = _GetSelectionPhiPsiCategoriesResiduesInfo(SelectionCmd)

    return GeneralPhiPsiInfo, GlyPhiPsiInfo, ProPhiPsiInfo, PreProPhiPsiInfo

def GetSurfaceAndBuriedResiduesInfo(MoleculeName, ChainName, SASACutoff = 2.5):
    """Get information for surafce and buried residues present in a chain of a
    molecule. The surface residues correspond to residues with Solvent Accessible
    Surface Area (SASA) greater than or equal to the cutoff value. Otherwise, these
    residues are considered as buried residues.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        SASACutoff (float): SASA cutoff for heavy atoms corresponding to
            surface residues in chain. Units: Angstroms ** 2

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.
        dict2: Buried residues in a chain.

    Examples:

        SurfaceResiduesInfo, BurriedResiduesInfo =
            GetSurfaceResiduesInfo(MolName, ChainName, 2.5)
        for ResName in SurfaceResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None
    
    # Get ready to calculate solvent accessible surface area..
    CurrentDotSolvent = cmd.get("dot_solvent")
    cmd.set("dot_solvent", 1)

    # Setup tmp object name...
    TmpChainName = "Tmp_%s" % (MoleculeName)
    TmpChainSelection = "(%s and chain %s and polymer and not hydrogens)" % (MoleculeName, ChainName)
    cmd.create(TmpChainName, "(%s)" % (TmpChainSelection))

    # Calculate SASA for chain and load it into b...
    cmd.get_area(TmpChainName, load_b = 1)
    
    # Retrieve SASA for residues as B values...
    ResiduesBValuesInfoMap = _GetSelectionResiduesBValuesInfo(TmpChainName)

    # Retrieve information for surface and buried residues...
    SurfaceResiduesInfoMap = _GetResiduesInfoFromResiduesBValues(ResiduesBValuesInfoMap, ">=", SASACutoff)
    BuriedResiduesInfoMap = _GetResiduesInfoFromResiduesBValues(ResiduesBValuesInfoMap, "<", SASACutoff)
    
    # Delete tmp objects...
    cmd.delete(TmpChainName)
    
    # Restore current dot solvent...
    cmd.set("dot_solvent", CurrentDotSolvent)

    return SurfaceResiduesInfoMap, BuriedResiduesInfoMap

def GetInterfaceChainsResiduesByCAlphaAtomsDistance(MoleculeName1, ChainNames1, MoleculeName2, ChainNames2, DistanceCutoff = 8.0):
    """Get information for interface residues between chains in two
    molecules based on the distance between CAlpha atoms. The chain
    specification for molecules may contain multiple chain names delimited
    by commas.
    
    The interface residues are identified using PyMOL 'bycalpha' selection
    operator with in a specified distance.
    
    Arguments:
        MoleculeName1 (str): Name of a PyMOL molecule object.
        ChainNames1 (str): A chain name or comma delimited list of chain
            names in a molecule.
        MoleculeName2 (str): Name of a PyMOL molecule object.
        ChainNames2 (str): A chain name or comma delimited list of chain
            names in a molecule.
        DistanceCutoff (float): Distance cutoff for distance between
            any two CAlpha atoms in interface residues in different chains.

    Returns:
        dict1: Interface residues in a chain for first molecule. It is a
            dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.
        dict2: Interface residues in the chain for second molecule.

    Examples:

        ResiduesInfo1, ResiduesInfo2 = 
            GetInterfaceResiduesByHeavyAtomsDistance(MolName1,
            ChainName1, MolName2, ChainName2, DistanceCutoff)
        for ChainID in ResiduesInfo1["ChainIDs"]:
            for ResName in ResiduesInfo1["ResNames"][ChainID]:
                ResCount = ResiduesInfo1["ResCount"][ChainID][ResName]
                ResNums = ResiduesInfo1["ResNum"][ChainID][ResName]
                MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums:
                %s" % (ResName, ResCount, ResNums))

    """
    ChainNames1Selection, ChainNames2Selection = _ValidateAndSetupChainSelectionsForInterfaceResidues(MoleculeName1, ChainNames1, MoleculeName2, ChainNames2)
    
    if ChainNames1Selection is None or ChainNames2Selection is None:
        return None

    MoleculeSelection1 = "%s and chain %s" % (MoleculeName1, ChainNames1Selection)
    MoleculeSelection2 = "%s and chain %s" % (MoleculeName2, ChainNames2Selection)

    SelectionCmd1 = "((bycalpha (%s) within %.1f of  (%s)) and polymer)" % (MoleculeSelection1, DistanceCutoff, MoleculeSelection2)
    ResiduesInfoMap1= _GetSelectionChainsAndResiduesInfo(SelectionCmd1)

    SelectionCmd2 = "((bycalpha (%s) within %.1f of  (%s)) and polymer)" % (MoleculeSelection2, DistanceCutoff, MoleculeSelection1)
    ResiduesInfoMap2= _GetSelectionChainsAndResiduesInfo(SelectionCmd2)

    return ResiduesInfoMap1, ResiduesInfoMap2

def GetnterfaceChainsResiduesByHeavyAtomsDistance(MoleculeName1, ChainNames1, MoleculeName2, ChainNames2, DistanceCutoff = 5.0):
    """Get information for interface residues between chains in two
    molecules based on the distance between heavy atoms. The chain
    specification for molecules may contain multiple chain names delimited
    by commas.
    
    The interface residues are identified using PyMOL 'byresidue' selection
    operator with in a specified distance.
    
    Arguments:
        MoleculeName1 (str): Name of a PyMOL molecule object.
        ChainNames1 (str): A chain name or comma delimited list of chain
            names in a molecule.
        MoleculeName2 (str): Name of a PyMOL molecule object.
        ChainNames2 (str): A chain name or comma delimited list of chain
            names in a molecule.
        DistanceCutoff (float): Distance cutoff for distance between
            any two heavy atoms in interface residues in different chains.

    Returns:
        dict1: Interface residues in a chain for first molecule. It is a
            dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.
        dict2: Interface residues in the chain for second molecule.

    Examples:

        ResiduesInfo1, ResiduesInfo2 = 
            GetInterfaceResiduesByHeavyAtomsDistance(MolName1,
            ChainName1, MolName2, ChainName2, DistanceCutoff)
        for ChainID in ResiduesInfo1["ChainIDs"]:
            for ResName in ResiduesInfo1["ResNames"][ChainID]:
                ResCount = ResiduesInfo1["ResCount"][ChainID][ResName]
                ResNums = ResiduesInfo1["ResNum"][ChainID][ResName]
                MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums:
                %s" % (ResName, ResCount, ResNums))

    """
    ChainNames1Selection, ChainNames2Selection = _ValidateAndSetupChainSelectionsForInterfaceResidues(MoleculeName1, ChainNames1, MoleculeName2, ChainNames2)
    
    if ChainNames1Selection is None or ChainNames2Selection is None:
        return None

    MoleculeSelection1 = "%s and chain %s and (not hydrogens)" % (MoleculeName1, ChainNames1Selection)
    MoleculeSelection2 = "%s and chain %s and (not hydrogens)" % (MoleculeName2, ChainNames2Selection)

    SelectionCmd1 = "((byresidue (%s) within %.1f of  (%s)) and polymer)" % (MoleculeSelection1, DistanceCutoff, MoleculeSelection2)
    ResiduesInfoMap1= _GetSelectionChainsAndResiduesInfo(SelectionCmd1)

    SelectionCmd2 = "((byresidue (%s) within %.1f of  (%s)) and polymer)" % (MoleculeSelection2, DistanceCutoff, MoleculeSelection1)
    ResiduesInfoMap2= _GetSelectionChainsAndResiduesInfo(SelectionCmd2)

    return ResiduesInfoMap1, ResiduesInfoMap2

def GetInterfaceChainsResiduesBySASAChange(MoleculeName1, ChainNames1, MoleculeName2, ChainNames2, ChangeCutoff = 0.75):
    """Get information for interface residues between chains in two
    molecules based on the change in solvent accessible surface area
    (SASA) of a residue in chains in a molecule and chain complex
    containing specified chains across both molecules. The chain 
    specification for molecules may contain multiple chain names
    delimited by commas.

    Arguments:
        MoleculeName1 (str): Name of a PyMOL molecule object.
        ChainNames1 (str): A chain name or comma delimited list of chain
            names in a molecule.
        MoleculeName2 (str): Name of a PyMOL molecule object.
        ChainNames2 (str): A chain name or comma delimited list of chain
            names in a molecule.
        ChangeCutoff (float): SASA change cutoff for heavy atoms in
            a interface residue between an individual chain and a complex
            containing both chains.  Units: Angstroms ** 2

    Returns:
        dict1: Interface residues in the chain for first molecule. It is a
            dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.
        dict2: Interface residues in the chain for second molecule.

    Examples:

        ResiduesInfo1, ResiduesInfo2 = 
            GetInterfaceResiduesBySASAChange(MolName1,
            ChainName1, MolName2, ChainName2, DistanceCutoff)
        for ResName in ResiduesInfo1["ResNames"]:
            ResCount = ResiduesInfo1["ResCount"][ResName]
            ResNums = ResiduesInfo1["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    ChainNames1Selection, ChainNames2Selection = _ValidateAndSetupChainSelectionsForInterfaceResidues(MoleculeName1, ChainNames1, MoleculeName2, ChainNames2)
    
    if ChainNames1Selection is None or ChainNames2Selection is None:
        return None

    # Get ready to calculate solvent accessible surface area..
    CurrentDotSolvent = cmd.get("dot_solvent")
    cmd.set("dot_solvent", 1)

    # Setup tmp object names...
    TmpComplexName = "Tmp1_%s_%s" % (MoleculeName1, MoleculeName2)
    TmpChainName1 = "Tmp2_%s" % (MoleculeName1)
    TmpChainName2 = "Tmp3_%s" % (MoleculeName2)
    
    # Setup complex...
    TmpComplexSelection = "(%s and chain %s and polymer and not hydrogens) or (%s and chain %s and polymer and not hydrogens)" % (MoleculeName1, ChainNames1Selection, MoleculeName2, ChainNames2Selection)
    cmd.create(TmpComplexName, "(%s)" % (TmpComplexSelection))

    # Calculate SASA for complex and load it into b...
    cmd.get_area(TmpComplexName, load_b = 1)
    cmd.alter(TmpComplexName, "q=b")
    
    # Setup individual chains and calculate SASA...
    TmpChainNameSelection1 = "%s and (chain %s)" % (TmpComplexName, ChainNames1Selection)
    cmd.extract(TmpChainName1, "(%s)" % TmpChainNameSelection1)
    cmd.get_area(TmpChainName1, load_b = 1)

    TmpChainNameSelection2 = "%s and (chain %s)" % (TmpComplexName, ChainNames2Selection)
    cmd.extract(TmpChainName2, "(%s)" % TmpChainNameSelection2)
    cmd.get_area(TmpChainName2, load_b = 1)

    # Calculate SASA difference for individual chains...
    TmpChainNamesSelection = "(%s or %s)" % (TmpChainName1, TmpChainName2)
    cmd.alter(TmpChainNamesSelection, "b=b-q")

    ResiduesInfoMap1 = _GetChainsResiduesInfoForSASAChange(TmpChainNamesSelection, TmpChainName1, ChangeCutoff)
    ResiduesInfoMap2 = _GetChainsResiduesInfoForSASAChange(TmpChainNamesSelection, TmpChainName2, ChangeCutoff)

    # Delete tmp objects...
    cmd.delete(TmpComplexName)
    cmd.delete(TmpChainName1)
    cmd.delete(TmpChainName2)
    
    # Restore current dot solvent...
    cmd.set("dot_solvent", CurrentDotSolvent)
    
    return ResiduesInfoMap1, ResiduesInfoMap2

def _GetChainsResiduesInfoForSASAChange(SelectionCmd, MoleculeName, ChangeCutoff):
    """Get residue info for SASA change. """

    # Retrieve data...
    stored.ResiduesInfo = []
    cmd.iterate(SelectionCmd, "stored.ResiduesInfo.append([model, chain, resi, resn, b])")

    # Calculate SASA change for each residue...
    SASAChangeMap = {}
    SASAChangeMap["ChainIDs"] = []
    SASAChangeMap["ResIDs"] = {}
    SASAChangeMap["ResSASAChange"] = {}
    
    for Model, ChainID, ResNum, ResName, AtomSASAChange in stored.ResiduesInfo:
        if not re.match(Model, MoleculeName, re.I):
            continue
        ResID = "%s_%s" % (ResName, ResNum)
        
        if not ChainID in SASAChangeMap["ChainIDs"]:
            # Track new chain ID and intialize SASA change information...
            SASAChangeMap["ChainIDs"].append(ChainID)
            SASAChangeMap["ResIDs"][ChainID] = []
            SASAChangeMap["ResSASAChange"][ChainID]= {}

        if ResID in SASAChangeMap["ResSASAChange"][ChainID]:
            SASAChangeMap["ResSASAChange"][ChainID][ResID] += AtomSASAChange
        else:
            # Track new residue ID...
            SASAChangeMap["ResIDs"][ChainID].append(ResID)
            SASAChangeMap["ResSASAChange"][ChainID][ResID] = AtomSASAChange

    # Identify residues with SASA change greater than the specified cutoff value...
    SelectionInfoMap = {}
    SelectionInfoMap["ChainIDs"] = []
    
    SelectionInfoMap["ResNames"] = {}
    SelectionInfoMap["ResNum"] = {}
    SelectionInfoMap["ResCount"] = {}

    for ChainID in SASAChangeMap["ChainIDs"]:
        for ResID in SASAChangeMap["ResIDs"][ChainID]:
            ResSASAChange = SASAChangeMap["ResSASAChange"][ChainID][ResID]
            if abs(ResSASAChange) < ChangeCutoff:
                continue
            ResName, ResNum = ResID.split("_")
            
            if not ChainID in SelectionInfoMap["ChainIDs"]:
                # Track new chain ID and intialize residues information...
                SelectionInfoMap["ChainIDs"].append(ChainID)
                SelectionInfoMap["ResNames"][ChainID] = []
                SelectionInfoMap["ResNum"][ChainID]= {}
                SelectionInfoMap["ResCount"][ChainID] = {}
            
            if ResName in SelectionInfoMap["ResNames"][ChainID]:
                if not ResNum in SelectionInfoMap["ResNum"][ChainID][ResName]:
                    # Same residue name but different residue number
                    SelectionInfoMap["ResNum"][ChainID][ResName].append(ResNum)
                    SelectionInfoMap["ResCount"][ChainID][ResName] += 1
            else:
                # Track new residue information...
                SelectionInfoMap["ResNames"][ChainID].append(ResName)
                
                SelectionInfoMap["ResNum"][ChainID][ResName] = []
                SelectionInfoMap["ResNum"][ChainID][ResName].append(ResNum)
                
                SelectionInfoMap["ResCount"][ChainID][ResName] = 1
    
    return SelectionInfoMap

def _ValidateAndSetupChainSelectionsForInterfaceResidues(MoleculeName1, ChainNames1, MoleculeName2, ChainNames2):
    """Validate and setup chains name selectons for idenfitying interface residues."""
    
    ChainNames1 = re.sub(" ", "", ChainNames1)
    ChainNames2 = re.sub(" ", "", ChainNames2)
    
    if not (len(MoleculeName1) and len(ChainNames1) and len(MoleculeName2) and len(ChainNames2)):
        return None, None

    # Setup chain names for PyMOL selections...
    ChainNames1Selection = "+".join(ChainNames1.split(","))
    ChainNames2Selection = "+".join(ChainNames2.split(","))
    
    return ChainNames1Selection, ChainNames2Selection
    
def _GetSelectionChainsAndResiduesInfo(SelectionCmd):
    """Get chain names, residue names and count information for a selection. """
    
    # Retrieve atoms...
    stored.SelectionInfo = []
    cmd.iterate(SelectionCmd, "stored.SelectionInfo.append([chain, resi, resn])")

    # Retrieve chains and residues...
    SelectionInfoMap = {}
    SelectionInfoMap["ChainIDs"] = []
    
    SelectionInfoMap["ResNames"] = {}
    SelectionInfoMap["ResNum"] = {}
    SelectionInfoMap["ResCount"] = {}

    for ChainID, ResNum, ResName in stored.SelectionInfo:
        if not ChainID in SelectionInfoMap["ChainIDs"]:
            # Track new chain ID and intialize residues information...
            SelectionInfoMap["ChainIDs"].append(ChainID)
            SelectionInfoMap["ResNames"][ChainID] = []
            SelectionInfoMap["ResNum"][ChainID]= {}
            SelectionInfoMap["ResCount"][ChainID] = {}
            
        if ResName in SelectionInfoMap["ResNames"][ChainID]:
            if not ResNum in SelectionInfoMap["ResNum"][ChainID][ResName]:
                # Same residue name but different residue number
                SelectionInfoMap["ResNum"][ChainID][ResName].append(ResNum)
                SelectionInfoMap["ResCount"][ChainID][ResName] += 1
        else:
            # Track new residue information...
            SelectionInfoMap["ResNames"][ChainID].append(ResName)
            
            SelectionInfoMap["ResNum"][ChainID][ResName] = []
            SelectionInfoMap["ResNum"][ChainID][ResName].append(ResNum)
            
            SelectionInfoMap["ResCount"][ChainID][ResName] = 1
                
    return SelectionInfoMap

def _GetSelectionResiduesInfo(SelectionCmd):
    """Get residue names and count information for a selection. """
    
    # Retrieve atoms...
    stored.SelectionInfo = []
    cmd.iterate(SelectionCmd, "stored.SelectionInfo.append([resi, resn])")

    # Retrieve residues...
    SelectionInfoMap = {}
    SelectionInfoMap["ResNames"] = []
    SelectionInfoMap["ResNum"] = {}
    SelectionInfoMap["ResCount"] = {}

    for ResNum, ResName in stored.SelectionInfo:
        if ResName in SelectionInfoMap["ResNames"]:
            if not ResNum in SelectionInfoMap["ResNum"][ResName]:
                # Same residue name but different residue number
                SelectionInfoMap["ResNum"][ResName].append(ResNum)
                SelectionInfoMap["ResCount"][ResName] += 1
        else:
            SelectionInfoMap["ResNames"].append(ResName)
            
            SelectionInfoMap["ResNum"][ResName] = []
            SelectionInfoMap["ResNum"][ResName].append(ResNum)
            
            SelectionInfoMap["ResCount"][ResName] = 1
                
    return SelectionInfoMap

def _GetSelectionPhiPsiResiduesInfo(SelectionCmd, Categorize):
    """Get phi and psi torsion angle information for residues in a  selection
    with in a chain.
     
    The phi and psi angles are optionally categorized into the following groups
    corresponding to four types of  Ramachandran plots:
    
    General: All residues except glycine, proline, or pre-proline
    Glycine: Only glycine residues
    Proline: Only proline residues
    Pre-Proline: Only residues before proline not including glycine or proline
    
    """

    # Initialize...
    SelectionInfoMap = {}
    SelectionInfoMap["ResNums"] = []
    SelectionInfoMap["ResName"] = {}
    SelectionInfoMap["Phi"] = {}
    SelectionInfoMap["Psi"] = {}
    SelectionInfoMap["Category"] = {}

    # Retrieve phi and psi info...
    PhiPsiInfo = _GetSelectionPhiPsiInfo(SelectionCmd, Categorize)
    if PhiPsiInfo is None:
        return SelectionInfoMap

    # Track...
    for Model, ChainID, ResNum, ResName, Phi, Psi, Category in PhiPsiInfo:
        if ResNum in SelectionInfoMap["ResName"]:
            continue
        
        SelectionInfoMap["ResNums"].append(ResNum)
        SelectionInfoMap["ResName"][ResNum] = ResName
        SelectionInfoMap["Phi"][ResNum] = Phi
        SelectionInfoMap["Psi"][ResNum] = Psi
        SelectionInfoMap["Category"][ResNum] = Category
    
    return SelectionInfoMap

def _GetSelectionPhiPsiChainsAndResiduesInfo(SelectionCmd, Categorize):
    """Get phi and psi torsion angle information for residues in a  selection
    across chains in a molecule.
    
    The phi and psi angles are optionally categorized into the following groups
    corresponding to four types of  Ramachandran plots:
    
    General: All residues except glycine, proline, or pre-proline
    Glycine: Only glycine residues
    Proline: Only proline residues
    Pre-Proline: Only residues before proline not including glycine or proline
    
    """

    # Initialize...
    SelectionInfoMap = {}
    SelectionInfoMap["ChainIDs"] = []
    SelectionInfoMap["ResNums"] = {}
    SelectionInfoMap["ResName"] = {}
    SelectionInfoMap["Phi"] = {}
    SelectionInfoMap["Psi"] = {}
    SelectionInfoMap["Category"] = {}

    # Retrieve phi and psi info...
    PhiPsiInfo = _GetSelectionPhiPsiInfo(SelectionCmd, Categorize)
    if PhiPsiInfo is None:
        return SelectionInfoMap

    # Track...
    for Model, ChainID, ResNum, ResName, Phi, Psi, Category in PhiPsiInfo:
        if not ChainID in SelectionInfoMap["ResNums"]:
            # Track new chain ID and initialize residues information...
            SelectionInfoMap["ChainIDs"].append(ChainID)
            
            SelectionInfoMap["ResNums"][ChainID] = []
            SelectionInfoMap["ResName"][ChainID] = {}
            SelectionInfoMap["Phi"][ChainID] = {}
            SelectionInfoMap["Psi"][ChainID] = {}
            SelectionInfoMap["Category"][ChainID] = {}
            
        # Check for duplicates...
        if ResNum in SelectionInfoMap["ResName"][ChainID]:
            continue

        # Track new residues information...
        SelectionInfoMap["ResNums"][ChainID].append(ResNum)
        SelectionInfoMap["ResName"][ChainID][ResNum] = ResName
        SelectionInfoMap["Phi"][ChainID][ResNum] = Phi
        SelectionInfoMap["Psi"][ChainID][ResNum] = Psi
        SelectionInfoMap["Category"][ChainID][ResNum] = Category
    
    return SelectionInfoMap

def _GetSelectionPhiPsiCategoriesResiduesInfo(SelectionCmd):
    """Get phi and psi torsion angle information for residues in a  selection
    with in a chain.

    The phi and psi angles are optionally categorized into the following groups
    corresponding to four types of  Ramachandran plots:
    
    General: All residues except glycine, proline, or pre-proline
    Glycine: Only glycine residues
    Proline: Only proline residues
    Pre-Proline: Only residues before proline not including glycine or proline
    
    """
    
    # Initialize...
    PhiPsiInfoMaps = {}
    PhiPsiInfoMaps["Maps"] = []
    for Category in ["General", "Glycine", "Proline", "Pre-Proline"]:
        PhiPsiInfoMap = {}
        PhiPsiInfoMap["ResNums"] = []
        PhiPsiInfoMap["ResName"] = {}
        PhiPsiInfoMap["Phi"] = {}
        PhiPsiInfoMap["Psi"] = {}
        PhiPsiInfoMap["Category"] = {}
        
        PhiPsiInfoMaps[Category] = PhiPsiInfoMap
        PhiPsiInfoMaps["Maps"].append(PhiPsiInfoMap)
        
    # Retrieve phi and psi info...
    PhiPsiInfo = _GetSelectionPhiPsiInfo(SelectionCmd, True)
    if PhiPsiInfo is None:
        return  PhiPsiInfoMaps["Maps"]
    
    Index = -1
    for Model, ChainID, ResNum, ResName, Phi, Psi, Category in PhiPsiInfo:
        Index += 1
        
        if ResNum in PhiPsiInfoMaps[Category]["ResName"]:
            continue
        
        # Track...
        PhiPsiInfoMaps[Category]["ResNums"].append(ResNum)
        PhiPsiInfoMaps[Category]["ResName"][ResNum] = ResName
        PhiPsiInfoMaps[Category]["Phi"][ResNum] = Phi
        PhiPsiInfoMaps[Category]["Psi"][ResNum] = Psi
        PhiPsiInfoMaps[Category]["Category"][ResNum] = Category
    
    return  PhiPsiInfoMaps["Maps"]

def _GetSelectionPhiPsiInfo(SelectionCmd, Categorize = False):
    """Get phi and psi angles for a selection along with information regarding
    model, chain, residue number, and residue name.
    
    The phi and psi angles are optionally categorized into the following groups
    corresponding to four types of  Ramachandran plots:
    
    General: All residues except glycine, proline, or pre-proline
    Glycine: Only glycine residues
    Proline: Only proline residues
    Pre-Proline: Only residues before proline not including glycine or proline
    
    """

    # Retrieve phi and psi values...
    PhiPsiInfoMap = cmd.get_phipsi("(%s)" % SelectionCmd)
    if PhiPsiInfoMap is None:
        return None

    # Sort keys corresponding to model name and residue numbers...
    PhiPsiKeys = sorted(PhiPsiInfoMap.keys())
    
    # Retrieve related information...
    PhiPsiInfoList = []

    Category = None
    for Key in PhiPsiKeys:
        Model, Index = Key
        Phi, Psi = PhiPsiInfoMap[Key]
        
        stored.PhiPsiInfoList = []
        cmd.iterate("(%s`%d)" % (Model, Index), "stored.PhiPsiInfoList.append([model, chain, resi, resn])")
        for Model, Chain, ResNum, ResName in stored.PhiPsiInfoList:
            PhiPsiInfoList.append([Model, Chain, ResNum, ResName, Phi, Psi, Category])

    if Categorize:
        _CategorizePhiPsiAnglesInfo(PhiPsiInfoList)
    
    return PhiPsiInfoList

def _CategorizePhiPsiAnglesInfo(PhiPsiInfo):
    """The phi and psi angles are optionally categorized into the following groups
    corresponding to four types of  Ramachandran plots:
    
    General: All residues except glycine, proline, or pre-proline
    Glycine: Only glycine residues
    Proline: Only proline residues
    Pre-Proline: Only residues before proline not including glycine or proline
    
    """
    
    Index = -1
    LastIndex = len(PhiPsiInfo) - 1
    PhiPsiCategories = []
    for Model, ChainID, ResNum, ResName, Phi, Psi, CategoryPlaceHolder in PhiPsiInfo:
        Index += 1
        if re.match("^Gly$", ResName, re.I):
            # Gly: Only glycine residues
            Category = "Glycine"
        elif re.match("^Pro$", ResName, re.I):
            # Pro: Only proline residues
            Category = "Proline"
        elif _IsPreProlinePhiPsiAngle(PhiPsiInfo, Index, LastIndex):
            #  Pre-Proline: Only residues before proline not including glycine or proline
            Category = "Pre-Proline"
        else:
            # General: All residues except Gly, Pro, or pre-Pro
            Category = "General"
        # Track categories...
        PhiPsiCategories.append(Category)

    # Update category in PhiPsiInfo...
    CategoryIndex = 6
    for Index, Category in enumerate(PhiPsiCategories):
         PhiPsiInfo[Index][CategoryIndex] = Category

def _IsPreProlinePhiPsiAngle(PhiPsiInfo, Index, LastIndex):
    """Check for Pre-Proline phi and psi angles."""
    
    #  Pre-Proline: Only residues before proline not including glycine or proline

    # Not a last residue...
    NextIndex = Index + 1
    if NextIndex  >= LastIndex:
        return False

    # Next residue is proline...
    if not re.match("^Pro$", PhiPsiInfo[NextIndex][3], re.I):
        return False

    # Current residue is not Gly or Pro...
    NextResName = PhiPsiInfo[Index][3]
    if re.match("^(Gly|Pro)$", NextResName, re.I):
        return False

    # Next chain ID is same as current chain ID...
    CurrentChainID =  PhiPsiInfo[Index][1]
    NextChainID = PhiPsiInfo[NextIndex][1]
    if not re.match("^%s$" % CurrentChainID, NextChainID):
        return False

    # Current and next residue numbers are sequential...
    CurrentResNum = int(PhiPsiInfo[Index][2])
    NextResNum = int(PhiPsiInfo[NextIndex][2])
    if not (NextResNum - CurrentResNum == 1):
        return False
    
    return True
    
def _GetSelectionResiduesBValuesInfo(SelectionCmd):
    """Get B values info for residues in a chain. """

    # Retrieve data...
    stored.ResiduesInfo = []
    cmd.iterate(SelectionCmd, "stored.ResiduesInfo.append([resi, resn, b])")

    # Setup B-values for each residue...
    BValuesInfoMap = {}
    BValuesInfoMap["ResIDs"] = []
    BValuesInfoMap["BValue"] = {}
    
    for ResNum, ResName, AtomBValue in stored.ResiduesInfo:
        ResID = "%s_%s" % (ResName, ResNum)
        
        if ResID in BValuesInfoMap["BValue"]:
            BValuesInfoMap["BValue"][ResID] += AtomBValue
        else:
            # Track new residue ID...
            BValuesInfoMap["ResIDs"].append(ResID)
            BValuesInfoMap["BValue"][ResID] = AtomBValue

    return BValuesInfoMap

def _GetResiduesInfoFromResiduesBValues(ResiduesBValuesInfoMap, Mode, Cutoff):
    """Get residues info map using residues B values. """

    GreaterThanEquals, GreaterThan, LessThanEquals, LessThan, Equals = [False] * 5
    if re.match("^>=$", Mode, re.I):
        GreaterThanEquals = True
    elif re.match("^>$", Mode, re.I):
        GreaterThan = True
    elif re.match("^<=$", Mode, re.I):
        LessThanEquals = True
    elif re.match("^<$", Mode, re.I):
        LessThan = True
    elif re.match("^=$", Mode, re.I):
        Equals = True
    else:
        MiscUtil.PrintError("PyMOLUtil._GetResiduesInfoFromResiduesBValues: The value, %s, specified for Mode is not valid. Supported values: >=, >, <=, <, =")
    
    # Setup residues info map...
    ResiduesInfoMap = {}
    ResiduesInfoMap["ResNames"] = []
    ResiduesInfoMap["ResNum"] = {}
    ResiduesInfoMap["ResCount"] = {}

    for ResID in ResiduesBValuesInfoMap["ResIDs"]:
        ResName, ResNum = ResID.split("_")
        ResBValue = ResiduesBValuesInfoMap["BValue"][ResID]

        if GreaterThanEquals:
            if ResBValue < Cutoff:
                continue
        elif GreaterThan:
            if ResBValue <= Cutoff:
                continue
        elif LessThanEquals:
            if ResBValue > Cutoff:
                continue
        elif LessThan:
            if ResBValue >= Cutoff:
                continue
        elif Equals:
            if ResBValue != Cutoff:
                continue
        
        if ResName in ResiduesInfoMap["ResNames"]:
            if not ResNum in ResiduesInfoMap["ResNum"][ResName]:
                # Same residue name but different residue number
                ResiduesInfoMap["ResNum"][ResName].append(ResNum)
                ResiduesInfoMap["ResCount"][ResName] += 1
        else:
            ResiduesInfoMap["ResNames"].append(ResName)
            
            ResiduesInfoMap["ResNum"][ResName] = []
            ResiduesInfoMap["ResNum"][ResName].append(ResNum)
            
            ResiduesInfoMap["ResCount"][ResName] = 1
    
    return ResiduesInfoMap

def ProcessChainsAndLigandsOptionsInfo(ChainsAndLigandsInfo, ChainsOptionName, ChainsOptionValue, LigandsOptionName = None, LigandsOptionValue = None):
    """Process specified chain and ligand IDs using command line options.

    Arguments:
        ChainsAndLigandsInfo (dict): A dictionary containing information
            existing chains and ligands. 
        ChainsOptionName (str): Name of command line chains option.
        ChainsOptionValue (str): Value for command line chains option.
        LigandsOptionName (str): Name of command line ligands option.
        LigandsOptionValue (str): Value for command line ligands option.

    Returns:
        dict: A dictionary containing list of chain identifiers and dictionaries
            of chains containing lists of ligand names for each chain.

    Examples:

        ChainsAndLigandsInfo = ProcessChainsAndLigandsOptionsInfo(
            ChainsAndLigandsInfo, "-c, --chainIDs", OptionsInfo["ChainIDs"],
            "-l, --ligandIDs", OptionsInfo["LigandIDs"])
        for ChainID in ChainsAndLigandsInfo["ChainIDs"]:
            for LigandID in ChainsAndLigandsInfo["LigandIDs"][ChainID]:
                MiscUtil.PrintInfo("ChainID: %s; LigandID: %s" % (ChainID,
                    LigandID))

    """
    SpecifiedChainsAndLigandsInfo = {}
    SpecifiedChainsAndLigandsInfo["ChainIDs"] = []
    SpecifiedChainsAndLigandsInfo["LigandIDs"] = {}

    if ChainsOptionValue is None:
        return SpecifiedChainsAndLigandsInfo
        
    _ProcessChainIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo, ChainsOptionName, ChainsOptionValue)

    if LigandsOptionValue is None:
        return SpecifiedChainsAndLigandsInfo
    
    _ProcessLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo, LigandsOptionName, LigandsOptionValue)
    
    return SpecifiedChainsAndLigandsInfo

def _ProcessChainIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo, ChainsOptionName, ChainsOptionValue):
    """Process chain IDs"""

    MiscUtil.PrintInfo("Processing chain IDs...")
    
    if re.match("^All$", ChainsOptionValue, re.I):
        SpecifiedChainsAndLigandsInfo["ChainIDs"] = ChainsAndLigandsInfo["ChainIDs"]
        return
    elif re.match("^(First|Auto)$", ChainsOptionValue, re.I):
        FirstChainID = ChainsAndLigandsInfo["ChainIDs"][0] if (len(ChainsAndLigandsInfo["ChainIDs"])) else None
        if FirstChainID is not None:
            SpecifiedChainsAndLigandsInfo["ChainIDs"].append(FirstChainID)
        return
    
    ChainIDs = re.sub(" ", "", ChainsOptionValue)
    if not ChainIDs:
        MiscUtil.PrintError("No valid value specified using \"%s\" option." % ChainsOptionName)

    ChainIDsList = ChainsAndLigandsInfo["ChainIDs"]
    SpecifiedChainIDsList = []
    
    ChainIDsWords = ChainIDs.split(",")
    for ChainID in ChainIDsWords:
        if not ChainID in ChainIDsList:
            MiscUtil.PrintWarning("The chain ID, %s, specified using \"%s\" option is not valid. It'll be ignored. Valid chain IDs: %s" % (ChainID, ChainsOptionName, ", ".join(ChainIDsList)))
            continue
        if ChainID in SpecifiedChainIDsList:
            MiscUtil.PrintWarning("The chain ID, %s, has already been specified using \"%s\" option. It'll be ignored." % (ChainID, ChainsOptionName))
            continue
        SpecifiedChainIDsList.append(ChainID)
    
    if not len(SpecifiedChainIDsList):
        MiscUtil.PrintError("No valid chain IDs \"%s\"  specified using \"%s\" option." % (ChainsOptionValue, ChainsOptionName))
    
    SpecifiedChainsAndLigandsInfo["ChainIDs"] = SpecifiedChainIDsList
    
def _ProcessLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo, LigandsOptionName, LigandsOptionValue):
    """Process ligand IDs"""

    MiscUtil.PrintInfo("Processing ligand IDs...")
    
    # Intialize ligand IDs...
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
        SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID] = []
        
    if re.match("^All$", LigandsOptionValue, re.I):
        for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"] :
            SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID] = ChainsAndLigandsInfo["LigandIDs"][ChainID]
        return
    elif re.match("^(Largest|Auto)$", LigandsOptionValue, re.I):
        # Setup largest ligand ID for each chain...
        for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"] :
            LargestLigandID = ChainsAndLigandsInfo["LigandIDs"][ChainID][0] if (len(ChainsAndLigandsInfo["LigandIDs"][ChainID])) else None
            if LargestLigandID is not None:
                SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID].append(LargestLigandID)
        return
    
    LigandIDs = re.sub(" ", "", LigandsOptionValue)
    if not LigandIDs:
        MiscUtil.PrintError("No valid value specified using \"%s\" option." % LigandsOptionName)
    
    LigandIDsWords = LigandIDs.split(",")
    
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
        LigandIDsList = ChainsAndLigandsInfo["LigandIDs"][ChainID]
        SpecifiedLigandIDsList = []

        for LigandID in LigandIDsWords:
            if not LigandID in LigandIDsList:
                LigandIDsListNames = ",".join(LigandIDsList) if len(LigandIDsList) else "None"
                MiscUtil.PrintWarning("The ligand ID, %s, specified using \"%s\" option is not valid for chain, %s. It'll be ignored. Valid ligand IDs: %s" % (LigandID, LigandsOptionName, ChainID, LigandIDsListNames))
                continue
            if LigandID in SpecifiedLigandIDsList:
                MiscUtil.PrintWarning("The ligand ID, %s, has already been specified using \"%s\" option. It'll be ignored." % (LigandID, LigandsOptionName))
                continue
            SpecifiedLigandIDsList.append(LigandID)
            
        if not len(SpecifiedLigandIDsList):
            MiscUtil.PrintWarning("No valid ligand IDs \"%s\" specified using \"%s\" option for chain ID, %s." % (LigandsOptionValue, LigandsOptionName, ChainID))
        
        SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID] = SpecifiedLigandIDsList

def ProcessResidueTypesOptionsInfo(ResidueTypesOptionName, ResidueTypesOptionValue):
    """Process specified residue types using command line option.

    Arguments:
        ResidueTypesOptionName (str): Name of command line  option.
        ResidueTypesOptionValue (str): Value for command line option.

    Returns:
        list: A list containing names of valid residue types.
        dict: A dictionary containing residue types pointing to dictionaries of
            color names and list of residues for a residue type.

    Examples:

        ResidueTypesNamesInfo, ResidueTypesParamsInfo =
            ProcessChainsAndLigandsOptionsInfo("-r, --residueTypes",
            OptionsInfo["ResidueTypes"])
        for ResidueTypeName in ResidueTypesNamesInfo:
            MiscUtil.PrintInfo("ResidueType: %s; Color: %s; Residues: %s" %
            (ResidueTypeName, ResidueTypeName[ResidueTypeName]["Color"],
            " ".join(ResidueTypeName[ResidueTypeName]["Residues"]))

    """
    
    # Set up  default values for residue types, colors, and names.
    ResidueTypesNamesInfo = ["Aromatic", "Hydrophobic", "Polar", "Positively_Charged", "Negatively_Charged"]
    
    ResidueTypesParamsInfo = {}
    ResidueTypesParamsInfo["Aromatic"] = {"Color": "brightorange", "Residues" : ["HIS", "PHE", "TRP", "TYR"]}
    ResidueTypesParamsInfo["Hydrophobic"] = {"Color": "orange", "Residues" : ["ALA", "GLY", "VAL", "LEU", "ILE", "PRO", "MET"]}
    ResidueTypesParamsInfo["Polar"] = {"Color": "palegreen", "Residues" : ["ASN", "GLN", "SER", "THR", "CYS"]}
    ResidueTypesParamsInfo["Positively_Charged"] = {"Color": "marine", "Residues" : ["ARG", "LYS"]}
    ResidueTypesParamsInfo["Negatively_Charged"] = {"Color": "red", "Residues" : ["ASP", "GLU"]}
    
    ResidueTypes = ResidueTypesOptionValue
    if re.match("^auto$", ResidueTypesOptionValue, re.I):
        _SetupOtherResidueTypes(ResidueTypesParamsInfo)
        return ResidueTypesNamesInfo, ResidueTypesParamsInfo
    
    # Parse specified residue types...
    ResidueTypesWords = ResidueTypes.split(",")
    if len(ResidueTypesWords) % 3:
        MiscUtil.PrintError("The number of comma delimited residue type, color, and name triplets, %d, specified using \"%s\" option must be a multple of 3." % (len(ResidueTypesWords), ResidueTypesOptionName))

    # Set up canonical residue type names...
    ValidResidueTypeNames = []
    CanonicalResidueTypeNamesMap = {}
    for Name in sorted(ResidueTypesNamesInfo):
        ValidResidueTypeNames.append(Name)
        CanonicalResidueTypeNamesMap[Name.lower()] = Name

    # Validate and set residue types, colors, and names..
    for Index in range(0, len(ResidueTypesWords), 3):
        TypeName = ResidueTypesWords[Index].strip()
        ResidueTypeColor = ResidueTypesWords[Index + 1].strip()
        ResidueNames = ResidueTypesWords[Index + 2].strip()
        
        ResidueNames = re.sub("[ ]+", " ", ResidueNames)
        ResidueNamesWords = ResidueNames.split(" ")

        CanonicalTypeName = TypeName.lower()
        if not CanonicalTypeName in CanonicalResidueTypeNamesMap:
            MiscUtil.PrintError("The residue type, %s, specified using \"%s\" option is not a valid type. Supported residue types: %s" % (TypeName, ResidueTypesOptionName, ", ".join(ValidResidueTypeNames)))
        ResidueTypeName = CanonicalResidueTypeNamesMap[CanonicalTypeName]

        if not ResidueTypeColor:
            MiscUtil.PrintError("No color name specified for residue type, %s, using \"%s\" option, %s" % (TypeName, ResidueTypesOptionName, ResidueTypes))
        
        if not ResidueNames:
            MiscUtil.PrintError("No residue names specified for residue type, %s, using \"%s\" option, %s" % (TypeName, ResidueTypesOptionName, ResidueTypes))
            
        ResidueTypesParamsInfo[ResidueTypeName]["Color"] = ResidueTypeColor
        ResidueTypesParamsInfo[ResidueTypeName]["Residues"] = ResidueNamesWords
        
        SetupOtherResidueTypes()
        return ResidueTypesNamesInfo, ResidueTypesParamsInfo

def _SetupOtherResidueTypes(ResidueTypesParamsInfo):
    """Setup other residue types."""

    # Set other residues to all specified residues. The other residues are selected by
    # performing a negation operation on this list during the creation of PyMOL objects.
    
    ResidueNames = []
    for ResiduesType in ResidueTypesParamsInfo:
        ResidueNames.extend(ResidueTypesParamsInfo[ResiduesType]["Residues"] )
    
    ResidueTypesParamsInfo["Other"]= {}
    ResidueTypesParamsInfo["Other"]["Color"] = None
    ResidueTypesParamsInfo["Other"]["Residues"] = ResidueNames
    
def ProcessSurfaceAtomTypesColorsOptionsInfo(ColorOptionName, ColorOptionValue):
    """Process specified surafce atom types colors using command line option.

    Arguments:
        ColorOptionName (str): Name of command line  option.
        ColorOptionValue (str): Value for command line option.

    Returns:
        dict: A dictionary containing atom types and colors.

    Examples:

        AtomTypesColorNamesInfo =
            PyMOLUtil.ProcessSurfaceAtomTypesColorsOptionsInfo(
            "--surfaceAtomTypesColors", OptionsInfo["SurfaceAtomTypesColors"])

    """
    
    # Set up  default values for atom type colors...
    AtomTypesColorNamesInfo = {"HydrophobicAtomsColor": "yellow", "NegativelyChargedAtomsColor": "red", "PositivelyChargedAtomsColor": "blue", "OtherAtomsColor": "gray90"}
    
    AtomTypesColors = ColorOptionValue
    if re.match("^auto$", AtomTypesColors, re.I):
        return AtomTypesColorNamesInfo
    
    # Parse atom type colors and values...
    AtomTypesColorsWords = AtomTypesColors.split(",")
    if len(AtomTypesColorsWords) % 2:
        MiscUtil.PrintError("The number of comma delimited surface atom color types and values, %d, specified using \"%s\" option must be a multple of 2." % (len(AtomTypesColorsWords), ColorOptionName))

    # Set up canonical atom type colors...
    ValidAtomTypesColors = []
    CanonicalAtomTypesColorsMap = {}
    for Type in sorted(AtomTypesColorNamesInfo):
        ValidAtomTypesColors.append(Type)
        CanonicalAtomTypesColorsMap[Type.lower()] = Type
    
    # Validate and process specified values...
    for Index in range(0, len(AtomTypesColorsWords), 2):
        Type = AtomTypesColorsWords[Index].strip()
        Color = AtomTypesColorsWords[Index + 1].strip()
        
        Color = re.sub("[ ]+", " ", Color)
        ColorWords = Color.split(" ")
    
        CanonicalType = Type.lower()
        if not CanonicalType in CanonicalAtomTypesColorsMap:
            MiscUtil.PrintError("The surface atom color type, %s, specified using \"%s\" option  is not a valid type: Supported atom color types: %s" % (Type, ColorOptionName, ", ".join(ValidAtomTypesColors)))
        Type = CanonicalAtomTypesColorsMap[CanonicalType]

        if not Color:
            MiscUtil.PrintError("No color type specified for atom color type, %s, using \"%s\" option." % (Type, ColorOptionName))

        ColorWordsLen = len(ColorWords)
        if not (ColorWordsLen == 1 or ColorWordsLen == 3):
            MiscUtil.PrintError("The number, %s, of color name or space delimited RGB values, %s, specified for atom color type, %s, using \"%s\" option must be 1 or 3." % (ColorWordsLen, Color, ColorOptionName, Type))

        AtomTypesColorNamesInfo[Type] = " ".join(ColorWords)

        return AtomTypesColorNamesInfo

def ProcessChainSelectionsOptionsInfo(SelectionOptionName, SelectionOptionValue):
    """Process names and selections specified  using command line option. It
    is a pairwise list comma delimited values corresponding to PyMOL object
    names and selection specification

    Arguments:
        SelectionOptionName (str): Name of command line  option.
        SelectionOptionValue (str): Value for command line option.

    Returns:
        dict: A dictionary containing lists f names and selection commands.

    Examples:

        ChainSelectionsInfo =
            PyMOLUtil.ProcessChainSelectionsOptionsInfo("--selectionsChain",
            OptionsInfo["SelectionChains"])

    """
    
    # Initialize...
    ChainSelectionsInfo = {}
    ChainSelectionsInfo["Names"] = []
    ChainSelectionsInfo["Selections"] = []

    if re.match("^None$", SelectionOptionValue, re.I):
        return ChainSelectionsInfo
    
    # Parse selection chains names and selections...
    ChainSelectionsWords = SelectionOptionValue.split(",")
    if len(ChainSelectionsWords) % 2:
        MiscUtil.PrintError("The number of comma delimited selection chains names and selections, %d, specified using \"%s\" option must be a multple of 2." % (len(ChainSelectionsWords), SelectionOptionName))

    CanonicalNamesMap = {}
    for Index in range(0, len(ChainSelectionsWords), 2):
        Name = ChainSelectionsWords[Index].strip()
        Selection = ChainSelectionsWords[Index + 1].strip()

        if not len(Name):
            MiscUtil.PrintError("A name specified, \"%s\",  using \"%s\" option is empty." % (SelectionOptionValue, SelectionOptionName))
        if not len(Selection):
            MiscUtil.PrintError("A selection specified, \"%s\", using \"%s\" option is empty." % (SelectionOptionValue, SelectionOptionName))

        CanonicalName = Name.lower()
        if CanonicalName in CanonicalNamesMap:
            MiscUtil.PrintError("The name %s  specified using \"%s\" option is a duplicate name." % (Name, SelectionOptionName))
        CanonicalNamesMap[CanonicalName] = Name

        if re.search("[^a-zA-Z0-9 _\-]", Name, re.I):
            MiscUtil.PrintError("The name %s  specified using \"%s\" option contains invalid charactors. Supportted characters: alphanumeric,  space, hyphen and underscore.." )
        
        ChainSelectionsInfo["Names"].append(Name)
        ChainSelectionsInfo["Selections"].append(Selection)
    
    return ChainSelectionsInfo

def CalculateCenterOfMass(Selection = "all", Quiet = 0):
    """Calculate center of mass for a selection.

    Arguments:
        Selection (str): A PyMOL selection.
        Quiet (int): Print information.

    Returns:
        list: X, Y, Z coordinates for center of mass.

    """
    MassTotal = 0.0
    X, Y, Z = [0.0, 0.0, 0.0]
    
    Atoms = cmd.get_model(Selection)
    for Atom in Atoms.atom:
        Mass = Atom.get_mass()
        MassTotal += Mass

        X += Atom.coord[0] * Mass
        Y += Atom.coord[1] * Mass
        Z += Atom.coord[2] * Mass

    XCOM = X/MassTotal
    YCOM = Y/MassTotal
    ZCOM = Z/MassTotal

    if not Quiet:
        MiscUtil.PrintInfo("PyMOLUtil.CalculateCenterOfMass: %f, %f, %f" % (XCOM, YCOM, ZCOM))
    
    return [XCOM, YCOM, ZCOM]

def ConvertFileFormat(Infile, Outfile, Reinitialize = True, OutputFeedback = True):
    """Convert infile to outfile by automatically detecting their formats
    from the file extensions.

    The formats of both input and output files must be a valid format supported
    by PyMOL.

    Arguments:
        Infile (str): Name of input file.
        Outfile (str): Name of outfile file.
        Reinitialize (bool): Reinitialize PyMOL before loading input file.
        OutputFeedback (bool): Control output feedback.

    """
    
    if not os.path.exists(Infile):
        MiscUtil.PrintWarning("The input file, %s, doesn't exists.%s..." % (Infile))

    if Reinitialize:
        cmd.reinitialize()

    if not OutputFeedback:
        # Turn off output feedback...
        MiscUtil.PrintInfo("Disabling output feedback for PyMOL...")
        cmd.feedback("disable", "all", "output")

    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Infile)
    MolName = FileName
    
    cmd.load(Infile, MolName)
    cmd.save(Outfile, MolName)
    cmd.delete(MolName)
    
    if not OutputFeedback:
        # Turn it back on...
        MiscUtil.PrintInfo("\nEnabling output feedback for PyMOL...")
        cmd.feedback("enable", "all", "output")

def ConvertPMLFileToPSEFile(PMLFile, PSEFile, Reinitialize = True, OutputFeedback = True):
    """Convert PML file to PME file.

    Arguments:
        PMLFile (str): Name of PML file.
        PSEFile (str): Name of PSE file.
        Reinitialize (bool): Reinitialize PyMOL before loading PML file.
        OutputFeedback (bool): Control output feedback.

    """
    
    if not os.path.exists(PMLFile):
        MiscUtil.PrintWarning("The PML file, %s, doesn't exists.%s..." % (PMLFile))

    if Reinitialize:
        cmd.reinitialize()

    if not OutputFeedback:
        # Turn off output feedback...
        MiscUtil.PrintInfo("Disabling output feedback for PyMOL...")
        cmd.feedback("disable", "all", "output")

    cmd.do("@%s" % PMLFile)
    cmd.save(PSEFile)
    
    if not OutputFeedback:
        # Turn it back on...
        MiscUtil.PrintInfo("\nEnabling output feedback for PyMOL...")
        cmd.feedback("enable", "all", "output")

def SetupPMLHeaderInfo(ScriptName = None, IncludeLocalPython = True):
    """Setup header information for generating PML files. The local Python
    functions are optionally embedded in the header information for their
    use in PML files.
 
    Arguments:
        ScriptName (str): Name of script calling the function.
        IncludeLocalPython (bool): Include local Python functions.

    Returns:
        str: Text containing header information for generating PML files.

    """
    if ScriptName is None:
        HeaderInfo = """\
#
# This file is automatically generated by a script available in MayaChemTools.
#
cmd.reinitialize()"""
    else:
        HeaderInfo = """\
#
# This file is automatically generated by the following PyMOL script available in
# MayaChemTools: %s
#
cmd.reinitialize() """ % (ScriptName)

    if IncludeLocalPython:
        PMLForLocalPython = _SetupPMLForLocalPython()
        HeaderInfo = "%s\n\n%s" % (HeaderInfo, PMLForLocalPython)
    
    return HeaderInfo

def _SetupPMLForLocalPython():
    """Setup local Python functions for PML file.
    """
    
    PMLForPython = """\
""
"Setting up local Python functions  for PML script..."
""
python

from __future__ import print_function
import re

def ColorByHydrophobicity(Selection, ColorPalette = "RedToWhite"):
    \"""Color by hydrophobicity using hydrophobic values for amino acid
    residues corresponding to the Eisenberg hydrophobicity scale.
    
    Possible values for ColorPalette: RedToWhite or WhiteToGreen from most
    hydrophobic amino acid to least hydrophobic.
    
    The colors values for amino acids are taken from color_h script avaiable
    as part of the Script Library at PyMOL Wiki. 
        
    \"""
        
    if not re.match("^(RedToWhite|WhiteToGreen)$", ColorPalette, re.I):
        print("Invalid ColorPalette value: %s. Valid values: RedToWhite, WhiteToGreen" % ColorPalette)
    
    ResColors = {}
    ColorType = ""
    
    if re.match("^WhiteToGreen$", ColorPalette, re.I):
        ColorType = "H2"
        ResColors = {"ile" : [0.938,1,0.938], "phe" : [0.891,1,0.891], "val" : [0.844,1,0.844], "leu" : [0.793,1,0.793], "trp" : [0.746,1,0.746], "met" : [0.699,1,0.699], "ala" : [0.652,1,0.652], "gly" : [0.606,1,0.606], "cys" : [0.555,1,0.555], "tyr" : [0.508,1,0.508], "pro" : [0.461,1,0.461], "thr" : [0.414,1,0.414], "ser" : [0.363,1,0.363], "his" : [0.316,1,0.316], "glu" : [0.27,1,0.27], "asn" : [0.223,1,0.223], "gln" : [0.176,1,0.176], "asp" : [0.125,1,0.125], "lys" : [0.078,1,0.078], "arg" : [0.031,1,0.031]}
    else:
        ColorType = "H1"
        ResColors = {"ile" : [0.996,0.062,0.062], "phe" : [0.996,0.109,0.109], "val" : [0.992,0.156,0.156], "leu" : [0.992,0.207,0.207], "trp" : [0.992,0.254,0.254], "met" : [0.988,0.301,0.301], "ala" : [0.988,0.348,0.348], "gly" : [0.984,0.394,0.394], "cys" : [0.984,0.445,0.445], "tyr" : [0.984,0.492,0.492], "pro" : [0.980,0.539,0.539], "thr" : [0.980,0.586,0.586], "ser" : [0.980,0.637,0.637], "his" : [0.977,0.684,0.684], "glu" : [0.977,0.730,0.730], "asn" : [0.973,0.777,0.777], "gln" : [0.973,0.824,0.824], "asp" : [0.973,0.875,0.875], "lys" : [0.899,0.922,0.922], "arg" : [0.899,0.969,0.969]}
        
    # Set up colors...
    for ResName in ResColors:
        ColorName = "color_%s_%s" % (ResName, ColorType)
        cmd.set_color(ColorName, ResColors[ResName])
        
        ResSelection = "(%s and resn %s*)" % (Selection, ResName)
        cmd.color(ColorName, ResSelection)
    
cmd.extend("ColorByHydrophobicity", ColorByHydrophobicity)

def ColorAtomsByHydrophobicityAndCharge(Selection, HydrophobicAtomsColor = "yellow", NegativelyChargedAtomsColor = "red", PositivelyChargedAtomsColor = "blue", OtherAtomsColor = "gray90"):
    \"""Color atoms in amino acids by their propensity to make hydrophobic and
    charge interactions [REF 140].  The atom names in standard amino acid
    residues are used to identify atom types as shown below:
    
        Hydrophobic: C atoms not bound to N or O atoms
        NegativelyCharged: Side chain O atoms in ASP and GLU
        PositivelyCharged: Side chain N atoms in ARG and LYS
        Others: Remaining atoms in polar and other residues
         
    The following color scheme is used by default:

        Hydrophobic: yellow
        NegativelyCharged: red
        PositivelyCharged: blue
        Others: gray90
    
    The  amino acid atom names to color specific atoms  are taken from YRB.py
    script [ REF 140]. The color values may also be specified as comma delimited
    RGB triplets. For example: HydrophobicAtomsColor = "0.950 0.78 0.0",
    NegativelyChargedAtomsColor = "1.0 0.4 0.4", PositivelyChargedAtomsColor
    = "0.2 0.5 0.8", OtherAtomsColor = "0.95 0.95 0.95"
    
    \"""

    # Process colors...
    Colors = {"HydrophobicAtomsColor": HydrophobicAtomsColor, "NegativelyChargedAtomsColor": NegativelyChargedAtomsColor, "PositivelyChargedAtomsColor": PositivelyChargedAtomsColor, "OtherAtomsColor": OtherAtomsColor}
    ColorsRGB = {}
    for ColorKey, ColorValue in Colors.items():
        if re.search(" ", ColorValue):
            ColorRGB = ColorValue.split(" ")
        else:
            ColorRGB = cmd.get_color_tuple(cmd.get_color_index(ColorValue))
        ColorsRGB[ColorKey] = ColorRGB
    
    # Set up colors...
    for ColorKey, ColorValue in ColorsRGB.items():
        cmd.set_color(ColorKey, ColorValue)

    #  Set up colors for atom names across all resiudes...
    AtomColors = {"OtherAtomsColor": "N,C,CA,O", "HydrophobicAtomsColor": "CB"}
    
    #  Set up colors for atom names across specific resiudes...
    ResAtomColors = {"arg" : {"HydrophobicAtomsColor": "CG", "PositivelyChargedAtomsColor": "NE,NH2,NH1", "OtherAtomsColor": "CD,CZ"}, "asn" : {"OtherAtomsColor": "CG,OD1,ND2"}, "asp" : {"NegativelyChargedAtomsColor": "OD2,OD1", "OtherAtomsColor": "CG"}, "cys" : {"OtherAtomsColor": "SG"}, "gln" : {"HydrophobicAtomsColor": "CG", "OtherAtomsColor": "CD,OE1,NE2"}, "glu" : {"HydrophobicAtomsColor": "CG", "NegativelyChargedAtomsColor": "OE1,OE2", "OtherAtomsColor": "CD"}, "his" : {"OtherAtomsColor": "CG,CD2,ND1,NE2,CE1"}, "ile" : {"HydrophobicAtomsColor": "CG1,CG2,CD1"}, "leu" : {"HydrophobicAtomsColor": "CG,CD1,CD2"}, "lys" : {"HydrophobicAtomsColor": "CG,CD", "PositivelyChargedAtomsColor": "NZ", "OtherAtomsColor": "CE"}, "met" : {"HydrophobicAtomsColor": "CG,CE", "OtherAtomsColor": "SD"}, "phe" : {"HydrophobicAtomsColor": "CG,CD1,CE1,CZ,CE2,CD2"}, "pro" : {"HydrophobicAtomsColor": "CG", "OtherAtomsColor": "CD"}, "ser" : {"OtherAtomsColor": "CB,OG"}, "thr" : {"HydrophobicAtomsColor": "CG2", "OtherAtomsColor": "CB,OG1"}, "trp" : {"HydrophobicAtomsColor": "CG,CD2,CZ2,CH2,CZ3,CE3", "OtherAtomsColor": "CD1,NE1,CE2"}, "tyr" : {"HydrophobicAtomsColor": "CG,CE1,CD1,CE2,CD2", "OtherAtomsColor": "CZ,OH"}, "val" : {"HydrophobicAtomsColor": "CG1,CG2"}}
    
    #  Color atom names across all resiudes...
    for ColorName in AtomColors:
        AtomNames = AtomColors[ColorName]
        AtomsSelection = "(%s and name %s)" % (Selection, AtomNames)
        cmd.color(ColorName, AtomsSelection)
    
    #  Color hydrogen atoms across all residues to other color...
    AtomsSelection = "(%s and hydro)" % (Selection)
    cmd.color("OtherAtomsColor", AtomsSelection)
    
    #  Color atom names across specific resiudes...
    for ResName in ResAtomColors:
        for ColorName in ResAtomColors[ResName]:
            AtomNames = ResAtomColors[ResName][ColorName]
            AtomsSelection = "(%s and resn %s and name %s)" % (Selection, ResName, AtomNames)
            cmd.color(ColorName, AtomsSelection)
    
cmd.extend("ColorAtomsByHydrophobicityAndCharge", ColorAtomsByHydrophobicityAndCharge)

def CheckAndDeleteEmptyObjects(ObjectNames, ParentObjectName = None):
    \"""Delete an empty objects along with optionally deleting their parent.
        
    \"""
    
    ObjectNamesList = ObjectNames.split(",")

    AllObjectsEmpty = True
    for ObjectName in ObjectNamesList:
        ObjectName = ObjectName.strip()
        if  cmd.count_atoms("(%s)" % ObjectName):
            AllObjectsEmpty = False
        else:
            cmd.delete("%s" % ObjectName)
    
    if AllObjectsEmpty and ParentObjectName is not None:
        cmd.delete("%s" % ParentObjectName)
    
cmd.extend("CheckAndDeleteEmptyObjects", CheckAndDeleteEmptyObjects)

python end"""
    
    return PMLForPython

def SetupPMLForEnableDisable(Name, Enable = True):
    """Setup PML command for enabling or disabling display of a PyMOL object.

    Arguments:
        Name (str): Name of a PyMOL object.
        Enable (bool): Display status.

    Returns:
        str: PML command for enabling or disabling display of an object.

    """
    
    if Enable:
        PML = """cmd.enable("%s")""" % Name
    else:
        PML = """cmd.disable("%s")""" % Name
        
    return PML

def SetupPMLForGroup(GroupName, GroupMembersList, Enable = None, Action = None):
    """Setup PML commands for creating a group from a list of group members. The
    display and open status of the group may be optionally set. The 'None' values
    for Enable and Action imply usage of PyMOL defaults for the creation of group.

    Arguments:
        GroupName (str): Name of a PyMOL group.
        GroupMembersList (list): List of group member names.
        Enable (bool): Display status of group.
        Action (str): Open or close status of group object.

    Returns:
        str: PML commands for creating a group object.

    """

    PMLCmds = []
    
    GroupMembers = " ".join(GroupMembersList)
    PMLCmds.append("""cmd.group("%s", "%s")""" % (GroupName, GroupMembers))
    
    if Enable is not None:
        if Enable:
            PMLCmds.append("""cmd.enable("%s")""" % GroupName)
        else:
            PMLCmds.append("""cmd.disable("%s")""" % GroupName)
    
    if Action is not None:
        PMLCmds.append("""cmd.group("%s", action="%s")""" % (GroupName, Action))

    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForLigandView(Name, Selection, LigandResName, Enable = True):
    """Setup PML commands for creating a ligand view corresponding to a ligand 
    present in a selection. The ligand is identified using organic selection
    operator available in PyMOL in conjunction with the specified ligand ID.
    The ligand is colored by atom types and displayed as 'sticks'.

    Arguments:
        Name (str): Name of a new PyMOL ligand object.
        Selection (str): PyMOL selection containing ligand.
        LigandResName (str): Ligand ID.
        Enable (bool): Display status of ligand object.

    Returns:
        str: PML commands for a ligand view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((%s) and organic and (resn %s))")""" % (Name, Selection, LigandResName))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""util.cbag("%s", _self = cmd)""" % (Name))
    PMLCmds.append("""cmd.show("sticks", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForLigandPocketView(Name, Selection, LigandSelection, DistanceCutoff, Enable = True):
    """Setup PML commands for creating a ligand binding pocket view
    corresponding all residues present in a selection within a specified
    distance from a ligand selection. The solvent and inorganic portions of
    the selection are not included in the binding pocket. The pocket residues
    are shown as 'lines'. The hydrogen atoms are not displayed.

    Arguments:
        Name (str): Name of a new PyMOL binding pocket object.
        Selection (str): PyMOL selection containing binding pocket residues.
        LigandSelection (str): PyMOL selection containing ligand.
        DistanceCutoff (float): Distance cutoff from ligand for selecting
            binding pockect residues.
        Enable (bool): Display status of binding pocket object.

    Returns:
        str: PML commands for a ligand binding pocket view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((byresidue (%s) within %.1f of (%s)) and (not solvent) and (not inorganic) and (not organic))")""" % (Name, Selection, DistanceCutoff, LigandSelection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("lines", "(%s)")""" % (Name))
    PMLCmds.append("""cmd.hide("(%s and hydro)")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForLigandPocketSolventView(Name, Selection, LigandSelection, DistanceCutoff, Enable = True):
    """Setup PML commands for creating a ligand binding pocket view
    corresponding to only solvent residues present in a selection within a
    specified distance from a ligand selection. The solvent pocket residues
    are shown as 'lines' and 'nonbonded'.

    Arguments:
        Name (str): Name of a new PyMOL solvent binding pocket object.
        Selection (str): PyMOL selection containing binding pocket residues.
        LigandSelection (str): PyMOL selection containing ligand.
        DistanceCutoff (float): Distance cutoff from ligand for selecting
            binding pocket solvent residues.
        Enable (bool): Display status of binding pocket object.

    Returns:
        str: PML commands for a ligand binding pocket view only showing solvent
            residues.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((byresidue (%s) within %.1f of (%s)) and solvent)")""" % (Name, Selection, DistanceCutoff, LigandSelection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("nonbonded", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("lines", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForLigandPocketInorganicView(Name, Selection, LigandSelection, DistanceCutoff, Enable = True):
    """Setup PML commands for creating a ligand binding pocket view
    corresponding to only inorganic residues present in a selection within a
    specified distance from a ligand selection. The inorganic pocket residues
    are shown as 'lines' and 'nonbonded'.

    Arguments:
        Name (str): Name of a new PyMOL solvent binding pocket object.
        Selection (str): PyMOL selection containing binding pocket residues.
        LigandSelection (str): PyMOL selection containing ligand.
        DistanceCutoff (float): Distance cutoff from ligand for selecting
            binding pocket inorganic residues.
        Enable (bool): Display status of binding pocket object.

    Returns:
        str: PML commands for a ligand binding pocket view only showing inorganic
            residues.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((byresidue (%s) within %.1f of (%s)) and inorganic)")""" % (Name, Selection, DistanceCutoff, LigandSelection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("nonbonded", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("lines", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForPolarContactsView(Name, Selection1, Selection2, Enable = True, Color = "yellow", Cutoff = None):
    """Setup PML commands for creating polar contacts view between a pair of
    selections. The polar contact view is generated using 'util.dist' command. The
    distance labels are shown by default.

    Arguments:
        Name (str): Name of a new PyMOL polar contacts object.
        Selection1 (str): First PyMOL selection.
        Selection2 (str): Second PyMOL selection.
        Enable (bool): Display status of polar contacts object.
        Color (str): Color for polar contact lines and labels.
        DistanceCutoff (float): None or distance cutoff for polar contacts.

    Returns:
        str: PML commands for polar contacts view between a pair of selections.

    """
    
    PMLCmds = []
    if Cutoff is None:
        PMLCmds.append("""cmd.dist("%s","(%s)","(%s)", quiet = 1, mode = 2, label = 1, reset = 1)""" % (Name, Selection1, Selection2))
    else:
        PMLCmds.append("""cmd.dist("%s","(%s)","(%s)", cutoff = %.1f, quiet = 1, mode = 2, label = 1, reset = 1)""" % (Name, Selection1, Selection2, Cutoff))
        
    PMLCmds.append("""util.color_deep("%s", "%s")""" % (Color, Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForHydrophobicContactsView(Name, Selection1, Selection2, Enable = True, Color = "yellow", Cutoff = None):
    """Setup PML commands for creating hydrophobic contacts view between a pair of
    selections. The hydrophobic contacts are shown between pairs of carbon atoms not
    connected to hydrogen bond donor or acceptors atoms as identified by PyMOL. The
    distance labels are shown by default.

    Arguments:
        Name (str): Name of a new PyMOL polar contacts object.
        Selection1 (str): First PyMOL selection.
        Selection2 (str): Second PyMOL selection.
        Enable (bool): Display status of polar contacts object.
        Color (str): Color for polar contact lines and labels.
        Cutoff (float): None or distance cutoff for hydrophobic contacts.

    Returns:
        str: PML commands for polar contacts view between a pair of selections.

    """

    PMLCmds = []
    
    HydrophobicSelectionAtoms1 = "((%s) and (elem C) and (not bound_to (donors or acceptors)))" % (Selection1)
    HydrophobicSelectionAtoms2 = "((%s) and (elem C) and (not bound_to (donors or acceptors)))" % (Selection2)
    
    if Cutoff is None:
        PMLCmds.append("""cmd.dist("%s","(%s)","(%s)", quiet = 1, label = 1, reset = 1)""" % (Name, HydrophobicSelectionAtoms1, HydrophobicSelectionAtoms2))
    else:
        PMLCmds.append("""cmd.dist("%s","(%s)","(%s)", cutoff = %.1f, quiet = 1,  label = 1, reset = 1)""" % (Name, HydrophobicSelectionAtoms1, HydrophobicSelectionAtoms2, Cutoff))
        
    PMLCmds.append("""util.color_deep("%s", "%s")""" % (Color, Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForAlignment(Method, RefSelection, FitSelection):
    """Setup PML commands for aligning a pair of selection using  a specified
    alignment method.

    Arguments:
        Method (str): Alignment method. Possible values: align, cealign, super.
        RefSelection (str): Name of reference selection which stays stationary.
        FitSelection (str): Name of selection to align to reference selection.

    Returns:
        str: PML commands for aligning  a pair of selections.

    """

    PMLCmds = []
    if re.match("^align$", Method, re.I):
        PMLCmds.append("""cmd.align("(%s)", "(%s)")""" % (FitSelection, RefSelection))
    elif re.match("^cealign$", Method, re.I):
        PMLCmds.append("""cmd.cealign("(%s)", "(%s)")""" % (RefSelection, FitSelection))
    elif re.match("^super$", Method, re.I):
        PMLCmds.append("""cmd.super("(%s)", "(%s)")""" % (FitSelection, RefSelection))
    else:
        MiscUtil.PrintWarning("PyMOLUtil.SetupPMLForAlignment: Invalid method name: %s" % Method)

    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForBFactorPuttyView(Name, Selection, ColorPalette = "blue_white_red", Enable = True):
    """Setup PML commands for creating a B factor putty view for a specified
    selection. The B factor values must be available for the atoms. The atoms
    are colored using a color spectrum corresponding to a specified color
    palette. Any valid PyMOL color palette name may be used.

    Arguments:
        Name (str): Name of a new PyMOL B factor putty object.
        Selection (str): Name of PyMOL selection.
        ColorPalette (str): Name of color palette to use for color spectrum.
        Enable (bool): Display status of B factor putty object.

    Returns:
        str: PML commands for B factor putty view.

    """
    
    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "(%s)")""" % (Name, Selection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.spectrum("b", "%s", "(%s)")""" % (ColorPalette, Name))
    PMLCmds.append("""cmd.show("cartoon", "%s")""" % (Name))
    PMLCmds.append("""cmd.cartoon("putty", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForHydrophobicSurfaceView(Name, Selection, ColorPalette = "RedToWhite", Enable = True, DisplayAs = "cartoon"):
    """Setup PML commands for creating a hydrophobic surface view for a specified
    selection. The surfaces are colored using a specified color palette. This is only valid
    for amino acids.

    Arguments:
        Name (str): Name of a new PyMOL hydrophobic surface object.
        Selection (str): Name of PyMOL selection.
        ColorPalette (str): Name of color palette to use for coloring surfaces.
            Possible values: RedToWhite or WhiteToGreen for most hydrophobic
            to least hydrophobic amino acids.
        Enable (bool): Display status of surface object.
        DisplayAs (str): Any additional valid display type such as lines,
            sticks, ribbon, cartoon, or None. 

    Returns:
        str: PML commands for hydrophobic surface view.

    """

    PMLCmds = _GetPMLCmdsForSurfaceView(Name, Selection, DisplayAs)
    PMLCmds.append("""ColorByHydrophobicity("%s", "%s")""" % (Name, ColorPalette))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForHydrophobicAndChargeSurfaceView(Name, Selection, HydrophobicAtomsColor = "yellow", NegativelyChargedAtomsColor = "red", PositivelyChargedAtomsColor = "blue", OtherAtomsColor = "gray90", Enable = True, DisplayAs = "cartoon"):
    """Setup PML commands for creating a surface colored by hydrophobic and
    charge [ REF 140] properties of atoms in amino acids. The atom names in
    standard amino acid residues are used to identify atom types as shown below:

    Hydrophobic: C atoms not bound to N or O atoms; NegativelyCharged: Side
    chain O atoms in ASP and GLU; PositivelyCharged: Side chain N atoms in
    ARG and LYS; Others: Remaining atoms in polar and other residues

    The  amino acid atom names to color specific atoms  are taken from YRB.py
    script [ REF 140]. The color values may also be specified as comma delimited
    RGB triplets. For example: HydrophobicAtomsColor = "0.950 0.78 0.0",
    NegativelyChargedAtomsColor = "1.0 0.4 0.4", PositivelyChargedAtomsColor
    = "0.2 0.5 0.8", OtherAtomsColor = "0.95 0.95 0.95"

    Arguments:
        Name (str): Name of a new PyMOL hydrophobic surface object.
        Selection (str): Name of PyMOL selection.
        HydrophobicAtomsColor (str): Color name or space delimited RGB values
        NegativelyChargedAtomsColor (str): Color name or space delimited RGB values
        PositivelyChargedAtomsColor (str): Color name or space delimited RGB values
        OtherAtomsColor (str): Color name or space delimited RGB values
        Enable (bool): Display status of surface object.
        DisplayAs (str): Any additional valid display type such as lines,
            sticks, ribbon, cartoon, or None. 

    Returns:
        str: PML commands for hydrophobic and charge surface view.

    """

    PMLCmds = _GetPMLCmdsForSurfaceView(Name, Selection, DisplayAs)
    PMLCmds.append("""ColorAtomsByHydrophobicityAndCharge("%s", "%s", "%s", "%s", "%s")""" % (Name, HydrophobicAtomsColor, NegativelyChargedAtomsColor, PositivelyChargedAtomsColor, OtherAtomsColor))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForSurfaceView(Name, Selection, Enable = True, DisplayAs = "cartoon", Color = "None"):
    """Setup PML commands for creating a molecular surface view for a specified
    selection.

    Arguments:
        Name (str): Name of a new PyMOL molecular surface object.
        Selection (str): Name of PyMOL selection.
        Enable (bool): Display status of surface object.
        DisplayAs (str): Any additional valid display type such as lines,
            sticks, ribbon, cartoon, or None. 
        Color (str): Surafce color.

    Returns:
        str: PML commands for molecular surface view.

    """

    PMLCmds = _GetPMLCmdsForSurfaceView(Name, Selection, DisplayAs)
    if Color is not None:
        PMLCmds.append("""util.color_deep("%s", "%s")""" % (Color, Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def _GetPMLCmdsForSurfaceView(Name, Selection, DisplayAs = "cartoon"):
    """Setup PML command for surface view."""

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "(%s)")""" % (Name, Selection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    if DisplayAs is not None:
        PMLCmds.append("""cmd.show("%s", "%s")""" % (DisplayAs, Name))
    PMLCmds.append("""cmd.show("surface", "%s")""" % (Name))

    return PMLCmds
    
def SetupPMLForSelectionDisplayView(Name, Selection, DisplayAs, Color = None, Enable = True):
    """Setup PML commands for creating a specific molecular display view for a
     selection.

    Arguments:
        Name (str): Name of a new PyMOL object.
        Selection (str): Name of PyMOL selection.
        DisplayAs (str): Any valid display type such as lines, sticks, ribbon,
            cartoon, or surface
        Color (str): Color name or use default color.
        Enable (bool): Display status of object.

    Returns:
        str: PML commands for molecular selection view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "(%s)")""" % (Name, Selection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("%s", "%s")""" % (DisplayAs, Name))
    if Color is not None:
        PMLCmds.append("""util.color_deep("%s", "%s")""" % (Color, Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForBallAndStickView(Name, Selection, Enable = True, SphereScale = 0.3, StickRadius = 0.2):
    """Setup PML commands for creating a ball and stick view for a specified
    selection.

    Arguments:
        Name (str): Name of a new PyMOL ball and stick object.
        Selection (str): Name of PyMOL selection.
        Enable (bool): Display status of ball and stick object.
        SphereScale (float): Scaling factor for sphere radii.
        StickScale (float): Scaling factor for stick radii.

    Returns:
        str: PML commands for ball and stick view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "(%s)")""" % (Name, Selection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("sphere", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("sticks", "%s")""" % (Name))
    PMLCmds.append("""cmd.set("sphere_scale", %.1f, "%s")""" % (SphereScale, Name))
    PMLCmds.append("""cmd.set("stick_radius", %.1f, "%s")""" % (StickRadius, Name))
    
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForInorganicView(Name, Selection, Enable = True):
    """Setup PML commands for creating a inorganic view corresponding to
    inorganic residues present in a selection. The inorganic residues are
    identified using inorganic selection operator available in PyMOL. The
    inorganic residues are displayed as 'lines' and 'nonbonded'.

    Arguments:
        Name (str): Name of a new PyMOL inorganic object.
        Selection (str): Name of PyMOL selection.
        Enable (bool): Display status of inorganic object.

    Returns:
        str: PML commands for inorganic view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((%s) and inorganic)")""" % (Name, Selection))
    PMLCmds.append("""cmd.show("nonbonded", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("lines", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForSolventView(Name, Selection, Enable = True):
    """Setup PML commands for creating a solvent view corresponding to
    solvent residues present in a selection. The solvent residues are
    identified using solvent selection operator available in PyMOL. The
    solvent residues are displayed as 'nonbonded'.

    Arguments:
        Name (str): Name of a new PyMOL solvent object.
        Selection (str): Name of PyMOL selection.
        Enable (bool): Display status of inorganic object.

    Returns:
        str: PML commands for solvent view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((%s) and solvent)")""" % (Name, Selection))
    PMLCmds.append("""cmd.show("nonbonded", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForPolymerChainView(Name, Selection, Enable = True):
    """Setup PML commands for creating a polymer chain view corresponding
    to backbone and sidechain residues in a selection. The polymer chain is
    displayed as 'cartoon'.

    Arguments:
        Name (str): Name of a new PyMOL polymer chain object.
        Selection (str): Name of PyMOL selection.
        Enable (bool): Display status of chain object.

    Returns:
        str: PML commands for polymer chain view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((%s) and (backbone or sidechain))")""" % (Name, Selection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""util.cbag("%s", _self = cmd)""" % (Name))
    PMLCmds.append("""cmd.show("cartoon", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    
    PML = "\n".join(PMLCmds)
    
    return PML


def SetupPMLForPolymerComplexView(MoleculeName, PDBFile, Enable = True, ShowSolvent = True, ShowInorganic = True, ShowLines = True):
    """Setup PML commands for creating a polymer complex view for all chains
    in a PDB file. The solvent and inorganic residues are also shown by default.
     The polymer chains are displayed as 'cartoon'. The 'line' display for the
    polymer chains is also shown and may be turned off. The organic residues are
    displayed as 'sticks'. The solvent and inorganic residues are displayed as
    'nonbonded' and 'lines'.
    
    Arguments:
        MoleculeName (str): Name of a new PyMOL molecule object.
        PDBFile (str): Name of PDB file.
        Enable (bool): Display status of chain object.
        ShowSolvent (bool): Display solvent residues.
        ShowInorganic (bool): Display inorganic residues.
        ShowLines (bool): Display lines for polymer chains.

    Returns:
        str: PML commands for polymer complex view.

    """

    PMLCmds = []
    
    PMLCmds.append("""cmd.load("%s", "%s")""" % (PDBFile, MoleculeName))
    PML = _SetupPMLForPolymerComplexView(MoleculeName, Enable, ShowSolvent, ShowInorganic, ShowLines)
    PMLCmds.append(PML)
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForPolymerChainComplexView(ChainComplexName, Selection, ChainName, Enable = True, ShowSolvent = True, ShowInorganic = True, ShowLines = True):
    """Setup PML commands for creating a polymer chain complex view for a specified
    chain in a selection. The solvent and inorganic residues are also shown by
    default. The polymer chain is displayed as 'cartoon'. The 'line' display for the
    polymer chain is also shown and may be turned off. The organic residues are
    displayed as 'sticks'. The solvent and inorganic residues are displayed as
    'nonbonded' and 'lines'.

    Arguments:
        ChainComplexName (str): Name of a new PyMOL polymer chain complex.
        Selection (str): Name of PyMOL selection.
        ChainName (str): Name of a chain.
        Enable (bool): Display status of chain object.
        ShowSolvent (bool): Display solvent residues.
        ShowInorganic (bool): Display inorganic residues.
        ShowLines (bool): Display lines for polymer chain.

    Returns:
        str: PML commands for polymer chain complex view.

    """

    PMLCmds = []
    
    PMLCmds.append("""cmd.create("%s", "(%s and chain %s)")""" % (ChainComplexName, Selection, ChainName))
    PML = _SetupPMLForPolymerComplexView(ChainComplexName, Enable, ShowSolvent, ShowInorganic, ShowLines)
    PMLCmds.append(PML)

    PML = "\n".join(PMLCmds)
    
    return PML

def _SetupPMLForPolymerComplexView(Name, Enable = True, ShowSolvent = True, ShowInorganic = True,  ShowLines = False):
    """Setup PML for creating a polymer complex view."""

    PMLCmds = []
    
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("cartoon", "%s")""" % (Name))
    PMLCmds.append("""util.cba(33, "%s", _self = cmd)""" % (Name))
    PMLCmds.append("""cmd.show("sticks", "(organic and (%s))")""" % (Name))
    if ShowSolvent:
        PMLCmds.append("""cmd.show("nonbonded", "(solvent and (%s))")""" % (Name))
    if ShowInorganic:
        PMLCmds.append("""cmd.show("nonbonded", "(inorganic and (%s))")""" % (Name))
    
    if ShowLines:
        PMLCmds.append("""cmd.show("lines", "%s")""" % (Name))
    else:
        if ShowInorganic:
            PMLCmds.append("""cmd.show("lines", "(inorganic and (%s))")""" % (Name))
    
    PMLCmds.append("""cmd.set_bond("valence", "1", "%s", quiet = 1)""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))

    PML = "\n".join(PMLCmds)
    
    return PML
