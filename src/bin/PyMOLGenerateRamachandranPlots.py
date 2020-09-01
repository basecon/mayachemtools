#!/usr/bin/env python
#
# File: PyMOLGenerateRamachandranPlots.py
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
import csv
import matplotlib.pyplot as plt
import numpy as np

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
    GenerateRamachandranPlots()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def GenerateRamachandranPlots():
    """Calculate phi and psi angles for macromolecules containing amino acids
    and generate Ramachandran plots.
    """

    # Calculate phi and psi angles..
    CalculatePhiPsiAngles()

    # Read phi and psi densities...
    ReadPhiAndPsiDensities()

    # Setup contour info...
    SetupContoursInfo()

    # Generate plots...
    if OptionsInfo["MultipleOutFiles"]:
        GenerateMultiplePlotFiles()
    else:
        GenerateSinglePlotFile()

def GenerateSinglePlotFile():
    """Generate a single plot file containg all four types of Ramachandran plots."""
    
    Outfile = OptionsInfo["Outfile"]
    MiscUtil.PrintInfo("\nGenerating output file %s..." % (Outfile))

    SetupFontFamily()

    # Setup figure...
    PlotFigure, Axes = plt.subplots(2, 2, figsize = (OptionsInfo["FigWidth"], OptionsInfo["FigHeight"]), dpi = OptionsInfo["FigDPI"])
    PlotAxes = [Axis for RowAxes in Axes for Axis in RowAxes]

    # Adjust space between subplots...
    plt.subplots_adjust(left  = 0.1, right = 0.9,  bottom = 0.1, top = 0.9, wspace = 0.3, hspace = 0.25)

    for  PlotIndex, PlotType in enumerate(OptionsInfo["PlotTypesInfo"]["Types"]):
        DrawPlot(PlotAxes[PlotIndex], PlotType)

    # Save figure...
    plt.savefig(Outfile)

def GenerateMultiplePlotFiles():
    """Generate multiple plot files corresponding to four types of Ramachandran plots."""
    
    MiscUtil.PrintInfo("\nGenerating multiple output files...")
    
    SetupFontFamily()

    for  PlotType in OptionsInfo["PlotTypesInfo"]["Types"]:
        Outfile = OptionsInfo["PlotTypesInfo"]["Outfiles"][PlotType]
        MiscUtil.PrintInfo("Generating  output file %s..." % Outfile)
        
        PlotFigure, PlotAxis = plt.subplots(1, 1, figsize = (OptionsInfo["FigWidth"], OptionsInfo["FigHeight"]), dpi = OptionsInfo["FigDPI"])
        DrawPlot(PlotAxis, PlotType)

        # Save figure...
        plt.savefig(Outfile)

        # Get ready for the next figure...
        plt.clf()

def DrawPlot(PlotAxis, PlotType):
    """Draw contour and scatter plot """
    
    PlotTypesInfo = OptionsInfo["PlotTypesInfo"]
    
    # Draw filled contours...
    PhiPsiContourInfo = PlotTypesInfo["PhiPsiContourInfo"][PlotType]
    PlotAxis.contourf(PhiPsiContourInfo["X"], PhiPsiContourInfo["Y"], PhiPsiContourInfo["Z"], levels = PlotTypesInfo["Levels"][PlotType], colors = PlotTypesInfo["Colors"][PlotType])
    
    # Draw scatter plot for phi and psi angles...
    if PlotTypesInfo["ResCount"][PlotType]:
        PlotAxis.scatter(PlotTypesInfo["PhiAngles"][PlotType], PlotTypesInfo["PsiAngles"][PlotType], s = OptionsInfo["ScatterMarkerSize"], c = OptionsInfo["ScatterMarkerColor"], marker = OptionsInfo["ScatterMarkerStyle"])
    
    # Setup limits...
    PlotAxis.set_xlim(PlotTypesInfo["Limits"][PlotType])
    PlotAxis.set_ylim(PlotTypesInfo["Limits"][PlotType])
    
    # Setup major tick marks...
    PlotAxis.set_xticklabels(PlotTypesInfo["MajorTickLabels"][PlotType], fontdict = {"fontsize": OptionsInfo["FontTicksSize"], "fontweight": OptionsInfo["FontTicksWeight"]})
    PlotAxis.set_xticks(PlotTypesInfo["MajorTickPositions"][PlotType])
    PlotAxis.set_yticklabels(PlotTypesInfo["MajorTickLabels"][PlotType], fontdict = {"fontsize": OptionsInfo["FontTicksSize"], "fontweight": OptionsInfo["FontTicksWeight"]})
    PlotAxis.set_yticks(PlotTypesInfo["MajorTickPositions"][PlotType])

    # Set up minor ticks...
    if OptionsInfo["TicksMinor"]:
        PlotAxis.set_xticks(PlotTypesInfo["MinorTickPositions"][PlotType], minor = True)
        PlotAxis.set_yticks(PlotTypesInfo["MinorTickPositions"][PlotType], minor = True)

    # Setup grid...
    if OptionsInfo["Grid"]:
        PlotAxis.grid(True, color = OptionsInfo["GridLineColor"],  linestyle = OptionsInfo["GridLineStyle"], linewidth = OptionsInfo["GridLineWidth"])
    
    # Setup title...
    PlotAxis.set_title(PlotTypesInfo["Titles"][PlotType], fontsize = OptionsInfo["FontTitleSize"], fontweight = OptionsInfo["FontTitleWeight"])
    
    # Setup axes labels...
    if PlotTypesInfo["DrawXLabel"][PlotType]:
        XLabel = r"$\Phi$" if OptionsInfo["Greek"] else "Phi"
        PlotAxis.set_xlabel(XLabel, fontsize = OptionsInfo["FontAxesSize"], fontweight = OptionsInfo["FontAxesWeight"])
    if PlotTypesInfo["DrawYLabel"][PlotType]:
        # Setup a horizontal ylabel close to the axis...
        YLabel = r"$\Psi$" if OptionsInfo["Greek"] else "Psi"
        YRotation = 0 if OptionsInfo["Greek"] else 90
        PlotAxis.set_ylabel(YLabel, fontsize = OptionsInfo["FontAxesSize"], fontweight = OptionsInfo["FontAxesWeight"], rotation = YRotation, labelpad = 0)

def SetupFontFamily():
    """Setuo global font family. """
    
    if re.match("^auto$", OptionsInfo["FontFamily"], re.I):
        return
    plt.rcParams["font.family"] = OptionsInfo["FontFamily"]

def SetupContoursInfo():
    """Setup contour info for generating contour plots."""

    MiscUtil.PrintInfo("\nProcessing phi and psi densities for contour plots...")
    
    OptionsInfo["PlotTypesInfo"]["PhiPsiContourInfo"] = {}
    for  PlotType in OptionsInfo["PlotTypesInfo"]["Types"]:
        PhiPsiContourInfo = SetupPhiAndPsiContourInfo(PlotType)
        OptionsInfo["PlotTypesInfo"]["PhiPsiContourInfo"][PlotType] = PhiPsiContourInfo

def SetupPhiAndPsiContourInfo(PlotType):
    """Setup X, Y and Z contour arrays for generating contour plots. """

    DensityInfo = OptionsInfo["PlotTypesInfo"]["PhiPsiDensityInfo"][PlotType]
    
    X, Y = np.meshgrid(DensityInfo["PhiValues"], DensityInfo["PsiValues"])
    Z = np.zeros((len(DensityInfo["PhiValues"]), len(DensityInfo["PsiValues"])))
    
    # Initialize X, Y, and Z arrays for contour plots...
    for ZRowIndex, PsiID in enumerate(DensityInfo["PhiIDs"]):
        for ZColIndex, PhiID in enumerate(DensityInfo["PsiIDs"]):
            Z[ZRowIndex][ZColIndex] = DensityInfo["Density"][PhiID][PsiID]
            
    # Track contour data...
    ContourInfo = {}
    ContourInfo["X"] = X
    ContourInfo["Y"] = Y
    ContourInfo["Z"] = Z
    
    return ContourInfo
    
def ReadPhiAndPsiDensities():
    """Read phi and psi densities for generating filled contours. """

    OptionsInfo["PlotTypesInfo"]["PhiPsiDensityInfo"] = {}
    for  PlotType in OptionsInfo["PlotTypesInfo"]["Types"]:
        DensityFile = OptionsInfo["PlotTypesInfo"]["PhiPsiDensityFiles"][PlotType]
        PhiPsiDensityInfo = ReadPhiAndPsiDensityFile(DensityFile)
        OptionsInfo["PlotTypesInfo"]["PhiPsiDensityInfo"][PlotType] = PhiPsiDensityInfo
    
def ReadPhiAndPsiDensityFile(DensityFile):
    """Read phi and psi desnsity file.

    Format: 
    Phi,Psi,Density
    -179.0,-179.0,0.00782923406455425
    -179.0,-177.0,0.00641357067237856
    ... ... ...
    """

    MiscUtil.PrintInfo("\nReading psi and psi density grid file %s..." % DensityFile)
    
    DensityFH = open(DensityFile, "r")
    if DensityFH is None:
        MiscUtil.PrintError("Couldn't open phi and psi density file: %s.\n" % (DensityFile))
    
    HeaderLine = True
    DensityLines = []
    for Line in DensityFH:
        Line = Line.rstrip()
        # Ignore comments...
        if re.match("^#", Line, re.I):
            continue
        # Ignore header line...
        if HeaderLine:
            HeaderLine = False
            continue
        DensityLines.append(Line)

    DensityInfo = {}
    DensityInfo["Density"] = {}
    
    DensityInfo["PhiIDs"] = []
    DensityInfo["PhiValues"] = []
    
    DensityInfo["PsiIDs"] = []
    DensityInfo["PsiValues"] = []
    
    PhiValuesMap = {}
    PsiValuesMap = {}
    
    Count = 0
    MinDensity = 99999.0
    MaxDensity = - MinDensity
    
    DensityReader = csv.reader(DensityLines, delimiter=',', quotechar='"')
    for LineWords in DensityReader:
        Count += 1
        
        Phi = LineWords[0]
        Psi = LineWords[1]
        Density = LineWords[2]

        # Track unique phi and psi value...
        if not Phi in PhiValuesMap:
            PhiValuesMap[Phi] = float(Phi)
        if not Psi in PsiValuesMap:
            PsiValuesMap[Psi] = float(Psi)
        
        # Track density data...
        if not Phi in DensityInfo["Density"]:
            DensityInfo["Density"][Phi] = {}

        Density = float(Density)
        DensityInfo["Density"][Phi][Psi] = Density
        if Density < MinDensity:
            MinDensity = Density
        if Density > MaxDensity:
            MaxDensity = Density

    # Sort and track values for phi and psi angles...
    DensityInfo["PhiIDs"] = sorted(PhiValuesMap.keys(), key = lambda Phi: PhiValuesMap[Phi])
    DensityInfo["PhiValues"] = [PhiValuesMap[Phi] for Phi in DensityInfo["PhiIDs"]]
    
    DensityInfo["PsiIDs"] = sorted(PsiValuesMap.keys(), key = lambda Psi: PsiValuesMap[Psi])
    DensityInfo["PsiValues"] = [PsiValuesMap[Psi] for Psi in DensityInfo["PsiIDs"]]
    
    MiscUtil.PrintInfo("Minimum density: %.4f; Maximum density: %.4f" % (MinDensity, MaxDensity))
    MiscUtil.PrintInfo("Number of phi and psi angles: %s" % Count)
    
    MiscUtil.PrintInfo("\nDimensions of phi and psi grid angles:")
    MiscUtil.PrintInfo("Phi - Min: %s; Max: %s; Bin size: %s; Count: %s" % (DensityInfo["PhiValues"][0], DensityInfo["PhiValues"][-1], abs(DensityInfo["PhiValues"][1] - DensityInfo["PhiValues"][0]), len(DensityInfo["PhiValues"])))
    MiscUtil.PrintInfo("Psi - Min: %s; Max: %s; Bin size: %s; Count: %s" % (DensityInfo["PsiValues"][0], DensityInfo["PsiValues"][-1], abs(DensityInfo["PsiValues"][1] - DensityInfo["PsiValues"][0]), len(DensityInfo["PsiValues"])))

    return DensityInfo
    
def CalculatePhiPsiAngles():
    """Calculate phi and psi angles for scatter plots."""

    Infile = OptionsInfo["Infile"]
    MolName = OptionsInfo["InfileRoot"]

    # Load molecule...
    pymol.cmd.reinitialize()
    pymol.cmd.load(Infile, MolName)
    
    MiscUtil.PrintInfo("\nCalculating phi and psi torsion angles for input file %s..." % Infile)

    # Initialize...
    OptionsInfo["PlotTypesInfo"]["PhiAngles"] = {}
    OptionsInfo["PlotTypesInfo"]["PsiAngles"] = {}
    OptionsInfo["PlotTypesInfo"]["ResCount"] = {}
    for PlotType in  OptionsInfo["PlotTypesInfo"]["Types"]:
        OptionsInfo["PlotTypesInfo"]["PhiAngles"][PlotType] = []
        OptionsInfo["PlotTypesInfo"]["PsiAngles"][PlotType] = []
        OptionsInfo["PlotTypesInfo"]["ResCount"][PlotType] = 0

    Precision = OptionsInfo["Precision"]
    
    TotalResCount = 0
    # Go over specified chain IDs..
    for ChainID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainIDs"]:
        PhiPsiInfoList = []
        GeneralPhiPsiInfo, GlycinePhiPsiInfo, ProlinePhiPsiInfo, PreProlinePhiPsiInfo = PyMOLUtil.GetPhiPsiCategoriesResiduesInfo(MolName, ChainID)
        PhiPsiInfoList.extend([GeneralPhiPsiInfo, GlycinePhiPsiInfo, ProlinePhiPsiInfo, PreProlinePhiPsiInfo])
        
        for Index, PlotType in enumerate(OptionsInfo["PlotTypesInfo"]["Types"]):
            PhiPsiInfo = PhiPsiInfoList[Index]
            ResCount = len(PhiPsiInfo["ResNums"])
            if not ResCount:
                continue

            TotalResCount += ResCount
            OptionsInfo["PlotTypesInfo"]["ResCount"][PlotType] += ResCount
            
            PhiAngles, PsiAngles = ProcessPsiInfo(PhiPsiInfo, Precision)
            OptionsInfo["PlotTypesInfo"]["PhiAngles"][PlotType].extend(PhiAngles)
            OptionsInfo["PlotTypesInfo"]["PsiAngles"][PlotType].extend(PsiAngles)
        
    # Delete MolName object
    pymol.cmd.delete(MolName)

    MiscUtil.PrintInfo("\nTotal number of phi and psi angles: %d" % TotalResCount)

    MiscUtil.PrintInfo("")
    for PlotType in  OptionsInfo["PlotTypesInfo"]["Types"]:
        MiscUtil.PrintInfo("Number of \"%s\" phi and psi angles: %s" % (PlotType, OptionsInfo["PlotTypesInfo"]["ResCount"][PlotType]))
    
    if not TotalResCount:
        MiscUtil.PrintInfo("")
        MiscUtil.PrintWarning("No valid phi and psi angles found in input file. Ramachandran plots will be generated without phi and psi scatter plots...")
        
def ProcessPsiInfo(PhiPsiInfo, Precision):
    """Process phi and psi angels for scatter plots. """

    PhiAngles = []
    PsiAngles = []
    for ResNum in PhiPsiInfo["ResNums"]:
        Phi = "%.*f" % (Precision, PhiPsiInfo["Phi"][ResNum])
        Psi = "%.*f" % (Precision, PhiPsiInfo["Psi"][ResNum])
        PhiAngles.append(float(Phi))
        PsiAngles.append(float(Psi))

    return PhiAngles, PsiAngles

def RetrieveInfileInfo():
    """Retrieve information for input file."""
    
    Infile = OptionsInfo["Infile"]
    InfileRoot = OptionsInfo["InfileRoot"]
    
    ChainsAndLigandsInfo = PyMOLUtil.GetChainsAndLigandsInfo(Infile, InfileRoot)
    OptionsInfo["ChainsAndLigandsInfo"] = ChainsAndLigandsInfo

def ProcessChainIDs():
    """Process specified chain IDs for infile."""
    
    MiscUtil.PrintInfo("\nProcessing specified chain IDs for input file %s..." % OptionsInfo["Infile"])        
    ChainsAndLigandsInfo = OptionsInfo["ChainsAndLigandsInfo"]
    SpecifiedChainsAndLigandsInfo = PyMOLUtil.ProcessChainsAndLigandsOptionsInfo(ChainsAndLigandsInfo, "-c, --chainIDs", OptionsInfo["ChainIDs"], None, None)
    
    OptionsInfo["SpecifiedChainsAndLigandsInfo"] = SpecifiedChainsAndLigandsInfo
    
    MiscUtil.PrintInfo("Specified chain IDs: %s" % (", ".join(SpecifiedChainsAndLigandsInfo["ChainIDs"])))

def SetupPlotsInfo():
    """Setup information for generating plots. """

    InitializePlotTypesInfo()
    SetupPlotTypesOutfiles()
    
    ProessContourLevelsAndColors()

def ProessContourLevelsAndColors():
    """Process specified contour levels and colors. """

    if re.match("^auto$", OptionsInfo["LevelsAndColors"], re.I):
        return
    
    # Setup canonical plot types for validation...
    CanonicalPlotTypes = {}
    for PlotType in OptionsInfo["PlotTypesInfo"]["Types"]:
        CanonicalPlotTypes[PlotType.lower()] = PlotType
    
    LevelsAndColors = re.sub(" ", "", OptionsInfo["LevelsAndColors"])
    if not len(LevelsAndColors):
        MiscUtil.PrintError("The levels and colors specified using \"-l, --levelsAndColors\" option are empty.")
        
    for TypeLevelsColorsWord in LevelsAndColors.split(";"):
        if not len(TypeLevelsColorsWord):
            MiscUtil.PrintError("The plot types, levels, and colors, \"%s\" specified using \"-l, --levelsAndColors\" option in, \"%s\", is empty." % (TypeLevelsColorsWord, LevelsAndColors))
            
        TypeLevelsColorsWords = TypeLevelsColorsWord.split(":")
        if len(TypeLevelsColorsWords) !=2:
            MiscUtil.PrintError("The format of plot type, levels, and colors specification, \"%s\", specified using \"-l, --levelsAndColors\" option, in \"%s\",  is not valid: Supported format: <PlotType>: <Level>, <Color>, <Level>,..." % (TypeLevelsColorsWord, OptionsInfo["LevelsAndColors"]))
        
        PlotType = TypeLevelsColorsWords[0]
        if not len(PlotType):
            MiscUtil.PrintError("The plot type, \"%s\" specified using \"-l, --levelsAndColors\" option in, \"%s\", is empty." % (PlotType, TypeLevelsColorsWord))
        CanonicalPlotType = PlotType.lower()
        
        if not CanonicalPlotType in CanonicalPlotTypes:
            MiscUtil.PrintError("The plot type, \"%s\" specified using \"-l, --levelsAndColors\" option in, \"%s\", is not valid. Supported valus: %s" % (PlotType, TypeLevelsColorsWord, ", ".join(OptionsInfo["PlotTypesInfo"]["Types"])))
        PlotType = CanonicalPlotTypes[CanonicalPlotType]
            
        LevelsColorsWords = TypeLevelsColorsWords[1].split(",")
        if not (len(LevelsColorsWords) % 2):
            MiscUtil.PrintError("The format of levels and colors specifification, \"%s\", specified using \"-l, --levelsAndColors\" option in, \"%s\",  is not valid. It must contain odd number of values. Supported format: <PlotType>: <Level>, <Color>, <Level>,..." % (", ".join(LevelsColorsWords), TypeLevelsColorsWord))

        # Retrieve levels and colors...
        Levels = []
        Colors = []
        for Index, SpecWord in enumerate(LevelsColorsWords):
            if not len(SpecWord):
                MiscUtil.PrintError("The level or color, \"%s\" specified using \"-l, --levelsAndColors\" option in, \"%s\", is empty." % (SpecWord, TypeLevelsColorsWord))
            
            if Index % 2:
                Colors.append(SpecWord)
                continue
            
            # Process level...
            if not MiscUtil.IsFloat(SpecWord):
                MiscUtil.PrintError("The level, \"%s\" specified using \"-l, --levelsAndColors\" option in, \"%s\", must be a number." % (SpecWord, TypeLevelsColorsWord))
            
            Level = float(SpecWord)
            if len(Levels):
                # The current level must be greater than all previous levels..
                for PreviousLevel in Levels:
                    if Level <= PreviousLevel:
                        MiscUtil.PrintError("The level, \"%s\" specified using \"-l, --levelsAndColors\" option in, \"%s\", must be greater than all previous levels." % (SpecWord, TypeLevelsColorsWord))

            Levels.append(Level)
        
        OptionsInfo["PlotTypesInfo"]["Levels"][PlotType] = Levels
        OptionsInfo["PlotTypesInfo"]["Colors"][PlotType] = Colors
        
def InitializePlotTypesInfo():
    """Initialize information for generating plots. """

    PlotTypesInfo = {}
    PlotTypesInfo["Types"] = []
    PlotTypesInfo["PhiPsiDensityFiles"] = {}
    PlotTypesInfo["PhiPsiDensityInfo"] = {}
    PlotTypesInfo["PhiPsiContourInfo"] = {}
    
    PlotTypesInfo["Titles"] = {}
    PlotTypesInfo["DrawXLabel"] = {}
    PlotTypesInfo["DrawYLabel"] = {}
    
    PlotTypesInfo["Limits"] = {}
    PlotTypesInfo["MajorTickPositions"] = {}
    PlotTypesInfo["MajorTickLabels"] = {}
    PlotTypesInfo["MinorTickPositions"] = {}
    
    PlotTypesInfo["Levels"] = {}
    PlotTypesInfo["Colors"] = {}
    PlotTypesInfo["Outfiles"] = {}
    PlotTypesInfo["PhiAngles"] = {}
    PlotTypesInfo["PsiAngles"] = {}

    MayaChemToolsDataDir = MiscUtil.GetMayaChemToolsLibDataPath()

    # Setup contour colors for supported default schemes...
    ContourColorSchemes = {}
    ContourColorSchemes["General"] = {"MuttedColorShades1": ["#FFFFFF", "#EBF1DE", "#C3D69B"],
                                      "MuttedColorShades2": ["#FFFFFF", "#EBF1DE", "#D7E4BD"],
                                      "BrightColorShades": ["#FFFFFF", "#B3E8FF", "#7FD9FF"]}
    ContourColorSchemes["Glycine"] = {"MuttedColorShades1": ["#FFFFFF", "#FDEADA", "#FAC090"],
                                      "MuttedColorShades2": ["#FFFFFF", "#FDEADA", "#FCD5B5"],
                                      "BrightColorShades": ["#FFFFFF", "#FFE8C5", "#FFCC7F"]}
    ContourColorSchemes["Proline"] = {"MuttedColorShades1": ["#FFFFFF", "#E6E0EC", "#B3A2C7"],
                                      "MuttedColorShades2": ["#FFFFFF", "#E6E0EC", "#CCC1DA"],
                                      "BrightColorShades": ["#FFFFFF", "#D0FFC5", "#7FFF8C"]}
    ContourColorSchemes["PreProline"] = {"MuttedColorShades1": ["#FFFFFF", "#DCE6F2", "#95B3D7"],
                                         "MuttedColorShades2": ["#FFFFFF", "#DCE6F2", "#B9CDE5"],
                                         "BrightColorShades": ["#FFFFFF", "#B3E8FF", "#7FD9FF"]}

    if re.match("^MuttedColorShades1$", OptionsInfo["LevelsAndColorsScheme"], re.I):
        DefaultColorScheme = "MuttedColorShades1"
    elif re.match("^MuttedColorShades2$", OptionsInfo["LevelsAndColorsScheme"], re.I):
        DefaultColorScheme = "MuttedColorShades2"
    elif re.match("^BrightColorShades$", OptionsInfo["LevelsAndColorsScheme"], re.I):
        DefaultColorScheme = "BrightColorShades"
    else:
        MiscUtil.PrintError("The color scheme, %s, specified using \"--levelsAndColorsScheme\" option is not supported." % (OptionsInfo["LevelsAndColorsScheme"]))
    
    for Type in ["General", "Glycine", "Proline", "PreProline"]:
        PlotTypesInfo["Types"].append(Type)
        
        # Setup phi and psi density file...
        DensityFile = os.path.join(MayaChemToolsDataDir, "PhiPsiDensity%s.csv" % (Type))
        if not os.path.exists(DensityFile):
            MiscUtil.PrintError("The phi and psi density file file, %s, doesn't exist. This is required for generating contour plots.\n" % (DensityFile))
        PlotTypesInfo["PhiPsiDensityFiles"][Type] = DensityFile
        
        # Setup plot title...
        Title = Type
        if re.match("^PreProline$", Type, re.I):
            Title = "pre-Proline"
        PlotTypesInfo["Titles"][Type] = Title

        # Setup flags for drawing axis labels...
        DrawXLabel, DrawYLabel = [True] * 2
        if not OptionsInfo["MultipleOutFiles"]:
            # Turn off XLabel for plots in first row...
            DrawXLabel = False if re.match("^(General|Glycine)$", Type, re.I) else True
            
            # Turn off YLabel for plots in second column...
            DrawYLabel = False if re.match("^(Glycine|PreProline)$", Type, re.I) else True
        PlotTypesInfo["DrawXLabel"][Type] = DrawXLabel
        PlotTypesInfo["DrawYLabel"][Type] = DrawYLabel

        # Setup limits...
        (MinLimit, MaxLimit) = [-180, 180]
        PlotTypesInfo["Limits"][Type] = [MinLimit, MaxLimit]

        # Setup major tick labels and positions...
        MajorTickPositions = list(range(MinLimit, MaxLimit, OptionsInfo["TicksMajorInterval"]))
        MajorTickPositions.append(MaxLimit)
        MajorTickLabels = ["%s" % Position for Position in MajorTickPositions]
        PlotTypesInfo["MajorTickPositions"][Type] = MajorTickPositions
        PlotTypesInfo["MajorTickLabels"][Type] = MajorTickLabels

        # Setup minor tick positions without any labels...
        MinorTickPositions = list(range(MinLimit, MaxLimit, OptionsInfo["TicksMinorInterval"]))
        MinorTickPositions.append(MaxLimit)
        PlotTypesInfo["MinorTickPositions"][Type] = MinorTickPositions

        # Setup contour levels and colors...
        Levels = []
        Colors = []
        if re.match("^General$", Type, re.I):
            Levels = [0.0, 0.0005, 0.02, 1.0]
            Colors = ContourColorSchemes[Type][DefaultColorScheme]
        elif re.match("^Glycine$", Type, re.I):
            Levels = [0.0, 0.002, 0.02, 1.0]
            Colors = ContourColorSchemes[Type][DefaultColorScheme]
        elif re.match("^Proline$", Type, re.I):
            Levels = [0.0, 0.002,  0.02, 1.0]
            Colors = ContourColorSchemes[Type][DefaultColorScheme]
        elif re.match("^PreProline$", Type, re.I):
            Levels = [0.0, 0.002, 0.02, 1.0]
            Colors = ContourColorSchemes[Type][DefaultColorScheme]
        
        PlotTypesInfo["Levels"][Type] = Levels
        PlotTypesInfo["Colors"][Type] = Colors

    OptionsInfo["PlotTypesInfo"] = PlotTypesInfo

def SetupPlotTypesOutfiles():
    """Setup output file names for plot types """

    OptionsInfo["OutfilesList"] = []
    OptionsInfo["OutfilesList"].append(OptionsInfo["Outfile"])

    OptionsInfo["PlotTypesInfo"]["Outfiles"] = {}
    for PlotType in  OptionsInfo["PlotTypesInfo"]["Types"]:
        OptionsInfo["PlotTypesInfo"]["Outfiles"][PlotType] = None

    if not OptionsInfo["MultipleOutFiles"]:
        return

    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Outfile"])
    OutfileRoot = FileName
    OutfileExt = FileExt
    
    for PlotType in  OptionsInfo["PlotTypesInfo"]["Types"]:
        PlotOutfile = "%s_%s.%s" % (OutfileRoot, PlotType, OutfileExt)
        if os.path.exists(PlotOutfile):
            if not OptionsInfo["Overwrite"]:
                MiscUtil.PrintError("The plot output file, %s, already exist. Use option \"--ov\" or \"--overwrite\" and try again.\n" % (PlotOutfile))
                
        OptionsInfo["PlotTypesInfo"]["Outfiles"][PlotType] = PlotOutfile
        OptionsInfo["OutfilesList"].append(PlotOutfile)
        
def ProcessOptions():
    """Process and validate command line arguments and options"""

    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["OutMode"] = Options["--outMode"]
    OptionsInfo["MultipleOutFiles"] = True if re.match("^MultipleFiles$", OptionsInfo["OutMode"], re.I) else False
    MultipleOutFiles = OptionsInfo["MultipleOutFiles"]

    OptionsInfo["FigDPI"] = int(Options["--figDPI"])
    
    FigSize = Options["--figSize"]
    Width = 6.4
    Height = 4.8 if MultipleOutFiles else 6.4
    
    if not re.match("^auto$", FigSize, re.I):
        FigSizeWords = FigSize.split(",")
        Width = float(FigSizeWords[0])
        Height = float(FigSizeWords[1])
    OptionsInfo["FigSize"] = FigSize
    OptionsInfo["FigWidth"] = Width
    OptionsInfo["FigHeight"] = Height
    
    OptionsInfo["FontFamily"] = Options["--fontFamily"]
    OptionsInfo["FontAxesSize"] = Options["--fontAxesSize"]
    OptionsInfo["FontAxesWeight"] = Options["--fontAxesWeight"]
    OptionsInfo["FontTicksSize"] = Options["--fontTicksSize"]
    OptionsInfo["FontTicksWeight"] = Options["--fontTicksWeight"]
    OptionsInfo["FontTitleSize"] = Options["--fontTitleSize"]
    OptionsInfo["FontTitleWeight"] = Options["--fontTitleWeight"]
    
    OptionsInfo["Greek"] = True if re.match("^Yes$", Options["--greek"], re.I) else False
    
    OptionsInfo["Grid"] = True if re.match("^Yes$", Options["--grid"], re.I) else False
    OptionsInfo["GridLineColor"] = Options["--gridLineColor"]
    OptionsInfo["GridLineStyle"] = Options["--gridLineStyle"]
    OptionsInfo["GridLineWidth"] = float(Options["--gridLineWidth"])
    
    OptionsInfo["Infile"] = Options["--infile"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Infile"])
    OptionsInfo["InfileRoot"] = FileName
    
    OptionsInfo["LevelsAndColorsScheme"] = Options["--levelsAndColorsScheme"]
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Outfile"])
    OptionsInfo["OutfileRoot"] = FileName
    
    OptionsInfo["Overwrite"] = Options["--overwrite"]
    OptionsInfo["Precision"] = int(Options["--precision"])
    
    OptionsInfo["ScatterMarkerColor"] = Options["--scatterMarkerColor"]
    OptionsInfo["ScatterMarkerSize"] = float(Options["--scatterMarkerSize"])
    OptionsInfo["ScatterMarkerStyle"] = Options["--scatterMarkerStyle"]

    TicksMajorInterval = 90 if MultipleOutFiles else 180
    if not re.match("^auto$", Options["--ticksMajorInterval"], re.I):
        TicksMajorInterval = int(Options["--ticksMajorInterval"])
    OptionsInfo["TicksMajorInterval"] = TicksMajorInterval
        
    OptionsInfo["TicksMinor"] = True if re.match("^Yes$", Options["--ticksMinor"], re.I) else False
    TicksMinorInterval = 10 if MultipleOutFiles else 45
    if not re.match("^auto$", Options["--ticksMinorInterval"], re.I):
        TicksMinorInterval = int(Options["--ticksMinorInterval"])
    OptionsInfo["TicksMinorInterval"] = TicksMinorInterval
    
    RetrieveInfileInfo()
    OptionsInfo["ChainIDs"] = Options["--chainIDs"]
    
    ProcessChainIDs()
    
    OptionsInfo["LevelsAndColors"] = Options["--levelsAndColors"]
    SetupPlotsInfo()

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

    MiscUtil.ValidateOptionIntegerValue("--figDPI", Options["--figDPI"], {">": 0})
    if not re.match("^auto$", Options["--figSize"], re.I):
        MiscUtil.ValidateOptionNumberValues("--figSize", Options["--figSize"], 2, ",", "float", {">": 0})
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "pdb cif")

    MiscUtil.ValidateOptionTextValue("-g, --greek", Options["--greek"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--grid", Options["--grid"], "yes no")
    MiscUtil.ValidateOptionFloatValue("--gridLineWidth", Options["--gridLineWidth"], {">": 0})
    
    MiscUtil.ValidateOptionTextValue("--levelsAndColorsScheme", Options["--levelsAndColorsScheme"], "MuttedColorShades1 MuttedColorShades2 BrightColorShades")
    
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    
    MiscUtil.ValidateOptionTextValue("--outMode", Options["--outMode"], "SingleFile MultipleFiles")
    MiscUtil.ValidateOptionIntegerValue("-p, --precision", Options["--precision"], {">": 0})
    
    MiscUtil.ValidateOptionFloatValue("--scatterMarkerSize", Options["--scatterMarkerSize"], {">": 0})

    if not re.match("^auto$", Options["--ticksMajorInterval"], re.I):
        MiscUtil.ValidateOptionIntegerValue("--ticksMajorInterval", Options["--ticksMajorInterval"], {">": 0, "<": 360})
    
    MiscUtil.ValidateOptionTextValue("--ticksMinor", Options["--ticksMinor"], "yes no")
    if not re.match("^auto$", Options["--ticksMinorInterval"], re.I):
        MiscUtil.ValidateOptionIntegerValue("--ticksMinorInterval", Options["--ticksMinorInterval"], {">": 0, "<": 360})

# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLGenerateRamachandranPlots.py - Generate Ramachandran plots

Usage:
    PyMOLGenerateRamachandranPlots.py [--chainIDs <First, All or ID1,ID2...>]
                                      [--figDPI <number>] [--figSize <width, height>] [--fontFamily <text>]
                                      [--fontAxesSize <number or text>] [--fontAxesWeight <number or text>]
                                      [--fontTicksSize <number or text>] [--fontTicksWeight <number or text>]
                                      [--fontTitleSize <number or text>] [--fontTitleWeight <number or text>] [--greek <yes or no>]
                                      [--grid <yes or no>] [--gridLineColor <text>] [--gridLineStyle <text>] [--gridLineWidth <number>]
                                      [--levelsAndColors <PlotType:Level,color,Level,...;...>] [--levelsAndColorsScheme <text>]
                                      [--outMode <SingleFile or MultipleFiles>] [--overwrite]  [--precision <number>]
                                      [--scatterMarkerColor <text>] [--scatterMarkerSize <number>] [--scatterMarkerStyle <text>]
                                      [--ticksMajorInterval <number>] [--ticksMinor <yes or no>] [--ticksMinorInterval <number>]
                                      [-w <dir>] -i <infile> -o <outfile>
    PyMOLGenerateRamachandranPlots.py -h | --help | -e | --examples

Description:
    Generate Ramachandran plots for amino acid residues present in macromolecules.
    
    The Ramachandran plots are generated by plotting phi and psi backbone angles
    corresponding to the following four categories of amino acids:
    
        General: All residues except glycine, proline, or pre-proline
        Glycine: Only glycine residues
        Proline: Only proline residues
        PreProline: Only residues before proline not including glycine or
            proline
    
    In addition to the scatter plots for phi and psi angles, the filled contours
    are generated for the density of phi and psi angles [ Ref 144 ] for the
    Ramachandran plots. The contours are generated for "favored" and "allowed"
    regions. The phi and psi density is retrieved from the following density files
    available in MAYACHEMTOOLS/lib/data/ directory:
    
        General: PhiPsiDensityGeneral.csv 
        Glycine: PhiPsiDensityGlycine.csv
        Proline: PhiPsiDensityProline.csv
        PreProline: PhiPsiDensityPreProline.csv
    
    The supported input  file format are: PDB (.pdb), mmCIF (.cif)
    
    The output image file can be saved in any format supported by Python
    module Matplotlib. The image format is automatically detected from the
    output file extension. 
    
    Some of the most common output image file formats are: EPS (.eps), PDF (.pdf),
    PNG (.png), PS (.ps), SVG (.svg).

Options:
    -c, --chainIDs <First, All or ID1,ID2...>  [default: All]
        List of chain IDs to use for calculating phi and psi angles for residues
        in chains. Possible values: First, All, or a comma delimited list of chain
        IDs. The default is to use all chain IDs in input file.
    -e, --examples
        Print examples.
    --figDPI <number>  [default: 300]
        Figure resolution in dots per inches. The DPI value must be supported
        by Matplotlib during generation of an image of a specific format. No
        validation is performed.
    --figSize <width, height>  [default: auto]
        Figure dimensions in inches. The default values are dependent on the
        the value of '--outMode' option as shown below:
        
            SingleFile: 6.4, 6.4
            MultipleFiles: 6.4, 4.8
        
    --fontFamily <text>  [default: auto]
        Font family to use for title, axes labels, and tick marks. It must be a
        valid Matplotlib value. The default value corresponds to the value 
        plt.rcParams["font.family"] in your environment. For example: serif,
        sans-serif, cursive, etc. 
    --fontAxesSize <number or text>  [default: 10]
        Font size for labels on axes. It must be valid Matplotlib font size. For
        example: size in points, xx-small, x-small, small, medium, etc.
    --fontAxesWeight <number or text>  [default: regular]
        Font weight for labels on axes. It must be valid Matplotlib value. For
        example: a numeric value in range 0-1000, ultralight, light, normal,
        regular, book, medium, etc.
    --fontTicksSize <number or text>  [default: 8]
        Font size for tick labels. It must be a valid Matplotlib font size. For
        example: size in points, xx-small, x-small, small, medium, etc.
    --fontTicksWeight <number or text>  [default: regular]
        Font weight for tick labels. It must be valid Matplotlib value. For
        example: a numeric value in range 0-1000, ultralight, light,
        normal, regular, book, medium, etc.
    --fontTitleSize <number or text>  [default: 10]
        Font size for title. It must be a valid Matplotlib font size. For example:
        size in points, xx-small, x-small, small, medium, etc.
    --fontTitleWeight <number or text>  [default: bold]
        Font weight for title. It must be a valid Matplotlib value. For example: a
        numeric value in range 0-1000, ultralight, light, normal, regular, book,
        medium, etc.
    -g, --greek <yes or no>  [default: yes]
        Show phi and psi labels as greek characters.
    --grid <yes or no>  [default: yes]
        Display grid lines at major tick marks.
    --gridLineColor <text>  [default: #b0b0b0]
        Grid line color. It must be a valid Matplotlib value. The default color
        is light gray.
    --gridLineStyle <text>  [default: dotted]
        Grid line style. It must be a valid Matplotlib value. For example:
        '-' or 'solid', --' or 'dashed', '-.' or 'dashdot', ':' or 'dotted' etc.
    --gridLineWidth <number>  [default: 0.8]
        Grid line width. It must be a valid Matplotlib value.
    -h, --help
        Print this help message.
    -i, --infile <infile>
        Input file name.
    -l, --levelsAndColors <PlotType:Level,color,Level,...;...>  [default: auto]
        Semicolon delimited list of contour levels and colors for four types
        of Ramachandran plots.
        
        Three default contour level and color scheme may be specified by
        '--levelsAndColorsScheme' option. By default, the 'MuttedColorShades1'
        scheme is used. The default contour levels correspond to 'favored' and
        'allowed' regions [ Ref 144 ] for phi and psi angles.
        
        The colors are used to fill spaces between contour levels. The values
        for contour levels must be ascending order. The number of colors
        must be one less than the number contour levels.
        
        The format of contour level and color specification is as follows:
        
            PlotType:Level,Color,Level,...;PlotType:Level,Color,Level,...
        
        The valid values for plot type are:
        
            General, Glycine, Proline, or PreProline
        
        The contour level must be a number. The color value must be a valid color
        name or a hexadecimal color string supported by Matplotlib. No validation
        is performed.
        
        For example:
        
            General: 0.0, #FFFFFF, 0.0005, #EBF1DE, 0.02, #C3D69B, 1.0
        
    --levelsAndColorsScheme <text>  [default: MuttedColorShades1]
        Default contour levels and colors scheme.  Possible values:
        MuttedColorShades1, MuttedColorShades2, or BrightColorShades.
        
        This option is only used during 'auto' value of '--levelsAndColors' option.
        The default contour levels correspond to 'favored' and 'allowed' regions
        [ Ref 144 ] for phi and psi angles.
        
        The default contour and color values for different default schemes are
        shown below:
        
        MuttedColorShades1:
        
            General: 0.0, #FFFFFF, 0.0005, #EBF1DE, 0.02, #C3D69B, 1.0
            Glycine: 0.0, #FFFFFF, 0.002, #7FD9FF, 0.02, #FAC090, 1.0
            Proline: 0.0, #FFFFFF, 0.002, #E6E0EC, 0.02, #B3A2C7, 1.0
            PreProline: 0.0, #FFFFFF, 0.002, #DCE6F2, 0.02, #95B3D7, 1.0
        
        MuttedColorShades2:
        
            General: 0.0, #FFFFFF, 0.0005, #EBF1DE, 0.02, #D7E4BD, 1.0
            Glycine: 0.0, #FFFFFF, 0.002, #FDEADA, 0.02, #FCD5B5, 1.0
            Proline: 0.0, #FFFFFF, 0.002, #E6E0EC, 0.02, #CCC1DA, 1.0
            PreProline: 0.0, #FFFFFF, 0.002, #DCE6F2, 0.02, #B9CDE5, 1.0
        
        BrightColorShades: [ Ref 145 ]
        
            General: 0.0, #FFFFFF, 0.0005, #B3E8FF, 0.02, #7FD9FF, 1.0
            Glycine: 0.0, #FFFFFF, 0.002, #FFE8C5, 0.02, #FFCC7F, 1.0
            Proline: 0.0, #FFFFFF, 0.002, #D0FFC5, 0.02, #7FFF8C, 1.0
            PreProline: 0.0, #FFFFFF, 0.002, #B3E8FF, 0.02, #7FD9FF, 1.0
        
    -o, --outfile <outfile>
        Output image file name.
        
        A set of output files is optionally generated for 'MultipleFiles' value of
        '--outMode' option. The names of these output files are automatically
        generated from the the name of the specified output file as shown
        below:
        
            General: <OutfileRoot>_General.<OutfileExt>
            Glycine: <OutfileRoot>_Glycine.<OutfileExt>
            Proline: <OutfileRoot>_Proline.<OutfileExt>
            PreProline: <OutfileRoot>_PreProline.<OutfileExt>
        
    --outMode <Single or Multiple>  [default: SingleFile]
        A single output file containing all four Ramachandran plots or multiple
        output files corresponding to different types of Ramachandran plots.
        
        The phi and psi angles are categorized into the following groups
        corresponding to four types of Ramachandran plots:
        
            General: All residues except glycine, proline, or pre-proline
            Glycine: Only glycine residues
            Proline: Only proline residues
            PreProline: Only residues before proline not including glycine or
                proline
        
    --overwrite
        Overwrite existing files.
    -p, --precision <number>  [default: 2]
        Floating point precision for plotting the calculated phi and psi angles.
    --scatterMarkerColor <text>  [default: #1f77b4]
        Scatter marker color for plotting to phi and psi angles. It must be a
        valid Matplotlib value. The default color is dark blue.
    --scatterMarkerSize <number>  [default: 1.0]
        Scatter marker size for piloting phi and psi angles. It must be a valid
        Matplotlib value.
    --scatterMarkerStyle <text>  [default: .]
        Scatter marker style for piloting phi and psi angles. It must be a valid
        Matplotlib value. For example: '.' (point), ',' (pixel), 'o' (circle), etc.
    --ticksMajorInterval <number>  [default: auto]
        Display major marks on axes at intervals specified in degrees for phi and
        psi angles. The default value is dependent on the the value of '--outMode'
        option: SingleFile: 180; MultipleFiles: 90
        
        The grid lines are drawn at the locations of major tick marks.
    --ticksMinor <yes or no>  [default: yes]
        Display minor tick marks. The major tick mark are always displayed.
    --ticksMinorInterval <number>  [default: auto]
        Display minor marks on axes at intervals specified in degrees for phi and
        psi angles. The default value is dependent on the the value of '--outMode'
        option: SingleFile:  45; MultipleFiles: 10
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To generate Ramachandran plot for all residues across all chains in input
    file and write out a single SVG file containing all four types of plots, type:

        % PyMOLGenerateRamachandranPlots.py -i Sample3.pdb -o Sample3Out.svg

    To generate Ramachandran plot for all residues across all chains in input
    file and write out four SVG files corresponding to four types of plots, type:

        % PyMOLGenerateRamachandranPlots.py --outMode MultipleFiles
          -i Sample3.pdb -o Sample3Out.svg

    To generate Ramachandran plot for all residues in a specific chain in input
    file and write out a single PDF file containing all four types of plots, type:

        % PyMOLGenerateRamachandranPlots.py -c E -i Sample3.pdb
          -o Sample3Out.pdf

    To generate Ramachandran plot for all residues across all chains in input
    file using specific options and write out four PNG files containing all four
    types of plots, type:

        % PyMOLGenerateRamachandranPlots.py --outMode MultipleFiles
          --figSize "6,4" --figDPI 600 --fontTitleSize 10 --fontTitleWeight
          normal --greek no --grid no --levelsAndColors
          "General: 0.0, #FFFFFF, 0.0005, #B3E8FF, 0.02, #7FD9FF, 1.0"
          -i Sample3.pdb -o Sample3Out.png

Author:
    Manish Sud(msud@san.rr.com)

See also:
    DownloadPDBFiles.pl, PyMOLCalculatePhiPsiAngles.py, PyMOLCalculateRMSD.py,
    PyMOLCalculateProperties.py

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
