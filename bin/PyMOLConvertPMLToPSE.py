#!/bin/env python
#
# File: PyMOLConvertPMLToPSE.py
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
    GeneratePSEFile()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def GeneratePSEFile():
    """Comvert PML file to PSE file."""

    PMLFile = OptionsInfo["Infile"]
    PSEFile = OptionsInfo["Outfile"]
    
    MiscUtil.PrintInfo("\nGenerating file %s..." % PSEFile)

    PyMOLUtil.ConvertPMLFileToPSEFile(PMLFile, PSEFile, OutputFeedback = OptionsInfo["Feedback"])
    
    if not os.path.exists(PSEFile):
        MiscUtil.PrintWarning("Failed to generate PSE file, %s..." % (PSEFile))

def ProcessOptions():
    """Process and validate command line arguments and options"""

    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["Feedback"] = True if re.match("^Yes$", Options["--feedback"], re.I) else False
    
    OptionsInfo["Infile"] = Options["--infile"]
    OptionsInfo["Outfile"] = Options["--outfile"]

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

    MiscUtil.ValidateOptionTextValue("-f, --feedback", Options["--feedback"], "yes no")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "pml")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "pse")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])

# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLConvertPMLToPSE.py - Convert PML to PSE

Usage:
    PyMOLConvertPMLToPSE.py [--feedback <yes or no>] [--overwrite]
                            [-w <dir>] -i <infile> -o <outfile>
    PyMOLConvertPMLToPSE.py -h | --help | -e | --examples

Description:
    Convert PyMOL script language (PML) file to PyMOL session (PSE) file.

    The supported input and output file formats are PML (.pml) and PSE (.pse).

Options:
    -f, --feedback <yes or no>  [default: yes]
        PyMOL output feedback during loading of PML file. This option may not
        work in all versions of PyMOL across various platforms.
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
    To convert a PML file to a PSE file, type:

        % PyMOLConvertPMLToPSE.py -i Sample.pml -o Sample.pse

    To convert a PML file to a PSE file along with turning off PyMOL feedback
    during loading of PML file, type:

        % PyMOLConvertPMLToPSE.py -f no -i Sample.pml -o Sample.pse

Author:
    Manish Sud(msud@san.rr.com)

See also:
    PyMOLConvertLigandFileFormat.py, PyMOLSplitChainsAndLigands.py,
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
