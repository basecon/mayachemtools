#!/usr/bin/perl -w
#
# File: DownloadPDBFiles.pl
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2020 Manish Sud. All rights reserved.
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

use strict;
use FindBin; use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use File::Basename;
use File::Fetch;
use File::Copy;
use Text::ParseWords;
use Benchmark;
use FileUtil;
use TextUtil;
use PDBFileUtil;

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

# Starting message...
$ScriptName = basename($0);
print "\n$ScriptName: Starting...\n\n";
$StartTime = new Benchmark;

# Get the options and setup script...
SetupScriptUsage();
if ($Options{help} || @ARGV < 1) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

# Process options...
print "Processing options...\n";
my(%OptionsInfo, %PDBIDsFileInfo);
ProcessOptions();

# Collect PDB IDs and download corresponding files...
my(%PDBFilesInfo);
SetupPDBFilesInfo();
DownloadPDBFiles();

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Download appropriate PDB fies...
sub DownloadPDBFiles {
  my($PDBIDCount, $PDBIDOkayCount, $PDBIDFailedCount, $PDBIDsIgnoredCount, $PDBID, $PDBStatus, $CIFStatus, $DownloadDensityMap, $CryoEMDataStatus, $EMDBID, $TotalEDMapOkayCount , $TotalEDMapFailedCount, $EDMapOkayCount, $EDMapFailedCount, $TotalCryoEMMapOkayCount, $TotalCryoEMMapFailedCount, $CryoEMMapOkayCount, $CryoEMMapFailedCount);

  print "\nDownloading PDB files...\n";

  ($PDBIDCount, $PDBIDOkayCount, $PDBIDFailedCount, $PDBIDsIgnoredCount, $TotalEDMapOkayCount , $TotalEDMapFailedCount, $TotalCryoEMMapOkayCount, $TotalCryoEMMapFailedCount) = (0) x 8;

  $DownloadDensityMap = $OptionsInfo{DensityMap};

  # Turn off warnings from File::Fetch
  $File::Fetch::WARN = 0;

  PDBID: for $PDBID (@{$PDBFilesInfo{PDBIDs}}) {
    $PDBIDCount++;

    print "\nProcessing PDB ID $PDBID...\n";

    if ($PDBID =~ /\./) {
      $PDBIDsIgnoredCount++;
      warn "Warning: Ignoring invalid PDB ID $PDBID\n";
      next PDBID;
    }

    $PDBStatus = 0;
    # Download PDB format file...
    if ($PDBFilesInfo{DownloadPDB}{$PDBID}) {
      print "Downloading PDB file: $PDBFilesInfo{RemoteFilePDBFormat}{$PDBID}\n";
      $PDBStatus = DownloadFile($PDBFilesInfo{RemoteFilePDBFormat}{$PDBID}, $PDBFilesInfo{LocalFilePDBFormat}{$PDBID});
    }

    # Try downloading CIF format file...
    $CIFStatus = 0;
    if ($PDBFilesInfo{DownloadCIF}{$PDBID} && !$PDBStatus) {
      print "Downloading PDB file: $PDBFilesInfo{RemoteFileCIFFormat}{$PDBID}\n";
      $CIFStatus = DownloadFile($PDBFilesInfo{RemoteFileCIFFormat}{$PDBID}, $PDBFilesInfo{LocalFileCIFFormat}{$PDBID});
    }

    # Check status of download...
    if (!($PDBStatus || $CIFStatus)) {
      $PDBIDFailedCount++;
      next PDBID;
    }
    $PDBIDOkayCount++;

    # Any need to download density files...
    if (!$DownloadDensityMap) {
      next PDBID;
    }

    # Check whether it's cryo-EM data file...
    #
    $CryoEMDataStatus = 0;
    $EMDBID = 0;
    if ($PDBStatus) {
      ($CryoEMDataStatus, $EMDBID) = RetrieveEMDBIDFromPDBFile($PDBFilesInfo{LocalFilePDBFormat}{$PDBID});
    }
    elsif ($CIFStatus) {
      ($CryoEMDataStatus, $EMDBID) = RetrieveEMDBIDFromCIFFile($PDBFilesInfo{LocalFileCIFFormat}{$PDBID});
    }

    ($EDMapOkayCount, $EDMapFailedCount) = DownloadEDMapFiles($PDBID, $CryoEMDataStatus);
    $TotalEDMapOkayCount += $EDMapOkayCount;
    $TotalEDMapFailedCount += $EDMapFailedCount,

    ($CryoEMMapOkayCount, $CryoEMMapFailedCount) = DownloadCryoEMMapFiles($PDBID, $CryoEMDataStatus, $EMDBID);
    $TotalCryoEMMapOkayCount += $CryoEMMapOkayCount;
    $TotalCryoEMMapFailedCount += $CryoEMMapFailedCount,
  }

  print "\nTotal number of PDB IDs:  $PDBIDCount\n";
  print "Number of PDB IDs ignored:  $PDBIDsIgnoredCount\n";

  print "\nNumber of successful downloads:  $PDBIDOkayCount\n";
  print "Number of failed downloads:  $PDBIDFailedCount\n";

  if ($DownloadDensityMap) {
    print "\nNumber of successful ED map downloads:  $TotalEDMapOkayCount\n";
    print "Number of failed ED map downloads:  $TotalEDMapFailedCount\n";

    print "\nNumber of successful cryo-EM map downloads:  $TotalCryoEMMapOkayCount\n";
    print "Number of failed cryo-EM map downloads:  $TotalCryoEMMapFailedCount\n";
  }
}

# Download x-ray electron density files...
sub DownloadEDMapFiles {
  my($PDBID, $CryoEMDataStatus) = @_;
  my($Index, $RemoteEDMapFile, $LocalEDMapFile, $TmpLocalEDMapFile, $FinalLocalEDMapFile, $EDMapFailedCount, $EDMapOkayCount, $Status);

  ($EDMapOkayCount, $EDMapFailedCount) = (0) x 2;

  if ($CryoEMDataStatus) {
    print "Skipping download of x-ray electron density files for cryo-EM PDB data...\n";
    return ($EDMapOkayCount, $EDMapFailedCount);
  }

 EDMAPFILE: for $Index (0 .. $#{$PDBFilesInfo{RemoteEDMapFiles}{$PDBID}}) {
    $RemoteEDMapFile = $PDBFilesInfo{RemoteEDMapFiles}{$PDBID}[$Index];
    $LocalEDMapFile = $PDBFilesInfo{LocalEDMapFiles}{$PDBID}[$Index];
    $TmpLocalEDMapFile = $PDBFilesInfo{TmpLocalEDMapFiles}{$PDBID}[$Index];
    $FinalLocalEDMapFile = $PDBFilesInfo{FinalLocalEDMapFiles}{$PDBID}[$Index];

    print "Downloading x-ray electron density map file: $RemoteEDMapFile\n";
    $Status = DownloadFile($RemoteEDMapFile, $LocalEDMapFile);

    if (!$Status) {
      $EDMapFailedCount++;
      next EDMAPFILE;
    }
    # Rename downloaded ED file...
    print "Moving file from $LocalEDMapFile to $FinalLocalEDMapFile\n";
    move $LocalEDMapFile, $TmpLocalEDMapFile or warn "Warning: Couldn't move file $LocalEDMapFile to $TmpLocalEDMapFile\n";
    move $TmpLocalEDMapFile, $FinalLocalEDMapFile or warn "Warning: Couldn't move file $TmpLocalEDMapFile to $FinalLocalEDMapFile\n";
    $EDMapOkayCount++;
  }

  return ($EDMapOkayCount, $EDMapFailedCount);
}

# Download cryo-EM density files...
sub DownloadCryoEMMapFiles {
  my($PDBID, $CryoEMDataStatus, $EMDBID) = @_;
  my($Index, $RemoteCryoEMMapFile, $LocalCryoEMMapFile, $CryoEMMapFailedCount, $CryoEMMapOkayCount, $Status, $FileType);

  ($CryoEMMapOkayCount, $CryoEMMapFailedCount) = (0) x 2;

  if (!$CryoEMDataStatus) {
    print "Skipping download of cryo-EM density files for non cryo-EM PDB data...\n";
    return ($CryoEMMapOkayCount, $CryoEMMapFailedCount);
  }

  CRYOEMMAPFILE: for $Index (0 .. $#{$PDBFilesInfo{RemoteCyroEMMapFiles}{$PDBID}}) {
    $FileType = $PDBFilesInfo{RemoteCyroEMMapFileTypes}{$PDBID}[$Index];
    $RemoteCryoEMMapFile = $PDBFilesInfo{RemoteCyroEMMapFiles}{$PDBID}[$Index];
    $LocalCryoEMMapFile = $PDBFilesInfo{LocalCyroEMMapFiles}{$PDBID}[$Index];

    # Update file names with actual EMDBID...
    $RemoteCryoEMMapFile =~ s/EMDBIDPlaceHolder/$EMDBID/ig;
    $LocalCryoEMMapFile =~ s/EMDBIDPlaceHolder/$EMDBID/ig;

    print "Downloading cryo-EM density $FileType file: $RemoteCryoEMMapFile\n";
    $Status = DownloadFile($RemoteCryoEMMapFile, $LocalCryoEMMapFile);

    if (!$Status) {
      $CryoEMMapFailedCount++;
      next CRYOEMMAPFILE;
    }
    $CryoEMMapOkayCount++;
  }

  return ($CryoEMMapOkayCount, $CryoEMMapFailedCount);
}


# Download specified file...
sub DownloadFile {
  my($RemoteFileURL, $LocalFileName) = @_;
  my($Status, $FileFetch, $FetchedFilePath);

  $Status = 1;

  # Setup a fetch object...
  $FileFetch = File::Fetch->new(uri => $RemoteFileURL);

  # Fetch file to the CWD...
  $FetchedFilePath = $FileFetch->fetch();

  if (IsEmpty($FetchedFilePath)) {
    warn "Warning: Download failed for file $RemoteFileURL: " . $FileFetch->error() . "\n";
    if (-e $LocalFileName) {
      warn "Warning: Deleting empty file $LocalFileName\n";
      unlink $LocalFileName or warn "Warning: Couldn't delete file $LocalFileName\n";
    }
    $Status = 0;
  }
  return $Status;
}

# Collect specified PDB IDs along with settting up PDB and ED file names for all
# specified PDB IDs...
#
sub SetupPDBFilesInfo {

  %PDBFilesInfo = ();
  RetrievePDBIDs();
  SetupPDBandEDFileNames();
}

# Retrieve EMDB ID from PDB file...
sub RetrieveEMDBIDFromPDBFile {
  my($PDBFile) = @_;
  my($EMDBID, $CryoEMDataType, $Line);

  $EMDBID = 0;

  if (!-e $PDBFile) {
    return $EMDBID;
  }

  $CryoEMDataType = 0;

  open PDBFILE, "$PDBFile" or die "Couldn't open $PDBFile: $! \n";
  LINE: while ($Line = GetTextLine(\*PDBFILE)) {
    if ($Line =~ /^EXPDTA/i) {
      if ($Line =~ /ELECTRON MICROSCOPY/i) {
         $CryoEMDataType = 1;
      }
    }
    elsif ($Line =~ /^REMARK/i) {
      if ($Line =~ /DB: EMDB/i) {
         (undef, $EMDBID, undef) = ($Line =~ /^(.*?) EMD-([0-9]+) (.*?)$/);
         if (!defined $EMDBID) {
           $EMDBID = 0;
         }
         last LINE;
      }
    }
  }
  close PDBFILE;

  return ($CryoEMDataType, $EMDBID);
}

# Retrieve EMDB ID from CIF file...
sub RetrieveEMDBIDFromCIFFile {
  my($CIFFile) = @_;
  my($EMDBID, $CryoEMDataType, $Line);

  $EMDBID = 0;

  if (!-e $CIFFile) {
    return $EMDBID;
  }

  $CryoEMDataType = 0;

  open CIFFILE, "$CIFFile" or die "Couldn't open $CIFFile: $! \n";
  LINE: while ($Line = GetTextLine(\*CIFFILE)) {
    if ($Line =~ /^_exptl.method/i) {
      if ($Line =~ /ELECTRON MICROSCOPY/i) {
         $CryoEMDataType = 1;
         last LINE;
      }
    }
    elsif ($Line =~ /^EMDB  EMD-/i) {
      (undef, $EMDBID, undef) = ($Line =~ /^(.*?) EMD-([0-9]+)(.*?)$/);
       if (!defined $EMDBID) {
         $EMDBID = 0;
       }
    }
  }
  close CIFFILE;

  return ($CryoEMDataType, $EMDBID);
}


# Set up PDB and ED file names for downloading....
sub SetupPDBandEDFileNames {
  my($PDBID, $DownloadDensityMap, $PDBDataLocationURL, $EDMapDataLocaltionURL, $EDMapType, $EDMapPDBID, $EDMapSuffix, $EDMapFileExt, $CryoEMDataLocationURL, $EMDBID);

  @{$PDBFilesInfo{PDBIDs}} = ();

  %{$PDBFilesInfo{DownloadPDB}} = ();
  %{$PDBFilesInfo{RemoteFilePDBFormat}} = ();
  %{$PDBFilesInfo{LocalFilePDBFormat}} = ();

  %{$PDBFilesInfo{DownloadCIF}} = ();
  %{$PDBFilesInfo{RemoteFileCIFFormat}} = ();
  %{$PDBFilesInfo{LocalFileCIFFormat}} = ();

  # Initilaize X-ray electron density file names...
  %{$PDBFilesInfo{DownloadEDMap}} = ();
  %{$PDBFilesInfo{RemoteEDMapFiles}} = ();
  %{$PDBFilesInfo{LocalEDMapFiles}} = ();
  %{$PDBFilesInfo{TmpLocalEDMapFiles}} = ();
  %{$PDBFilesInfo{FinalLocalEDMapFiles}} = ();

  # Initilaize cryo-EM  density file names...
  %{$PDBFilesInfo{DownloadCryoEMMap}} = ();
  %{$PDBFilesInfo{RemoteCyroEMMapFileTypes}} = ();
  %{$PDBFilesInfo{RemoteCyroEMMapFiles}} = ();
  %{$PDBFilesInfo{LocalCyroEMMapFiles}} = ();

  $DownloadDensityMap = $OptionsInfo{DensityMap};

  $PDBDataLocationURL = $OptionsInfo{DataLocationURL};
  if ($PDBDataLocationURL !~ /\/$/) {
    $PDBDataLocationURL .= "/";
  }

  $EDMapDataLocaltionURL = $OptionsInfo{DenistyMapLocationURLXRay};
  if ($EDMapDataLocaltionURL !~ /\/$/) {
    $EDMapDataLocaltionURL .= "/";
  }

  $CryoEMDataLocationURL = $OptionsInfo{DensityMapLocationURLCryoEM};
  if ($CryoEMDataLocationURL !~ /\/$/) {
    $CryoEMDataLocationURL .= "/";
  }

  PDBID: for $PDBID (@{$OptionsInfo{PDBIDs}}) {
    # Track PDB IDs..
    push @{$PDBFilesInfo{PDBIDs}}, $PDBID;

    # Intialize PDB file names...
    $PDBFilesInfo{DownloadPDB}{$PDBID} = 0;
    $PDBFilesInfo{RemoteFilePDBFormat}{$PDBID} = "";
    $PDBFilesInfo{LocalFilePDBFormat}{$PDBID} = "";

    $PDBFilesInfo{DownloadCIF}{$PDBID} = 0;
    $PDBFilesInfo{RemoteFileCIFFormat}{$PDBID} = "";
    $PDBFilesInfo{LocalFileCIFFormat}{$PDBID} = "";

    if ($OptionsInfo{PDBFormat} =~ /^(PDB|Auto)$/i) {
      $PDBFilesInfo{DownloadPDB}{$PDBID} = 1;
      $PDBFilesInfo{RemoteFilePDBFormat}{$PDBID} = "${PDBDataLocationURL}${PDBID}.pdb";
      $PDBFilesInfo{LocalFilePDBFormat}{$PDBID} = "${PDBID}.pdb";
    }
    if ($OptionsInfo{PDBFormat} =~ /^(CIF|Auto)$/i) {
      $PDBFilesInfo{DownloadCIF}{$PDBID} = 1;
      $PDBFilesInfo{RemoteFileCIFFormat}{$PDBID} = "${PDBDataLocationURL}${PDBID}.cif";
      $PDBFilesInfo{LocalFileCIFFormat}{$PDBID} = "${PDBID}.cif";
    }

    # Initialize x-ray ED map file names...
    $PDBFilesInfo{DownloadEDMap}{$PDBID} = 0;
    @{$PDBFilesInfo{RemoteEDMapFiles}{$PDBID}} = ();
    @{$PDBFilesInfo{LocalEDMapFiles}{$PDBID}} = ();
    @{$PDBFilesInfo{TmpLocalEDMapFiles}{$PDBID}} = ();
    @{$PDBFilesInfo{FinalLocalEDMapFiles}{$PDBID}} = ();

    # Initialize cryo-EM map file names...
    $PDBFilesInfo{DownloadCryoEMMap}{$PDBID} = 0;
    @{$PDBFilesInfo{RemoteCyroEMMapFileTypes}{$PDBID}} = ();
    @{$PDBFilesInfo{RemoteCyroEMMapFiles}{$PDBID}} = ();
    @{$PDBFilesInfo{LocalCyroEMMapFiles}{$PDBID}} = ();

    if (!$DownloadDensityMap) {
      next PDBID;
    }

    # Setup x-ray ED file names...
    if ($OptionsInfo{DensityMapMode} =~ /^(XRayElectronDensity|Auto)$/i) {
      $PDBFilesInfo{DownloadEDMap}{$PDBID} = 1;

      $EDMapPDBID = lc $PDBID;

      for $EDMapType (@{$OptionsInfo{EDMapTypesList}}) {
        $EDMapSuffix = $OptionsInfo{EDMapLocationSuffixesMap}{$EDMapType};
        $EDMapFileExt = $OptionsInfo{EDMapLocationFileExtMap}{$EDMapType};

        push @{$PDBFilesInfo{RemoteEDMapFiles}{$PDBID}}, "${EDMapDataLocaltionURL}${EDMapPDBID}${EDMapSuffix}.${EDMapFileExt}";
        push @{$PDBFilesInfo{LocalEDMapFiles}{$PDBID}}, "${EDMapPDBID}${EDMapSuffix}.${EDMapFileExt}";

        push @{$PDBFilesInfo{TmpLocalEDMapFiles}{$PDBID}}, "${EDMapPDBID}${EDMapSuffix}Tmp.${EDMapFileExt}";
        push @{$PDBFilesInfo{FinalLocalEDMapFiles}{$PDBID}}, "${PDBID}${EDMapSuffix}.${EDMapFileExt}";
      }
    }
    if ($OptionsInfo{DensityMapMode} =~ /^(CryoEMDensity|Auto)$/i) {
      # Set up cryo-EM map file names using "EMDBIDPlaceHolder" to be replaced later by
      # a valid ID retrieved from PDB or CIF file...
      #
      $EMDBID = "EMDBIDPlaceHolder";
      $PDBFilesInfo{DownloadCryoEMMap}{$PDBID} = 1;

      # Map files...
      push @{$PDBFilesInfo{RemoteCyroEMMapFileTypes}{$PDBID}}, "map";
      push @{$PDBFilesInfo{RemoteCyroEMMapFiles}{$PDBID}}, "${CryoEMDataLocationURL}EMD-${EMDBID}/map/emd_${EMDBID}.map.gz";
      push @{$PDBFilesInfo{LocalCyroEMMapFiles}{$PDBID}}, "emd_${EMDBID}.map.gz";

      # Metadata files...
      push @{$PDBFilesInfo{RemoteCyroEMMapFileTypes}{$PDBID}}, "header";
      push @{$PDBFilesInfo{RemoteCyroEMMapFiles}{$PDBID}}, "${CryoEMDataLocationURL}EMD-${EMDBID}/header/emd-${EMDBID}.xml";
      push @{$PDBFilesInfo{LocalCyroEMMapFiles}{$PDBID}}, "emd-${EMDBID}.xml";

    }
  }
}

# Collect PDB IDs...
#
sub RetrievePDBIDs {
  @{$OptionsInfo{PDBIDs}} = ();

  if ($OptionsInfo{Mode} =~ /^IDsOnCmdLine$/i) {
    RetriveCommandLinePDBIDs();
  }
  elsif ($OptionsInfo{Mode} =~ /^IDsInFile$/i) {
    RetriveTextFilePDBIDs();
  }
}

# Collect PDB IDs specified on the command line...
#
sub RetriveCommandLinePDBIDs {
  my($SpecifiedPDBID, @PDBIDs, @ProcessedPDBIDs);

  print "\nProcessing PDB ID(s) from command line...\n";

  @PDBIDs = ();
  for $SpecifiedPDBID (@{$OptionsInfo{CmdLinePDBIDs}}) {
    @ProcessedPDBIDs = ProcessSpecifiedPDBIDs($SpecifiedPDBID);
    if (@ProcessedPDBIDs) {
      push @PDBIDs, @ProcessedPDBIDs;
    }
  }
  @{$OptionsInfo{PDBIDs}} = @PDBIDs;

  RetrieveUniquePDBIDs();
}

# Collect PDB IDs specified in the text file...
#
sub RetriveTextFilePDBIDs {
  my($TextFile, $InDelim, $IDsColIndex, $LineCount, $ProcessedLineCount, $IgnoredLineCount, $PDBID, $Line, @PDBIDs, @ProcessedPDBIDs, @LineWords);

  $TextFile = $PDBIDsFileInfo{Name};

  $IDsColIndex = $PDBIDsFileInfo{IDsColIndex};
  $InDelim = $PDBIDsFileInfo{InDelim} ;

  ($LineCount, $ProcessedLineCount, $IgnoredLineCount) = (0) x 3;

  print "\nProcessing PDB ID(s) from PDB IDs  file $TextFile...\n";

  open TEXTFILE, "$TextFile" or die "Couldn't open $TextFile: $! \n";
  # Skip label line...
  $_ = <TEXTFILE>;

  @PDBIDs = ();
  LINE: while ($Line = GetTextLine(\*TEXTFILE)) {
    $LineCount++;
    @LineWords = quotewords($InDelim, 0, $Line);

    if ($IDsColIndex >= scalar @LineWords) {
      $IgnoredLineCount++;
      warn "Warning: Ignoring line number $LineCount: PDB IDs column number, ". ($IDsColIndex + 1) . ", doesn't exist in the line containing, " . (scalar @LineWords) .  ", columns.\nLine: $Line\n";
      next LINE;
    }
    $PDBID = $LineWords[$IDsColIndex];
    if (IsEmpty($PDBID )) {
      $IgnoredLineCount++;
      warn "Warning: Ignoring line number $LineCount: PDB ID value is empty.\nLine: $Line\n";
      next LINE;
    }
    $ProcessedLineCount++;

    @ProcessedPDBIDs = ProcessSpecifiedPDBIDs($PDBID);
    if (@ProcessedPDBIDs) {
      push @PDBIDs, @ProcessedPDBIDs;
    }
  }
  close TEXTFILE;

  @{$OptionsInfo{PDBIDs}} = @PDBIDs;

  print "\nTotal number of lines in PDB IDs text file: $LineCount\n";
  print "Total number of lines processed: $ProcessedLineCount\n";
  print "Total number of lines ignored: $IgnoredLineCount\n";

  RetrieveUniquePDBIDs();
}

# Process specified PDB IDs...
#
# Notes:
#   . Commas and spaces in the specification of PBD IDs are allowed.
#   . All PDB IDs are turned into uppercase letters.
#
sub ProcessSpecifiedPDBIDs {
  my($SpecifiedPDBID) = @_;
  my($PDBID, @PDBIDWords, @PDBIDs);

  $SpecifiedPDBID = RemoveLeadingAndTrailingWhiteSpaces($SpecifiedPDBID);
   if ($SpecifiedPDBID =~ / /) {
    @PDBIDWords = split " ",  $SpecifiedPDBID;
  }
  elsif ($SpecifiedPDBID =~ /,/) {
    @PDBIDWords = split ",",  $SpecifiedPDBID;
  }
  else {
    push @PDBIDWords, $SpecifiedPDBID;
  }

  @PDBIDs = ();
  for $PDBID (@PDBIDWords) {
    $PDBID =~ s/( |,)//g;
    push @PDBIDs, uc $PDBID;
  }
  return @PDBIDs;
}

# Collect unique PDB IDs...
sub RetrieveUniquePDBIDs {
  my($PDBID, $PDBIDsCount, $PDBIDsIgnoredCount, @UniquePDBIDs, %PDBIDsMap);

  %PDBIDsMap = ();
  @UniquePDBIDs = ();

  $PDBIDsCount = 0;
  $PDBIDsIgnoredCount = 0;

  PDBID: for $PDBID (@{$OptionsInfo{PDBIDs}}) {
   $PDBIDsCount++;
    if (exists $PDBIDsMap{$PDBID}) {
      $PDBIDsIgnoredCount++;
      warn "Warning: Ignoring duplicate PDB ID $PDBID\n";
      next PDBID;
    }
    $PDBIDsMap{$PDBID} = $PDBID;
    push @UniquePDBIDs, $PDBID;
  }
  @{$OptionsInfo{PDBIDs}} = @UniquePDBIDs;
  print "\nTotal number of PDB IDs:  $PDBIDsCount\n";
  print "Number of duplicate PDB IDs ignored:  $PDBIDsIgnoredCount\n";
}

# Process option values...
sub ProcessOptions {
  my($EDMapTypes, $EDMapType, $EDLocationSuffixes, $EDMapSuffix, $Index, @EDMapTypesList, @EDLocationSuffixesList, %EDLocationSuffixesMap);

  %OptionsInfo = ();
  %PDBIDsFileInfo = ();

  $OptionsInfo{Mode} = $Options{mode};
  $OptionsInfo{ColMode} = $Options{colmode};

  $OptionsInfo{DataLocationURL} = $Options{datalocationurl};
  if (IsEmpty($OptionsInfo{DataLocationURL} )) {
    die "Error: PDB data location URL specified using \"-d, --dataLocationURL\" is empty. Allowed value: Non empty string\n";
  }

  $OptionsInfo{DensityMap} = $Options{densitymap} =~ /^Yes$/i ? 1 : 0;
  $OptionsInfo{DensityMapMode} = $Options{densitymapmode};

  $OptionsInfo{DensityMapLocationURLCryoEM} = $Options{densitymaplocationurlcryoem};
  $OptionsInfo{DenistyMapLocationURLXRay} = $Options{denistymaplocationurlxray};

  # Process x-ray ED map location file suffixes...
  $EDLocationSuffixes = $Options{edmaplocationsuffixes};
  $EDLocationSuffixes =~ s/ //g;
  %EDLocationSuffixesMap = ();

  @EDLocationSuffixesList = split(",", $EDLocationSuffixes);
  if (@EDLocationSuffixesList % 2) {
    die "Invalid number  of values specified using \"--EDMapLocationSuffixes\" option: It must contain even number of valid values.\n";
  }
  for ($Index = 0; $Index < @EDLocationSuffixesList; $Index += 2) {
    $EDMapType = $EDLocationSuffixesList[$Index];
    $EDMapSuffix = $EDLocationSuffixesList[$Index + 1];

    if ($EDMapType !~ /^(CompositeMap|DifferenceMap|ReflectionMap)$/i) {
      die "Error: The value specified, $EDMapType, for option \"--EDMapLocationSuffixes\" is not valid. Allowed values: CompositeMap, DifferenceMap, ReflectionMap\n";
    }
    if (exists $EDLocationSuffixesMap{$EDMapType}) {
      die "Error: Duplicate ED map type, $EDMapType, specified for option \"--EDMapLocationSuffixes\"\n";
    }

    # Track suffixes...
    if ($EDMapSuffix =~ /^None$/i) {
      $EDMapSuffix = "";
    }
    $EDLocationSuffixesMap{$EDMapType} = $EDMapSuffix;
  }
  $OptionsInfo{EDMapLocationSuffixes} = $EDLocationSuffixes;
  @{$OptionsInfo{EDMapLocationSuffixesList}} = ();
  @{$OptionsInfo{EDMapLocationSuffixesList}} = @EDLocationSuffixesList;

  %{$OptionsInfo{EDMapLocationSuffixesMap}} = ("CompositeMap" => "", "DifferenceMap" => "_diff", "ReflectionMap" => "_map");
  %{$OptionsInfo{EDMapLocationFileExtMap}} = ("CompositeMap" => "ccp4", "DifferenceMap" => "ccp4", "ReflectionMap" => "mtz");

  for $EDMapType (keys %EDLocationSuffixesMap) {
    $EDMapSuffix = $EDLocationSuffixesMap{$EDMapType};
    $OptionsInfo{EDMapLocationSuffixesMap}{$EDMapType} = $EDMapSuffix;
  }

  # Process x-ray ED map types...
  $EDMapTypes = $Options{edmaptypes};
  $EDMapTypes =~ s/ //g;
  @EDMapTypesList = ();
  if ($EDMapTypes =~ /^All$/i) {
    push @EDMapTypesList, ("CompositeMap", "DifferenceMap", "ReflectionMap");
  }
  else {
    @EDMapTypesList = split(",", $EDMapTypes);
    for $EDMapType (@EDMapTypesList) {
      if ($EDMapType !~ /^(CompositeMap|DifferenceMap|ReflectionMap|All)$/i) {
        die "Error: The value specified, $EDMapType, for option \"--EDMapTypes\" is not valid. Allowed values: CompositeMap, DifferenceMap, ReflectionMap, All\n";
      }
      if ($EDMapType =~ /^All$/i) {
        die "Error: The value specified, $EDMapType, for option \"--EDMapTypes\" must be specified alone. It can't be specified with other values.\n";
      }
    }
  }

  $OptionsInfo{EDMapTypes} = $EDMapTypes;
  @{$OptionsInfo{EDMapTypesList}} = ();
  push @{$OptionsInfo{EDMapTypesList}}, @EDMapTypesList;


  $OptionsInfo{InDelim} = $Options{indelim};

  $OptionsInfo{PDBIDsCol } = defined $Options{pdbidscol} ? $Options{pdbidscol} : '';

  $OptionsInfo{PDBFormat} = $Options{pdbformat};

  @{$OptionsInfo{CmdLinePDBIDs}} = ();
  $OptionsInfo{PDBIDsFile} = "";

  if ($OptionsInfo{Mode} =~ /^IDsOnCmdLine$/i) {
    push @{$OptionsInfo{CmdLinePDBIDs}}, @ARGV;
  }
  elsif ($OptionsInfo{Mode} =~ /^IDsInFile$/i) {
    if (@ARGV != 1) {
      die "Error: Invalid number of PDB IDs text files, ". (scalar @ARGV) . ",specified on the command line for \"IDsInFile\" value of  $Options{mode}, for option \"-m --mode\". Allowed value: Only one text file\n";
    }
    $OptionsInfo{PDBIDsFile} = $ARGV[0];

    RetrievePDBIDsTextFileInfo();
  }
  else {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: IDsOnCmdLine or IDsInFile\n";
  }
}

# Retrieve information for PDB IDs text file...
#
sub RetrievePDBIDsTextFileInfo {
  my($TextFile, $FileDir, $FileName, $FileExt, $InDelim, $Line, $ColNum, $ColLabel, $PDBIDsColIndex, $ColMode, $PDBIDsCol , $ColCount, $ColIndex, @ColLabels);

  $TextFile = $OptionsInfo{PDBIDsFile};

  %PDBIDsFileInfo = ();
  $PDBIDsFileInfo{Name} = $TextFile;

  $PDBIDsFileInfo{ColCount} = 0;
  @{$PDBIDsFileInfo{ColLabels}} = ();
  %{$PDBIDsFileInfo{ColLabelToNumMap}} = ();
  $PDBIDsFileInfo{InDelim} = "";

  $PDBIDsFileInfo{IDsColIndex} = "";

  if (!-e $TextFile) {
    die "Error: PDBIDs text file, $TextFile, doesn't exist\n";
  }

  if (!CheckFileType($TextFile, "csv tsv")) {
    die "Error: Ignoring file $TextFile: It's not a csv or tsv file\n";
  }

  ($FileDir, $FileName, $FileExt) = ParseFileName($TextFile);
  if ($FileExt =~ /^tsv$/i) {
    $InDelim = "\t";
  }
  else {
    $InDelim = "\,";
    if (!($OptionsInfo{InDelim} =~ /^(comma|semicolon)$/i)) {
      die "Error: Ignoring file $TextFile: The value specified, $OptionsInfo{InDelim}, for option \"--indelim\" is not valid for textfile\n";
    }
    if ($OptionsInfo{InDelim} =~ /^semicolon$/i) {
      $InDelim = "\;";
    }
  }

  if (!open TEXTFILE, "$TextFile") {
    die "Error: Ignoring file $TextFile: Couldn't open it: $! \n";
  }

  $Line = GetTextLine(\*TEXTFILE);
  @ColLabels = quotewords($InDelim, 0, $Line);
  close TEXTFILE;
  $ColCount = scalar @ColLabels;

  push @{$PDBIDsFileInfo{ColLabels}}, @ColLabels;
  $PDBIDsFileInfo{ColCount} = $ColCount ;
  $PDBIDsFileInfo{InDelim} = $InDelim;

  # Setup collabel to colnum map...
  %{$PDBIDsFileInfo{ColLabelToNumMap}} = ();
  for $ColNum (0 .. $#ColLabels) {
    $ColLabel = $ColLabels[$ColNum];
    $PDBIDsFileInfo{ColLabelToNumMap}{lc $ColLabel} = $ColNum;
  }

  # Identify column containing PDB IDs...
  $PDBIDsColIndex = "";

  $ColMode = $OptionsInfo{ColMode};
  $PDBIDsCol = $OptionsInfo{PDBIDsCol };

  if (IsNotEmpty($PDBIDsCol )) {
    if ($ColMode =~ /^collabel$/i) {
      $ColLabel = lc $PDBIDsCol;
      if (!exists $PDBIDsFileInfo{ColLabelToNumMap}{$ColLabel} ) {
	die "Error: Ignoring file $TextFile: The column name, $PDBIDsCol, specified for option \"-p, --PDBIDsCol \" is not valid for text file\n";
      }
      $PDBIDsColIndex = $PDBIDsFileInfo{ColLabelToNumMap}{$ColLabel};
    }
    else {
      $ColNum = $PDBIDsCol;
      if ($ColNum <= 0 || $ColNum > $ColCount) {
	die "Error: Ignoring file $TextFile: The column number, $PDBIDsCol, specified for option \"-p, --PDBIDsCol \" is not valid for text file. It must be > 0 and <= $ColCount\n";
      }
      $PDBIDsColIndex = $ColNum - 1;
    }
  }
  else {
    # Look for column name containing PDB_ID or PDBID text string...
    $PDBIDsCol = "";
    $ColIndex = 0;
    COLLABEL: for $ColLabel (@ColLabels) {
      if ($ColLabel =~ /(PDB_ID|PDBID)/i) {
	$PDBIDsCol = $ColLabel;
	$PDBIDsColIndex = $ColIndex;
	last COLLABEL;
      }
      $ColIndex++;
    }
    if (IsEmpty($PDBIDsCol)) {
      die "Error: Ignoring file $TextFile: Couldn't find PDB IDs default column containing text string PDB_ID or PDBID in its name\n";
    }
  }
  $PDBIDsFileInfo{IDsColIndex} = $PDBIDsColIndex;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();

  $Options{colmode} = "colnum";
  $Options{datalocationurl} = "http://www.rcsb.org/pdb/files/";

  $Options{densitymap} = "no";
  $Options{densitymapmode} = "auto";

  $Options{densitymaplocationurlcryoem} = "ftp://ftp.wwpdb.org/pub/emdb/structures/";
  $Options{denistymaplocationurlxray} = "http://www.ebi.ac.uk/pdbe/coordinates/files/";

  $Options{edmaplocationsuffixes} = "CompositeMap,None,DifferenceMap,_diff,ReflectionMap, _map";
  $Options{edmaptypes} = "CompositeMap,DifferenceMap";

  $Options{indelim} = "comma";
  $Options{mode} = "IDsOnCmdLine";

  $Options{pdbformat} = "Auto";

  if (!GetOptions(\%Options, "colmode|c=s", "datalocationurl|d=s",  "densitymap=s", "densitymapmode=s", "densitymaplocationurlcryoem=s", "denistymaplocationurlxray=s", "edmaplocationsuffixes=s", "edmaptypes=s",  "help|h",  "indelim=s", "mode|m=s", "pdbidscol|p=s", "pdbformat=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{colmode} !~ /^(colnum|collabel)$/i) {
    die "Error: The value specified, $Options{colmode}, for option \"-c, --colmode\" is not valid. Allowed values: colnum or collabel\n";
  }
  if ($Options{densitymap} !~ /^(yes|no)$/i) {
    die "Error: The value specified, $Options{densitymap}, for option \"--DensityMap\" is not valid. Allowed values: yes or no\n";
  }
  if ($Options{densitymapmode} !~ /^(XRayElectronDensity|CryoEMDensity|Auto)$/i) {
    die "Error: The value specified, $Options{densitymapmode}, for option \"--DensityMapMode\" is not valid. Allowed values: XRayElectronDensity, CryoEMDensity, Auto\n";
  }
  if ($Options{indelim} !~ /^(comma|semicolon)$/i) {
    die "Error: The value specified, $Options{indelim}, for option \"--indelim\" is not valid. Allowed values: comma or semicolon\n";
  }
  if ($Options{mode} !~ /^(IDsOnCmdLine|IDsInFile)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: IDsOnCmdLine or IDsInFile\n";
  }
  if ($Options{pdbformat} !~ /^(PDB|CIF|Auto)$/i) {
    die "Error: The value specified, $Options{pdbformat}, for option \"--indelim\" is not valid. Allowed values: PDB, CIF or Auto\n";
  }
}

__END__

=head1 NAME

DownloadPDBFiles.pl - Download PDB files for PDB ID(s)

=head1 SYNOPSIS

DownloadPDBFiles.pl PDBID(s) or PDBIDsTextFile

DownloadPDBFiles.pl [B<-c, --colmode> I<colnum | collabel>]
[B<-d, --dataLocationURL> I<PDB URL>] [B<--DensityMap> I<yes | no>]
[B<--DensityMapMode> I<XRayElectronDensity, CryoEMDensity, Auto>]
[B<--DensityMapLocationURLCryoEM> I<Map URL>] [B<--DenistyMapLocationURLXRay> I<Map URL>]
[B<--EDMapLocationSuffixes> I<CompositeMap,None,...>] [B<--EDMapTypes>][B<-h, --help>]
[B<--indelim> I<comma | semicolon>] [B<-m, --mode> <IDsOnCmdLine | IDsInFile>]
[B<--PDBIDsCol > I<number | string>] [B<-p, --PDBFormat> I<PDB, CIF or Auto>]
[B<-w, --WorkingDir> dirname] PDBID(s) or PDBIDsTextFile

=head1 DESCRIPTION

Download PDB files corresponding to PDB IDs specified in a column in a CSV/TSV text file or
on the command line as space delimited parameters.

It is also possible to download x-ray electron density and cryo-EM density maps for the
specified PDB IDs.

=head1 OPTIONS

=over 4

=item B<-c, --colmode> I<colnum | collabel>

Specify how columns are identified in a I<TextFile> containing PDB IDs: using column number
or column label. Possible values: I<colnum or collabel>. Default value: I<colnum>.

=item B<-d, --DataLocationURL> I<PDB URL>

Specify location of PDB URL where data files are available for download. Default value:
I<http://www.rcsb.org/pdb/files/>.

=item B<--DensityMap> I<yes | no>

Download x-ray electron density and cryo-EM density map file(s). Possible values:
I<Yes or No>. Default value: I<no>.

=item B<--DensityMapMode> I<XRayElectronDensity, CryoEMDensity, Auto>

Specify type of density map files to download. Possible values: I<XRayElectronDensity,
CryoEMDensity, Auto>. Default value: I<Auto>. The x-ray electron density and cryo-EM
density map files are downloaded for I<XRayElectronDensity> and I<CryoEMDensity>
values. The availability of both B<XRayElectronDensity> and B<XRayElectronDensity>
is checked for I<Auto> value by attempting to download x-ray map files followed by
cryo-EM map files.

X-ray Electron Density (ED ) map file(s) are downloaded in CCP4 and MTZ format. Three
different types of ED map files may be downloaded using option B<--EDMapTypes>:
CompositeMap (2Fobs - Fcalc), DifferenceMap (Fobs - Fcalc), ReflectionMap. The format of
ED data in first two file types is CCP4. The third file type contains ED data in MTZ format.

The names of the downloaded ED files are derived from input PDB IDs as shown below:

    CompositeMap (2Fobs - Fcalc):  <PDBID>.ccp4
    DifferenceMap (Fobs - Fcalc): <PDBID>_Diff.ccp4
    ReflectionMap:  <PDBID>_Map.mtz

CryoEM density map file(s) are also downloaded in CCP4 format. The names of the cyroEM
density map files is derived from EMDB ID in downloaded PDB or CIF file:

    CryoEMFile:  emd_<EMDBID>.map.gz
    Path: <CryoEMMapLocationURL>/EMD-<EMDBID>/map/emd_<EMDBID>.map.gz

=item B<--DensityMapLocationURLCryoEM> I<Map URL>

Specify location of cryoEM map URL where data files are available for download. Default
value: I<ftp://ftp.wwpdb.org/pub/emdb/structures/>.

The cryo-EM map files are also availabe at the following FTP server:

    ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/

The cryoEM database ID is automatically retrieved from the downloaded PDB or CIF file.
It is used to generate the complete path name of the cryoEM map files:

    ftp://ftp.wwpdb.org/pub/emdb/structures/EMD-<ID>/map/emd_<ID>.map.gz

In addition to map file, the following metadata file is automatically downloaded from
FTP server:

    ftp://ftp.wwpdb.org/pub/emdb/structures/EMD-<ID>/header/emd_<ID>.xml

=item B<--DenistyMapLocationURLXRay> I<Map URL>

Specify location of x-ray electron density map URL where data files are available for
download. Default value: I<http://www.ebi.ac.uk/pdbe/coordinates/files/>.

=item B<--EDMapLocationSuffixes> I<CompositeMap,None,...>

Specify file root suffixes for generating file names for x-ray electron density map files on
a remote server. It is a pariwise comma delimited list of B<EDMapTypes> and remote file
suffixes. Default value: I<CompositeMap, None, DifferenceMap, _diff, ReflectionMap, _map>.

The default names of the x-ray ED map files available on the server are shown below:

    CompositeMap (2Fobs - Fcalc): <LowercasePDBID>.ccp4
    DifferenceMap (Fobs - Fcalc): <LowercasePDBID>_diff.ccp4
    ReflectionMap: <LowercasePDBID>_map.mtz

=item B<--EDMapTypes> I<CompositeMap,DifferenceMap,ReflectionMap,All>

Specify types of x-ray Electron Density (ED) map file(s) to download. It is either a comma
delimited list of valid file types or All available file types. Possible values: I<CompositeMap,
DifferenceMap, ReflectionMap, All>. Default value: I<CompositeMap,DifferenceMap>.

The CompositeMap (2Fobs - Fcalc) and DifferenceMap (Fobs - Fcalc) correspond to ED
map data in CCP4 format. The ReflectionMap corresponds to ED map data in MTZ format.

=item B<-h, --help>

Print this help message.

=item B<--indelim> I<comma | semicolon>

Input delimiter for CSV I<TextFile> containing PDB IDs. Possible values: I<comma or semicolon>.
Default value: I<comma>. For TSV files, this option is ignored and I<tab> is used as a delimiter.

=item B<-m, --mode> <IDsOnCmdLine | IDsInFile>

Indicate how PDB IDs are specified: PDB IDs are either present as space delimited command line
parameters or in a specific column in a CSV/TSV text file. Possible values: I<IDsOnCmdLine or  IDsInFile>.
Default: I<IDsOnCmdLine>.

=item B<-p, --PDBIDsCol > I<number | string>

Column used to identify PDB ID(s) in a text file. Default value: First column containing text
string B<PDB_ID> or <PDBID>.

For I<colnum> value of B<-c, --colmode> option, input value is a column number.
Example: I<1>.

For I<collabel> value of B<-c, --colmode> option, input value is a column label.
Example: I<PDB_ID>.

This option is ignored during I<IDsOnCmdLine> value of B<m, --mode> option.

=item B<--PDBFormat> I<PDB, CIF or Auto>

Specify file format for downloading PDB files. Possible values: I<PDB, CIF, auto>. Default
value: I<Auto>. The B<PDBID>.pdb and B<PDBID>.cif files are downloaded for I<PDB> and
I<CIF> option values. The availability of PDB fies in both I<PDB> and I<CIF> format is checked
for I<Auto> option by attempting to download B<PDB>.pdb file followed by B<PDBID>.cif
file.

The I<PDB> format files are usually not available for structures determined using
cryo-EM methodology.

=item B<-w, --WorkingDir> I<dirname>

Location of working directory. Default: current directory.

=back

=head1 EXAMPLES

To retrieve a PDB file for PDB ID 2HYY and generate a local 2HYY.pdb file, type:

    % DownloadPDBFiles.pl 2HYY

To retrieve a PDB file for PDB ID 2HYY along with electron density files and generate a
local 2HYY.pdb and electron density map files 2HYY.ccp4 and 2HYY_diff.ccp4 corresponding
to composit (2Fo - Fc) and difference maps (Fo - Fc), type:

    % DownloadPDBFiles.pl --densityMap yes 2HYY

To retrieve PDB file for 5K12 in CIF format along with cryo-EM density file and generate
a local 5K12.cif and density map file emd_8194.map.gz, type:

    % DownloadPDBFiles.pl --densityMap yes --pdbFormat CIF 5K12

To retrieve PDB files for multiple PDB IDs 2HYY and 1KV2 and generate corresponding
local PDB files, type:

    % DownloadPDBFiles.pl 2HYY 1KV2

To retrieve PDB files for multiple PDB IDs 2HYY and 1KV2 and generate corresponding
local PDB files along with appropriate x-ray electron density and cryo-EM density files,
type:

    % DownloadPDBFiles.pl --densityMap yes 2HYY 5K12

To download PDB files for PDB IDs present in column name PDB_ID or PDBID in
SamplePDBIDs.csv file and generate correponding PDB files, type

    % DownloadPDBFiles.pl -m IDsInFile SamplePDBIDs.csv

To download PDB files for PDB IDs present in a specific column name in
SamplePDBIDs.csv file and generate correponding PDB files, type

    % DownloadPDBFiles.pl -m IDsInFile -c collabel -p PDB_ID SamplePDBIDs.csv

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

ExtractFromPDBFiles.pl,  InfoPDBFiles.pl, ModifyPDBFiles.pl

=head1 COPYRIGHT

Copyright (C) 2020 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
