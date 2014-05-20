#!/usr/bin/perl -w
#
# seqConverter.pl v1.1
# Last modified March 2, 2006 15:33
# (c) Olaf R.P. Bininda-Emonds
#
# Input:
#   Sequence data in any of fasta, nexus, (classic or extended) phylip, or Se-Al formats
#
# Output:
#   Sequence data in any of fasta, nexus, (classic or extended) phylip, and/or Se-Al formats.
#
# Usage: seqConverter.pl -d<filename> -o<n|pc|pe|s> [-i<f|n|p|s>] [-r<a|i>] [-u] [-v] [-h]
#	options: -d<filename> = file containing raw sequence information; * = batch convert all specified file types in working directory
#							(suffixes must be .fasta, .nex, .phylip, or .seal as appropriate)
#            -i<f|n|p|s> = format of sequence file (fasta (f), nexus (n), phylip (p), or Se-Al (s)); default = autodetect
#            -o<n|pc|pe|s> = output results additionally in fasta (f), nexus (n), classic or extended phylip (pc or pe), and/or Se-Al (s) formats
#			 -r<a|i> = order sequences in final output alphabetically by name (a; default)or in input order from file (i);
#            -u = interactive user-input mode
#            -h = print this message and quit
#            -v = verbose output

use strict;

# Set user-input defaults and associated parameters
	# Data set variables
		my $inputType = "";	# Options are "fasta", "nexus", "phylip", and "Se-Al"
		my ($readFile, @seqFiles);
			my $dataSource;
 		my $globalGenCode = 1;

		my (@accNum, %nameLabel, %sequence, %geneticCode, %accPresent);
		my (%deletedSeq, %finalSeq);
		my $seqCount;
		my $ntax;

	# User input variables
 		my $seqOrder = "alphabetical";	# Options are "alphabetical" (default) and "input"

	# Output variables
		my $maxLength;
		my $fastaPrint = 0;
			my $fastaOut;
		my $nexusPrint = 0;
			my $nexusOut;
		my ($phylipTradPrint, $phylipExtPrint) = (0, 0);
			my $phylipOut;
		my $sealPrint = 0;
			my $sealOut;

	# Miscellaneous variables
		my $verbose = 0;
		my $debug = 0;
		my $version = "1.1";

# Read in user input
	if (not @ARGV or join(' ', @ARGV) =~ /\s-u/ or $ARGV[0] =~ /^-u/)	# Enter interactive user-input mode
		{
		print "Entering interactive user-input mode. Type \"q\" at any prompt to exit program.\n";

		# Get datafile
			until ($readFile)
				{
				print "\tEnter name of data file (* = batch convert) [$readFile]: ";
				$readFile = <stdin>;
					chomp ($readFile);
					exit(0) if ($readFile eq "q");
				unless (-e $readFile or $readFile eq "*")
					{
					print "\t\tFile '$readFile' does not exist\n";
					$readFile = "";
					}
				}
			push @seqFiles, $readFile;
		
		# Get format of datafile
			my $defaultInput = "autodetect";
				undef $inputType;
			until (defined $inputType)
				{
				print "\tEnter format of file $readFile (fasta|nexus|phylip|Se-Al) [$defaultInput]:";
				$inputType = <stdin>;
					chomp ($inputType);
					exit(0) if ($inputType =~ /^q/i);
					if (substr($inputType, 0, 1) =~ /^a/i or $inputType eq "")
						{
						$inputType = "autodetect";
						}
					elsif (substr($inputType, 0, 1) =~ /^f/i)
						{
						$inputType = "fasta";
						}
					elsif (substr($inputType, 0, 1) =~ /^n/i)
						{
						$inputType = "nexus";
						}
					elsif (substr($inputType, 0, 1) =~ /^p/i)
						{
						$inputType = "phylip";
						}
					elsif (substr($inputType, 0, 1) =~ /^s/i)
						{
						$inputType = "Se-Al";
						}
					else
						{
						print "\t\tInvalid input ($inputType)\n";
						undef $inputType;
						}
				}
			$inputType = "" if ($inputType eq "autodetect");
		
		# Get output order of sequences
			my $defaultOrder = $seqOrder;
				undef $seqOrder;
			until (defined $seqOrder)
				{
				print "\tEnter output order for sequences (alphabetical|clustal|input file) [$defaultOrder]: ";
				$seqOrder = <stdin>;
					chomp ($seqOrder);
					exit(0) if ($seqOrder =~ /^q/i);
					if (substr($seqOrder, 0, 1) =~ /^i/i)
						{
						$seqOrder = "input";
						}
					elsif (substr($seqOrder, 0, 1) =~ /^a/i or $seqOrder eq "")
						{
						$seqOrder = "alphabetical";
						}
					else
						{
						print "\t\tInvalid input ($seqOrder)\n";
						undef $seqOrder;
						}
				}

		# Get output formats
			my $defaultFasta = "y";
				undef $fastaPrint;
			until (defined $fastaPrint)
				{
				print "\tOutput results in fasta format (y|n) [$defaultFasta]: ";
				$fastaPrint = <stdin>;
					chomp ($fastaPrint);
					exit(0) if ($fastaPrint =~ /^q/i);
					if (substr($fastaPrint, 0, 1) =~ /^y/i)
						{
						$fastaPrint = 1;
						}
					elsif (substr($fastaPrint, 0, 1) =~ /^n/i or $fastaPrint eq "")
						{
						$fastaPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($fastaPrint)\n";
						undef $fastaPrint;
						}
				}

			my $defaultNexus = "y";
				undef $nexusPrint;
			until (defined $nexusPrint)
				{
				print "\tOutput results in nexus format (y|n) [$defaultNexus]: ";
				$nexusPrint = <stdin>;
					chomp ($nexusPrint);
					exit(0) if ($nexusPrint =~ /^q/i);
					if (substr($nexusPrint, 0, 1) =~ /^y/i)
						{
						$nexusPrint = 1;
						}
					elsif (substr($nexusPrint, 0, 1) =~ /^n/i or $nexusPrint eq "")
						{
						$nexusPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($nexusPrint)\n";
						undef $nexusPrint;
						}
				}

			my $defaultPhylip = "n";
				undef $phylipTradPrint;
			until (defined $phylipTradPrint or $phylipExtPrint)
				{
				print "\tOutput results in traditional phylip format (y|n) [$defaultPhylip]: ";
				$phylipTradPrint = <stdin>;
					chomp ($phylipTradPrint);
					exit(0) if ($phylipTradPrint =~ /^q/i);
					if (substr($phylipTradPrint, 0, 1) =~ /^y/i)
						{
						$phylipTradPrint = 1;
						}
					elsif (substr($phylipTradPrint, 0, 1) =~ /^n/i or $phylipTradPrint eq "")
						{
						$phylipTradPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($phylipTradPrint)\n";
						undef $phylipTradPrint;
						}
				}
				
				if ($phylipTradPrint == 0)	# Check for extended format
					{
					my $defaultPhylip = "n";
						undef $phylipExtPrint;
					until (defined $phylipExtPrint or $phylipExtPrint)
						{
						print "\tOutput results in extended phylip format (y|n) [$defaultPhylip]: ";
						$phylipExtPrint = <stdin>;
							chomp ($phylipExtPrint);
							exit(0) if ($phylipExtPrint =~ /^q/i);
							if (substr($phylipExtPrint, 0, 1) =~ /^y/i)
								{
								$phylipExtPrint = 1;
								}
							elsif (substr($phylipExtPrint, 0, 1) =~ /^n/i or $phylipExtPrint eq "")
								{
								$phylipExtPrint = 0;
								}
							else
								{
								print "\t\tInvalid input ($phylipExtPrint)\n";
								undef $phylipExtPrint;
								}
						}
					}

			my $defaultSeal = "n";
				undef $sealPrint;
			until (defined $sealPrint)
				{
				print "\tOutput results in Se-Al format (y|n) [$defaultSeal]: ";
				$sealPrint = <stdin>;
					chomp ($sealPrint);
					exit(0) if ($sealPrint =~ /^q/i);
					if (substr($sealPrint, 0, 1) =~ /^y/i)
						{
						$sealPrint = 1;
						}
					elsif (substr($sealPrint, 0, 1) =~ /^n/i or $sealPrint eq "")
						{
						$sealPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($sealPrint)\n";
						undef $sealPrint;
						}
				}
		
		# Get verbose output mode
			my $defaultVerbose = "n";
				undef $verbose;
			until (defined $verbose)
				{
				print "\tOutput verbose results to screen (y|n) [$defaultVerbose]: ";
				$verbose = <stdin>;
					chomp ($verbose);
					exit(0) if ($verbose =~ /^q/i);
					if (substr($verbose, 0, 1) =~ /^y/i)
						{
						$verbose = 1;
						print "\n";
						}
					elsif (substr($verbose, 0, 1) =~ /^n/i or $verbose eq "")
						{
						$verbose = 0;
						}
					elsif (substr($verbose, 0, 1) =~ /^x/i or $verbose eq "")
						{
						$verbose = $debug = 1;
						}
					else
						{
						print "\t\tInvalid input ($verbose)\n";
						undef $verbose;
						}
				}
		}

	elsif (join(' ', @ARGV) =~ /\s-h/ or $ARGV[0] =~ /^-h/)	# Print help screen
		{
		print "Usage: seqConverter.pl -d<filename> -o<n|pc|pe|s> [-i<f|n|p|s>] [-r<a|i>] [-u] [-v] [-h]\n";
		print "Version: $version\n";
		print "\nOptions: -d<filename> = file containing raw sequence information; * = batch convert all specified file types in working directory\n";
		print "                        (suffixes must be .fasta, .nex, .phylip, or .seal as appropriate)\n";
		print "         -i<f|n|p|s> = format of sequence file (fasta (f), nexus (n), phylip (p), or Se-Al (s))\n";
		print "         -o<n|pc|pe|s> = output results in fasta (f), nexus (n), classic or extended phylip (pc or pe), and/or Se-Al (s) formats\n";
		print "         -r<a|i> = order sequences in final output alphabetically by name (a; default) or in input order from file (i)\n";
		print "         -u = interactive user-input mode\n";
		print "         -h = print this message and quit\n";
		print "         -v = verbose output\n";
		exit(0);
		}
	else	# Process switches
		{
		for (my $i = 0; $i <= $#ARGV; $i++)
			{
			if ($ARGV[$i] =~ /^-d(.*)/)
				{
				$readFile = $1;
				unless ($readFile eq "*")
					{
					die "ERROR: Data file $readFile does not exist.\n" unless (-e $readFile);
					}
				push @seqFiles, $readFile;
				}
			elsif ($ARGV[$i] eq "-if")
				{
				$inputType = "fasta";
				}
			elsif ($ARGV[$i] eq "-in")
				{
				$inputType = "nexus";
				}
			elsif ($ARGV[$i] eq "-ip")
				{
				$inputType = "phylip";
				}
			elsif ($ARGV[$i] eq "-is")
				{
				$inputType = "Se-Al";
				}
			elsif ($ARGV[$i] eq "-on")
				{
				$nexusPrint = 1;
				}
			elsif ($ARGV[$i] eq "-opc")
				{
				$phylipTradPrint = 1;
				}
			elsif ($ARGV[$i] eq "-ope")
				{
				$phylipExtPrint = 1;
				}
			elsif ($ARGV[$i] eq "-os")
				{
				$sealPrint = 1;
				}
			elsif ($ARGV[$i] eq "-ra")
				{
				$seqOrder = "alphabetical";
				}
			elsif ($ARGV[$i] eq "-ri")
				{
				$seqOrder = "input";
				}
			elsif ($ARGV[$i] eq "-v")
				{
				$verbose = 1;
				}
			elsif ($ARGV[$i] eq "-x")
				{
				$debug = 1;
				$verbose = 1;
				}
			else
				{
				print "Don't understand argument: $ARGV[$i]\n";
				print "Usage: seqConverter.pl -d<filename> -o<n|pc|pe|s> [-i<f|n|p|s>] [-r<a|i>] [-u] [-v] [-h]\n";
				print "Version: $version\n";
				exit(1); 
				}
			}
		}

die "ERROR: Must supply name of file containing sequence data.\n" if (not @seqFiles);
die "ERROR: Must supply at least one output format.\n" unless ($fastaPrint or $nexusPrint or $phylipTradPrint or $phylipExtPrint or $sealPrint);
	
# Read in sequence data
	if ($seqFiles[0] eq "*")	# Batch convert all appropriate files in working directory
		{
		if ($inputType)
			{
			my $suffix = $inputType;
				$suffix = "nex" if ($inputType eq "nexus");
				$suffix = "seal" if ($inputType eq "Se-Al");
	
			system ("ls *.$suffix > convertList.txt");
				undef @seqFiles;	
			setLineBreak("convertList.txt");
			open (LIST, "<convertList.txt") or die "Cannot open file containing names of all sequence files, convertList.txt\n";
				while (<LIST>)
					{
					chomp;
					next unless ($_);
					push @seqFiles, $_;
					}
			close LIST;
	
			unlink ("convertList.txt") unless ($debug);
			die "ERROR: No files of file type $inputType found for batch conversion.\n" if (not @seqFiles);
			}
		else
			{
			die "ERROR: Must specify input file type for batch conversion\n";
			}
		}

	foreach my $seqFile (@seqFiles)
		{
		print "\nConverting file $seqFile ...\n";

		# Set output file names
			$dataSource = $seqFile;
				$dataSource =~ s/\.\w+$//;
		
			if ($fastaPrint)
				{
				$fastaOut = $dataSource.".fasta";
					$fastaOut =~ s/\.fasta$/_new.fasta/ if ($fastaOut eq $seqFile);
				}
			if ($nexusPrint)
				{
				$nexusOut = $dataSource.".nex";
					$nexusOut =~ s/\.nex$/_new.nex/ if ($nexusOut eq $seqFile);
				}
			if ($phylipTradPrint or $phylipExtPrint)
				{
				$phylipOut = $dataSource.".phylip";
					$phylipOut =~ s/\.phylip$/_new.phylip/ if ($phylipOut eq $seqFile);
				}
			if ($sealPrint)
				{
				$sealOut = $dataSource.".seal";
					$sealOut =~ s/\.seal$/_new.seal/ if ($sealOut eq $seqFile);
				}
			
			$phylipExtPrint = 0 if ($phylipTradPrint);
		
		# Read in sequence data
			# Clear variables
				undef @accNum;
				undef %nameLabel;
				undef %sequence;
				undef %geneticCode;
				undef %accPresent;
				undef %deletedSeq;
				undef %finalSeq;
				undef $seqCount;
				undef $ntax;
				$maxLength = 0;

			seqRead($seqFile);
			if (not @accNum)
				{
				print "\tERROR: Could not read in sequences from file $seqFile; skipping to next file\n";
				next;
				}
		
		# Process for printing
			foreach my $seq (@accNum)
				{
				$finalSeq{$seq} = $sequence{$seq};
				$maxLength = length($sequence{$seq}) if (length($sequence{$seq}) > $maxLength);
				}
	
		# Add gaps to end of any sequence less than maximum length
			foreach my $seq (@accNum)
				{
				$finalSeq{$seq} .= "-" x ($maxLength - length($sequence{$seq}));
				}
	
		# Print results!
			$ntax = scalar(@accNum);
			@accNum = sort { $nameLabel{$a} cmp $nameLabel{$b} } keys %nameLabel if ($seqOrder eq "alphabetical");
		
			print "\nPrinting results ...\n";
				seqPrint($seqFile);
		}

exit(0);

### Subroutines used in the program

sub setLineBreak	# Check line breaks of input files and set input record separator accordingly
	{
	my $inFile = shift;
	$/ ="\n";
	open (IN, "<$inFile") or die "Cannot open $inFile to check form of line breaks.\n";
		while (<IN>)
			{
			if ($_ =~ /\r\n/)
				{
				print "\tDOS line breaks detected ...\n" if ($verbose);
				$/ ="\r\n";
				last;
				}
			elsif ($_ =~ /\r/)
				{
				print "\tMac line breaks detected ...\n" if ($verbose);
				$/ ="\r";
				last;
				}
			else
				{
				print "\tUnix line breaks detected ...\n" if ($verbose);
				$/ ="\n";
				last;
				}
			}
	close IN;
	}

sub seqRead
	{
	my $seqFile = shift;

	print "\nReading in sequence data from file $seqFile (type is $inputType) ...\n" if ($inputType);
	setLineBreak($seqFile);
	open (SEQ, "<$seqFile") or die "Cannot open file containing sequences, $seqFile\n";
		my ($header, $tempAcc, $tempName, $tempSeq);
		my $fastaAcc;
		my (%nexusSpecies, %nexusAcc, $nexusRead);
		my ($phylipLineCount, $phylipTaxa, $phylipChars, %phylipSeq);
		my $sealCode;
		my ($sealDelFlag, $owner) = (0, 0);

		while (<SEQ>)
			{
			chomp;
			my $lineRead = $_;
			next unless ($lineRead);
			
			# Autodetect sequence format
				if (not $inputType)
					{
					$inputType = "fasta" if ($lineRead =~ /^>/);
					$inputType = "nexus" if ($lineRead =~ /\#nexus/i);
					$inputType = "phylip" if ($lineRead =~ /^\s*\d+\s+\d+/);
					$inputType = "Se-Al" if ($lineRead =~ /^\s*Database=\{/i);
					print "\nReading in sequence data from file $seqFile (type determined to be $inputType) ...\n" if ($inputType);
					}
			
			if ($inputType eq "nexus")
				{
				# Only read in data lines
					if ($lineRead =~ /^\s*matrix/i)
						{
						$nexusRead = 1;
						next;
						}
					$nexusRead = 0 if ($lineRead =~ /;\s*$/);
					next unless ($nexusRead);
					next unless ($lineRead =~ /a/i or $lineRead =~ /c/i or $lineRead =~ /g/i or $lineRead =~ /t/i);
				# Clean up input line
					$lineRead =~ s/^\s+//;
					$lineRead =~ s/\'//g;
				my @nexusLine = split(/\s+/, $lineRead);
					my $species = shift(@nexusLine);
						$species =~ s/\s+/_/g;
					my $seq = join('', @nexusLine);
						$seq =~ s/\s+//g;
				if (not defined $nexusSpecies{$species})
					{
					$nexusSpecies{$species} = 1;
					$seqCount++;
					$nexusAcc{$species} = "tAlign_".$seqCount;
					push @accNum, $nexusAcc{$species};
						$nameLabel{$nexusAcc{$species}} = $species;
						$sequence{$nexusAcc{$species}} = uc($seq);
						$geneticCode{$nexusAcc{$species}} = $globalGenCode;
					}
				else	# Sequences are in interleaved format; append sequence
					{
					$sequence{$nexusAcc{$species}} .= uc($seq);
					}
				}

			if ($inputType eq "fasta")
				{
				if ($lineRead =~/^\s*>/)
					{
					my $species;
					$seqCount++;
					(my $tempSpecies = $lineRead) =~ s/^\s*>//;
					
						if ($tempSpecies =~ /^Mit\.\s+/)	# Entry comes from European RNA project
							{
							$tempSpecies =~ s/^Mit\.\s+//i;	# To fix entries from European RNA project
							my @speciesInfo = split(/\s+/, $tempSpecies);
								$species = join('_', $speciesInfo[0], $speciesInfo[1]);
							if (defined $speciesInfo[2])
								{
								$fastaAcc = $speciesInfo[2];
								}
							else
								{
								$fastaAcc = "tAlign_".$seqCount;
								}
							}
						else
							{
							my @speciesLine = split(/\s+/, $tempSpecies);
							if ($speciesLine[$#speciesLine] =~ /^\(?[A-Z]+\d+\)?$/ and scalar(@speciesLine) > 1)	# Check whether last entry is an accession number
								{
								$fastaAcc = pop (@speciesLine);
								$fastaAcc =~ s/^\(//g;
								$fastaAcc =~ s/\)$//g;
								}
							else
								{
								$fastaAcc = "tAlign_".$seqCount;
								}
							$species = join('_', @speciesLine);
								$species = "Sequence_".$seqCount if ($species eq "");
							}
					push @accNum, $fastaAcc;
						$geneticCode{$fastaAcc} = $globalGenCode;
					$nameLabel{$fastaAcc} = $species;
					}
				else
					{
					$lineRead =~ s/\s+//g;
					$sequence{$fastaAcc} .= uc($lineRead);
					}
				}

			if ($inputType eq "Se-Al")
				{
				my $header;
				$sealDelFlag = 1 if ($lineRead =~/MCoL/);	# Se-Al sometimes places deleted species at end of file; do not read in remainder of file
					next if ($sealDelFlag == 1);
				next unless ($lineRead =~/NumSites/i or $lineRead =~/Owner/i or $lineRead =~/Name/i or $lineRead =~/Accession/i or $lineRead =~/Sequence/i or $lineRead =~/GeneticCode/i);
				if ($lineRead =~/Owner\s*\=\s*(\d+)/i)
					{
					$owner = $1;
					}
				if ($lineRead =~/Accession/i and $owner == 2)
					{
					$seqCount++;
					if ($lineRead =~ /null/ or $lineRead =~ /\"\"/)
						{
						$tempAcc = "tAlign_".$seqCount;
						}
					else
						{
						($header, $tempAcc) = split (/=/, $lineRead);
							$tempAcc =~ s/\"//g;
							$tempAcc =~ s/;//g;
						}
					push @accNum, $tempAcc;
					}
				if ($lineRead =~/Name/i and $owner == 2)
					{
					($header, $tempName) = split (/=/, $lineRead);
						$tempName =~ s/\"//g;
						$tempName =~ s/\s*;//g;
					}
				if ($lineRead =~/GeneticCode/i and $owner == 2)
					{
					($header, $sealCode) = split (/=/, $lineRead);
						$sealCode =~ s/\"//g;
						$sealCode =~ s/\s*;//g;
						$geneticCode{$tempAcc} = $sealCode + 1;
					}
				if ($lineRead =~/Sequence/i and $owner == 2)
					{
					($header, $tempSeq) = split (/=/, $lineRead);
						$tempSeq =~ s/\"//g;
						$tempSeq =~ s/;//g;
						$tempSeq =~ s/\s+//g;
					$nameLabel{$tempAcc} = $tempName;
					$sequence{$tempAcc} = uc($tempSeq);
					}
				}

			if ($inputType eq "phylip")
				{
				if ($lineRead =~ /^\s*(\d+)\s+(\d+)/)
					{
					$phylipTaxa = $1;
					$phylipChars = $2;
					$phylipLineCount = 0;
					}
				else
					{
					$phylipLineCount++;
					
					$lineRead =~ s/\s//g;
					
					$phylipSeq{$phylipLineCount} .= $lineRead;
					
					$phylipLineCount = 0 if ($phylipLineCount == $phylipTaxa);
					}
				}
			}
	close SEQ;
	
	if ($inputType eq "phylip")	# Postprocess input to derive taxon names and sequence; accounts for both sequential and extended formatting
		{
		for (my $i = 1; $i <= $phylipTaxa; $i++)
			{
			my $phylipAcc = "tAlign_" . $i;
			
			push @accNum, $phylipAcc;
			$geneticCode{$phylipAcc} = $globalGenCode;
			
			# Derive taxon name and sequence
				$sequence{$phylipAcc} = uc(substr($phylipSeq{$i}, 0 - $phylipChars));
				$nameLabel{$phylipAcc} = substr($phylipSeq{$i}, 0, length($phylipSeq{$i}) - $phylipChars);
					$nameLabel{$phylipAcc} =~ s/\s+//g;
			}
		}
	}
	
sub seqPrint
	{
	my $seqFile = shift;
	# Print fasta-formatted file (always)
		if ($fastaPrint)
			{
			print "\tWriting to fasta-formatted file $fastaOut ...\n";
			open (FASTA, ">$fastaOut") or die "Cannot open fasta file for aligned DNA sequences, $fastaOut";
				foreach my $entry (@accNum)
					{
					next if ($deletedSeq{$entry});
					my $fastaSeq = $finalSeq{$entry};
						my $breakPoint = 79;
						until ($breakPoint > length($fastaSeq))
							{
							my $replaceString = "\n" . substr($fastaSeq, $breakPoint, 1);
							substr($fastaSeq, $breakPoint, 1) = $replaceString;
							$breakPoint += 80;
							}
					print FASTA ">$nameLabel{$entry}";
						print FASTA "\t($entry)" unless ($entry =~ /^tAlign/);
					print FASTA "\n$fastaSeq\n";
					}
			close FASTA;
			}

	# Print nexus-formatted file (on demand)
		if ($nexusPrint)
			{
			print "\tWriting to nexus file $nexusOut ...\n";
			open (NEX, ">$nexusOut") or die "Cannot open nexus file for aligned DNA sequences, $nexusOut";
				print NEX "#nexus\n\n";
				print NEX "[File created from $seqFile using seqConverter.pl v$version on ".localtime()."]\n\n";
				print NEX "begin data;\n";
				print NEX "\tdimensions ntax = $ntax nchar = $maxLength;\n";
				print NEX "\tformat datatype = DNA gap = - missing = ?;\n\n";
				print NEX "\tmatrix\n\n";
				foreach my $entry (@accNum)
					{
					next if ($deletedSeq{$entry});
					if ($nameLabel{$entry} =~ /\W/)
						{
						print NEX "'$nameLabel{$entry}'";
						}
					else
						{
						print NEX "$nameLabel{$entry}";
						}
					print NEX "\t$finalSeq{$entry}\n";
					}
				print NEX "\t;\nend;\n";
			close NEX;
			}

	# Print phylip-formatted file (on demand)
		if ($phylipTradPrint or $phylipExtPrint)
			{
			my $maxTaxLength = 50;
				$maxTaxLength = 10 if ($phylipTradPrint);
			my %shortNameCount;	
				
			print "\tWriting to phylip file $phylipOut ...\n";
			open (PHYLIP, ">$phylipOut") or die "Cannot open phylip file for aligned DNA sequences, $phylipOut";
				print PHYLIP "\t$ntax\t$maxLength\n";
				foreach my $entry (@accNum)
					{
					next if ($deletedSeq{$entry});
					
					my $phylipName = $nameLabel{$entry};

					# Check name label and adjust to proper length if needed
						if (length($phylipName) < $maxTaxLength)
							{
							$shortNameCount{$phylipName}++;
							$phylipName .= " " x ($maxTaxLength - length($phylipName));	# Pad end of name with spaces as needed
							}
						else
							{
							my $trimmedName = substr($phylipName, 0 , $maxTaxLength);
							$shortNameCount{$trimmedName}++;
							if ($shortNameCount{$trimmedName} > 1)	# Check for duplicates among shortened names and make unique by adding numbers
								{
								$phylipName = substr($phylipName, 0, $maxTaxLength - length($shortNameCount{$trimmedName}));
									$phylipName .= $shortNameCount{$trimmedName};
									$phylipName .= " " x ($maxTaxLength - length($phylipName));	# Pad end of name with spaces as needed
								}
							else
								{
								$phylipName = $trimmedName;
								}
							}
						
					print PHYLIP "$phylipName";
						print PHYLIP " " if ($phylipExtPrint);
					print PHYLIP "$finalSeq{$entry}\n";
					}
			close PHYLIP;
			}

	# Print Se-Al-formatted file (on demand)
		if ($sealPrint)
			{
			print "\tWriting to Se_Al file $sealOut ...\n";
			open (SEAL, ">$sealOut") or die "Cannot open Se-Al file for aligned DNA sequences, $sealOut\n";
				print SEAL "Database={\n";
				print SEAL "\tID='MLst';\n";
				print SEAL "\tOwner=null;\n";
				print SEAL "\tName=null;\n";
				print SEAL "\tDescription=null;\n";
				print SEAL "\tFlags=0;\n";
				print SEAL "\tCount=2;\n";
				print SEAL "\t{\n\t\t{\n";
				
				print SEAL "\t\t\tID='PAli';\n";
				print SEAL "\t\t\tOwner=1;\n";
				print SEAL "\t\t\tName=\"$seqFile\";\n";
				print SEAL "\t\t\tDescription=null;\n";
				print SEAL "\t\t\tFlags=0;\n";
				print SEAL "\t\t\tNumSites=$maxLength;\n";
				print SEAL "\t\t\tType=\"Nucleotide\";\n";
				print SEAL "\t\t\tFeatures=null;\n";
				print SEAL "\t\t\tColourMode=1;\n";
				print SEAL "\t\t\tLabelMode=0;\n";
				print SEAL "\t\t\ttriplets=false;\n";
				print SEAL "\t\t\tinverse=true;\n";
				print SEAL "\t\t\tCount=$ntax;\n";
				print SEAL "\t\t\t{\n";
				
				my $i = 0;
				foreach my $sequence (@accNum)
					{
					next if ($deletedSeq{$sequence});
					$i++;
					print SEAL "\t\t\t\t{\n";
					print SEAL "\t\t\t\t\tID='PSeq';\n";
					print SEAL "\t\t\t\t\tOwner=2;\n";
					print SEAL "\t\t\t\t\tName=\"$nameLabel{$sequence}\";\n";
					print SEAL "\t\t\t\t\tDescription=null;\n";
					print SEAL "\t\t\t\t\tFlags=0;\n";
					print SEAL "\t\t\t\t\tAccession=";
						if ($sequence =~/^tAlign_/)
							{
							print SEAL "null;\n";
							}
						else
							{
							print SEAL "$sequence;\n";
							}
							
					print SEAL "\t\t\t\t\tType=\"DNA\";\n";
					print SEAL "\t\t\t\t\tLength=".length($finalSeq{$sequence}).";\n";
					print SEAL "\t\t\t\t\tSequence=\"$finalSeq{$sequence}\";\n";
					my $sealCode = $geneticCode{$sequence} - 1;
					print SEAL "\t\t\t\t\tGeneticCode=$sealCode;\n";
					print SEAL "\t\t\t\t\tCodeTable=null;\n";
					print SEAL "\t\t\t\t\tFrame=1;\n";
					print SEAL "\t\t\t\t\tFeatures=null;\n";
					print SEAL "\t\t\t\t\tParent=null;\n";
					print SEAL "\t\t\t\t\tComplemented=false;\n";
					print SEAL "\t\t\t\t\tReversed=false;\n";
					print SEAL "\t\t\t\t}";
					print SEAL "," unless ($i == $ntax);
					print SEAL "\n";
					}
				
				print SEAL "\t\t\t};\n";
				print SEAL "\t\t},\n";
				print SEAL "\t\t{\n";
				print SEAL "\t\t\tID='MCoL';\n";
				print SEAL "\t\t\tOwner=1;\n";
				print SEAL "\t\t\tName=\"Genetic Codes\";\n";
				print SEAL "\t\t\tDescription=\"Custom Genetic Codes\";\n";
				print SEAL "\t\t\tFlags=0;\n";
				print SEAL "\t\t\tCount=0;\n";
				print SEAL "\t\t}\n";
				print SEAL "\t};\n";
				print SEAL "};\n";
			close SEAL;

			if (-e "/Developer/Tools/SetFile")
				{
				system ("/Developer/Tools/SetFile -t 'TEXT' $sealOut");
				system ("/Developer/Tools/SetFile -c 'SEAL' $sealOut");
				}
			}
	}

# Version history
#
#	v1.1 (March 2, 2006)
#		- added ability to batch convert all specified file types in working directory (use -d*)
#		- updated to seqRead module 1.1.1 (includes autodetection of sequence format)
#		- checks that necessary input file(s) exists before proceeding
#		- added GNU GPL statement
#		- sets TYPE and CREATOR codes for Se-Al files on Mac systems when SetFile is present
#		- minor bug fixes
#
#	v1.0 (May 30, 2005)
#		- initial release
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is available at
# http://www.gnu.org/copyleft/gpl.html or by writing to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Fifth Floor, Boston, MA, 02110-1301, USA.

