#!/usr/bin/perl -w
#compareclades.pl
#Brian O'Meara, 24 August 2011
#http://www.brianomeara.info
#released under GPL2
#version 1.1
use diagnostics;
use strict;
my $assignmentfile="a";
my $observedfile="o";
my $simulatedfile="s";
my $hasSimulated=0;
my $countsim=0;
my $allowUnresolvedObservedTrees=0;
my $doSNPS=0;
#the how many trees array is indexed at 0, so to find out how many rooted, bifurcating, labeled trees there are for a polytomy with four edges you just do $howmanytrees[4]
my @howmanytrees = (1, 1, 1, 3, 15, 105, 945, 10395, 135135, 2027025, 34459425, 654729075, 13749310575, 316234143225, 7905853580625, 213458046676875, 6190283353629375, 191898783962510624, 6332659870762850304, 2.216430954767e+20, 8.20079453263789e+21, 3.19830986772878e+23, 1.3113070457688e+25, 5.63862029680584e+26, 2.53737913356263e+28, 1.19256819277443e+30, 5.84358414459473e+31, 2.98022791374331e+33, 1.57952079428395e+35, 8.68736436856175e+36, 4.9517976900802e+38, 2.92156063714732e+40, 1.78215198865986e+42, 1.12275575285571e+44, 7.29791239356214e+45, 4.88960130368663e+47, 3.37382489954378e+49, 2.39541567867608e+51, 1.74865344543354e+53, 1.31149008407515e+55, 1.00984736473787e+57, 7.97779418142917e+58, 6.46201328695763e+60, 5.36347102817483e+62, 4.5589503739486e+64, 3.96628682533529e+66, 3.5299952745484e+68, 3.21229569983905e+70, 2.98743500085031e+72, 2.8380632508078e+74, 2.75292135328357e+76, 2.72539213975073e+78, 2.75264606114824e+80, 2.83522544298268e+82, 2.97698671513182e+84, 3.18537578519105e+86, 3.47205960585824e+88, 3.85398616250265e+90, 4.35500436362799e+92, 5.00825501817219e+94, 5.85965837126146e+96, 6.97299346180114e+98, 8.43732208877938e+100, 1.03779061691986e+103, 1.29723827114983e+105, 1.64749260436028e+107, 2.12526545962477e+109, 2.78409775210844e+111, 3.70285001030423e+113, 4.99884751391071e+115, 6.84842109405767e+117, 9.51930532074016e+119, 1.34222205022436e+122, 1.91937753182084e+124, 2.78309742114022e+126, 4.09115320907612e+128, 6.09581828152342e+130, 9.20468560510036e+132, 1.40831689758035e+135, 2.18289119124955e+137, 3.42713917026179e+139, 5.44915128071625e+141, 8.77313356195316e+143, 1.43002077059837e+146, 2.3595342714873e+148, 3.9404222333838e+150, 6.65931357441862e+152, 1.13874262122558e+155, 1.97002473472026e+157, 3.44754328576045e+159, 6.102151615796e+161, 1.09228513922748e+164, 1.97703610200175e+166, 3.6179760666632e+168, 6.69325572332692e+170, 1.25163882026213e+173, 2.36559737029543e+175, 4.51829097726427e+177, 8.72030158612005e+179, 1.70045880929341e+182, 3.34990385430802e+184);
if ($#ARGV<1) {
	print "Usage: perl compareclades.pl -aAssignment.txt -oObserved.txt [-sSimulated.txt] [-c] [-u] [-d]\n\nThis script counts how many times each tree in the observed file occurs in the simulated file, taking into account the many-to-one mapping of samples to species.\n\nAssignments are given in a tab-delimited file with a column of species letters (A, B,....) and a column of corresponding sample numbers (1, 2, ...). The observed and simulated trees are lists of newick trees with sample numbers, not names. The optional -c argument just gives the count of each tree in the simulated trees file (after making samples from the same species equivalent: ((1,2),3) and ((3,1),2) are the same if samples 1 and 3 are from species A and sample 2 is from species B: the clades in each tree are thus A-B and A-A-B\nPassing it -u as an argument looks for matches for unresolved gene trees.\n\nImportant: if you omit -s you can pipe a ms input (what would normally go to sim.tre) directly into here";
}
else {
	for (my $i = 0; $i <= $#ARGV; $i++) {
		if ($ARGV[$i] =~ /-a(\S+)/) {
			$assignmentfile=$1;
		}
		elsif ($ARGV[$i] =~ /-o(\S+)/) {
			$observedfile=$1;
		}
		elsif ($ARGV[$i] =~ /-s(S+)/) {
			$simulatedfile=$1;
			$hasSimulated=1;
		}
		elsif ($ARGV[$i] =~ /-c/) {
			$countsim=1;
		}
		elsif ($ARGV[$i] =~ /-u/) {
			$allowUnresolvedObservedTrees=1;
		}
		elsif ($ARGV[$i] =~ /-d/) {
			$doSNPS=1;
		}

	}
	
	open(OBS,"$observedfile");
	open(ASSIGN,"$assignmentfile");
	my %assignments=();
	foreach my $inline (<ASSIGN>) {
		chomp $inline;
		my @args=split(/\s+/,$inline);
		$assignments{ $args[1] } = $args[0];
	}
	close ASSIGN;
	
	my %simtreescount;

	if ($hasSimulated==1) {
		open(SIM,"$simulatedfile");	
		foreach my $inline (<SIM>) {
			chomp $inline;
			my $clades=returnclades($inline, %assignments);
			$simtreescount{$clades}++;
		}
		close SIM;
	}
	else {
		while(<STDIN>) {
			my $inline=$_;
			chomp $inline;
			my $clades=returnclades($inline, %assignments);
			$simtreescount{$clades}++;
		}
	}
	
	if ($countsim==1) {
		print "Count\tTree\n";
		foreach my $clade (sort { $simtreescount{$b} <=> $simtreescount{$a} } keys %simtreescount) {
			print "$simtreescount{$clade}\t$clade\n";
		}
	}
	else {
		my $obstreecount=0;
		foreach my $inline (<OBS>) {
			chomp $inline;
			$obstreecount++;
			#$print "\n\nNow working on observed tree $obstreecount\n";
			my $clades=returnclades($inline, %assignments);
			#print "\t$clades\n";
			my $matchcount=0;
			if (exists($simtreescount{$clades})) {
				$matchcount=$simtreescount{$clades};
			}
			elsif ($allowUnresolvedObservedTrees==1 && $doSNPs==0) { #so, we know that the observed tree did not match exactly any of the simulated trees, which means it's either partially unresolved or just doesn't match. If we allow unresolved matches, let's do that
				foreach my $simtree (sort { $simtreescount{$b} <=> $simtreescount{$a} } keys %simtreescount) {
					my $potentialMatchcount=$simtreescount{$simtree};
					my @splitSimClades=split(/\|/,$simtree); #idea here is that if observed tree has some polytomies, we still want each of its clades to be paired with EXACTLY ONE of the clades in the simulated tree, leaving only those clades that are resolved in the simulated tree but not the observed tree. One concern is that you don't want two of the clades from the observed tree to match the same clade in the sim tree if that clade is only present once on the sim tree (this can happen due to the redundancy of taxon names once the many-to-one assignment to species is done). This code takes care of that problem
					my @splitObsClades=split(/\|/,$clades);
					my $allmatch=1;
					#print "\n\nsplitObsClades = @splitObsClades\nsplitSimClades = @splitSimClades";
					foreach my $obsClade (@splitObsClades) {
						my $thisCladeMatches=0;
						#print "size of splitSimClades = ".scalar(@splitSimClades)."\n";
						if ($allmatch==1) {
							for (my $simCladeIndex = 0; $simCladeIndex<scalar(@splitSimClades); $simCladeIndex++) {
								#print "\n$splitSimClades[$simCladeIndex]\t@splitSimClades\t";
								if ($obsClade eq $splitSimClades[$simCladeIndex]) {
									$thisCladeMatches=1;
									#print "\nmatch: $obsClade eq $splitSimClades[$simCladeIndex]";
									splice(@splitSimClades,$simCladeIndex,1);
								}
							}
							if ($thisCladeMatches!=1) {
								$allmatch=0;
							}
						}
						#print "\nsplitSimClades = @splitSimClades";

					}
					if ($allmatch==1) {
						#each missing edge means that one of three possible resolutions was true but all were possible
						my @polytomyarray = returnpolytomydegree($inline, %assignments);
						my $correction = 1;
						foreach my $polytomydegree (@polytomyarray) {
							$correction *= 1/$howmanytrees[$polytomydegree];
						}
						$matchcount+=$potentialMatchcount*$correction;
					}
					#print "$simtreescount{$clade}\t$clade\n";
				}
			}
			elsif ($doSNPs==1) { #all we care about is whether any clades on observed tree match (we assume there is one edge on the observed tree, which induces a bipartition of taxa, at least one of which is a clade
				foreach my $simtree (sort { $simtreescount{$b} <=> $simtreescount{$a} } keys %simtreescount) {
					my $potentialMatchcount=$simtreescount{$simtree};
					my @splitSimClades=split(/\|/,$simtree); 
					my @splitObsClades=split(/\|/,$clades);
					my $anymatch=0;
					foreach my $obsClade (@splitObsClades) {
						for (my $simCladeIndex = 0; $simCladeIndex<scalar(@splitSimClades); $simCladeIndex++) {
							#print "\n$splitSimClades[$simCladeIndex]\t@splitSimClades\t";
							if ($obsClade eq $splitSimClades[$simCladeIndex]) {
								$anymatch=1;
							}
						}
					}
					$matchcount += $anymatch;
				}				
			}
			print "$matchcount\n";
		}
	}
	
	sub returnclades {
		my ($intree, %localassignments) = @_;
		my @cladearray=();
		$intree=~s/:\d+\.?\d*e?\-?\d*//ig; #remove branch length (regex from Olaf Bininda-Emonds' partitionmetric.pl script)
		while ($intree=~m/(\([\w\,]+\))/) {
			my $matchstring=$1;
			my $innermatch="";
			if ($matchstring=~m/\(([\w\,]+)\)/) {
				$innermatch=$1;
			}
			my @splitinner=split(/\,/,$innermatch);
			my @renamedinner=();
			foreach my $taxon (@splitinner) {
				push(@renamedinner,$assignments{ $taxon });
			}
			@renamedinner = sort(@renamedinner);
			push(@cladearray,join('-',@renamedinner));
			$intree=~s/\($innermatch\)/$innermatch/;
		}
		@cladearray=sort(@cladearray);
		my $clades=join('|',@cladearray);
		return $clades;
	}
	
	sub returnpolytomydegree {
		my ($intree, %localassignments) = @_;
		my @polytomyarray=();
		$intree=~s/:\d+\.?\d*e?\-?\d*//ig; #remove branch length (regex from Olaf Bininda-Emonds' partitionmetric.pl script)
		while ($intree=~m/(\([\w\,]+\))/) {
			my $matchstring=$1;
			my $innermatch="";
			if ($matchstring=~m/\(([\w\,]+)\)/) {
				$innermatch=$1;
			}
			my @splitinner=split(/\,/,$innermatch);
			my @renamedinner=();
			if(scalar @splitinner > 2) {
				push(@polytomyarray, scalar @splitinner)
			}
			foreach my $taxon (@splitinner) {
				push(@renamedinner,$assignments{ $taxon });
			}
			@renamedinner = sort(@renamedinner);
			$intree=~s/\($innermatch\)/$innermatch/;
		}
		return @polytomyarray;
	}

	close OBS;
}
