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
if ($#ARGV<1) {
	print "Usage: perl compareclades.pl -aAssignment.txt -oObserved.txt [-sSimulated.txt] [-c] [-u]\n\nThis script counts how many times each tree in the observed file occurs in the simulated file, taking into account the many-to-one mapping of samples to species.\n\nAssignments are given in a tab-delimited file with a column of species letters (A, B,....) and a column of corresponding sample numbers (1, 2, ...). The observed and simulated trees are lists of newick trees with sample numbers, not names. The optional -c argument just gives the count of each tree in the simulated trees file (after making samples from the same species equivalent: ((1,2),3) and ((3,1),2) are the same if samples 1 and 3 are from species A and sample 2 is from species B: the clades in each tree are thus A-B and A-A-B\nPassing it -u as an argument looks for matches for unresolved gene trees.\n\nImportant: if you omit -s you can pipe a ms input (what would normally go to sim.tre) directly into here";
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
			elsif ($allowUnresolvedObservedTrees==1) { #so, we know that the observed tree did not match exactly any of the simulated trees, which means it's either partially unresolved or just doesn't match. If we allow unresolved matches, let's do that
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
						my $numberCollapsedEdges=scalar(@splitSimClades);
						#$print "\nnumberCollapsedEdges = $numberCollapsedEdges\ninitial match count = $matchcount\nnew count = $matchcount + $potentialMatchcount * 1.0 / (3** $numberCollapsedEdges) = ";
						$matchcount+=$potentialMatchcount*1.0/(3**$numberCollapsedEdges); #=number of sim trees / 3^number of unresolved edges. Since more than one simtree topology can work, have to do += rather than =
						#print "$matchcount\n";
					}
					#print "$simtreescount{$clade}\t$clade\n";
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
	close OBS;
}
