#!/usr/bin/env perl
# script to convert pathseq scores.txt file to the MOTHUR-type .taxsummary file
 
use warnings;
use strict;

=pod

=head1 NAME

Pathseq2taxsummary is a Perl script to convert a slightly modified and concatenated PathSeq (http://software.broadinstitute.org/pathseq/) scores.txt file to a  MOTHUR style tax.summary file (https://mothur.org/wiki/summary.tax/), which can then be used to make various plots and to compute various statistics using R. Another useful feature of this script is that it enables combined analysis of multiple samples while the PathSeq scores.txt file only includes read mapping results from a single sample.

=head1 LICENSE

Copy (C) 2020-2021  The J. Craig Venter Institute (JCVI).  All rights reserved
Written by Derrick E. Fouts and Thomas "Toby" Clarke.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 CITATION

The following publications used this program to convert a slightly modified and concatenated PathSeq (http://software.broadinstitute.org/pathseq/) scores.txt file to a  MOTHUR style tax.summary file (https://mothur.org/wiki/summary.tax/), which can then be used to make various plots and to compute various statistics using R.:

Lang S, Demir M, Martin A, Jiang L, Zhang X, Duan Y, Gao B, Wisplinghoff H, Kasper P, Roderburg C, Tacke F, Steffen H-M, Goeser T, Abraldes J G, Tu X M, Loomba R, Pride D, Fouts D E, Schnabl B. Intestinal virome signature associated with non-alcoholic steatohepatitis and fibrosis. Gastroenterology. in press.

Jiang L, Lang S, Duan Y, Zhang X, Gao B, Chopyk J, Schwanemann L K, Ventura-Cots M, Bataller R, Bosques-Padilla F, Verna E C, Abraldes J G, Brown Jr R S, Vargas V, Altamirano J, Caballería J, Shawcross D L, Ho S B, Louvet A, Lucey M R, Mathurin P, Garcia-Tsao G, Kisseleva T, Brenner D A, Tu X M, Stärkel P, Pride D, Fouts D E and Schnabl B. Intestinal virome in patients with alcoholic hepatitis. Hepatology. in press.

=head1 SYNOPSIS

  USAGE:  pathseq2taxsummary.pl -s <modified PathSeq scores.txt file> [options] > output_taxsummary.txt

=head1 OPTIONS

B<-s>            : [s]cores.txt file. This is a slightly modified (i.e., by appending "sample_prefix/scores.txt" to the end of each line) of the PathSeq scores.txt file
                 : ([0]tax_id\t[1]taxonomy\t[2]type\t[3]name\t[4]kingdom\t[5]score\t[6]score_normalized\t[7]reads\t[8]unambiguous\t[9]reference_length\t[10]sample_id)

B<-u>            : [u]nambiguous. This flags the use of unambiguous read counts in the tax.summary output file [optional].
                 : default is to use the ambiguous read counts.

B<-h>            : [h]elp menu

=head1 DESCRIPTON

This program generates a MOTHUR style tax.summary file (https://mothur.org/wiki/summary.tax/) from a slightly modified and concatenated PathSeq (http://software.broadinstitute.org/pathseq/) scores.txt file, which can then be used to make various plots and to compute various statistics using R or other software packages.

=head1 INPUT

Required input consists of a modified PathSeq scores.txt file (i.e., by appending "sample_prefix/scores.txt" to the end of each line).

=head1 CONTACT

    Derick E. Fouts
    dfouts@jcvi.org

=cut

use Getopt::Std;
use Cwd;
getopts ('hs:u');

our ($opt_h,$opt_s,$opt_u); # define imported variables
my ($scores_file,$unambiguous);

if ($opt_h)  { & option_help; }
if ($opt_s)  {$scores_file = $opt_s;} else { &option_help; }
if ($opt_u)  {$unambiguous = 1;}

############## Declare variables #################
my $cwd = getcwd();
my @a = ();
my @b = ();
my $line = "";
my $key = "";
my $Child_hashr = "";
my $node_id = "";
my %Count_hash = (); # keys are tax_id and sample_id, value is counts
my %Taxa_hash = (); # key is tax_id, values are taxonomy and level
my %Sample_hash = (); # key is the sample_id, value is 1 if the key is defined
my %Child_hash = (); # key is each taxonomic level (e.g., root|Bacteria, root|Bacteria|Firmicutes, root|Bacteria|Firmicutes|Bacilli, ...) and the value of the next level added (level+1)
my $RankID_hash = (); # key is taxonomy, value is rankID
 
sub option_help  {
    print "u = <u>nambibious results [default is 0, parse all results]\n";
    print "s = <s>cores.txt output file from PathSeq ([0]tax_id\t [1]taxonomy\t [2]type\t [3]name\t [4]kingdom\t [5]score\t [6]score_normalized\t [7]reads\t [8]unambiguous\t [9]reference_length\t [10]sample_id)\n";
    exit;
}

sub populate_taxa { # code written by Thomas Clarke to populate the "Child_hash"
    my ($taxaref,$Childref) = @_;
    my $curr = $taxaref->[0];
    my $i;
    my $next;
    for ($i = 1; $i < scalar(@{$taxaref}); $i++)
    {
	$next = $curr . "|" . $taxaref->[$i];
	my $search = "<$next>"; # because the regex below is treating "|" as an "or", not as a character in a string, adding the "<>" and the \Q$search\E fixed the problem
	if (!$Childref->{$curr})
	{
	    $Childref->{$curr}->{'children'} = $search . ";";
	}
	else  {
	    
	    if ($Childref->{$curr}->{'children'} !~ /\Q$search\E/) {
		$Childref->{$curr}->{'children'} .= "$search;";
	    }
	}
	$curr = $next;
    }
}

sub get_number  { # written by Toby Clarke to generate the MOTHUR-style taxsummary file rankIDs
    my ($hash,$id,$curr_val) = @_; # key is the taxonomy, value is the rankID
    my $i = "";
    my $list = $hash->{$id}->{children};
    $list =~ s/[\<\>]//g; # remove the "<>" flanking each taxonomic level
    my @list = split(/;/,$list);
    my @sort = sort { $a cmp $b } @list;
    my $c = "";
    for ($i = 0; $i < scalar(@sort); $i++)  {
	if (exists($hash->{$sort[$i]}->{children})) { # needed to check for this or else empty values spring up.
	    $c = $i + 1;
	    &get_number($hash, $sort[$i], $curr_val . "." . ($c));
	}
	$hash->{$sort[$i]}->{'id'} = $curr_val . "." . ($c);
    }
}

sub print_taxsummary { # written by Derrick Fouts to print out the MOTHUR-stype taxsummary file
    my ($Sampleref,$Countref,$Taxaref,$Childref) = @_; # read in hash references (pointers)
    my $key = ""; # first key is NCBI taxon_id
    my $yek = ""; # second key is sample_id
    
    print "taxon_id\ttaxlevel\trankID\ttaxon\tfull_taxonomy\ttotal\t";
    foreach $key (sort {$a cmp $b} keys %{$Sampleref}) { # print out the column headers one time
	print "$key\t";
    }
    print "\n";
    $key = ""; # clear out the value stored in the $key var
    foreach $key (sort {$a <=> $b} keys %{$Countref}) { # first key is the NCBI taxon_id
	print "$key\t$Taxaref->{$key}->{'level'}\t$Childref->{$Taxaref->{$key}->{'taxonomy'}}->{'id'}\t$Taxaref->{$key}->{'clean_Taxa'}\t$Taxaref->{$key}->{'taxonomy'}\t$Taxaref->{$key}->{'total'}\t";
	foreach $yek (sort {$a cmp $b} keys %{$Sampleref}) { # and sort by our defined sample_id hash
	    if ($Countref->{$key}{$yek}->{'cnts'}) { # if we have data for this taxon for this sample, then print the counts
		print "$Countref->{$key}{$yek}->{'cnts'}\t";
	    }
	    else {
		print "0\t"; # else no matches for this sample, print zero
	    }
	}
	print "\n"; # print carriage return to get to new line for the next taxon match
    }
}

############## M A I N #################

open (SCORESFILE, "<$cwd/$scores_file");
while (<SCORESFILE>)  {
    #bacterial example (remove those levels denoted by ""):
    #root [0]|"cellular_organisms"|Bacteria [1]|"Terrabacteria_group"|Firmicutes [2]|Bacilli [3]|Lactobacillales [4]|Streptococcaceae [5]|Streptococcus [6]|Streptococcus_sp._HMSC11C05 [7]
    #phage example:
    #root [0]|kingdom [1](Viruses)|group [2] (dsDNA_viruses)|order [3](Caudovirales)|family [4](Siphoviridae)|subfamily [5](Tunaviranae)|genus [6](Kp36virus)|species [7](Klebsiella_virus_KP36)
    chomp;
    if (/^tax_id/) { next; } # skip if header line
    $line = $_;
    $line =~ s/\/scores.txt$//; # remove the text "/scores.txt" from the sample_id
    $line =~ s/cellular_organisms\|//; # remove cellular organisms from all annotations
    if ($line =~ /Bacteria/) {
	$line =~ s/\|\w+\/.*group.*\|/\|/; # remove the extra "group-containing" level between kingdom/phylum and family/genus (for bacteria).
	$line =~ s/\|\w+_group//;
	$line =~ s/\|\w+\/\w+subdivisions//; # remove stuff like "delta/epsilon_subdivisions" that do not fit into 7 levels of taxonomy
    }
    @a=split(/\s+/,$line); #split the cleaned-up line on whitespace
    @b=split(/\|/,$a[1]); #split taxonomy on "|"
    
    &populate_taxa(\@b,\%Child_hash); # pass the taxonomy array and the Child_hash
    
    if ($unambiguous) { #if the user only wants to parse "unambiguous" hits
	$Count_hash{$a[0]}{$a[10]}->{'cnts'} = $a[8];
    }
    else {
	$Count_hash{$a[0]}{$a[10]}->{'cnts'} = $a[7];
    }
    if (!defined $Taxa_hash{$a[0]}) { #if no taxonomy info stored, then record it
	$Taxa_hash{$a[0]}->{'level'} = $#b;
	$Taxa_hash{$a[0]}->{'taxonomy'} = $a[1];
	$Taxa_hash{$a[0]}->{'clean_Taxa'} = $b[$#b]; # get the hightest level taxonomy
    }
    $Taxa_hash{$a[0]}->{'total'} += $Count_hash{$a[0]}{$a[10]}->{'cnts'}; # keep a running total of all hits to this taxon_id
    if (!defined $Sample_hash{$a[10]}) { #if no sample_id has been recorded, then do so to keep track of how many unique sample_ids we have in the data. We will sort by this later to make columns in the .taxsummary file
	$Sample_hash{$a[10]} = 1;
    }
    @a = (); #clear out @a var
    @b = (); #clear out @b var
}
close (SCORESFILE);
$Child_hash{'root'}->{'id'} = "0"; #initialize root
&get_number(\%Child_hash, "root", "0"); # fetch the rankIDs using Toby's subroutine
&print_taxsummary(\%Sample_hash,\%Count_hash,\%Taxa_hash,\%Child_hash);
exit(0);
