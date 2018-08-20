#!/usr/bin/perl -w
######coords-nucmer-parser-circos.pl######################
###### It parses .coords files and creates .txt files that can be used as input in circos######

use strict;
use warnings;
my $file=shift; #Takes ref_qry.coords file output from nucmer run
my $closeRel=shift; #Takes the name of the query organism (for file naming purposes)
my $outfile="$closeRel-nucmer-circos";
open (FILE, $file); ##from here: http://perl.about.com/od/filesystem/a/perl_parse_tabs.htm
open OUTPUT, ">", "$outfile.txt" or die "cannot open > $outfile: $!";
open OUTPUT2, ">", "$outfile-highlight.txt" or die "cannot open > $outfile: $!";
open OUTPUT3, ">", "$outfile-links.txt" or die "cannot open > $outfile: $!";
<<<<<<< HEAD
while (<FILE>) { 
#chomp;
next if $. < 6;# https://stackoverflow.com/a/14407167
#print $_;
my @F = split ' ';
#print "@F\n";
print "This $F[9] is 9 and percent ID \n";
print "This $F[18] is 18 and query ID \n";
print "This $F[17] is 17 and the reference ID \n";
#my ($ref,$query,$perc_id,$contigid) = @F[1,2,3,8];  # FROM: chrome-extension://gbkeegbaiigmenfmjfclcdgdpimamgkj/views/app.html
my $query = $F[19];
my $perc_id = $F[10];
=======
while (<FILE>) {
next unless /contig/;
my @F = split(/\t/);
#my ($ref,$query,$perc_id,$contigid) = @F[1,2,3,8];
my $query = $F[1];
my $perc_id = $F[6];
>>>>>>> 3fa82f1a65ae115d445f0343786dc4de916f65ac
$perc_id=~ s/\s+//g;
my $contigid = $F[18];
$contigid=~ s/\s+//g;
my $refID = $F[17];
$refID =~ s/\s+//g;

my $qstart = $F[3];
my $qend = $F[4];
$qend= $qend-1; #because the tracks start at 0 and not at 1 as alignments, alignment will always start at 1

my $rstart = $F[0];
my $rend = $F[1];
$rend= $rend-1;	#because the tracks start at 0 and not at 1 as alignments, alignment will always start at 1

<<<<<<< HEAD
#print "sm$contigid $qstart $qend $perc_id\n";
print OUTPUT "sm$contigid $qstart $qend $perc_id\n";
print OUTPUT2 "sm$contigid $qstart $qend\n";
print OUTPUT3 "sm$contigid $qstart $qend $refID $rstart $rend\n";
=======
print "$contigid $qstart $qend $perc_id\n";
print OUTPUT "$contigid $qstart $qend $perc_id\n";
print OUTPUT2 "$contigid $qstart $qend\n";
print OUTPUT3 "$contigid $qstart $qend $refID $rstart $rend\n";
>>>>>>> 3fa82f1a65ae115d445f0343786dc4de916f65ac
#print OUTPUT "sm$contigid $start $end 100\n";
# print "query:\t$query\tstart:\t$start\t0:\t$F[0]\t1:\t$F[1]\t2:\t$F[2]\t3:\t$F[3]\t4:\t$F[4]\t5:\t$F[5]\t6:\t$F[6]\t7:\t$F[7]\t8:\t$F[8]\n";
}
