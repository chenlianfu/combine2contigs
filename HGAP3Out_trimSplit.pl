#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 HGAP3_out.fasta > HGAP3.trim-split.fasta

USAGE
if (@ARGV==0){die $usage}

my ($seq_id, %seq, %seq_len);
while (<>) {
    chomp;
    if (/^>(\w+)/) { $seq_id = $1; }
    else { $seq{$seq_id} .= $_; $seq_len{$seq_id} += length; }
}

my (%contig, %contig_len);
foreach (sort {$seq_len{$b} <=> $seq_len{$a}} keys %seq) {
    my $seq = $seq{$_};
    $seq =~ s/^[atcg]*//;
    $seq =~ s/[atcg]*$//;
    my @seq = split /[atcg]+/, $seq;
    foreach (@seq) { $contig{$_} = 1; $contig_len{$_} = length; }
}

my $number = 0;
foreach (sort {$contig_len{$b} <=> $contig_len{$a}} keys %contig) {
    $number ++;
    print ">contig_$number\n$_\n";
}
