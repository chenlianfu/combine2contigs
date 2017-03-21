#!/usr/bin/perl
use strict;
use Getopt::Long;
use FileHandle;

my $usage = <<USAGE;
Usage:
    perl $0 master.contigs.fasta slave.contigs.fasta > out.fasta

    --nucmer <string>    default: None
    设置nucmer命令的参数。

    --show_coords <string>    default: "-I 98 -L 200"
    设置show-coords命令的参数。仅用于设置show-coords命令的 -I 和 -L 参素。-I 参数表示给出Identity大于此百分比的比对结果；-L 参数表示给出比对长度高于此值的比对结果。-L 参数要低于master和slave contigs序列的最短序列长度。

    --cpu <int>    default: 2
    设置并行化运行nucmer的并发数。程序将输入的master.contigs.fasta文件中的序列按从长到短进行排列，并逐条按顺序写入8个fasta文件中，然后将这8个fasta文件与slave.contigs.fasta使用nucmer进行序列比对；同样对slave.contigs.fasta进行分割后，再使用nucmer比对到master.contigs.fasta。

    --overlap_distance <int>    default: 50
    通过nucmer的结果对master和slave序列的首尾部分进行overlap分析时，允许在末端有该长度序列的不匹配。

    --min_overlap_length <int>    default: 8000
    若overlap长度低于此值，则需要对overlap部分的序列使用blast进行分析该overlap是否为重复序列。若overlap序列不是重复序列，则可以用于进行contig连接。

    --min_overlap_percentage <float>    default: 10
    若overlap对序列的覆盖率低于此百分数，则需要对overlap部分的序列使用blast进行分析该overlap是否为重复序列。若overlap序列不是重复序列，则可以用于进行contig连接。

    --blast_identity_percentage <float>    default: 95
    设置Blast结果Identity阈值。Identity不满足此阈值的blast结果是无效的比对结果。

    --blast_coverage_ratio <float>    default: 0.95
    设置Blast结果覆盖率阈值。比对的HSP结果的覆盖率不满足此阈值的blast结果是无效的比对结果。

USAGE
if (@ARGV==0){die $usage}

my ($prefix, $nucmer, $show_coords, $cpu, $overlap_distance, $min_overlap_length, $min_overlap_percentage, $blast_identity_percentage, $blast_coverage_ratio);
GetOptions(
    "prefix:s" => \$prefix,
    "nucmer:s" => \$nucmer,
    "show_coords:s" => \$show_coords,
    "cpu:i" => \$cpu,
    "overlap_distance:i" => \$overlap_distance,
    "min_overlap_length:i" => \$min_overlap_length,
    "min_overlap_percentage:f" => \$min_overlap_percentage,
    "blast_identity_percentage:f" => \$blast_identity_percentage,
    "blast_coverage_ratio:f" => \$blast_coverage_ratio,
);
$prefix ||= "out";
$show_coords ||= "-I 98 -L 200";
$cpu ||= 1;
$overlap_distance ||= 50;
$min_overlap_length ||= 8000;
$min_overlap_percentage ||= 10;
$blast_identity_percentage ||= 95;
$blast_coverage_ratio ||= 0.95;

mkdir "$prefix.tmp" unless -e "$prefix.tmp";
open CMD, ">", "command.nucmer.list" or die $!;
my (@nucmer_out_files_master, @nucmer_out_files_slave);

open IN, $ARGV[0] or die $!;
my ($seq_id, %seq, %seq_len, %master_id);
while (<IN>) {
    chomp;
    if (/^>(\S+)/) { $seq_id = $1; $master_id{$seq_id} = 1; }
    else { $seq{$seq_id} .= $_; $seq_len{$seq_id} += length; }
}
close IN;
my %masterID = %master_id;

my %fh;
foreach (0 .. $cpu - 1) {
    open $fh{$_}, ">", "$prefix.tmp/master.$_.fasta" or die $!;
    print CMD "nucmer -p $prefix.tmp/master.$_ $nucmer $ARGV[1] $prefix.tmp/master.$_.fasta; show-coords -c -d -l -q -T $show_coords $prefix.tmp/master.$_.delta > $prefix.tmp/master.$_.show\n";
    push @nucmer_out_files_master, "$prefix.tmp/master.$_.show";
}
my $number = 0;
foreach (sort {$seq_len{$b} <=> $seq_len{$a}} keys %master_id) {
    $fh{$number}->print(">$_\n$seq{$_}\n");
    $number ++;
    $number = 0 if $number == $cpu;
}

open IN, $ARGV[1] or die $!;
my ($seq_id, %slave_id);
while (<IN>) {
    chomp;
    if (/^>(\S+)/) { $seq_id = $1; $slave_id{$seq_id} = 1; }
    else { $seq{$seq_id} .= $_; $seq_len{$seq_id} += length; }
}
close IN;

my %fh;
foreach (0 .. $cpu - 1) {
    open $fh{$_}, ">", "$prefix.tmp/slave.$_.fasta" or die $!;
    print CMD "nucmer -p $prefix.tmp/slave.$_ $nucmer $ARGV[0] $prefix.tmp/slave.$_.fasta; show-coords -c -d -l -q -T $show_coords $prefix.tmp/slave.$_.delta > $prefix.tmp/slave.$_.show\n";
    push @nucmer_out_files_slave, "$prefix.tmp/slave.$_.show";
}
my $number = 0;
foreach (sort {$seq_len{$b} <=> $seq_len{$a}} keys %slave_id) {
    $fh{$number}->print(">$_\n$seq{$_}\n");
    $number ++;
    $number = 0 if $number == $cpu;
}

close CMD;
my $cmdString = "ParaFly -c command.nucmer.list -CPU $cpu &> /dev/null";
(system $cmdString) == 0 or die "Failed to excute: $cmdString\n";

my %aligment;
foreach (@nucmer_out_files_master) {
    open IN, $_ or die $!;
    while (<IN>) {
        next unless m/^\d/;
        @_ = split /\s+/;
        my $aligment;
        if ($_[12] eq "1") {
            $aligment = "$_[2]\t$_[3]\t$_[0]\t$_[1]\t$_[5]\t$_[4]\t$_[6]\t$_[8]\t$_[7]\t$_[10]\t$_[9]\t$_[11]\t$_[12]\t$_[14]\t$_[13]";
        }
        else {
            $aligment = "$_[3]\t$_[2]\t$_[1]\t$_[0]\t$_[5]\t$_[4]\t$_[6]\t$_[8]\t$_[7]\t$_[10]\t$_[9]\t$_[11]\t$_[12]\t$_[14]\t$_[13]";
        }
        $aligment{$aligment} ++;
    }
    close IN;
}
foreach (@nucmer_out_files_slave) {
    open IN, $_ or die $!;
    while (<IN>) {
        next unless m/^\d/;
        chomp;
        $aligment{$_} ++;
    }
    close IN;
}

$cmdString = "cat $ARGV[0] $ARGV[1] > for_blastn_db.fasta";
(system $cmdString) == 0 or die "Failed to excute: $cmdString!\n";
$cmdString = "makeblastdb -in for_blastn_db.fasta -dbtype nucl -out for_blastn_db -logfile for_blastn_db.log";
(system $cmdString) == 0 or die "Failed to excute: $cmdString!\n";
open OUT, ">", "relationship.overlap_need_validation.fasta" or die $!;
open OUT1, ">", "relationship.overlap_noNedd_validataion.txt" or die $!;
open OUT2, ">", "relationship.overlap_need_validation.txt" or die $!;
open OUT3, ">", "relationship.cross_info.txt" or die $!;
open OUT4, ">", "relationship.master_contain_slave.txt" or die $!;
open OUT5, ">", "relationship.slave_contain_master.txt" or die $!;
open OUT6, ">", "relationship.nucmer_out.txt" or die $!;

my (%overlap, %overlap_validation, %slave_contain_master, %slave_cover);
foreach (keys %aligment) {
    next if $aligment{$_} != 2;
    print OUT6 "$_\n";
    @_ = split /\t/;
    my $validation = 0;
    if ($_[4] <= $min_overlap_length or $_[5] <= $min_overlap_length or $_[9] <= $min_overlap_percentage or $_[10] <= $min_overlap_percentage) {
        $validation = 1;
    }

    if ($_[10] == 100) {
        print OUT4 "$_[13]\_$_[7]\t$_[0]\t$_[1]\t$_[14]\_$_[8]\t$_[2]\t$_[3]\n";
        delete $slave_id{$_[14]};
    }
    elsif ($_[9] == 100) {
        print OUT5 "$_[13]\_$_[7]\t$_[0]\t$_[1]\t$_[14]\_$_[8]\t$_[2]\t$_[3]\n";
        my $min_value = $_[2];
        $min_value = $_[3] if $_[3] < $_[2];
        $slave_contain_master{$_[14]}{$min_value} = "$_[13]\t$_[12]\t$_[0]\t$_[1]\t$_[2]\t$_[3]";
        delete $master_id{$_[13]};
    }
    elsif ($_[0] <= $overlap_distance) {
        if ($_[3] >= $_[8] - $overlap_distance) {
            if ($validation == 0) {
                print OUT1 "$_[13]\tstart\t$_[14]\tend\t$_[0]\t$_[1]\t$_[2]\t$_[3]\n";
                $overlap{$_[13]}{"start"}{$_[14]} = "end\t$_[0]\t$_[1]\t$_[2]\t$_[3]";
                $overlap{$_[14]}{"end"}{$_[13]} = "start\t$_[2]\t$_[3]\t$_[0]\t$_[1]";
            }
            else {
                print OUT2 "$_[13]\tstart\t$_[14]\tend\t$_[0]\t$_[1]\t$_[2]\t$_[3]\n";
                my $seq1 = substr($seq{$_[13]}, $_[0] - 1, $_[1] - $_[0] + 1);
                my $seq2 = substr($seq{$_[14]}, $_[2] - 1, $_[3] - $_[2] + 1);
                print OUT ">$_[13]_$_[0]_$_[1]\n$seq1\n>$_[14]_$_[2]_$_[3]\n$seq2\n";
                $overlap_validation{"$_[13]\tstart\t$_[14]\tend\t$_[0]\t$_[1]\t$_[2]\t$_[3]"} = "$_[13]_$_[0]_$_[1]\t$_[14]_$_[2]_$_[3]";
            }
        }
        elsif ($_[3] <= $overlap_distance) {
            if ($validation == 0) {
                print OUT1 "$_[13]\tstart\t$_[14]\tstart\t$_[0]\t$_[1]\t$_[3]\t$_[2]\n";
                $overlap{$_[13]}{"start"}{$_[14]} = "start\t$_[0]\t$_[1]\t$_[3]\t$_[2]";
                $overlap{$_[14]}{"start"}{$_[13]} = "start\t$_[3]\t$_[2]\t$_[0]\t$_[1]";
            }
            else {
                print OUT2 "$_[13]\tstart\t$_[14]\tstart\t$_[0]\t$_[1]\t$_[3]\t$_[2]\n";
                my $seq1 = substr($seq{$_[13]}, $_[0] - 1, $_[1] - $_[0] + 1);
                my $seq2 = substr($seq{$_[14]}, $_[3] - 1, $_[2] - $_[3] + 1);
                print OUT ">$_[13]_$_[0]_$_[1]\n$seq1\n>$_[14]_$_[3]_$_[2]\n$seq2\n";
                $overlap_validation{"$_[13]\tstart\t$_[14]\tstart\t$_[0]\t$_[1]\t$_[3]\t$_[2]"} = "$_[13]_$_[0]_$_[1]\t$_[14]_$_[3]_$_[2]";
            }
        }
        else {
            print OUT3 "$_\n";
            $slave_cover{$_[14]}{"$_[2]\t$_[3]"} = 1;
        }
    }
    elsif ($_[1] >= $_[7] - $overlap_distance) {
        if ($_[2] <= $overlap_distance) {
            if ($validation == 0) {
                print OUT1 "$_[13]\tend\t$_[14]\tstart\t$_[0]\t$_[1]\t$_[2]\t$_[3]\n";
                $overlap{$_[13]}{"end"}{$_[14]} = "start\t$_[0]\t$_[1]\t$_[2]\t$_[3]";
                $overlap{$_[14]}{"start"}{$_[13]} = "end\t$_[2]\t$_[3]\t$_[0]\t$_[1]";
            }
            else {
                print OUT2 "$_[13]\tend\t$_[14]\tstart\t$_[0]\t$_[1]\t$_[2]\t$_[3]\n";
                my $seq1 = substr($seq{$_[13]}, $_[0] - 1, $_[1] - $_[0] + 1);
                my $seq2 = substr($seq{$_[14]}, $_[2] - 1, $_[3] - $_[2] + 1);
                print OUT ">$_[13]_$_[0]_$_[1]\n$seq1\n>$_[14]_$_[2]_$_[3]\n$seq2\n";
                $overlap_validation{"$_[13]\tend\t$_[14]\tstart\t$_[0]\t$_[1]\t$_[2]\t$_[3]"} = "$_[13]_$_[0]_$_[1]\t$_[14]_$_[2]_$_[3]";
            }
        }
        elsif ($_[2] >= $_[8] - $overlap_distance) {
            if ($validation == 0) {
                print OUT1 "$_[13]\tend\t$_[14]\tend\t$_[0]\t$_[1]\t$_[3]\t$_[2]\n";
                $overlap{$_[13]}{"end"}{$_[14]} = "end\t$_[0]\t$_[1]\t$_[3]\t$_[2]";
                $overlap{$_[14]}{"end"}{$_[13]} = "end\t$_[3]\t$_[2]\t$_[0]\t$_[1]";
            }
            else {
                print OUT2 "$_[13]\tend\t$_[14]\tend\t$_[0]\t$_[1]\t$_[3]\t$_[2]\n";
                my $seq1 = substr($seq{$_[13]}, $_[0] - 1, $_[1] - $_[0] + 1);
                my $seq2 = substr($seq{$_[14]}, $_[3] - 1, $_[2] - $_[3] + 1);
                print OUT ">$_[13]_$_[0]_$_[1]\n$seq1\n>$_[14]_$_[3]_$_[2]\n$seq2\n";
                $overlap_validation{"$_[13]\tend\t$_[14]\tend\t$_[0]\t$_[1]\t$_[3]\t$_[2]"} = "$_[13]_$_[0]_$_[1]\t$_[14]_$_[3]_$_[2]";
            }
        }
        else {
            print OUT3 "$_\n";
            $slave_cover{$_[14]}{"$_[2]\t$_[3]"} = 1;
        }
    }
    else {
        print OUT3 "$_\n";
        $slave_cover{$_[14]}{"$_[2]\t$_[3]"} = 1;
    }
}
close OUT; close OUT1; close OUT2; close OUT3; close OUT4; close OUT5; close OUT6;

$cmdString = "blastn -query relationship.overlap_need_validation.fasta -db for_blastn_db -num_threads $cpu -outfmt 7 -out relationship.overlap_need_validation.blastout -max_target_seqs 10";
(system $cmdString) == 0 or die "Failed to excute: $cmdString!\n";

open IN, "relationship.overlap_need_validation.blastout" or die $!;
my %blast_hsp_number;
while (<IN>) {
    next if m/^#/;
    @_ = split /\t/;
    if ($_[2] >= $blast_identity_percentage) {
        my $query_len = $2 - $1 if $_[0] =~ m/(\d+)_(\d+)$/;
        my $coverage = abs($_[7] - $_[6]) / $query_len;
        $blast_hsp_number{$_[0]} ++ if $coverage >= $blast_coverage_ratio;
    }
}
close IN;

open OUT1, ">", "relationship.overlap_pass_validation.txt" or die $!;
open OUT2, ">", "relationship.overlap_fail_validation.txt" or die $!;
foreach (keys %overlap_validation) {
    my ($seq1, $seq2) = split /\t/, $overlap_validation{$_};
    #print STDERR "$blast_hsp_number{$seq1}\t$blast_hsp_number{$seq2}\n";
    if ($blast_hsp_number{$seq1} == 2 && $blast_hsp_number{$seq2} == 2) {
        print OUT1 "$_\n";
        @_ = split /\t/;
        $overlap{$_[0]}{$_[1]}{$_[2]} = "$_[3]\t$_[4]\t$_[5]\t$_[6]\t$_[7]";
        $overlap{$_[2]}{$_[3]}{$_[0]} = "$_[1]\t$_[6]\t$_[7]\t$_[4]\t$_[5]";
    }
    else {
        print OUT2 "$_\n";
    }
}
close OUT1; close OUT2;

my %overlap_f = %overlap ;
open OUT, ">", "relationship.overlap.effective.txt" or die $!;
foreach my $seqID1 (keys %overlap_f) {
    foreach my $SD (keys %{$overlap_f{$seqID1}}) {
        my @seqID2 = keys %{$overlap_f{$seqID1}{$SD}};
        if (@seqID2 >= 2) {
            foreach my $seqID2 (keys %{$overlap_f{$seqID1}{$SD}}) {
                #print STDERR "$seqID1\t$SD\t$seqID2\t$overlap{$seqID1}{$SD}{$seqID2}\n";
                @_ = split /\t/, $overlap_f{$seqID1}{$SD}{$seqID2};
                delete $overlap{$seqID2}{$_[0]};
            }
            delete $overlap{$seqID1}{$SD};
        }
        else {
            foreach my $seqID2 (keys %{$overlap_f{$seqID1}{$SD}}) {
                print OUT "$seqID1\t$SD\t$seqID2\t$overlap{$seqID1}{$SD}{$seqID2}\n";
            }
        }
    }
}
close OUT;

my %slave_contain_master_f = %slave_contain_master;
foreach my $slave_id (keys %slave_contain_master_f) {
    my (@region_overlap, @region_contain);
    foreach (keys %{$overlap{$slave_id}}) {
        my @seqID = keys %{$overlap{$slave_id}{$_}};
        @_ = split /\t/, $overlap{$slave_id}{$_}{$seqID[0]};
        push @region_overlap, "$_[1]\t$_[2]";
    }
    foreach (keys %{$slave_contain_master_f{$slave_id}}) {
        @_ = split /\t/, $slave_contain_master_f{$slave_id}{$_};
        push @region_contain, "$_[4]\t$_[5]";
    }
    my @delete = &cal_delete_region(\@region_contain, \@region_overlap);
    foreach (@delete) {
        delete $slave_contain_master{$slave_id}{$_};
    }
    my @number = keys %{$slave_contain_master{$slave_id}};
    delete $slave_contain_master{$slave_id} if @number < 1;
}

my %path;
while (%master_id) {
    #my @master_id = sort {$seq_len{$b} <=> $seq_len{$a}} keys %master_id;
    my @master_id = keys %master_id;
    my $master_id = $master_id[0];
    delete $master_id{$master_id};

    my ($path_right, $path_left);
    my @sides = ("start", "end");
    foreach my $side (@sides) {
        my $seq_id = $master_id;
        my $right_or_left = $side;
        my ($path_out, $end);
        if ($side eq "end") { $end = $seq_len{$seq_id}; } else { $end = 1; }

        while (exists $overlap{$seq_id}{$side}) {
            my @next_seq_id = keys %{$overlap{$seq_id}{$side}};
            my $next_seq_id = $next_seq_id[0];
            delete $master_id{$next_seq_id};
            @_ = split /\t/, $overlap{$seq_id}{$side}{$next_seq_id};
            my $another_side = "start";
            $another_side = "end" if $_[0] eq "start";

            my ($start1, $end1, $start2, $end2, $next_start, $last_end);
            if ($side eq "end") {
                if ($_[0] eq "start") {
                    ($start1, $end1, $start2, $end2) = ($_[1], $_[2], $_[3], $_[4]);
                    $next_start = $_[4] + 1;
                    $end = $seq_len{$next_seq_id};
                }
                else {
                    ($start1, $end1, $start2, $end2) = ($_[1], $_[2], $_[4], $_[3]);
                    $next_start = $_[3] - 1;
                    $end = 1;
                }
                $last_end = $_[1] - 1;
            }
            else {
                if ($_[0] eq "start") {
                    ($start1, $end1, $start2, $end2) = ($_[2], $_[1], $_[3], $_[4]);
                    $next_start = $_[4] + 1;
                    $end = $seq_len{$next_seq_id};
                }
                else {
                    ($start1, $end1, $start2, $end2) = ($_[2], $_[1], $_[4], $_[3]);
                    $next_start = $_[3] - 1;
                    $end = 1;
                }
                $last_end = $_[2] + 1;
            }
            $path_out .= "$seq_id:$last_end;$seq_id:$start1-$seq_id:$end1,$next_seq_id:$start2-$next_seq_id:$end2;$next_seq_id:$next_start-";

            if (exists $slave_id{$next_seq_id} && exists $slave_contain_master{$next_seq_id}) {
                if ($_[0] eq "start") {
                    foreach (sort {$a <=> $b} keys %{$slave_contain_master{$next_seq_id}}) {
                        my @align = split /\t/, $slave_contain_master{$next_seq_id}{$_};
                        if ($align[1] == 1) {
                            ($start1, $end1, $start2, $end2) = ($align[4], $align[5], $align[2], $align[3]);
                            $last_end = $align[4] - 1;
                            $next_start = $align[5] + 1;
                        }
                        else {
                            ($start1, $end1, $start2, $end2) = ($align[5], $align[4], $align[3], $align[2]);
                            $last_end = $align[5] - 1;
                            $next_start = $align[4] + 1;
                        }
                        $path_out .= "$next_seq_id:$last_end;$next_seq_id:$start1-$next_seq_id:$end1,$align[0]:$start2-$align[0]:$end2;$next_seq_id:$next_start-";
                    }
                }
                else {
                    foreach (sort {$b <=> $a} keys %{$slave_contain_master{$next_seq_id}}) {
                        my @align = split /\t/, $slave_contain_master{$next_seq_id}{$_};
                        if ($align[1] == 1) {
                            ($start1, $end1, $start2, $end2) = ($align[5], $align[4], $align[3], $align[2]);
                            $last_end = $align[5] + 1;
                            $next_start = $align[4] - 1;
                        }
                        else {
                            ($start1, $end1, $start2, $end2) = ($align[4], $align[5], $align[2], $align[3]);
                            $last_end = $align[4] + 1;
                            $next_start = $align[5] - 1;
                        }
                        $path_out .= "$next_seq_id:$last_end;$next_seq_id:$start1-$next_seq_id:$end1,$align[0]:$start2-$align[0]:$end2;$next_seq_id:$next_start-";
                    }
                }
            }

            delete $slave_id{$next_seq_id};
            $side = $another_side;
            $seq_id = $next_seq_id;
        }
        $path_out .= "$seq_id:$end";

        if ($right_or_left eq "end") {
            $path_right = $path_out;
        }
        else {
            $path_left = $path_out;
        }
    }

    #print STDERR "$master_id\t$path_left\n$master_id\t$path_right\n";
    my @path_left = reverse split /;/, $path_left;
    my @path_left_reverse;
    foreach (@path_left) {
        s/([^,;]+?):(\d+)-([^,;]+?):(\d+)/$3:$4-$1:$2/g;
        push @path_left_reverse, $_;
    }
    $path_left = join ";", @path_left_reverse;
    $path{"$path_left-$path_right"} = 1;
    #print STDERR "$master_id\t$path_left-$path_right\n";
}

open OUT1, ">", "remaining_slave_mapping_ratio.txt";
open OUT2, ">", "not_used_slave_seq.txt";
foreach my $slave_id (keys %slave_id) {
    my $path_out;
    my @cover = keys %{$slave_cover{$slave_id}};
    my @contain_region;
    my @site = sort {$a <=> $b} keys %{$slave_contain_master{$slave_id}};
    my $site_number = @site;
    foreach my $site (@site) {
        @_ = split /\t/, $slave_contain_master{$slave_id}{$site};
        push @contain_region, "$_[4]\t$_[5]";
    }
    my $coverage_length_exclude_contain = &cal_regionA_length_excludeB(\@cover, \@site);
    my @empty;
    my $coverage_length_contain = &cal_regionA_length_excludeB(\@site, \@empty);
    my $coverage_rate_repeat = $coverage_length_exclude_contain / $seq_len{$slave_id};
    my $coverage_rate_contain = $coverage_length_contain / $seq_len{$slave_id};
    my $rate = $coverage_rate_repeat + $coverage_rate_contain;
    print OUT1 "$slave_id\t$coverage_rate_repeat\t$coverage_rate_contain\t$rate\n";

    my $keep = 0;
    if ($coverage_rate_repeat <= 0) { $keep = 1; }
    elsif ($site_number >= 1) { $keep = 1; }

    if ($keep == 1) {
        if (exists $slave_contain_master{$slave_id}) {
            my ($slave_start, $slave_end);
            my $extend_start = &cal_extend("1\t$site[0]", \@cover);
            if ($extend_start == 1) {
                $path_out = "$slave_id:1-";
            }
            foreach (@site) {
                @_ = split /\t/, $slave_contain_master{$slave_id}{$_};
                my ($min, $max) = ($_[4], $_[5]);
                if ($_[4] > $_[5]) { $min = $_[5], $max = $_[4]; }
                $slave_end = $min - 1;
                $slave_start = $max + 1;
                if ($_[1] == 1) {
                    $path_out .= "$slave_id:$slave_end;$slave_id:$_[4]-$slave_id:$_[5],$_[0]:$_[2]-$_[0]:$_[3];$slave_id:$slave_start-";
                }
                else {
                    $path_out .= "$slave_id:$slave_end;$slave_id:$_[5]-$slave_id:$_[4],$_[0]:$_[3]-$_[0]:$_[2];$slave_id:$slave_start-";
                }
            }
            my $extend_end = &cal_extend("$site[-1]\t$seq_len{$slave_id}", \@cover);
            if ($extend_end == 1) {
                $path_out .= "$slave_id:$seq_len{$slave_id}";
            }
            if ($extend_start != 1) { $path_out =~ s/^.*?;//; }
            if ($extend_end != 1) { $path_out =~ s/(.*);.*/$1/; }
        }
        else {
            $path_out = "$slave_id:1-$slave_id:$seq_len{$slave_id}";
        }
        #print STDERR "$slave_id\t$path_out\n";
    }
    else {
        print OUT2 "$slave_id\t$coverage_rate_repeat\t$coverage_rate_contain\t$rate\n";
    }
    $path{$path_out} = 1;
}
close OUT1; close OUT2;

my (%path_seq, %path_len, %path_modify);
open OUT1, ">", "path.origin.txt";
open OUT2, ">", "path.modify.txt";
foreach my $path (keys %path) {
    next unless $path;
    print OUT1 "$path\n";
    my $path_m = &path_modify($path);
    $path_modify{$path_m} = 1;
    print OUT2 "$path_m\n";

    my @region = split /;/, $path_m;
    my $seq;
    foreach (@region) {
        if (m/(.*?):(\d+)-(\d+)/ && exists $seq{$1}) {
            if ($3 >= $2) {
                $seq .= substr($seq{$1}, $2 - 1, $3 - $2 + 1);
            }
            else {
                my $seq0 = substr($seq{$1}, $3 - 1, $2 - $3 + 1);
                $seq0 = reverse $seq0;
                $seq0 =~ tr/ATCGatcg/TAGCtagc/;
                $seq .= $seq0;
            }
        }
        else {
            $path_seq{"$_\t$path_m"} = $seq;
            $path_len{"$_\t$path_m"} = length $seq;
            $seq = "";
            print STDERR "PATH format is not right: $_\t$path_m\n";
        }
    }
    $path_seq{$path_m} = $seq;
    $path_len{$path_m} = length $seq;
}

my $combine_number = 1;
foreach my $path (sort {$path_len{$b} <=> $path_len{$a}} keys %path_seq) {
    my @region = split /;/, $path;
    if (@region > 1) {
        print ">combine_$combine_number\t$path\n$path_seq{$path}\n";
        $combine_number ++;
    }
    else {
        $region[0] =~ m/(.*?):\d+-\d+/;
        print ">$1\n$path_seq{$path}\n";
    }
}

%path = %path_modify;
foreach my $path (keys %path) {
    next unless $path;
    my @region = split /;/, $path;
    my %last_region;
    foreach (@region) {
        my %region;
        foreach (split /,/, $_) {
            $region{$1} = "$2\t$3" if m/(.*):(\d+)-.*?:(\d+)/;
        }
        foreach my $seqID (keys %region) {
            if (exists $last_region{$seqID}) {
                my ($start1, $end1) = split /\t/, $last_region{$seqID};
                my ($start2, $end2) = split /\t/, $region{$seqID};
                if ($end1 > $start1) {
                    print STDERR "$seqID\t$start1-$end1\t$start2-$end2\n" if $start2 != $end1 + 1;
                    print STDERR "$seqID\t$start1-$end1\t$start2-$end2\n" if $end2 < $start2;
                }
                elsif ($end1 < $start1) {
                    print STDERR "$seqID\t$start1-$end1\t$start2-$end2\n" if $start2 != $end1 - 1;
                    print STDERR "$seqID\t$start1-$end1\t$start2-$end2\n" if $end2 > $start2;
                }
            }
        }
        %last_region = %region;
    }
}

sub cal_extend {
    my $aa = $_[0];
    my ($start, $end) = split /\t/, $aa;
    my @bb = @{$_[1]};
    my %bb;
    foreach (@bb) {
        @_ = split /\t/;
        my $min = $_[0]; my $max = $_[1];
        if ($_[1] < $_[0]) { $min = $_[1]; $max = $_[0]; }
        $bb{"$min\t$max"} = $min;
    }
    my $return = 1;
    foreach (sort {$bb{$a} <=> $bb{$b}} keys %bb) {
        @_ = split /\t/;
        $return = 0 if ($start < $_[1] && $end > $_[0]);
    }
    return $return;
}

sub cal_regionA_length_excludeB {
    my @aa = @{$_[0]};
    my @bb = @{$_[1]};
    my (%aa, %bb);
    foreach (@aa) {
        @_ = split /\t/;
        my $min = $_[0]; my $max = $_[1];
        if ($_[1] < $_[0]) { $min = $_[1]; $max = $_[0]; }
        $aa{"$min\t$max"} = $min;
    }
    foreach (@bb) {
        @_ = split /\t/;
        my $min = $_[0]; my $max = $_[1];
        if ($_[1] < $_[0]) { $min = $_[1]; $max = $_[0]; }
        $bb{"$min\t$max"} = 1;
    }

    my $length = 0;
    my $last_end = 0;
    #my @region = sort {$aa{$a} <=> $aa{$b}} keys %aa;
    foreach my $region (sort {$aa{$a} <=> $aa{$b}} keys %aa) {
        my ($start, $end) = split /\t/, $region;
        if ($end <= $last_end) {
            next;
        }
        elsif ($start <= $last_end) {
            $start = $last_end + 1;
        }
        $last_end = $end;
        foreach (keys %bb) {
            @_ = split /\t/;
            delete $bb{$_} if $_[1] < $start;
            if ($start < $_[1] && $end > $_[0]) {
                if ($start > $_[0]) {
                    $start = $_[1] + 1;
                }
                else {
                    $end = $_[0] - 1;
                }
            }
        }
        if ($start < $end) {
            $length += $end - $start + 1;
        }
        else {
            next;
        }
    }

    return $length;
}

sub cal_delete_region {
    my @aa = @{$_[0]};
    my @bb = @{$_[1]};
    my (%aa, %bb);
    foreach (@aa) {
        @_ = split /\t/;
        my $min = $_[0]; my $max = $_[1];
        if ($_[1] < $_[0]) { $min = $_[1]; $max = $_[0]; }
        $aa{"$min\t$max"} = $min;
    }
    foreach (@bb) {
        @_ = split /\t/;
        my $min = $_[0]; my $max = $_[1];
        if ($_[1] < $_[0]) { $min = $_[1]; $max = $_[0]; }
        $bb{"$min\t$max"} = 1;
    }

    my %return;
    my ($last_start, $last_end, $last_length) = (0, 0, 0);
    my $last_region;
    foreach my $region (sort {$aa{$a} <=> $aa{$b}} keys %aa) {
        my ($start, $end) = split /\t/, $region;
        my $keep = 1;
        foreach (keys %bb) {
            @_ = split /\t/;
            if ($start <= $_[1] && $end >= $_[0]) {
                $return{$aa{$region}} = 1;
                $keep = 0;
                ($last_start, $last_end, $last_length) = (0, 0, 0);
            }
        }
        next if $keep == 0;

        my $length = $end - $start + 1;
        if ($start <= $last_end && $end >= $last_start) {
            if ($length < $last_length) {
                $return{$aa{$region}} = 1;
            }
            else {
                $return{$aa{$last_region}} = 1;
                $last_start = $start;
                $last_end = $end;
                $last_region = $region;
            }
        }
        else {
            $last_start = $start;
            $last_end = $end;
            $last_region = $region;
        }
    }
    my @return = keys %return;
    return @return;
}

sub path_modify {
    my $path = shift @_;
    my (%cut, %remove);
    my @region = split /;/, $path;
    my %last_region;
    foreach (@region) {
        my %region;
        foreach (split /,/, $_) {
            $region{$1} = "$2\t$3" if m/(.*):(\d+)-.*?:(\d+)/;
        }
        my $cut = 0;
        foreach my $seqID (keys %region) {
            if (exists $last_region{$seqID}) {
                my ($start1, $end1) = split /\t/, $last_region{$seqID};
                my ($start2, $end2) = split /\t/, $region{$seqID};
                if ($end1 > $start1) {
                    warn "$seqID\t$start1-$end1\t$start2-$end2\n" if $start2 != $end1 + 1;
                    if ($end2 < $start2) { $cut = 1; }
                }
                elsif ($end1 < $start1) {
                    warn "$seqID\t$start1-$end1\t$start2-$end2\n" if $start2 != $end1 - 1;
                    if ($end2 > $start2) { $cut =  1; }
                }
                if ($cut == 1) {
                    $cut{"$seqID:$start1-$seqID:$end1"} = abs($end2 - $end1);
                    $remove{"$seqID:$start2-$seqID:$end2"} = 1;
                }
            }
        }
        if ($cut == 1) {
            %last_region = ();
        }
        else {
            %last_region = %region;
        }
    }
    
    my ($path_out, @path_out);
    foreach (@region) {
        my @sub_region = split /,/, $_;
        my $cut_length = 0;
        foreach (@sub_region) {
            $cut_length = $cut{$_} if exists $cut{$_};
        }
        if (@sub_region == 2) {
            foreach (@sub_region) {
                m/(.*):(\d+)-.*?:(\d+)/;
                if (exists $masterID{$1}) {
                    my $end = $3 - $cut_length;
                    $end = $3 + $cut_length if $2 > $3;
                    push @path_out, "$1:$2-$end";
                }
            }
        }
        else {
            if (!exists $remove{$sub_region[0]}) {
                $sub_region[0] =~ m/(.*):(\d+)-.*?:(\d+)/;
                push @path_out, "$1:$2-$3";
            }
        }
    }
    $path_out = join ";", @path_out;
    return $path_out;
}
