#!/usr/bin/perl
# $Revision: 0.6 $
# $Date: 2016/05/27 $
# $Id: class_sequences.pl $
# $Author: Michaël Bekaert $
#
#
# This file is part of lamp-virus. It generates RT-LAMP primer sets for each group
# of genomes provided.
#
# Copyright (C) 2014-2016, Michaël Bekaert <michael.bekaert@stir.ac.uk>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw/uniq/;
use Bio::SeqIO;
use File::Temp qw/tempfile/;

#----------------------------------------------------------
our ($VERSION) = 0.6;
my $pathlava = '/opt/lava-dna';

#----------------------------------------------------------
my ($verbose, $alt, $combi, $realign, $strict, $loose, $extra, $fasta, $list, $prefix) = (0, 0, 0, 0, 0, 0);
GetOptions(
           'p|prefix=s'          => \$prefix,
           'a|align|alignment=s' => \$fasta,
           'g|groups=s'          => \$list,
           'r|realign!'          => \$realign,
           'c|combinatory!'      => \$combi,
           'e|extra:s'           => \$extra,
           's|strict+'           => \$strict,
           'l|loose+'            => \$loose,
           'alt!'                => \$alt,
           'v|verbose!'          => \$verbose
          );
if (defined $prefix && defined $list && -r $list && defined $fasta && -r $fasta)
{
    my (%pop, %sample);
    if (open my $infile, q{<}, $list)
    {
        while (<$infile>)
        {
            chomp;
            my @data = split m/\t/;
            if (scalar @data >= (2 + $alt) && defined $data[0] && length($data[0]) > 2 && defined $data[($alt ? 2 : 1)] && int($data[($alt ? 2 : 1)]) > 0)
            {
                $data[0] = $1 if ($data[0] =~ m/(^\w+\d+\|\d{4})/);
                $sample{$data[0]} = $data[($alt ? 2 : 1)];
                if   (exists $pop{$data[($alt ? 2 : 1)]}) { $pop{$data[($alt ? 2 : 1)]}++; }
                else                                      { $pop{$data[($alt ? 2 : 1)]} = 1; }
            }
        }
        close $infile;
    }
    my $seq_in = Bio::SeqIO->new(-file => q{<} . $fasta, -format => 'fasta');
    my %seq_out;
    foreach my $i (keys %pop) { $seq_out{$i} = Bio::SeqIO->new(-file => q{>} . $prefix . q{.} . $i . '.fasta', -format => 'fasta'); }
    while (my $inseq = $seq_in->next_seq)
    {
        my $id = $inseq->id();
        $id = $1 if ($id =~ m/(^\w+\d+\|\d{4})/);
        if   (exists $sample{$id}) { $seq_out{$sample{$id}}->write_seq($inseq); }
        else                       { print {*STDERR} '! Unknown sequence ', $inseq->id(), "\n"; }
    }
    foreach my $i (keys %pop) { $seq_out{$i}->close(); }
    $seq_in->close();
    if (scalar keys %pop > 0 && $combi)
    {
        #test combinaisons..
        my @bad;
        my @list = sort keys %pop;
        my $size = scalar @list;
        for (my $i = 0 ; $i < 2**$size ; $i++)
        {
            my $str = sprintf("%*.*b", $size, $size, $i);
            my %combination;
            for (my $j = 0 ; $j < $size ; $j++)
            {
                if (substr($str, $j, 1)) { $combination{$list[$j]} = $list[$j]; }
            }
            if (scalar keys %combination > 0 && scalar keys %combination <= $size)
            {
                my $flag = 0;
                foreach my $item (@bad)
                {
                    my $flag2 = scalar @{$item};
                    foreach my $key (@{$item})
                    {
                        if (exists $combination{$key}) { $flag2--; }
                    }
                    $flag++ if (!$flag2);
                }
                next if ($flag);
                if ($verbose)
                {
                    foreach my $item (@list) { print {*STDERR} q{ } . $item . q{:} . (exists $combination{$item} ? q{1} : q{0}); }
                    print {*STDERR} "\n";
                }
                my ($fd, $tmp_fasta) = tempfile(UNLINK => 1, OPEN => 0, SUFFIX => '.fasta');
                my $seq_in  = Bio::SeqIO->new(-file => q{<} . $fasta,     -format => 'fasta');
                my $seq_out = Bio::SeqIO->new(-file => q{>} . $tmp_fasta, -format => 'fasta');
                while (my $inseq = $seq_in->next_seq)
                {
                    my $id = $inseq->id();
                    $id = $1 if ($id =~ m/(^\w+\d+\|\d{4})/);
                    if (exists $sample{$id} && exists $combination{$sample{$id}}) { $seq_out->write_seq($inseq); }
                }
                $seq_out->close();
                $seq_in->close();
                if ($realign)
                {
                    my ($fd3, $tmp_align) = tempfile(UNLINK => 1, OPEN => 0, SUFFIX => '.fasta');
                    print {*STDERR} '>Realign ', join(q{ }, sort keys %combination), "\n" if ($verbose);
                    print {*STDERR} 'GramAlign -q -F 1 -f 2 -i ', $tmp_fasta, ' -o ', $tmp_align, "\n" if ($verbose);
                    system 'GramAlign -q -F 1 -f 2 -i ' . $tmp_fasta . ' -o ' . $tmp_align;
                    unlink '_ga_temp.page0' if (-e '_ga_temp.page0');
                    if (-e $tmp_align)
                    {
                        unlink $tmp_fasta;
                        $tmp_fasta = $tmp_align;
                    }
                }
                my $para_file = $pathlava . '/t_data/normal_parameters.xml';
                if (defined $extra && -r $extra) { $para_file = $extra; }
                elsif ($strict == 1 ) { $para_file = $pathlava . '/t_data/strict_parameters.xml'; }
                elsif ($strict == 2)  { $para_file = $pathlava . '/t_data/stricter_parameters.xml'; }
                elsif ($loose == 1)   { $para_file = $pathlava . '/t_data/loose_parameters.xml'; }
                elsif ($loose == 2)   { $para_file = $pathlava . '/t_data/looser_parameters.xml'; }
                my ($fd2, $tmp_output) = tempfile(UNLINK => 1, OPEN => 0);
                print {*STDERR} $pathlava, '/lava.pl --alignment_fasta ', $tmp_fasta, ' --output_file ', $tmp_output, ' --option_file ', $para_file, "\n" if ($verbose);
                system $pathlava . '/lava.pl --alignment_fasta ' . $tmp_fasta . ' --output_file ' . $tmp_output . ' --option_file ' . $para_file . ' >/dev/null 2>/dev/null';
                print {*STDOUT} '> Group ', join(q{ }, sort keys %combination), "\t";

                if (-r $tmp_output)
                {
                    if (open my $infile, q{<}, $tmp_output)
                    {
                        my $lines = 0;
                        $lines++ while (<$infile>);
                        close $infile;
                        print {*STDOUT} int($lines / 12);
                    }
                    print {*STDOUT} "\tPrimer(s) generated!\n";
                    unlink $tmp_output;
                    system 'cp ' . $tmp_output . '.dash ' . $prefix . q{.} . join(q{_}, sort keys %combination) . '.primers';
                    unlink $tmp_output . '.dash';
                }
                else
                {
                    my @tmp;
                    foreach my $item (@list) { push @tmp, $item if (exists $combination{$item}); }
                    push @bad, [@tmp];
                    print {*STDOUT} "0\tNo primer generated! Will skip " . join(q{ }, @tmp) . " in the next steps\n";
                }
                unlink $tmp_fasta;
            }
        }
    }
    exit 0;
}
else
{
    print {*STDOUT}
      "\nUsage: class_sequences.pl --prefix <output prefix> --align <aligned genomes> --groups <list of groups> [..]\nDescription: Generate RT-LAMP primer sets for each group of genomes provided.\n\n--prefix <file prefix>\n    Filename output prefix. [mandatory]\n--groups <path>\n    Path to a TSV file with the genome groups (e.g. output.group.tsv). [mandatory]\n        sequence_id<tab>grouping<tab>alt_grouping\n        sequence_1       1           1\n        sequence_2       1           2\n        sequence_3       2           3\n\n--align <aligned genome>\n    Alignment generated by collect_genomevirus.pl. [mandatory]\n--realign\n    Force each group or group combination to be realign prior running LAVA.\n--combinatory\n    Force to test group combinations (recommended).\n--extra <lava parameter file>\n    Provide LAVA XML parameter file.\n--strict\n    Force usage of strict LAMP parameters (see documentation).\n--strict --strict\n    Force usage of even stricter LAMP parameters (see documentation).\n--loose\n    Force usage of loose LAMP parameters (see documentation).\n--loose --loose\n    Force usage of even looser LAMP parameters (see documentation).\n--alt\n    Force usage alternative grouping column from the group file.\n--verbose\n    Becomes very chatty.\n\n";
    exit 1;
}