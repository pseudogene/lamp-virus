#!/usr/bin/perl
# $Revision: 0.5 $
# $Date: 2016/05/27 $
# $Id: collect_genomevirus.pl $
# $Author: Michaël Bekaert $
#
#
# This file is part of lamp-virus. It collects all complete genomes using the provided
# ENTREZ query and align them.
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
use Bio::SeqIO;
use Bio::DB::Query::GenBank;
use Bio::DB::GenBank;

#----------------------------------------------------------
our ($VERSION) = 0.5;
my $pathtemplate = '/opt';

#----------------------------------------------------------
my $query_string = 'txid64320[Organism:exp]';
my ($verbose, $after, $before, $run, $ok, $skip, $missing, $genome) = (0, 0, 1900 + (localtime(time))[5], 1, 0, 0);
GetOptions('p|prefix=s' => \$genome, 'q|query=s' => \$query_string, 'a|after:i' => \$after, 'b|before:i' => \$before, 'r|run!' => \$run, 'v|verbose!' => \$verbose);

if (defined $query_string && defined $genome && !-e ($genome . '.fasta') && defined $after && defined $before && open(my $seq_out, q{>}, $genome . '.fasta'))
{
    my $query  = Bio::DB::Query::GenBank->new(-db => 'nucleotide', -query => $query_string . ' AND (complete[All Fields] AND genome[All Fields])', -mindate => ($after > 0 ? "$after" : "0"), -maxdate => "$before");
    my $gb     = Bio::DB::GenBank->new();
    my $stream = $gb->get_Stream_by_query($query);
    while (my $seq = $stream->next_seq)
    {
        print {*STDERR} '>', $seq->id(), "\n" if ($verbose);
        my ($date, $country);
        for my $feat_object ($seq->get_SeqFeatures)
        {
            if ($feat_object->primary_tag eq 'source')
            {
                for my $tag ($feat_object->get_all_tags)
                {
                    if ($tag eq 'country')
                    {
                        for my $value ($feat_object->get_tag_values($tag)) { $country = $value; }
                    }
                    elsif ($tag eq 'collection_date')
                    {
                        for my $value ($feat_object->get_tag_values($tag)) { $date = $1 if ($value =~ m/(\d{4})/); }
                    }
                }
            }
        }
        if (defined $date && defined $country)
        {
            if ($date > $after && $date <= $before)
            {
                $ok++;
                print {$seq_out} q{>} . $seq->accession_number() . q{|} . $date . q{|} . $country . "\n" . $seq->seq() . "\n";
            }
            else { $skip++; }
        }
        else { print {*STDERR} '  !No metadata for accession: ', $seq->accession_number(), "\n"; $missing++; }
    }
    close $seq_out;
    print {*STDERR} "\nNumber of sequence retrieved         ", ($ok + $skip + $missing), ($after > 0 ? "\nNumber of sequence too old too young " . $skip : q{}), "\nNumber of sequence without data      ", $missing,
      "\nNumber of remaining sequences        ", $ok, "\n";
    if ($ok > 0 && -r ($genome . '.fasta'))
    {
        system 'GramAlign -q -F 1 -f 2 -i ' . $genome . '.fasta -o ' . $genome . '.align.fa';
        unlink '_ga_temp.page0' if (-e '_ga_temp.page0');
        if (-e ($genome . '.align.fa') && $verbose)
        {
            my $seqio = Bio::SeqIO->new('-format' => 'fasta', -file => $genome . '.align.fa');
            my $seqobj = $seqio->next_seq();
            print {*STDERR} 'Alignment length                 ', $seqobj->length, "\n";
        }
        if ($run && -e ($genome . '.align.fa') && -r ($pathtemplate . '/virus.R') && (open my $input, q{<}, $pathtemplate . '/virus.R') && (open my $output, q{>}, $genome . '.R'))
        {
            while (<$input>) { print {$output} $_; }
            close $input;
            print {$output} "\nprocess_virus('$genome','$genome.align.fa','$genome.results'", ($verbose ? ',TRUE' : q{}), ");\n";
            close $output;
            system 'Rscript --vanilla ' . $genome . '.R';
        }
        exit 0;
    }
    else { exit 2; }
}
else
{
    print {*STDOUT}
      "\nUsage: collect_genomevirus.pl --prefix <output prefix> --query <NCBI entrez query> [..]\nDescription: Collect all complete genomes using the provided ENTREZ query and align them.\n\n--prefix <file prefix>\n    Filename output prefix. [mandatory]\n--query <query string>\n    NCBI entrez query. [mandatory]\n--after <year>\n    minimum date to retrieve from (YYYY).\n--before <year>\n    maximum date to retrieve from (YYYY).\n--norun\n    Disable the automatic run of R and the PCA and cluster analysis\n--verbose\n    Becomes very chatty.\n\n";
    exit 1;
}