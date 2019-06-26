#!/usr/bin/perl

use strict;
use warnings;
use autodie;

use Carp;
use Getopt::Long;

=pod

    This script relates to the overall purpose of determining lengths of
    amplicon sequences from a database of reference (barcode) sequences.

    For now, this just parses output from EMBOSS fuzznuc, which specifies the
    locations, in each sequence in the database, of forward and reverse primer
    sequences, in order to determine:

=over

=item    which sequences have no primers
=item    which sequences have 1 (or more) forward primers
=item    which sequences have 1 (or more) reverse primers
=item    which sequences have at least 1 forward and at least 1 reverse primer
=item    the lengths of the resulting expected amplicons

=back

    Counts of the above are reported, and conditionally, lists of the sequence
    IDs. Sequences with 0 primers can be known only if either -ids_file or
    -seqs_file are used.

    If there is a forward primer but no reverse primer, then the assumption
    can be made (conditionally on using -trunc3) that the amplicon extends to
    the 3' end of the sequence (that was parsed by fuzznuc), i.e. the database
    sequence is truncated towards the 3' end, before the reverse primer
    sequence (or reverse complement thereof).

    Likewise, if there is a reverse primer but no forward, then the assumption
    can be made (conditionally on using -trunc5) that the amplicon extends to
    the 5' end of the sequence (that was parsed by fuzznuc), i.e. the database
    sequence is truncated towards the 5' end, 3' to the reverse primer sequence.

    In both the above cases, a notional amplicon length (albeit probably
    shorter than the actual amplicon would be, due to the truncation) can be
    calculated.

=cut


my $forward_file;
my $reverse_file;
my $sequence_id_file;
my $fasta_file;

my $assume_5pr_truncated;
my $assume_3pr_truncated;

my $max_mismatches;

my $verbosity = 1;


my $usage = qq{Usage:\n\n$0 [ options ] 
};

die $usage if !GetOptions(
	'forward=s'	=>	\$forward_file,
	'reverse=s'	=>	\$reverse_file,
	'ids_file|id_file=s'	=>	\$sequence_id_file,
	'seqs_file|seq_file=s'	=>	\$fasta_file,
	'trunc5'	=>	\$assume_5pr_truncated,
	'trunc3'	=>	\$assume_3pr_truncated,
	'verbosity=i'	=>	\$verbosity,
);

die $usage if !$forward_file || !$reverse_file;



sub parse_fuzznuc {

    my $file = shift or croak qq{supply file name};

    open my $fh, '<', $file;

    my $seq_id;
    my $seq_length;
    my $is_complement;
    my $n_hits;
    my $in_table;

=pod

    Keys of %match_locations are sequence IDs; values are refs to arrays of
    refs to hashes, whose keys are 'begin', 'end', 'mismatches', 'strand',
    'sequence'. The first 3 of those keys have integer values, while the
    latter 2 have string values.

=cut

    my %match_locations;

    LINE:
    while (defined (my $line = <$fh>)) {

        if (($line !~ m{ \S }xms) || ($line =~ m{ \A [#] \S }xms)) {
            undef $in_table if $in_table;
            next LINE;
        }

        chomp $line;

        ($line =~ m{ \A [#] \s (\w+) [:] (\S+) \s+ (.*) \z }xms) && do {

            my ($key, $value1, $remainder) = ($1, $2, $3);

            KEYWORD: {
                ($key eq 'Sequence') && do {

                    $seq_id = $value1;

                    croak qq{unexpected line format:\n$line}
                      if $remainder !~ m{ \A from[:] \s+ (\d+) \s+
                                               to[:] \s+ (\d+) \s* \z }xms;

                    my ($begin, $end) = ($1, $2);

                    print STDERR
                      qq{WARNING: 'from' value is $begin in line:\n$line\n}
                      if $begin != 1;

                    $seq_length = $end;

                    last KEYWORD;
                };

                ($key eq 'HitCount') && do {
                    $n_hits = $value1;
                    last KEYWORD;
                };

                ($key eq 'Complement') && do {
                    ($value1 eq 'Yes') && do {
                        $is_complement = 1;
                        last KEYWORD;
                    };
                    ($value1 eq 'No') && do {
                        $is_complement = 0;
                        last KEYWORD;
                    };
                    croak qq{unexpected line format:\n$line};
                };

                ($key eq 'Total_sequences') && do {
                    last KEYWORD;
                };

                ($key eq 'Total_length') && do {
                    last KEYWORD;
                };

                ($key eq 'Reported_sequences') && do {
                    last KEYWORD;
                };

                ($key eq 'Reported_hitcount') && do {
                    last KEYWORD;
                };

                ($key =~ m{ \A (Program|Rundate|Commandline|
                                   Report_format|Report_file) \z }xms) && do {
                    last KEYWORD;
                };

                croak qq{unexpected line format (unknown keyword '$key'):\n$line};

            } # END block KEYWORD

            next LINE;

        }; # end of do { } dealing with meaningful '#' lines

        ($line =~ m{ \A \s+ Start \s+ End \s+ Strand \s+ Pattern_name \s+
                            Mismatch \s+ Sequence \s+ \z }xms) && do {
            croak qq{unexpected header line: already in table:\n$line}
              if $in_table;

            $in_table++;

            next LINE;
        };

        ($line =~ m{ \A \s* (\d+) \s+ (\d+) \s+ (.) \s+ (\S+) \s+
                            (\d+|[.]) \s+ (\w+) \s* \z }xms) && do {
            croak qq{unexpected table row (no header encountered yet):\n$line}
              if !$in_table;

            my ($begin, $end, $strand, $pattern_name, $n_mismatches, $sequence)
              = ($1, $2, $3, $4, $5, $6);
            # .... do stuff
            next LINE;
        };

        croak qq{unexpected line format:\n$line};

    } # END block LINE (iterates over each line of input file)

    close $fh;

}
