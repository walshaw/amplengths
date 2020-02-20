#!/usr/bin/perl

use 5.016;
use strict;
use warnings;
use autodie;

use Carp;
use Getopt::Euclid;
use Data::Dumper;
use List::MoreUtils qw( uniq );

use Bio::AlignIO;

my ($alignment_file, $alignment_format, $gap_chars,
    $fwd_primer_seq, $rev_primer_seq,
    $verbosity,
)
  = @ARGV{'--sequences', '--format', '--gap_chars',
          '--forward_primer', '--reverse_primer',
          '--verbosity'};

print Dumper \%ARGV if $verbosity > 2;


die qq{'$alignment_file' does not exist or is not a file.\n}
  if !$alignment_file;

my $aln_obj =
  get_alignment($alignment_file, $alignment_format);

show_alignment($aln_obj) if $verbosity > 3;

# Implication from Bio::SimpleAlign is that '-' is expected to
# be the sole gap char; therefore, if the user has specified
# that only '-' be treated as a gap char, then nothing needs
# doing. Otherwise, any other gaps chars must be changed to
# '-'.

resolve_gap_chars($aln_obj, $gap_chars) if $gap_chars ne q{-};

show_alignment($aln_obj) if $verbosity > 3;

if ($fwd_primer_seq || $rev_primer_seq) {

    primer_locations_in_aln($aln_obj,
                     $fwd_primer_seq, $rev_primer_seq);

}


exit 0;

sub get_alignment {

    my $alignment_file = shift or croak qq{no alignment file};
    my $alignment_format = shift || q{fasta};

    my $alnio = Bio::AlignIO->new(-file   => $alignment_file,
                                -format => $alignment_format);

    # should be just the one alignment in the file
    my $alignment = $alnio->next_aln(); # Bio::SimpleAlign object

    return $alignment;
}

sub resolve_gap_chars {

    my $aln_obj   = shift; # Bio::SimpleAlign
    my $gap_char_string = shift; # string of concatenated characters

    # Note that a space is a valid char here.

    my @gap_chars = uniq(split q{}, $gap_char_string);

    for my $gchar (@gap_chars) {

        next if $gchar eq q{-};

        # Important, as the first arg to Bio::SimpleAlign->map_chars() is
        # a regexp:
        my $gchar_qm = quotemeta $gchar;

        $aln_obj->map_chars($gchar_qm, q{-});

    }
}

sub show_alignment {

    # Mainly for debugging.

    my $aln_obj   = shift; # Bio::SimpleAlign
    my $n_seqs    = shift || 0;

    # $seq_obj is a Bio::LocatableSeq object.

    for my $seq_obj ($aln_obj->each_seq()) {

        my $seq = $seq_obj->seq();

        (my $seq_gapless = $seq) =~ s{ [.\-]+ }{}gxms;

        print STDERR $seq_obj->display_id(),
                     qq{\tlength }, $seq_obj->length(),
                     qq{\tgapless length }, length($seq_gapless),
                     qq{\n$seq\n};
        last if ! --$n_seqs;
    }

}

sub primer_locations_in_aln {

    my $gapped_aln_obj = shift;
    my @query_seqs = @_;

    my @gapped_seq_objs  = $gapped_aln_obj->each_seq();

    for my $gapped_seq_obj (@gapped_seq_objs) {

        # $gapped_seq_obj->seq() is a Bio::LocatableSeq object

        (my $gapless_seq_string = $gapped_seq_obj->seq() )
          =~ s{ [-]+ }{}gxms;

        # Coordinate numbering here is base 0.
        my ($coord_begin_aref, $coord_end_aref) =

        primer_locations_in_seq_string(
            $gapless_seq_string,
            @query_seqs
        );

        print $gapped_seq_obj->display_id() if ($verbosity);

        for my $i (0 .. $#query_seqs) {
            my ($begin_coord, $end_coord, $begin_column, $end_column)
              = (q{NA}) x 4;
            # The original coord can be 0. It's expected that
            # either both or neither of the begin and end
            # coords can be defined.
            # Note that here the string index values are converted into
            # 'coordinates' (i.e. add 1)
            if (defined $coord_begin_aref->[$i]) {
                $begin_coord = $coord_begin_aref->[$i] + 1;
                $end_coord   = $coord_end_aref->[$i] + 1;

                # These refer to alignment positions (columns); numbered
                # in base 1
                ($begin_column, $end_column) = map {
                    $gapped_seq_obj->column_from_residue_number($_)
                } ($begin_coord, $end_coord);


            }
            print qq{\t$query_seqs[$i]\t}.
                  qq{seq coords: $begin_coord .. $end_coord\t}.
                  qq{alignment cols: $begin_column .. $end_column}
              if $verbosity;
        }
        print qq{\n} if $verbosity;
    }

}

sub primer_locations_in_seq { # not used

    my $gapless_seq_obj = shift;
    my @query_seqs = @_;

    my @begin; # location of each of @query_seqs
    my @end;   # end-location of each of @query_seqs

    for my $qseq (@query_seqs) {

        my $location = location_of_primer($gapless_seq_obj, $qseq)
          if defined $qseq;

        if (!defined $qseq || ($location == -1)) {
            push @begin, undef;
            push @end  , undef;
            next;
        }

        push @begin, $location;
        push @end,   $location + length($qseq) - 1;
    }

    return \@begin, \@end;

}

sub primer_locations_in_seq_string {

    my $gapless_seq_string = shift;
    my @query_seqs = @_;

    my @begin; # location of each of @query_seqs
    my @end;   # end-location of each of @query_seqs

    for my $qseq (@query_seqs) {

        my $location = index $gapless_seq_string, $qseq
          if defined $qseq;

        if (!defined $qseq || ($location == -1)) {
            push @begin, undef;
            push @end  , undef;
            next;
        }

        push @begin, $location;
        push @end,   $location + length($qseq) - 1;
    }

    # stringwise (i.e. base 0) results are returned, since this process
    # a string

    return \@begin, \@end;

}

sub location_of_primer {

    my $gapless_seq_obj = shift;
    my $query_seq       = shift;

    print STDERR qq{Looking for '$query_seq' in '},
      $gapless_seq_obj->seq(), qq{'\n}
        if $verbosity > 3;

    # It's base-0 numbering of course.
    my $location = index $gapless_seq_obj->seq(), $query_seq;
    return $location;

}


__END__

=head1 NAME

alnampli.pl

=head1 VERSION

This is version 0.1

=head1 REQUIRED ARGUMENTS

=over

=item --sequences <file> | -s[eq[uence][s]] [=] <file>

File containing sequences. 

=for Euclid: 
    file.type:    readable
    file.default: '-'

=back


=head1 OPTIONS

=over

=item --format <format> | -form[at] [=] <format>

Sequence file format. Many different sequence file formats are recognized
(anything understood by Bioperl, specifically by the Bio::SeqIO module).

=for Euclid:
    format.type:    string
    format.default: 'fasta'

=item --begin <begin_at_seq> | -b[egin[_at]] <begin_at_seq>

=for Euclid:
    begin_at_seq.type:    0+integer
    begin_at_seq.default: 0

=item --n_sequences <n_seqs> | -n[_seq[uence][s]] <n_seqs>

=for Euclid:
    n_seqs.type: +integer

=item --get_columns <first_column>-<last_column> | -[get_]col[umn][s] <first_column>-<last_column>

=for Euclid:
    first_column.type: +integer
    last_column.type:  +integer

=item --forward_primer <fwd_primer_seq> | -f[[or]w[ar]d][_primer] <fwd_primer_seq>

=for Euclid:
    fwd_primer_seq.type: string

=item --reverse_primer <rev_primer_seq> | -r[ev[erse]][_primer] <rev_primer_seq>

=for Euclid:
    rev_primer_seq.type: string

=item --gap_chars <gap_chars> | -gap[_char][s] <gap_chars>

=for Euclid:
    gap_chars.type: string
    gap_chars.default: '-.'

=item --verbosity <verbosity> | -verb[os[e|ity]] <verbosity>

=for Euclid:
    verbosity.type: 0+integer
    verbosity.default: 1

=back

=head1 AUTHOR

John Walshaw

=head1 BUGS

=head1 COPYRIGHT

Copyright (c) 2020, John Walshaw. All Rights Reserved.
This module is free software. It may be used, redistributed
and/or modified under the terms of the Perl Artistic License
(see http://www.perl.com/perl/misc/Artistic.html)


