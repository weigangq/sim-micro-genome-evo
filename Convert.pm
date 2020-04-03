package Convert;
# Gene Conversion algorithms:
# Wief & Heil 2000 (GenCon)
# Didelot et al 2007 & 2009 (SimMLST)
# Blocks of genes: conversion affecting a single block only
#
use Moose;
use Data::Dumper;
use Math::Random::MT::Auto qw /rand exponential/; # a much better RNG: no timestamp dependence
use strict;

has 'donor' => ( isa => 'Ref', is => 'rw', required => 1);
has 'acceptor' => ( isa => 'Ref', is => 'rw', required => 1);
has 'num_codons' => ( isa => 'Num', is => 'rw', required => 1);
has 'gen' => ( isa => 'Num', is => 'rw', required => 1);
has 'tract_length' => ( isa => 'Num', is => 'rw', required => 1);
has 'block_starts' => ( isa => 'Ref', is => 'rw', required => 1);

sub gene_con { 
    my $self=shift;
    my $donor=$self->donor();
    my $gam=$self->acceptor(); 
#    my $gam_id = $gam->{id};
#    my @gam_ancs = split /_/, $gam_id;
    my $newgam = $gam;
#    my $donor_id = $donor->{id};
#    my @donor_ancs = split /_/, $donor_id;
    my @donor_codons = @{$donor->{codons}};
    my @gam_codons = @{$gam->{codons}};

    my $bk = $self->pick_a_block(); # only one gene are converted (don't go between genes)
    my ($tract_start, $tract_end) = $self->get_tract($bk);

    my @imports = @donor_codons[$tract_start..$tract_end];
    my @to_be_replaced = @gam_codons[$tract_start..$tract_end];
    return $newgam if $self->same_codon(\@imports, \@to_be_replaced, $tract_end-$tract_start);
#    my @backend_ids = splice @donor_ancs, $self->gen-1;
#    my $accep_str; foreach (@to_be_replaced) { $accep_str .= $_->{codon} . ' ' }
#    warn "\tRecombine: $gam_id\n\t\trep:\t$accep_str\n";
#    warn "\tRecombine: $gam_id\n";
#    my $imp_str; foreach (@imports) { $imp_str .= $_->{codon} . ' ' }
#    warn "\t\twith:\t$imp_str between coords $tract_start and $tract_end\n";
#    warn "\t\tbetween coords $tract_start and $tract_end\n";
    splice @gam_codons, $tract_start, $tract_end-$tract_start+1, @imports;
#    splice @gam_ancs, $self->gen-1, $self->num_codons-$breakpoint, @backend_ids; 
#    my $new_str; foreach (@gam_codons) { $new_str .= $_->{codon} . ' ' }
#    warn "\t\tgot:\t$new_str\n";
    $newgam->{codons} = \@gam_codons;
#    $newgam->{id} = join "_", @gam_ancs;
    return $newgam;
}

sub get_tract { # Didelot et al 2007 & 2009 
    my $self = shift;
    my $block = shift;
    my @block_starts = @{$self->block_starts()}; # beginning to end coords
    my $cut = rand;
    my $size = $self->tract_length();

    my $p_begin = $size / ($self->num_codons-$size-$#block_starts*$size); 
    my $actual_size = int exponential($size);
    my ($bg, $ed);

    if ($cut <= $p_begin) {
	$bg = $block_starts[$block];
    } else {
	$bg = $block_starts[$block] + int rand ($block_starts[$block+1] - $block_starts[$block]+1);
    }

    if ($actual_size >= $block_starts[$block+1] - $bg ) {
	$ed = $block_starts[$block+1]-1; # tract doesn't go beyond a gene block
    } else {
	$ed = $bg + $actual_size-1;
    }
    return ($bg, $ed);
}

sub pick_a_block {
    my $self = shift;
    my @blocks = @{$self->block_starts()};
    my $pick = int rand($self->num_codons);
    my $track = 0;
    for (my $i=1; $i < $#blocks; $i++) { # start from 2nd, end at 2nd last
	if ($pick < $blocks[$i]) {
#	    warn "\t\twithin block $track\n";
	    return $track;
	} else {
	    $track++;
	}
    }
#    warn "\t\twithin block $track\n";
    return $track;
}

sub same_codon {
    my $self = shift;
    my ($ref1, $ref2, $len) = @_;
    for (my $i=0; $i<$len; $i++) {
	return 0 if $ref1->[$i]->{codon} ne $ref2->[$i]->{codon};
    }
    return 1;
}

1;
