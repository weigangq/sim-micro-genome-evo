package Recombine;
# prob=1-e(-r), where r is recom per gen per seq (ref. Padhukasahasram et al Genetics 2008)
# Limitation 1: only split between codons (otherwise fitness hard to reconstitute)

use Moose;
use Data::Dumper;
use Math::Random::MT::Auto qw /rand/; # a much better RNG: no timestamp dependence
use strict;

has 'donor' => ( isa => 'Ref', is => 'rw', required => 1);
has 'acceptor' => ( isa => 'Ref', is => 'rw', required => 1);
has 'num_codons' => ( isa => 'Num', is => 'rw', required => 1);
has 'gen' => ( isa => 'Num', is => 'rw', required => 1);

sub recombine { # replace back frag with that from a randomly chosen donor, with replacement
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
    my $breakpoint = int rand($self->num_codons); # possible to be zero
    my @backend = splice @donor_codons, $breakpoint;
#    my @backend_ids = splice @donor_ancs, $self->gen-1;
#    my $donor_str; foreach (@gam_codons) { $donor_str .= $_->{codon} . ' ' }
#    warn "\t\treplace:\t$donor_str\n";
#    my $back_str; foreach (@backend) { $back_str .= $_->{codon} . ' ' }
#    warn "\t\twith:\t$back_str at $breakpoint\n";
    splice @gam_codons, $breakpoint, $self->num_codons-$breakpoint, @backend;
#    splice @gam_ancs, $self->gen-1, $self->num_codons-$breakpoint, @backend_ids; 
#    my $new_str; foreach (@gam_codons) { $new_str .= $_->{codon} . ' ' }
#    warn "\t\tgot:\t$new_str\n";
    $newgam->{codons} = \@gam_codons;
#    $newgam->{id} = join "_", @gam_ancs;
    return $newgam;
}


1;
