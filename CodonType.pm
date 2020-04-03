package CodonType;
use Moose;
use Math::Random::MT::Auto qw /rand/;
use Data::Dumper;

has 'gene_start' => ( isa => 'Int', is => 'rw',  required => 1);
has 'gene_end' => ( isa => 'Int', is => 'rw',  required => 1);
has 'prop_neg' => ( isa => 'Num', is => 'rw', required=>1);
has 'prop_neutral' => ( isa => 'Num', is => 'rw', required=>1);
has 'prop_pos' => ( isa => 'Num', is => 'rw', required=>1);
has 'codons' => ( isa => 'Ref', is => 'rw', required=>1);

sub get_status {
    my $self=shift;
    my @constant_sites;
    my @codons = @{$self->codons};
    my $length=$self->gene_end - $self->gene_start + 1;
    for (my $i=0; $i< $length; $i++) {
	my $pos = $self->gene_start + $i;
	push @constant_sites, {
	    type => $self->_site_type(),
	    site => $pos,
	    ori=> $codons[$i]->{codon}}; 
    }    
    return \@constant_sites;
}

sub _site_type {
    my $self = shift;
    my $cut = rand;
    my $p1 = $self->prop_neg();
    my $p2 = $self->prop_neutral();
    my $p3 = $self->prop_pos();
    if ($cut <= $p1) {
	return -1; 
    } elsif ($cut <= $p1+$p2) {
	return 0;
    } else { return 1 }
}

1;
