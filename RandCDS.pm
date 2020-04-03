package RandCDS;
use Moose;
use Bio::CodonUsage::Table;
use Bio::Tools::CodonTable;
use Math::Random::MT::Auto qw /rand shuffle/;
use Data::Dumper;

has 'cutg' => ( isa => 'Ref', is => 'rw',  required => 1);
has 'gene_start' => ( isa => 'Int', is => 'rw',  required => 1);
has 'gene_end' => ( isa => 'Int', is => 'rw',  required => 1);

sub _codon_strings {
    my $self = shift;
    my @bases = qw (A T C G);
    my @codons;
    for (my $i=0; $i<=$#bases; $i++) {
	for (my $j=0; $j<=$#bases; $j++) {
	    for (my $k=0; $k<=$#bases; $k++) {
		my $codon = join '', ($bases[$i], $bases[$j], $bases[$k]);
		push @codons, {codon=>$codon, fit=>1, order=>undef};
	    }
	}
    }
    return \@codons;
}

sub get_rand_cds {
    my $self = shift;
    my $ref = $self->_codon_strings();
    my $cutg = $self->cutg;
    my $myCodonTable  = Bio::Tools::CodonTable->new(-id => 11);
    my $sum=0;
    my $length=$self->gene_end - $self->gene_start + 1;
    my (@seq, @codons);
    foreach (@$ref) {
	next if $myCodonTable->is_ter_codon($_->{codon});
	my $codon_count = int $cutg->codon_count($_->{codon});
	for (my $i=1; $i<= $codon_count; $i++) {
	    push @codons, $_;
	    $sum++;
	}
    }

    shuffle(\@codons);

    for (my $i=1; $i<= $length; $i++) {
	my $chosen=$codons[int rand($sum)];
	$chosen->{order} = $self->gene_start+$i-1;
	push @seq, $chosen;
    }
    return \@seq;
}


1;
