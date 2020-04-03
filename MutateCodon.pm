package MutateCodon;
use Moose;
use Bio::CodonUsage::Table;
use Data::Dumper;
use Math::Random::MT::Auto qw /rand shuffle/;
use Bio::Tools::CodonTable;

has 'source' => ( isa => 'Str', is => 'rw', required => 1);
has 'cutg' => ( isa => 'Ref', is => 'rw',  required => 1);
has 'titv' => ( isa => 'Int', is => 'rw',  required => 1, default => 1);
has 'neighbors' => ( isa => 'Ref', is => 'rw');
has 'dest' => ( isa => 'Str', is => 'rw');

sub _get_neighbors {
    my $self = shift; 
    my $source = $self->source; 
    my ($first, $second, $third) = split //, $source;
    my @neighbors;
    my @raw_bases = qw (A T C G);
    my %transition =   ( "A" => "G",
			 "G" => "A",
			 "C" => "T",
			 "T" => "C");

    foreach my $rawbase (@raw_bases) {
	next if $rawbase eq $first;
	my $ti = ($transition{$first} eq $rawbase) ? 1 : 0;
	push @neighbors, {
	    'codon' => $rawbase . $second . $third,
	    'ti' => $ti,
	}
    }

    foreach my $rawbase (@raw_bases) {
	next if $rawbase eq $second;
	my $ti = ($transition{$second} eq $rawbase) ? 1 : 0;
	push @neighbors, {
	    'codon' => $first . $rawbase . $third,
	    'ti' => $ti,
	}
    }

    foreach my $rawbase (@raw_bases) {
	next if $rawbase eq $third;
	my $ti = ($transition{$third} eq $rawbase) ? 1 : 0;
	push @neighbors, {
	    'codon' => $first . $second . $rawbase,
	    'ti' => $ti,
	}
    }
    $self->neighbors(\@neighbors);
}

sub pick_a_neighbor {
    my $self = shift;
    my $cutg = $self->cutg; 
    $self->_get_neighbors();
    my $ref_neighbor = $self->neighbors();
    my @codon_pool;
    my $sum=0;
    my $myCodonTable  = Bio::Tools::CodonTable->new(-id => 11);
    foreach (@$ref_neighbor) {
	next if $myCodonTable->is_ter_codon($_->{codon});
	my $codon_count = int $cutg->codon_count($_->{codon});
	$codon_count *= $self->titv if $_->{ti}; # adjusting titv
	$sum += $codon_count;
	for (my $i=1; $i<= $codon_count; $i++) {
	    push @codon_pool, $_->{codon};
	}
    }
    shuffle(\@codon_pool);
    my $pick = int rand($sum);
    $self->dest($codon_pool[$pick]);
}

1;
