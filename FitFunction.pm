package FitFunction;
use Moose;
use Data::Dumper;
use Math::Random::MT::Auto qw /rand/;
use strict;
use Bio::Variation::Allele;
use Bio::Variation::AAChange;

has 'aa_ori' => ( isa => 'Str', is => 'rw', required => 1);
has 'aa_des' => ( isa => 'Str', is => 'rw', required => 1);
has 'type' => ( isa => 'Num', is => 'rw', required => 1); 
has 's_uniform_neg' => ( isa => 'Num', is => 'rw');
has 's_uniform_pos' => ( isa => 'Num', is => 'rw');

sub _neg_mat { # blosum62->selection coefficient
    my $self=shift;
    my %score = (
	-4 => 0.5,
	-3 => 0.4,
	-2 => 0.3,
	-1 => 0.2,
	0 => 0.1,
	1 => 0.05,
	2 => 0.02,
	3 => 0.01
	);
    return \%score;
}


sub _pos_mat { # negative blosum are encouraged
    my $self=shift;
    my %score = (
	-4 => 1,
	-3 => 0.9,
	-2 => 0.8,
	-1 => 0.7,
	0 => 0.6,
	1 => 0.5,
	2 => 0.4,
	3 => 0.3
	);
    return \%score;
}

sub _fitness {
    my $self=shift;
    my $type = shift; 
#    my $aa_x=shift;
#    my $aa_y=shift;
#    my %neg_score = %{$self->_neg_mat()}; 
#    my %pos_score = %{$self->_pos_mat()}; 
    return 1 if $type == 0; # neutral, w=1 model
    return 1-$self->s_uniform_neg() if $type == -1; # uniform negative selection coefficeint
#    return 1-$neg_score{$self->_blosum($aa_x, $aa_y)} if $model == 2; # blosum-weighted s
    return 1+$self->s_uniform_pos() if $type == 1; # uniform positive selection coefficeint
#    return 1+$pos_score{$self->_blosum($aa_x, $aa_y)} if $model == 4; # blosum-weighted s
}

sub _blosum {
    my $self = shift;
    my $aa1 = shift;
    my $aa2 = shift;
    my $allele1 = Bio::Variation::Allele->new ( -seq => $aa1,
                                                -id  => 'dummy1',
                                                -alphabet => 'protein',
                                                -is_reference => 1
        );
    
    my $allele2 = Bio::Variation::Allele->new ( -seq => $aa2,
                                                -id  => 'dummy2',
                                                -alphabet => 'protein',
                                                -is_reference => 1
        );
    
    my  $aamut = Bio::Variation::AAChange->new
        ('-start'         => 1,
         '-end'           => 1,
         '-length'        => 1,
         '-isMutation'    => 1,
         '-mut_number'    => 1,
         allele_ori => $allele1,
         allele_mut=>$allele2,
        );
    return $aamut->similarity_score;    
}

sub find_fitness { 
    my $self = shift;  
    my $type = $self->type;
    my $aa_ori = $self->aa_ori;
    my $aa_des = $self->aa_des;
    my $fit;

    return 1  if ($aa_ori eq $aa_des) or $type == 0;
    return 1 - $self->s_uniform_neg if $type < 0;
    return 1 + $self->s_uniform_pos if $type > 0;
}

1;
