#!/usr/bin/env perl
# 
# Forward simulation of bact evolution
# 
#
use warnings;
use strict;
use lib '/home/weigang/SimPopGenome';
use Bio::PopGen::Population;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Algorithm::Numerical::Sample qw /sample/;
use Bio::CodonUsage::IO;
use Data::Dumper;
use RandCDS;
use Bio::PopGen::Utilities;
use Bio::PopGen::Statistics;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;
use Convert;
use FitFunction;
use Math::Random::MT::Auto qw(rand shuffle poisson exponential);
use Bio::Tools::CodonTable;
use Statistics::Basic qw(:all);
use CodonType;
use MutateCodon;
################################################################################
# Option parsing
################################################################################
my %opts;
GetOptions (\%opts, 
	    "help|h",
	    "verbose|V",
	    "popsize|N=i", 
	    "numcodon|l=i", 
	    "genomesize|w=i",
	    "mutrate|m=f",
	    "selection_cof|s=f",
	    "recomrate|r=f",
	    "gen|g=i",
#	    "composition|w") or pod2usage(2);
    );
pod2usage(1) if $opts{"help"};
################################################################################
# Initital parameters: CHECK CAREFULLY BEFORE RUN
################################################################################
my $pop_size = $opts{popsize} || 100;
my $gen = $opts{gen} || $pop_size * 6;
my $num_codons = $opts{numcodon} || 1000; # total number of codons
my $mut_per_seq = $opts{mutrate} || 0.1;
my $rec_per_seq = $opts{recomrate} || 0;
my $titv_ratio = 4;
my $sample_size = 50;
my $sample_freq = 100;
#my $s_uniform = $opts{selection_cof} || 0.02;
my $s_positive = 0.2;
my $s_negative = 0.2;
my $tract_len = 50; # unit: num of codons
my @genes;
&report_para() if $opts{verbose};
my $ospC = { # under diversifying/positive selection
    start => 1,
    end => 500,
    prop_neg => 0.3,
    prop_neutral => 0.2,
    prop_pos => 0.5,
};

push @genes, $ospC;

my $ospA = { # under purifying/negative selection
    start => 501,
    end => 1000,
    prop_neg => 0.8,
    prop_neutral => 0.2,
    prop_pos => 0,
};

push @genes, $ospA;

#my $dbpA = { # under slightly positive selection
#    start => 801,
#    end => 1200,
#    prop_neg => 0.6,
#    prop_neutral => 0.3,
#    prop_pos => 0.1,
#};
#push @genes, $dbpA;
my $blk_starts = [ (0, 500, 1000) ]; # start at 0 (not 1)
################################################################################
# END Initital parameters: CHECK ABOVE CAREFULLY BEFORE RUN
################################################################################
# Get CUTG
my $io=Bio::CodonUsage::IO->new(-file => shift @ARGV); 
my $cdtable=$io->next_data();
# Build founder 
my @founder_codons;
my @constrain;
my @founders;

foreach my $gene (@genes) {
    my $founder_partial = RandCDS->new(
	gene_start => $gene->{start},
	gene_end => $gene->{end},
	cutg => $cdtable
	);
    
    my $ran_cds = $founder_partial->get_rand_cds(); 
    push @founder_codons, @{$ran_cds}; 
    
    my $constr_part = CodonType->new(
	gene_start => $gene->{start},
	gene_end => $gene->{end},
	prop_neg => $gene->{prop_neg}, 
	prop_neutral => $gene->{prop_neutral}, 
	prop_pos => $gene->{prop_pos},
	codons => $ran_cds,
	);
    
    push @constrain, @{$constr_part->get_status()};
}

for (my $i=1; $i<=$pop_size; $i++) {
    push @founders, {
	codons => \@founder_codons,
	fitness => undef,
	id => "L".$i,
    };
}

warn "Founder pop initialized\n" if $opts{verbose};
if ($opts{verbose}) {
    for (my $i=0; $i<$num_codons; $i++) {
        print STDERR $constrain[$i]->{site}, "\t", $constrain[$i]->{ori}, "\t", $constrain[$i]->{type}, "\t", $founders[0]->{codons}->[$i]->{codon}, "\n";
    } 
}

################################################################################
# Stage I. Pre-Speciation: evolve to mutation-selection equilibrium
################################################################################
warn "\nEvolution for ".  $gen . " generations\n" if $opts{verbose};
&evolve($gen);
warn "Simulation finished\n" if $opts{verbose};
#&report_subs("Post-B",\@constrain) if $opts{verbose};
exit;

sub report_para {
    print STDERR "Pop size: $pop_size haploid individuals\n";
    print STDERR "Total generation: $gen\n";
    print STDERR "Num of Codons: $num_codons codons\n";
    print STDERR "Mutation rate: $mut_per_seq per gen per seq\n";
    print STDERR "Recom rate: $rec_per_seq per gen per seq\n";
    print STDERR "Ti:Tv=$titv_ratio\n";
    print STDERR "Sample size: $sample_size\n";
    print STDERR "Sample frequency: $sample_freq\n";
    print STDERR "Selective coefficient: Positive-$s_positive per nonsyn\n";
    print STDERR "Selective coefficient: Negative-$s_negative per nonsyn\n";
    print STDERR "Tract length for gene conversion: $tract_len codons\n";
}

sub report_subs {
    my $label=shift;
    my $ref=shift;
    foreach my $con (@$ref) {
	if ($con->{subs}) {
	    foreach (@{$con->{subs}}) {
		print STDERR join "\t", ($label, $_->{pos}, $_->{parent_id}, $_->{gam_id}, $_->{type});
		print STDERR "\n";
	    } 
	}
    }
}

sub report_aln {
    my $id = shift;
    my $ref= shift;
    my @sample = sample (-set =>$ref, -sample_size=> $sample_size);
    my $aln = Bio::SimpleAlign->new();
    foreach my $ind ( @sample) { 
	my $lineage=$ind->{id};
	my $seq_str;
	foreach (@{$ind->{codons}}) { $seq_str .= $_->{codon} }
	my $seq = Bio::LocatableSeq->new(-id=>$lineage, -seq=>$seq_str);
	$aln->add_seq($seq);
    }

    $aln->set_displayname_flat();
    my $out=Bio::AlignIO->new(-file=>">". "G_$id.fas", -format=>'fasta');
    $out->write_aln($aln);
}

sub report_fitness {
    my $id = shift;
    my $ref = shift;
    print STDERR "\t\tGeneration $id: ", "Fitness of gametes (original=1):\t";
    my @fits;
    foreach my $ind ( @$ref) { 
	push @fits, &gam_fitness($ind);	       
#	my $id=$ind->{id};
#	my $fit=sprintf "%.6f", &gam_fitness($ind);		
#	print STDERR "$fit,";	
    }
    my $meanfit = mean(@fits);
    my $vector = $meanfit->query_vector;
    printf STDERR "%.4f\t%.4f\n", $meanfit, stddev($vector);
}

sub evolve {
    my $generation = shift;
    shuffle(\@founders);
    for (my $g=1; $g<=$generation; $g++) {
	warn "\tGeneration\t$g...\n" if $opts{verbose};
	my $pa_ct = 1;
	my @gamete_pool;
	while ($pa_ct <= $pop_size) {
	    my ($pa1, $pa2) = sample (-set =>\@founders, -sample_size=> 2);
#	    &report_aln($label . "$g-$pa_ct", \@founders) if $opts{verbose};
	    my $gamete = { # Need to create a new hash: do NOT modify founders
		codons => $pa1->{codons},
		fitness => undef,
		id => $pa1->{id}
	    };
	    $gamete = &recombine($g, $gamete, $pa2);
	    $gamete = &mutate($gamete);
	    my $fit = &gam_fitness($gamete); 
	    my $cutoff_fit = rand;
	    if ($cutoff_fit <= 1-2.718**(-1*$fit)) {
		$gamete->{id} .= "_". $pa_ct; 
		push @gamete_pool, $gamete;
#		warn "\t\tGot parent No.". $pa_ct . ": ". $gamete->{id} . "\n" if $opts{verbose};
		$pa_ct++;
	    }
	}
	@founders = @gamete_pool;
	if ($g % $sample_freq == 0 and $opts{verbose}) {
	    &report_fitness($g, \@gamete_pool); 
	    &report_aln($g, \@founders);
	}
#	&report_stat($label, \@gamete_pool)  if $opts{verbose};
    }
#    &report_aln($label, \@founders);
}

sub mutate { # assign poisson-distributed number of mutations (more efficient than testing each codon)
    my $ind = shift;
    my $num_mut = poisson($mut_per_seq);
    return $ind unless $num_mut;
#    warn "\tMutate:\t".$ind->{id}."\n";
    my @codons = @{$ind->{codons}}; 
    for (my $i=1; $i<=$num_mut; $i++) {
	my $pick = int rand ($num_codons);
	my $source_codon = $codons[$pick]->{codon};
	my $mut = MutateCodon->new(source=>$source_codon, cutg=>$cdtable, titv=>$titv_ratio);
	$mut->pick_a_neighbor();
	my $des_codon = $mut->dest();
#	print STDERR "\t\tReplace " . $codons[$pick]->{codon} . "\t";
	my $myCodonTable  = Bio::Tools::CodonTable->new( -id => 11 );
	my $aa_ori = $myCodonTable->translate($constrain[$pick]->{ori}); # compare with original
	my $aa_des = $myCodonTable->translate($des_codon);
	my $fit = FitFunction->new(
	    aa_ori => $aa_ori,
	    aa_des => $aa_des,
	    type => $constrain[$pick]->{type},
	    s_uniform_neg => $s_negative,
	    s_uniform_pos => $s_positive,
	    );
	splice @codons, $pick, 1, {codon=>$des_codon, order=>$pick+1, fit=>$fit->find_fitness()};
#	warn "with " . $mut->dest() . " fitness ". $fit->find_fitness(), " at $pick\n";
    }
    $ind->{codons} = \@codons;
    return $ind;
}

sub recombine {
    my ($g, $acceptor, $donor) = @_;
    return $acceptor if &same_codons($acceptor->{codons}, $donor->{codons});
    my $cut = rand;   
    return $acceptor unless $cut <= 1-2.718**(-1*$rec_per_seq);
#    warn "\tRecombine:\t".$acceptor->{id}."\n";
    my $rec = Convert->new(
	acceptor=>$acceptor,	    
	donor=>$donor, 
	num_codons=>$num_codons,
	gen=>$g,
	tract_length => $tract_len,
	block_starts => $blk_starts,
	);
    return  $rec->gene_con();
}

sub same_codons {
    my ($ref1, $ref2) = @_;
    for (my $i=0; $i<$num_codons; $i++) {
	return 0 if $ref1->[$i]->{codon} ne $ref2->[$i]->{codon};
    }
    return 1;
}

sub gam_fitness {
    my $gamete = shift;
    my $fitness=1;
    my @codons = @{$gamete->{codons}};
    foreach (@codons) {
	$fitness *= $_->{fit};
    }
    return $fitness;
}
=cut
