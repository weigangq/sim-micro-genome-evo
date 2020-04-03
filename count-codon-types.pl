#!/usr/bin/perl -w

use Bio::SeqIO;
use Bio::Tools::CodonTable;

my $logfile = shift;
my $seqfile = shift;
my $myCodonTable  = Bio::Tools::CodonTable->new(-id => 11);
my (@ori, @seqs, %type_ct);
open LOG, "<".$logfile;
while (<LOG>) {
    next unless /^(\d+)\s+(\S\S\S)\s+(\S+)\s+\S\S\S$/;
    push @ori, {
	pos => $1,
	codon_ori => $2,
	constr_type => $3,
#	aa_ori => $myCodonTable->translate($2),
    };
    $type_ct{$3}++;
}

my $num_codons = scalar @ori;

my $in = Bio::SeqIO->new(-file=>$seqfile, -format=>'fasta');
while (my $seq = $in->next_seq) {
    push @seqs, {
	id => $seq->id(),
	codons => &split_seq($seq->seq),
    };
}

my @cts;
for (my $i=0; $i<$num_codons;$i++) {
    my %seen;
    foreach my $seq (@seqs) {
	$seen{$seq->{codons}->[$i]}++;
    }
    my @no_uniq = keys %seen;
    next if @no_uniq == 1 || @no_uniq > 2; # constant or more-than-2-states
    my $aa1 = $myCodonTable->translate($no_uniq[0]);
    my $aa2 = $myCodonTable->translate($no_uniq[1]);
    push @cts, {
	pos=>$i+1,
	con_type=>$ori[$i]->{constr_type},
	poly=> $aa1 eq $aa2 ? "syn": "non",
    }
}

#foreach my $type (keys %type_ct) {
#    print $type, "\t", $type_ct{$type}, "\n";
#}

foreach my $type (keys %type_ct) {
    my @non = grep {$_->{poly} eq "non" && $_->{con_type} eq $type} @cts;
    my @syn = grep {$_->{poly} eq "syn" && $_->{con_type} eq $type} @cts;
#    print $type, "\t", scalar @non, "\t", scalar @syn, "\n";
    printf "%d\t%.3f\t%.3f\n",  $type, scalar @non/$type_ct{$type}, scalar @syn/$type_ct{$type};
}

exit;

sub split_seq {
    my $dna = shift;
    my @codons;
    for (my $i=0; $i<length $dna; $i+=3) {
	push @codons, substr($dna, $i, 3);
    }
    return \@codons;
}
