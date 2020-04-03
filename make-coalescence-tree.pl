#!/usr/bin/env perl
# Produce a coalescent tree from simulated seqs in FASTA output
#  e.g., Post-A.fas
#
use Bio::SeqIO;
use strict;
use Data::Dumper;
use Getopt::Std;
use Bio::Tree::Tree;
use Bio::Tree::Node;
use Bio::TreeIO;
use warnings;
use 5.10.0;

my%opts;
&getopts('ABVt', \%opts);
my $file = shift;
my $in=Bio::SeqIO->new(-file=>$file);
my $allele=1;
my $last_gen=shift;
my (@lineages, @coales_events, @gens);
while (my $seq=$in->next_seq) {
    my $id=$seq->display_id;
    $id =~ s/^L//;
    my @anc = split /_/, $id;
    push @lineages, {
	allele=>$allele++,
	ancestors=>\@anc,
    }
}

if ($opts{V}) {
    for (my $i=0; $i<=$last_gen; $i++) {
	print $i, "\t";
	foreach (@lineages) {
	    print $_->{ancestors}->[$i], "\t";
	}
	print "\n";
    }
}

my $coal_ct=0;
my @coals;
my $root=Bio::Tree::Node->new(-id=>'root');

&split_lineage(0, \@lineages, $root);
my $tree=Bio::Tree::Tree->new(-id=>$file, -root=>$root, -nodelete=>1);
if ( $opts{t} ) {
    my $tr_str = $tree->as_text("newick");
    $tr_str =~ s/^\(//;
    $tr_str =~ s/\)root//;
    say $tr_str;
}

my $first_coal=shift @coals;
if ($opts{B}) { # before first coalescence
    foreach (@lineages) {
	my @anc = @{$_->{ancestors}};
	for (my $i=0; $i<$first_coal; $i++) {
	    my $pa = "L";
	    my @subs = @anc[0..$i];
	    my $name = join "_", @subs;
	    $pa .= $name;
	    my $ch = $pa."_".$anc[$i+1];
	    print $pa, "\t", $ch, "\n";
	}
    }
}

if ($opts{A}) { # after first coalescence
    foreach (@lineages) {
	my @anc = @{$_->{ancestors}};
	for (my $i=$first_coal; $i<$last_gen; $i++) {
	    my $pa = "L";
	    my @subs = @anc[0..$i];
	    my $name = join "_", @subs;
	    $pa .= $name;
	    my $ch = $pa."_".$anc[$i+1];
	    print $pa, "\t", $ch, "\n";
	}
    }
}

exit;

sub split_lineage {
    my $begin= shift;
    my $ref=shift;
    my $parent=shift;
    my @set = @$ref;
    for (my $i=$begin; $i<=$last_gen; $i++) {
	my %seen;
	foreach (@set) {   $seen{$_->{ancestors}->[$i]}++;  }
	my @mems = keys %seen;
	next if @mems == 1;
	$coal_ct++;
	push @coals, $i-1;
	print "Coalescent event $coal_ct at ", $i-1, "\n" if $opts{V}; 
	my $nd = Bio::Tree::Node->new(-id=>"co_".$coal_ct, -branch_length=>$i-$begin+1);
	$parent->add_Descendent($nd);
	
	foreach my $pid (@mems) {
	    my @subset =  grep {$_->{ancestors}->[$i] == $pid} @set;
	    if (@subset > 1) {
		&split_lineage($i+1, \@subset, $nd);
	    } else { # reached an OTU
		my $ln = shift @subset;
		my $otu = Bio::Tree::Node->new(-id=>"S" . sprintf "%03s", $ln->{allele});
		$otu->branch_length($last_gen-$i+1);
		print "$i\tadding OTU:\t", $otu->id(), "\t", $otu->branch_length(), "\n" if $opts{V}; 
		$nd->add_Descendent($otu);
	    }
	}

	last; # stop at $i; don't proceed further
    } 
}


exit;
