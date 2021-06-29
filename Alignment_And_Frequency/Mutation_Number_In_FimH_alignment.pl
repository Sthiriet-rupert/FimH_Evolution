#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;
use Smart::Comments;
use Tie::IxHash;
use File::Basename;
use Path::Class 'file';

unless (@ARGV == 2) {
	die <<"EOT";
Usage: $0 <RefSeq.txt> <msa.fasta>

This tool counts the number of single nucleotide or amino-acid changing (no Indels) in a 
multiple sequence alignment in fasta format and returns the number of such changings in 
each structural parts of the sequence (signal peptide, domaine1, domaine2...).
	The first argument is the part of the alignment corresponding to the sequence of reference 
without the defline. For exemple, in the following fasta formated alignment:

	>seq1
	ATTAGAGCCTAG--CG
	>seq2
	ATGAG-GCCTAGTGCG
	>seq3
	ATGAG-GCGTAGTGCG

The RefSeq.txt file will only containe the sequence 'ATTAGAGCCTAG--CG' on a single line.
	The second argument is the fasta formated alugnment without the sequence of reference. In 
the exemple alignment the msa.fasta file will containe the sequences seq2 and seq3 but not the 
seq1:

	>seq2
	ATGAG-GCCTAGTGCG
	>seq3
	ATGAG-GCGTAGTGCG

	Using these two files, the tool will identify the variants present in each seqence relatively
to the sequence of reference. The number of such variant sites is outputed in a file for all sequences
and for the 3 structural parts of the fimH protein.
EOT
}


my $signal_seq_end = 23;
my $Dom1_end = 254;
my $linker_end = 258;
my $Dom2_end = 435;


my $ref_seq = file(shift)->slurp;
my $alignment = shift;

tie my %msa, 'Tie::IxHash';
%msa = read_fasta($alignment);

my ($basename, $dir, $suffix) = fileparse($alignment, qr{\.[^.]*}xms);
my $outfile = file($dir, $basename . '_Frequences.txt');

my $Nb_mut_signal = 0;
my $Nb_mut_Dom1 = 0;
my $Nb_mut_link = 0;
my $Nb_mut_Dom2 = 0;

my $len_al = length $ref_seq;

open my $out, '>', $outfile;
say {$out} '#' . join "\t", qw(Seq_id Nb_mut_Signal Nb_mut_Dom1 Nb_mut_linker Nb_mut_Dom2);

while ( my ($seq_name, $seq) = each %msa ) {
	
	SCAN:
	for my $i ( 0..$len_al-1 ) {
		my $ref_pos = substr($ref_seq,$i,1);
		my $seq_pos = substr($seq,$i,1);
		
		if ( $seq_pos eq "-" || $ref_pos eq "-" ) {
			$i +=1;
			next SCAN;
		}
		
		if ( $seq_pos eq $ref_pos ) {
			$i +=1;
			next SCAN;
		}
		
		if ( $i+1 <= $signal_seq_end ) {
			$Nb_mut_signal += 1;
		}
		
		if ( $signal_seq_end < $i+1 && $i+1 <= $Dom1_end ) {
			$Nb_mut_Dom1 += 1;
		}
		
		if ( $Dom1_end < $i+1 && $i+1 <= $linker_end ) {
			$Nb_mut_link += 1;
		}
		
		if ( $linker_end < $i+1 && $i+1 <= $Dom2_end ) {
			$Nb_mut_Dom2 += 1;
		}
		$i +=1;
	}
	say {$out} join "\t", $seq_name, $Nb_mut_signal, $Nb_mut_Dom1, $Nb_mut_link, $Nb_mut_Dom2;
	
	$Nb_mut_signal = 0;
	$Nb_mut_Dom1 = 0;
	$Nb_mut_link = 0;
	$Nb_mut_Dom2 = 0;
}

close $out;
	


sub read_fasta {


	my $infile = shift;
	open my $in, '<', $infile;

	my $seq_id;
	my $seq;
	tie my %seq_for, 'Tie::IxHash';		#Preserve original seq order

	LINE:
	while (my $line = <$in>) {
		chomp $line;

		#At each '>' char...
		if (substr($line, 0, 1) eq '>') {

			#Add current seq to hash (if any)
			if ($seq) {
				$seq_for{$seq_id} = $seq;
				$seq = q{};
			}

			# Extract new seq_id
			$seq_id = substr($line, 1);
			next LINE;
		}

		#Elongate current seq (seqs can be broken ob several lines)
		$seq .= $line;
	}

	#Add last seq to hash (if any)
	$seq_for{$seq_id} = $seq if $seq;

	close $in;
	return %seq_for;

}
