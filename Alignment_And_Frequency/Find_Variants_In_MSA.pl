#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;
use Smart::Comments;
use Tie::IxHash;
use File::Basename;
use Path::Class 'file';

unless (@ARGV == 3) {
	die <<"EOT";
Usage: $0 <RefSeq.txt> <msa.fasta> <VcfHeader.txt>

This tool identifies variants in a multiple sequence aligment (msa) in fasta format relatively 
to a sequence of reference of this alignment.
	The first argument is the part of the alignment corresponding to the sequence of reference 
without the defline. For exemple, in the following fasta formated alignment:

	>seq1
	ATTAGAGCCTAG--CG
	>seq2
	ATGAG-GCCTAGTGCG
	>seq3
	ATGAG-GCGTAGTGCG

The RefSeq.txt file will only containe the sequence 'ATTAGAGCCTAG--CG' on a single line.

	The second argument is the fasta formated alignment without the sequence of reference. In 
the exemple alignment the msa.fasta file will containe the sequences seq2 and seq3 but not the 
seq1:

	>seq2
	ATGAG-GCCTAGTGCG
	>seq3
	ATGAG-GCGTAGTGCG

	Using these two files, the tool will identify the variants present in each seqence relatively
to the sequence of reference. These variants will be outputed in a vcf file (one file for each
sequence) using the file provided by the third argument as header. The minimal header should be as 
follows:

##fileformat=VCFv4.2
##fileDate=20190218
##source=FindVariantsInMSA.pl
##reference=E.coli_FimH.fasta
##contig=<ID=FimH,length=36,species="E. coli",taxonomy=x>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">

Just set the contig length to the number of alignment positions. The other informations are not
mendatory.
These vcf files could then be merged using 'bcftools merge' to obtain a vcf file containing all 
variants in all sequences.
EOT
}




my $ref_seq = file(shift)->slurp;
my $alignment = shift;
my $header = file(shift)->slurp;
chomp $header;

tie my %msa, 'Tie::IxHash';
%msa = read_fasta($alignment);



while ( my ($seq_name, $seq) = each %msa ) {

	my ($basename, $dir, $suffix) = fileparse($alignment, qr{\.[^.]*}xms);
				
	my @words = split ' ', $seq_name;
	my $short_seq_name = shift @words;
	my $outfile = file($dir, $short_seq_name . '_Variants.vcf');

	open my $out, '>', $outfile;
	
	say {$out} $header;
	say {$out} '#' . join "\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT), . "\t$short_seq_name";

	my $len_al = length $ref_seq;
	my $del_size = 0;
	my $del = 0;
	my $in_size = 0;
	my $in = 0;
	my $ref_pos;
	my $seq_pos;
	
	SEQ_SCAN:
	for my $i (0..$len_al-2) {
		$ref_pos = substr($ref_seq,$i,1);
		$seq_pos = substr($seq,$i,1);
		
		if ( $del ) {
			if ( $seq_pos eq "-" && $ref_pos ne "-" ) {
				$del_size+=1;
				
				
				if ( $i eq $len_al-2 ) {
			
					my $del_pos = $i - $del_size;
					my $del_REF = substr($ref_seq,($del_pos - 1),($del_size + 1));
					my $del_ALT = substr($seq,($del_pos - 1),1);
			
					say $out "FimH" . "\t$del_pos" . "\t." . "\t$del_REF" . "\t$del_ALT" . "\t." . "\tPASS" . "\t." . "\tGT" . "\t1";
				}
				$i+=1;
				next SEQ_SCAN;
				
			}
			my $del_pos = $i - $del_size;
			my $del_REF = substr($ref_seq,($del_pos - 1),($del_size + 1));
			my $del_ALT = substr($seq,($del_pos - 1),1);
			
			if ( $del_pos eq '0' ) {
				$del_pos = 1;
				$del_REF = substr($ref_seq,($del_pos - 1),$del_size);
				$del_ALT = substr($seq,($del_pos - 1),1);
			}			
			
			say $out "FimH" . "\t$del_pos" . "\t." . "\t$del_REF" . "\t$del_ALT" . "\t." . "\tPASS" . "\t." . "\tGT" . "\t1";
			
			$del_size = 0;
			$del = 0;

		}
		
		if ( $in ) {
			if ( $ref_pos eq "-" && $seq_pos ne "-" ) {
				$in_size+=1;
				
				
				if ( $i eq $len_al-2 ) {
					my $in_pos = $i - $in_size + 1;
					my $in_ALT = substr($seq,($in_pos - 1),($in_size + 1));
					my $in_REF = substr($ref_seq,($in_pos - 1),1);
					say $out "FimH" . "\t$in_pos" . "\t." . "\t$in_REF" . "\t$in_ALT" . "\t." . "\tPASS" . "\t." . "\tGT" . "\t1";

				}
				$i+=1;
				next SEQ_SCAN;
			}
			my $in_pos = $i - $in_size;
			my $in_ALT = substr($seq,($in_pos - 1),($in_size + 1));
			my $in_REF = substr($ref_seq,($in_pos - 1),1);
			
			if ( $in_pos eq '0' ) {
				$in_pos = 1;
				$in_ALT = substr($seq,($in_pos - 1),$in_size);
				$in_REF = substr($ref_seq,($in_pos - 1),1);
			}
			
			say $out "FimH" . "\t$in_pos" . "\t." . "\t$in_REF" . "\t$in_ALT" . "\t." . "\tPASS" . "\t." . "\tGT" . "\t1";
		
			$in_size = 0;
			$in = 0;

		}		
		
		unless ( $ref_pos eq $seq_pos ) {
			if ( $seq_pos eq "-" ) {
				$del_size+=1;
				$del = 1;
				$i+=1;

				next SEQ_SCAN;
			}
			
			if ( $ref_pos eq "-" ) {
				$in_size+=1;
				$in = 1;
				$i+=1;
				next SEQ_SCAN;
			}
			my $var_pos = $i + 1;
			say $out "FimH" . "\t$var_pos" . "\t." . "\t$ref_pos" . "\t$seq_pos" . "\t." . "\tPASS" . "\t." . "\tGT" . "\t1";
			
		}
		
		$i+=1;
	}
		
	close $out;
}



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
