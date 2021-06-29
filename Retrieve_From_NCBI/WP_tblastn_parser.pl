#!/usr/bin/env perl

use Modern::Perl '2011';
use Smart::Comments '###';
use Bio::DB::EUtilities;
use File::Basename;
use Path::Class 'file';
use Bio::SeqIO;
use List::AllUtils 'mesh';

unless (@ARGV == 2) {
	die <<"EOT";
Usage: $0 <report.fmt6> <sequences.fasta>

	This tool takes the report (outfmt 6/7) of a tblastn of NCBI RefSeq 
protein (accession numbers in WP_xxx) against the genomes in which the 
WP sequences are contained as a DB. The informations needed to run the
tblastn could be obtained using NCBI_FromWP-ProtAccToGenomesNucSeq.pl.
The second argument is the fasta file used as query for the tblastn.
	The fasta file is used to build a hash with ids as keys and seq
length as values. Queries and DB sequences should be included. The tblastn 
report is then parsed to keep only hits having 100% identity and an alignment 
length equal to the protein length. 
The coordinates of the hit are then used to retrieve the nucleotide sequence 
corresponding to the queried protein. All gene sequences are outputed in a fasta 
file named Retrieved_seqs.fasta.
	A second file is outputed which contains the correspondance between
the protein WP accession numbers and the related geneomes in which the
protein is comprised. This file is named based on the blast report name
and location: report_Ids_correspondance.txt.

#WP_protein_id	Genome_nuc_id	Start	End	Strand
WP_000832213.1	MOKN01000043.1	22545	21646	2
WP_000832213.1	NXMU01000046.1	5822	6721	1

EOT
}


my $report = shift;
my $parser = get_parser($report);
my $fasta = shift;
my $file = 'Retrieved_WP_nuc_seqs.fasta';


my $seqio = Bio::SeqIO->new(-file => $fasta);
my %len_for;
while( my $seq = $seqio->next_seq ) {
	my $id = $seq->id;
	my $len = $seq->length;
	$len_for{ $id } = $len;
 }

open(FILE,">$file") or die "cannot open $file....";

my $fetcher = Bio::DB::EUtilities->new(-eutil   => 'efetch',
                                       -db      => 'nucleotide',
                                       -email   => 'stanislas.thiriet-rupert@pasteur.fr',
                                       -rettype => 'fasta');
                                       
#Creating output file for ids correspondance based on input path
my ($basename, $dir, $suffix) = fileparse($report, qr{\.[^.]*}xms);
my $out_ids_corr = file($dir, $basename . '_Ids_correspondance.txt');

open my $out_corresp, '>', $out_ids_corr;
say {$out_corresp} "#" . join "\t", qw(WP_protein_id Genome_nuc_id Start End Strand);

#Parsing the blast report
HSP:
while ( my $hsp_ref = $parser->() ) {
	my ($qid, $hid, $identity, $al_len, $start, $end) = @{ $hsp_ref }{ qw(query_id hit_id identity length hit_start hit_end) };

	next HSP unless $identity == 100.00 && $al_len == $len_for{ $qid };		# skip unwanted hits
	if ( $al_len == $len_for{ $hid } ) {
		### Fetching nucleotide id: $hid
		my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
                                       -db      => 'nucleotide',
                                       -rettype => 'gb',
                                       -email   => 'stanislas.thiriet-rupert@pasteur.fr',
                                       -id      => $hid);


		my $file = 'WP_tblastn_GeneBank.gb';

		# dump <HTTP::Response> content to a file (not retained in memory)
		$factory->get_Response(-file => $file);

		#Converte the gbk file to fasta
		# create one SeqIO object to read in,and another to write out
		my $seq_in = Bio::SeqIO->new( -file   => $file,
                              -format => 'genbank',
                            );

		# write each entry in the input file to the output file
		# Scanning for corresponding CDS...
		SCAN:
		while (my $inseq = $seq_in->next_seq) {
			FEATURES:
			for my $feat_object ($inseq->get_SeqFeatures) {
				if ($feat_object->primary_tag eq "CDS") {
					my $cds_start = $feat_object->location->start;
					my $cds_stop = $feat_object->location->end;

					if ( $start > $end ) {
						next FEATURES unless ( $cds_start == $end-3 && $cds_stop == $start );
						my $strand = 2;
						### Retrieving nucleotide sequence...
						$fetcher->set_parameters(-id => $hid,
								-seq_start => $end,	#$end-3 if need to account for the stop codon
								-seq_stop  => $start,
								-strand    => $strand);

						print FILE $fetcher->get_Response->content;
						say {$out_corresp} join "\t", $qid, $hid, $start, $end, $strand;	#$end-3 if need to account for the stop codon
						next HSP;
					}
					my $strand = 1;
					next FEATURES unless ( $cds_start == $start && $cds_stop == $end+3 );
					### Retrieving nucleotide sequence...
					$fetcher->set_parameters(-id => $hid,
						-seq_start => $start,
						-seq_stop  => $end,  #$end+3 if need to account for the stop codon
						-strand    => $strand);

					print FILE $fetcher->get_Response->content;
					say {$out_corresp} join "\t", $qid, $hid, $start, $end, $strand;	#$end+3 if need to account for the stop codon
				}
			}
		}
	}
}

close FILE;




sub get_parser {
	my $infile = shift;
	
	open my $fh, '<', $infile;
	
	# define HSP hash keys
	# slightly different from column names
	my @ attrs = qw(
		query_id hit_id
		identity length mismatches gaps
		query_start query_end
		hit_start hit_end
		evalue bits
	);
	
	# setup parser
	my $closure = sub {
		
		LINE:
		while (my $line = <$fh>) {
			chomp $line;

			next LINE if $line =~ m/\A \s* \z/xms;
			next LINE if $line =~ m/\A \#/xms;
			
			# process HSP line
			my @fields = split /\t/xms, $line;
			
			# Fields for BLAST -outfmt 6 or 7
			#	 0. query id
			#	 1. subject id
			#	 2. % identity
			#	 3. alignment length
			#	 4. mismatches
			#	 5. gap opens
			#	 6. q. start
			#	 7. q. end
			#	 8. s. start
			#	 9. s. end
			#	10. evalue
			#	11. bit score
			
			# here, we may add or modify some fields if needed
			
			# build and return anonumous HSP hash ref
			return { mesh @attrs, @fields };
		}
		
		return;			# no more line to read
	};
	
	return $closure;
}
