#!/usr/bin/env perl

use Modern::Perl '2011';
use Smart::Comments;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use Template;
use IPC::System::Simple qw(system);


unless (@ARGV == 1) {
	die <<"EOT";
Usage: $0 <Prot_accs.txt>

	This tool takes as argument a text file with all protein accession numbers 
for which you want to get the genebank data.
Only one accession number per line:

WP_087711612.1
QBA19949.1

	Each accession number is used to retrieve the related protein sequence as well
as some features that are put in a file called "Correspondance_Accs_Features.txt". 
This file will store the following informations:

#protAcc	AltProtAcc	NucAcc		AssemblyAcc	Species					Strain		Serotype	Isolation_source	Host
WP_087711612.1	OVE55035.1	MVAG01000142.1	GCA_002177115.1	Chryseobacterium mucoviscidosis		VT16-26		na		sputum		Homo sapiens
WP_087711612.1	REC50469.1	QNVU01000015.1	GCA_003385535.1	Candidatus Chryseobacterium massiliae	CCUG 51329	na		nasal swab	Homo sapiens
QBA19949.1	QBA19949.1	CP035532.1	GCA_004137745.1	Chryseobacterium indologenes		StR 01		na		ponds		Nile tilapia

The protein sequences are witten in a fasta file called "Proteins_DATE.fasta"
where the DATE as the following format :

Proteins_Wed_18_Nov_2020_12-29-53.fasta

For RefSeq ids (WP_xxx), all realted proteins and informations are downloaded.
For classical ids, only the corresponding sequence and informations are downloaded.
EOT
}

my $file = shift;

open my $infile, '<', $file;
my @ids;

# Store protein accessions in an array
LINE:
while (my $line = <$infile>) {
	chomp $line;
	next LINE if $line =~ m/\A \s* \z/xms;
	next LINE if $line =~ m/\A \#/xms;
	
	push @ids, $line;
}
close $infile;

# Get the date
my ($wday, $mon, $day, $hrs, $year) = split " ", localtime();
my ($hr, $min, $sec) = split ':', $hrs;
my $date = join "_", $wday, $day, $mon, $year, join "-", $hr, $min, $sec;

# Create a a fasta file to store all protein sequences with the date to be sure that no error will be raised by an already existing directory
my $outfasta = "Proteins_$date.fasta";
system "touch Proteins_$date.fasta";


# Create output file for features
my $outfile = 'Correspondance_Accs_Features.txt';
open my $out, '>', $outfile;
say {$out} '#' . join "\t", qw(protAcc AltProtAcc NucAcc AssemblyAcc Species Strain Serotype Isolation_source Host);

open my $outf, '>', $outfasta;

# Create a counting hash to check if a protAcc or assembly accession was already encountered to avoid identical sequences from the same genome
my %Id_count;
my %Assembly_count;

# Loop over chunks of id to get features for all corresponding sequence (or sequenceS in case of RefSeq ids) and write them in the ouput file
my $chunk = natatime 10, @accs;

while ( my @values = $chunk->() ) {
	for my $id ( @values ) {
		my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
										   -db      => 'protein',
										   -rettype => 'ipg',
										   -email   => 'stanislas.thiriet-rupert@pasteur.fr',
										   -id      => $id);
		my $ipg_file = 'test.ipg';

		# dump <HTTP::Response> content to a file (not retained in memory)
		$factory->get_Response(-file => $ipg_file);

		open my $in, '<', $ipg_file;
		my $i = 1;
		IPG:
		while (my $line = <$in>) {
			chomp $line;
			
			# Do not read the header line
			next IPG if $line =~ m/\A Id /xms;
			
			
			my @infos = split "\t", $line;
			my $source = $infos[1];
			
			my $nucAcc = $infos[2];
			my $protAcc = $infos[6];
			my $sp = $infos[8];
			my $strain = $infos[9];
			my $AssemblyAcc = pop @infos;
			

			# check if the protein id was already encountered to avoid duplicate sequences
			next IPG if ( $Id_count{ $protAcc } && Assembly_count{ $AssemblyAcc } );
			
			$Id_count{ $protAcc }++;
			Assembly_count{ $AssemblyAcc }++;
			
			# download a genebank file for the protein to get informations from
			my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
											   -db      => 'protein',
											   -rettype => 'gb',
											   -email   => 'stanislas.thiriet-rupert@pasteur.fr',
											   -id      => $protAcc);
			my $gb_file = 'test.gb';
			$factory->get_Response(-file => $gb_file);


			# create one SeqIO object to read in the gbk file
			my $seq_in = Bio::SeqIO->new( -file   => $gb_file,
									-format => 'genbank',
									);
									
			
			### Writing informations for protein: $id
			# write the sequence in the output fasta file and the features in the text file
			my @tags = qw(serotype isolation_source host);
			FEATURES:
			while (my $inseq = $seq_in->next_seq) {
				say {$outf} ">$id" . "|$protAcc" . "|$sp" . "|$strain";
				say {$outf} $inseq->seq;
				
				for my $feat_object ($inseq->get_SeqFeatures) {

					if ( $feat_object->primary_tag eq 'source' ) {
						print {$out} $id . "\t" . $protAcc . "\t" . $nucAcc . "\t" . $AssemblyAcc . "\t" . $sp . "\t" . $strain;
						TAG:
						for my $tag ( @tags ) {
							if ( $feat_object->has_tag($tag) ) {
								for my $value ($feat_object->get_tag_values($tag)) {

									print {$out} "\t" . $value;
									next TAG;
								}
							}
							print {$out} "\tna";
						}
					print {$out} "\n";
					last FEATURES;
					}
				}
				
			}
			
			
		}
	}
}
