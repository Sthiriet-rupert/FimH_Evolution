#!/usr/bin/env perl

use Modern::Perl '2011';
use Smart::Comments '###';
use Bio::DB::EUtilities;
use List::MoreUtils qw( natatime );
use File::Basename;
use Path::Class 'file';
use Bio::SeqIO;

unless (@ARGV == 1) {
	die <<"EOT";
Usage: $0 <protAccs.txt> >summary.txt

	This tool takes NCBI protein accession numbers and outputs the related
nucleotide sequences in two separate fasta files. One for sequences contained
in a genome (retrieved-seqs-fromGenome.fasta) and one for sequences having
its own nucleotide accession number.
By redirecting the standard output in a summary.txt file, you can get informations
such as the different types of identifiers and their number at the end of the
proccess.
	A third file named based on the infile name 'infile_Ids_correspondance.txt'
is ouputed with the correspondance between the different types of ids:

#Prot_Accs Prot_uid Nuc_uid Nuc_acc Len_nuc_acc
ARE49919.1	1173487020	1173482665	CP020543	903
ARO85821.1	1189921537	1189921536	KY007009	429
EOT
}


#Get the prot Acc in an array
my $infile = shift;
my @accs;

open my $in, '<', $infile;
#Put accs in an array
while (my $line = <$in>) {
	chomp $line;
	push @accs, $line;
}
close $in;
my @prot_uids;

### Performing esearch to obtain prot UIDs
##Get correspondance beetwed prot Acc and prot UID, using chunks of 20 Accs
my $chunk = natatime 20, @accs;

while (my @values = $chunk->()) {

	
	my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch',
										   -email => 'stanislas.thiriet-rupert@pasteur.fr',
										   -db    => 'protein',
										   -term  => join(',',@values) );

	push @prot_uids, $factory->get_ids;
}


### Performing esummary to obtain prot summaries
#Get summary for each protein UID
my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary',
									   -email => 'stanislas.thiriet-rupert@pasteur.fr',
									   -db    => 'protein',
									   -id    => \@prot_uids);

while (my $ds = $factory->next_DocSum) {
    print "ID: ",$ds->get_id,"\n";

    # flattened mode
    while (my $item = $ds->next_Item('flattened'))  {
        # not all Items have content, so need to check...
        printf("%-20s:%s\n",$item->get_name,$item->get_content)
          if $item->get_content;
    }
    print "\n";
}

#Put prot acc and related prot uid in a hash
my %prot_uid_for_acc;
@prot_uid_for_acc{ @accs } = @prot_uids;

### Performing elink to obtain nuc UIDs
my @nuc_uids;
my %nuc_for;
for my $single_id ( @prot_uids ) {

	$factory->reset_parameters(-eutil 		   => 'elink',
                           -email          => 'stanislas.thiriet-rupert@pasteur.fr',
                           -db             => 'nuccore',
                           -dbfrom         => 'protein',

                           -id             => $single_id);




# iterate through the LinkSet objects
	while (my $ds = $factory->next_LinkSet) {
		my $link = $ds->get_link_name;

		if ( $link eq 'protein_nuccore' ) {
			my $nuc_uid = join(',',$ds->get_ids);
			my $prot = join(',',$ds->get_submitted_ids);

			#Put prot uid and related gene uid in a hash
			$nuc_for{ $prot } = $nuc_uid;
			push @nuc_uids, $nuc_uid;
    
		}
	}

}


my $outfile = 'NonWP_retrieved_nuc_seqs.fasta';

my $seq_out = Bio::SeqIO->new( -file   => ">$outfile",
                               -format => 'fasta',
                             );

my $outFromGenome = 'NonWP_retrieved_nuc_seqs_fromGenome.fasta';
open my $out, '>', $outFromGenome;


my $nuc_chunk = natatime 20, @nuc_uids;

my %len_seq_for;
my @nuc_accs;
my %nuc_acc_for_uid;
while (my @values = $nuc_chunk->()) {

	### Performing efetch to get nuc sequences
	$factory->reset_parameters(-eutil => 'efetch',
                                       -db      => 'nucleotide',
                                       -rettype => 'gb',
                                       -email   => 'stanislas.thiriet-rupert@pasteur.fr',
                                       -id      => \@values);


	my $file = 'GeneBank.gb';

	# dump <HTTP::Response> content to a file (not retained in memory)
	$factory->get_Response(-file => $file);

	#Converte the gbk file to fasta
	# create one SeqIO object to read in,and another to write out
	my $seq_in = Bio::SeqIO->new( -file   => $file,
                              -format => 'genbank',
                            );




	# write each entry in the input file to the output file
	SCAN:
	while (my $inseq = $seq_in->next_seq) {
		my $nuc_acc = $inseq->id;
		my $len = length $inseq->seq;
		next SCAN unless $len;
		if ( $len > 1000 ) {

			FEATURES:
			for my $feat_object ($inseq->get_SeqFeatures) {
				if ($feat_object->primary_tag eq "CDS") {
				
					my @protein_id = $feat_object->get_tag_values("protein_id") 
					if ($feat_object->has_tag("protein_id"));

					my $id = shift @protein_id;
				
					next FEATURES unless $id;
					next FEATURES unless $prot_uid_for_acc{ $id };

					my $sequence = $feat_object->spliced_seq->seq;
					$len = length $sequence;
					last FEATURES unless $len;
					say {$out} ">nuc_" . $nuc_acc . "_for_" . $id;
					say {$out} $sequence;
					$len_seq_for{ $nuc_acc } = $len;
					$nuc_acc_for_uid{ $id } = $nuc_acc;
					push @nuc_accs, $nuc_acc;
					next SCAN;
				}
			}
		}

		for my $feat_object ($inseq->get_SeqFeatures) {
			if ($feat_object->primary_tag eq "CDS") {
				
				my @protein_id = $feat_object->get_tag_values("protein_id") 
				if ($feat_object->has_tag("protein_id"));

				my $id = shift @protein_id;
		
				
				if ( $inseq->seq ) {
					$seq_out->write_seq($inseq);
					$len_seq_for{ $nuc_acc } = $len;
					$nuc_acc_for_uid{ $id } = $nuc_acc;
					push @nuc_accs, $nuc_acc;
				}
			}	
		}
	}
}


		
#### @accs
#### @prot_uids
#### @nuc_uids
#### @nuc_accs

my $Nb_prot_accs = @accs;
my $Nb_prot_uids = @prot_uids;
my $Nb_nuc_uids = @nuc_uids;
my $Nb_nuc_accs = @nuc_accs;
print "Nb Protein accs: ", $Nb_prot_accs,"\n";
print "Nb Protein UIDs: ", $Nb_prot_uids,"\n";
print "Nb Nucleotide UIDs: ", $Nb_nuc_uids,"\n";
print "Nb Nucleotide accs: ", $Nb_nuc_accs,"\n";


#Creating output file for ids correspondance based on input path
my ($basename, $dir, $suffix) = fileparse($infile, qr{\.[^.]*}xms);
my $out_ids_corr = file($dir, $basename . '_Ids_correspondance.txt');
	
#Outputing ids correspondace between prot acc, prot uid and gene uid
open my $out_corresp, '>', $out_ids_corr;
say {$out_corresp} join "\t", qw(Prot_Accs Prot_uid Nuc_uid Nuc_acc Len_nuc_acc);
for my $acc ( @accs ) {
	say {$out_corresp} join "\t", $acc, $prot_uid_for_acc{ $acc }, ( $nuc_for{ $prot_uid_for_acc{ $acc } } // 'na' ), ( $nuc_acc_for_uid{ $acc } // 'na' ), ( $len_seq_for{ $nuc_acc_for_uid{ $acc } } // 'na' );
}
