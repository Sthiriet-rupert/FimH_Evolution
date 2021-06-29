#!/usr/bin/env perl

use Modern::Perl '2011';
use Smart::Comments;
use Bio::DB::EUtilities;
use List::MoreUtils qw( natatime );
use File::Basename;
use Path::Class 'file';
use Bio::SeqIO;

unless (@ARGV == 1) {
	die <<"EOT";
Usage: $0 <protAccs.txt>

	This tool takes NCBI RefSeq protein accession numbers (WP_xxx) and 
outputes a fasta file containing the genomes the WP sequences are comprised 
in. This file could then be used as a databank for a tblastn search using 
the WP prot sequences as queries. The corresponding nuclotide sequences 
could then be extracted from the genomes using WP_tblastn_parser.pl.
	A second file named based on the infile name 'infile_Ids_correspondance.txt'
is ouputed with the correspondance between the different types of ids:

#Prot_Accs Prot_uid Nuc_uid Nuc_acc
WP_126715566.1	1546643859	1546686536	NZ_PNVT01000059
WP_126721538.1	1546682453	1546682596	NZ_PNVP01000059
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

my @prot_uids;


##Get correspondance beetwed prot Acc and prot UID, using chunks of 20 Accs
my $chunk = natatime 20, @accs;

while (my @values = $chunk->()) {

	
	my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch',
										   -email => 'stanislas.thiriet-rupert@pasteur.fr',
										   -db    => 'protein',
										   -term  => join(',',@values) );

	push @prot_uids, $factory->get_ids;
}



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



##Get correspondance between prot UID and nucl UID

$factory->reset_parameters(-eutil 		   => 'elink',
                           -email          => 'stanislas.thiriet-rupert@pasteur.fr',
                           -db             => 'nucleotide',
                           -dbfrom         => 'protein',
                           -correspondence => 1,
                           -id             => \@prot_uids);

my @nuc_uids;
my %nuc_for;

# iterate through the LinkSet objects
while (my $ds = $factory->next_LinkSet) {
	my $link = $ds->get_link_name;
    print "   Link name: ",$link,"\n";
    
    if ( $link eq 'protein_nuccore' ) {
		my $nuc_uid = join(',',$ds->get_ids);
		my $prot = join(',',$ds->get_submitted_ids);
		
		#Put prot uid and related gene uid in a hash
		$nuc_for{ $prot } = $nuc_uid;
		push @nuc_uids, $nuc_uid;
		print "Protein IDs: ",join(',',$ds->get_submitted_ids),"\n";
		print "    Nuc IDs: ",join(',',$ds->get_ids),"\n";
    
	}

}



$factory->reset_parameters(-eutil => 'efetch',
                                       -db      => 'nucleotide',
                                       -rettype => 'gb',
                                       -email   => 'stanislas.thiriet-rupert@pasteur.fr',
                                       -id      => \@nuc_uids);


my $file = 'WP_GeneBank.gb';


# dump <HTTP::Response> content to a file (not retained in memory)
$factory->get_Response(-file => $file);

#Converte the gbk file to fasta
# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new( -file   => $file,
                              -format => 'genbank',
                            );

my @nuc_accs;

# Put nuc accs in an array
while (my $inseq = $seq_in->next_seq) {
	my $nuc_acc = $inseq->id;
	
	say "The nucleotide accession number is:" . $nuc_acc;

	FEATURES:
	for my $feat_object ($inseq->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "CDS") {
			
			my @protein_id = $feat_object->get_tag_values("protein_id") 
			if ($feat_object->has_tag("protein_id"));

			my $id = shift @protein_id;
			
			next FEATURES unless $id;
			next FEATURES unless $prot_uid_for_acc{ $id };
			
			push @nuc_accs, $nuc_acc;
			last FEATURES;
		}
	}
}


my @WP_genomes;
for my $nuc_acc ( @nuc_accs ) {
	my (undef, $genome) = split "_", $nuc_acc;
	push @WP_genomes, $genome;
}

# Get genome or contig sequences to perform subsequent tblastn
$factory->reset_parameters(-eutil => 'efetch',
                                       -db      => 'nucleotide',
                                       -rettype => 'fasta',
                                       -email   => 'stanislas.thiriet-rupert@pasteur.fr',
                                       -id      => \@WP_genomes);

my $WP_genome = 'WP_genomes.fasta';

# dump <HTTP::Response> content to a file (not retained in memory)
$factory->get_Response(-file => $WP_genome);







		



### @accs
### @prot_uids
### @nuc_uids
### @nuc_accs
my $Nb_prot_accs = @accs;
my $Nb_prot_uids = @prot_uids;
my $Nb_nuc_uids = @nuc_uids;
my $Nb_nuc_accs = @nuc_accs;
print "Nb Protein accs: ", $Nb_prot_accs,"\n";
print "Nb Protein UIDs: ", $Nb_prot_uids,"\n";
print "Nb Nucleotide UIDs: ", $Nb_nuc_uids,"\n";
print "Nb Nucleotide accs: ", $Nb_nuc_accs,"\n";


#Put nuc uids and related accs in a hash
my %nuc_acc_for_uid;
@nuc_acc_for_uid{ @nuc_uids } = @nuc_accs;

#Creating output file for ids correspondance based on input path
my ($basename, $dir, $suffix) = fileparse($infile, qr{\.[^.]*}xms);
my $out_ids_corr = file($dir, $basename . '_Ids_correspondance.txt');
	
#Outputing ids correspondace between prot acc, prot uid and gene uid
open my $out_corresp, '>', $out_ids_corr;
say {$out_corresp} join "\t", qw(Prot_Accs Prot_uid Nuc_uid Nuc_acc);
for my $acc ( @accs ) {
	say {$out_corresp} join "\t", $acc, $prot_uid_for_acc{ $acc }, ( $nuc_for{ $prot_uid_for_acc{ $acc } } // 'na' ), ( $nuc_acc_for_uid{ $nuc_for{ $prot_uid_for_acc{ $acc } } } // 'na' );
}
