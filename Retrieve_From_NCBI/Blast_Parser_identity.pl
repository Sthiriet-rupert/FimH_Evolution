#!/usr/bin/env perl

use Modern::Perl '2011';
use Smart::Comments;
use File::Basename;
use Path::Class 'file';
use List::AllUtils 'mesh';

unless (@ARGV == 2) {
	die <<"EOT";
Usage: $0 <Blast_Report.fmt6> <identity_threshold> <evalue_threshold>

	This tool takes two arguments. The blast report in fmt6 and the identity
threshold to apply to all hits in the blast report.

EOT
}


my $report = shift;
my $parser = get_parser($report);

my $identity_threshold = shift;
my $evalue_threshold = shift;
my @fields = qw(
		query_id hit_id
		identity length mismatches gaps
		query_start query_end
		hit_start hit_end
		evalue bits
	);

#Creating output file based on blast report path
my ($basename, $dir, $suffix) = fileparse($report, qr{\.[^.]*}xms);
my $outfile = file($dir, $basename . '_Parsed(id' . "$identity_threshold" . ').fmt6');

open my $out, '>', $outfile;

#Parsing the blast report
HSP:
while ( my $hsp_ref = $parser->() ) {
	my ($identity, $evalue) = @{ $hsp_ref }{ qw(identity evalue) };
	
	next HSP if $identity < $identity_threshold;		# skip "weak" hits
	next HSP if $evalue > $evalue_threshold;
	
	say {$out} join "\t", @{ $hsp_ref }{ @fields };
}



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
