#!/usr/bin/perl

#usage = path_to/RetrieveFasta.pl IDFIle FastaLibrary > OUTPUT

if ($#ARGV !=1) {
        print "usage: RetrieveFasta.pl IDFile FastaLibrary > outputfile\n";
        exit;
}

use strict;
use Bio::DB::Fasta;

my $database;
my $fasta_library = $ARGV[1]; #opens second user-supplied variable as the
#library file
my %records;

open IDFILE, "<$ARGV[0]" or die $!; #opens first user-supplied variable as
#the ID file
open OUTPUT, <STDOUT>;

# creates the database of the library, based on the file
$database = Bio::DB::Fasta->new("$fasta_library") or die "Failed to create
Fasta DP object on fasta library\n";

# now, it parses the file with the fasta headers you want to get
while (<IDFILE>) {

      my ($id) = (/^>*(\S+)/);  # capture the id string (without the initial ">")
      my $header = $database->header($id);
      #print "$header\n";
      print ">$header\n", $database->seq( $id ), "\n";
      print OUTPUT  ">$header\n", $database->seq( $id ), "\n";
}

#remove the index file
unlink "$fasta_library.index";

#close filehandles
close IDFILE;
close OUTPUT;
exit;
