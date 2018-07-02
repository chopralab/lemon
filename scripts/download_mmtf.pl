#!/usr/bin/env perl

use warnings;
use strict;

use LWP::Simple;

my $url = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx';

my $entries = get $url;

open my $pdbs, '<', \$entries or die "$!\n";

# skip the first line
<$pdbs>;
<$pdbs>;

while ( my $line = <$pdbs> ) {

    my $pdb = substr( $line, 0, 4);
    next if -s "$pdb.mmtf.gz";

    print "$pdb\n";

    `wget http://mmtf.rcsb.org/v1.0/full/$pdb.mmtf.gz`;
}
