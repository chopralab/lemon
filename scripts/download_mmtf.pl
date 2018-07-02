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

    my $first = substr( $pdb, 0, 1);
    my $second= substr( $pdb, 1, 1);

    next if -s "$first/$second/$pdb.mmtf.gz";

    print "$pdb\n";

    `wget -q http://mmtf.rcsb.org/v1.0/full/$pdb.mmtf.gz -o $first/$second/$pdb.mmtf.gz`;
}
