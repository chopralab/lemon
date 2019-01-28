use strict;
use LWP::Simple qw( $ua );

my $request_file = shift;

if (!$request_file) {
    print STDERR "You must supply a request file.\n";
    exit 1;
}

open my $fh, '<', $request_file or die;
$/ = undef;
my $XML_query = <$fh>;
close $fh;

# you can configure a proxy...                                                                          
#$ua->proxy( http => 'http://yourproxy:8080' );

# Create a request                                                                                  

my $request = HTTP::Request->new( POST => 'http://www.rcsb.org/pdb/rest/search/');
$request->content_type( 'application/x-www-form-urlencoded' );
$request->content( $XML_query );

my $response = $ua->request( $request );

# Check to see if there is an error
unless( $response->is_success ) {
    print STDERR "An error occurred: ", $response->status_line, "\n";
    print STDERR "Please ensure output is correct\n";
}

# Print response content in either case
print $response->content;
