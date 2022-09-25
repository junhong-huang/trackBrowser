package bioUtils;

BEGIN {
    $VERSION = '0.1';
}

use strict;
use warnings;
use commonTools;

use base 'Exporter';
our @EXPORT = qw(getFastaSeq);

sub getFastaSeq
{
    my $seqFile = shift @_;
    my $name = "";
    my %fastaSeq = ();
    open(SEQ, "<$seqFile") || die "can't open the $seqFile\n";
    while (my $line = <SEQ>)
    {
        $line =~ s/\s+$//;
        if ($line =~ />(.*)/)
        {
            $name = $1;
        }
        else {
            $fastaSeq{$name} .= $line;
        }
    }
    close(SEQ);
    return \%fastaSeq;
}

1;
__END__
