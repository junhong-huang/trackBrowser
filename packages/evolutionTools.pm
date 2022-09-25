package evolutionTools;

BEGIN {
    $VERSION = '0.1';
}

use strict;
use warnings;
use commonTools;

use base 'Exporter';
our @EXPORT = qw(getEvoConVals);

sub getEvoConVals {    # get evolution conservation value
    my ( $program, $options, $conservationProfile, $targetFile ) = @_;
    if ( $program eq "" ) {
        $program = '/public/home/huangjh/bioTools/bin/bedConservation';
    }
    my %targetHash = ();
    my $annoLine = "";
    my $targetBed  = $targetFile . ".evolution.bed";
    open( TARGET, "<$targetFile" ) || die "can't open the $targetFile\n";
    open( OUTPUT, ">$targetBed" )
        || die "can't open the $targetBed\n";
    while ( my $line = <TARGET> ) {
        $line =~ s/\s+$//;
        if ($line =~ /\#/)
        {
            $annoLine = $line;
            next;
        }
        my @items = split /\t/, $line;
        $targetHash{ $items[3] } = $line;
        print OUTPUT $items[0], "\t", $items[1], "\t", $items[2], "\t",
            $items[3], "\t", $items[4], "\t", $items[5], "\n";
    }
    close(TARGET);
    close(OUTPUT);

    my %conservationHash = ();
    my $conservatedFile  = $targetFile . ".conservation";

    my $commandLine
        = $program . " "
        . $options . " " . " -o "
        . $conservatedFile . " "
        . $conservationProfile . " "
        . $targetBed;
    &executeCommand($commandLine);

    open( CONS, "<$conservatedFile" )
        || die "can't open $conservatedFile\n";
    while ( my $line = <CONS> ) {
        $line =~ s/\s+$//;
        my @items = split /\t/, $line;
        $conservationHash{ $items[3] } = $items[6];
    }
    close(CONS);

    my $outFile = $targetFile . ".bedConservation";
    open( OUTPUT, ">$outFile" )
        || die "can't open the $outFile\n";
    if (length($annoLine)>0)
    {
        print OUTPUT $annoLine, "\t", "conservation", "\n"; 
    }
    foreach my $name ( sort keys %targetHash ) {
        my $con = $conservationHash{$name};
        print OUTPUT $targetHash{$name}, "\t", $con, "\n";
    }
    close(OUTPUT);
    &executeCommand("rm -f $targetBed $conservatedFile");
    return $outFile;
}


1;
__END__
