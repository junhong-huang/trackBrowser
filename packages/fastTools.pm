package fastTools;

BEGIN {
    $VERSION = '0.1';
}

use strict;
use warnings;

use base 'Exporter';
our @EXPORT = qw(unzipFile runFastQC clipRead runFlash);

sub unzipFile {

    # options are "-d -c "
    # program is gzip
    my ( $program, $options, $seqFile, $outFile ) = @_;
    if ( $program eq "" ) {
        $program = 'gzip';
    }
    my $commandLine
        = $program . " " . $options . " " . $seqFile . " > " . $outFile;
    return $commandLine;
}

sub runFastQC {

    # options are -o xxx -f xxx
    my ( $program, $options, $seqFiles ) = @_;
    if ( $program eq "" ) {
        $program = '/home/jianhua/softwares/FastQC/fastqc';
    }
    my $commandLine = $program . " " . $options . " " . $seqFiles;
    return $commandLine;
}

sub runClipRead {

    # options are -l 20 -e 0.2 -o $outdir -f $prefix -a $adapter1 -b $adapter2
    # $seqFiles are -1 $read1File -2 $read2File
    # logFile
    my ( $program, $options, $seqFiles ) = @_;
    if ( $program eq "" ) {
        $program = '/home/jianhua/bioTools/bin/clipRead';
    }
    my $commandLine = $program . " " . $options . " " . $seqFiles;
    return $commandLine;
}

sub runFlash {

    # options are -m 10 -M 120 -x 0.25 -o $prefix
    # $seqFiles are $read1File $read2File
    # logFile
    my ( $program, $options, $seqFiles ) = @_;
    if ( $program eq "" ) {
        $program = '/home/jianhua/softwares/FLASH/flash';
    }
    my $commandLine = $program . " " . $options . " " . $seqFiles;
    return $commandLine;
}

sub runCollapser {

    # options are -Q 33
    # $seqFiles are  -i ". $filename." -o ". $outfile
    # logFile
    my ( $program, $options, $seqFile, $outFile ) = @_;
    if ( $program eq "" ) {
        $program = '/home/jianhua/softwares/fastx/fastx_collapser';
    }
    my $commandLine
        = $program . " " . $options . " -i " . $seqFile . " -o " . $outFile;
    return $commandLine;
}

1;
__END__
