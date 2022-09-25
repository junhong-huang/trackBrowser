package samTools;

BEGIN {
    $VERSION = '0.1';
}

use strict;
use warnings;
use commonTools;

use base 'Exporter';
our @EXPORT = qw(samToBam sortBam indexBam sortBamByName);

sub samToBam {
    my ( $program, $alignFile, $bamOutFile ) = @_;
    if ( $program eq "" ) {
        $program = '/public/home/huangjh/soft/miniconda3/envs/RNAseq/bin/samtools';
    }
    my $options = ' view -bS -@ 16 ';
    my $commandLine
        = $program . " " . $options . " " . $alignFile . " -o " . $bamOutFile;
    &executeCommand($commandLine);
    return $bamOutFile;
}

sub sortBam {
    my ( $program, $bamFile, $sortBamFile ) = @_;
    if ( $program eq "" ) {
        $program = '/public/home/huangjh/soft/miniconda3/envs/RNAseq/bin/samtools';
    }
    my $options = " sort  -@ 16 ";
    my $commandLine
        = $program . " "
        . $options . " "
        . $bamFile . " -o "
        . $sortBamFile;
    &executeCommand($commandLine);
    return $sortBamFile;
}

sub indexBam{
    my ( $program, $sortBamFile ) = @_;
    if ( $program eq "" ) {
        $program = '/public/home/huangjh/soft/miniconda3/envs/RNAseq/bin/samtools';
    }
    my $options = " index -@ 16 ";
    my $commandLine
        = $program . " "
        . $options . " "
        . $sortBamFile;
    &executeCommand($commandLine);
    return $sortBamFile;    
}

sub sortBamByName {
    my ( $program, $bamFile, $sortBamFile ) = @_;
    if ( $program eq "" ) {
        $program = '/public/home/huangjh/soft/miniconda3/envs/RNAseq/bin/samtools';
    }
    my $options = " sort -n -@ 16 ";
    my $commandLine
        = $program . " "
        . $options . " "
        . $bamFile . " -o "
        . $sortBamFile;
    &executeCommand($commandLine);    
    return $sortBamFile;
}


1;
__END__
