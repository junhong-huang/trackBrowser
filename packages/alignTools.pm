package alignTools;

BEGIN {
    $VERSION = '0.1';
}

use strict;
use warnings;
use commonTools;

use base 'Exporter';
our @EXPORT = qw(buildBowtieIndex buildBowtie2Index buildStarIndex alignReadBowtie alignReadBowtie2 alignReadStar);

sub buildBowtieIndex 
{
    my ( $program, $options, $baseIndexFile, $fastaFile ) = @_;
    if ( $program eq "" ) {
        $program = '/home/jianhua/softwares/bowtie1/bowtie-build';
    }
    my $commandLine
        = $program . " " . $options . " " . $fastaFile . " " . $baseIndexFile;
    return $commandLine;
}

sub buildBowtie2Index 
{
    my ( $program, $options, $baseIndexFile, $fastaFile ) = @_;
    if ( $program eq "" ) {
        $program = '/home/jianhua/softwares/bowtie2/bowtie2-build';
    }
    my $commandLine
        = $program . " " . $options . " " . $fastaFile . " " . $baseIndexFile;
    return $commandLine;
}

sub buildStarIndex 
{
    my ( $program, $options, $baseIndexFile, $fastaFile ) = @_;
    if ( $program eq "" ) {
        $program = '/public/apps/SOFTWARE/STAR-2.7.0f/bin/Linux_x86_64/STAR';
    }
    if ( $options eq "" ) {
        $options
            = "--runThreadN 20 --runMode genomeGenerate"
            . " --genomeDir "
            . $baseIndexFile . " "
            . " --genomeFastaFiles "
            . $fastaFile . " ";
    }
    my $commandLine = $program . " " . $options;
    return $commandLine;
}

sub alignReadBowtie 
{
    my ( $program, $options, $baseIndexFile, $fastqFile, $outFile ) = @_;
    if ( $program eq "" ) {
        $program = '/home/jianhua/softwares/bowtie2/bowtie';
    }
    if ( $options eq "" ) {
        $options
            = '-p 16 -t --no-unal -a --best --strata -v 2 -m 20 -l 15 --sam';
    }
    my $commandLine
        = $program . " "
        . $options . " "
        . $fastqFile . " " . " -x "
        . $baseIndexFile
        . " >$outFile";
    return $commandLine;
}

sub alignReadBowtie2 
{
    my ( $program, $options, $baseIndexFile, $fastqFile, $outFile ) = @_;
    if ( $program eq "" ) {
        $program = '/home/jianhua/softwares/bowtie2/bowtie2';
    }
    if ( $options eq "" ) {
        $options
            = '-p 16 -t --no-unal --reorder -k 21 -D 20 -R 3 -N 0 -L 15 -i S,1,0.5 -f --score-min L,-0.6,-0.6 -U';
    }
    my $commandLine
        = $program . " "
        . $options . " "
        . $fastqFile . " " . " -x "
        . $baseIndexFile . " "
        . " >$outFile";
    return $commandLine;
}

sub alignReadStar 
{
    my ( $program, $options, $baseIndexFile, $fastqFile, $outPrefix ) = @_;
    if ( $program eq "" ) {
        $program = '/public/apps/SOFTWARE/STAR-2.7.0f/bin/Linux_x86_64/STAR';
    }
    if ( $options eq "" ) {
        $options
            = "--runMode alignReads --outReadsUnmapped Fasta --outFilterMultimapNmax 100 "
            . " --outSAMtype BAM Unsorted --outSAMattributes All --alignIntronMax 1 "
            . " --outFilterMultimapScoreRange 50 --outFilterScoreMinOverLread 0 "
            . " --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 "
            . " --alignSJDBoverhangMin 1000 --chimJunctionOverhangMin 15 ";
    }
    my $commandLine
        = $program . " "
        . $options . " "
        . " --genomeDir "
        . $baseIndexFile . " "
        . " --readFilesIn "
        . $fastqFile . " "
        . " --outFileNamePrefix "
        . " $outPrefix ";
    return $commandLine;
}

1;
__END__
