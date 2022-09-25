#!/usr/bin/perl -w

BEGIN {
    unshift( @INC, "./packages" );
}

use Getopt::Long;
use strict;
use commonTools;
use bioUtils;
use fastTools;
use alignTools;
use samTools;
use annoTools;
use evolutionTools;
use XML::Simple;
use File::Basename;
use RNA;


my $usage="
    perl plotTrackBrowser.pl \
    -org hg38 \
    -i ./demo/plot_region.bed \
    -upExtend 20 \
    -downExtend 20 \
    -upOffset 0 \
    -downOffset 0 \
    -out . \
    -t1 ./configure/trackFile_human_plusStrand.json \
    -t2 ./configure/trackFile_human_minusStrand.json
";
print "usage:\n$usage\n";

my $upExtend = 0;
my $downExtend = 0;
my $upOffset = 0;
my $downOffset = 0;
my $prefix = "";
my $outDir  = "";
my $inputFile  = "";
my $track1 = "";
my $track2 = "";
my $org = "";#mm10;dm6;hg38;

my $optLong = GetOptions(
    "upOffset=s" => \$upOffset,
    "downOffset=s" => \$downOffset,
	"upExtend=s" => \$upExtend,
	"downExtend=s" => \$downExtend,
	"org=s" => \$org,
	"outDir=s" => \$outDir,
	"i=s" => \$inputFile,
    "t1=s" => \$track1,    # suffix Name of file
    "t2=s" => \$track2,     # program Directory
);


our $tmpFileDir="$outDir/results/tmpfiles";
executeCommand("mkdir -p $tmpFileDir") if (!(-e $tmpFileDir));
my $extendFile=$inputFile;
if( $upExtend > 0 ){
	$extendFile = &upExtend( $extendFile, $upExtend );
}
if( $downExtend > 0 ){
	$extendFile = &downExtend( $extendFile, $downExtend );
}
if( $upOffset > 0 ){
    $extendFile = &upOffset( $extendFile, $upOffset );
}
if( $downOffset > 0 ){
    $extendFile = &downOffset( $extendFile, $downOffset );
}

our $final_prefix=basename($extendFile);
$final_prefix=~s/\.txt//g;
# my $trackSVGdir = "$outDir/$final_prefix"."_trackSVG";
my $trackSVGdir = "$outDir/results";

&getTrackRange($extendFile,$outDir);
&trackBrowser($track1,$track2,$outDir);



sub downOffset {
    my ( $inputFile, $extendLength ) = @_;
    open( IN, "<$inputFile" ) || die "can't open the $inputFile\n";
    my $outPrefix = basename($inputFile);
    my $outFile = $tmpFileDir."/".$outPrefix.".downOffset" . $extendLength . "nt.txt";
    open( OUT, ">$outFile" ) || die "can't open the $outFile\n";

    while ( my $line = <IN> ) {
        $line =~ s/\s+$//;

        if ( $line =~ /\#/ ) {
            print OUT $line,"\n";
            next;
        }
        my @items  = split /\t/, $line;
        my $chrom      = $items[0];
        my $chromStart = $items[1];
        my $chromEnd   = $items[2];
        my $name       = $items[3];
        my $score      = $items[4];
        my $strand     = $items[5];
        my @lineInfo = @items[3..$#items];
        my $lineInfojoin  = join "\t", @lineInfo;
        my $extendStart = "";
        my $extendEnd   = "";
        if ( $strand =~ /\+/){
            $extendStart     = $chromStart+ $extendLength;
            $extendEnd   = $chromEnd + $extendLength;
            print OUT $chrom,"\t", $extendStart,"\t", $extendEnd,"\t",$lineInfojoin,"\n";
        }

        elsif ( $strand =~ /-/){
            $extendStart = $chromStart - $extendLength;
            $extendEnd   = $chromEnd - $extendLength;
            print OUT $chrom,"\t", $extendStart,"\t", $extendEnd,"\t",$lineInfojoin,"\n";
        }

    }

    close(IN);
    close(OUT);
    return $outFile;

}


sub upOffset {
    my ( $inputFile, $extendLength ) = @_;
    open( IN, "<$inputFile" ) || die "can't open the $inputFile\n";
    my $outPrefix = basename($inputFile);
    my $outFile = $tmpFileDir."/".$outPrefix.".upOffset" . $extendLength . "nt.txt";
    open( OUT, ">$outFile" ) || die "can't open the $outFile\n";
    while ( my $line = <IN> ) {
        $line =~ s/\s+$//;

        if ( $line =~ /\#/ ) {
            print OUT $line,"\n";
            next;
        }
        my @items  = split /\t/, $line;
        my $chrom      = $items[0];
        my $chromStart = $items[1];
        my $chromEnd   = $items[2];
        my $name       = $items[3];
        my $score      = $items[4];
        my $strand     = $items[5];
        my @lineInfo = @items[3..$#items];
        my $lineInfojoin  = join "\t", @lineInfo;
        my $extendStart = "";
        my $extendEnd   = "";
        if ( $strand =~ /\+/){
            $extendStart     = $chromStart - $extendLength;
            $extendEnd   = $chromEnd- $extendLength;
            print OUT $chrom,"\t", $extendStart,"\t", $extendEnd,"\t",$lineInfojoin,"\n";
        }

        elsif ( $strand =~ /-/){
            $extendStart = $chromStart + $extendLength;
            $extendEnd   = $chromEnd + $extendLength;
            print OUT $chrom,"\t", $extendStart,"\t", $extendEnd,"\t",$lineInfojoin,"\n";
        }

    }

    close(IN);
    close(OUT);
    return $outFile;

}


sub upExtend {
    my ( $inputFile, $extendLength ) = @_;
    open( IN, "<$inputFile" ) || die "can't open the $inputFile\n";
    my $outPrefix = basename($inputFile);
    my $outFile = $tmpFileDir."/".$outPrefix.".upExtend" . $extendLength . "nt.txt";
    open( OUT, ">$outFile" ) || die "can't open the $outFile\n";
    while ( my $line = <IN> ) {
        $line =~ s/\s+$//;

        if ( $line =~ /\#/ ) {
            print OUT $line,"\n";
            next;
        }
        my @items  = split /\t/, $line;
        my $chrom      = $items[0];
        my $chromStart = $items[1];
        my $chromEnd   = $items[2];
        my $name       = $items[3];
        my $score      = $items[4];
        my $strand     = $items[5];
        my @lineInfo = @items[3..$#items];
        my $lineInfojoin  = join "\t", @lineInfo;
        my $extendStart = "";
        my $extendEnd   = "";

        if ( $strand =~ /\+/){
            $extendStart     = $chromStart - $extendLength;
            $extendEnd   = $chromEnd;
            print OUT $chrom,"\t", $extendStart,"\t", $extendEnd,"\t",$lineInfojoin,"\n";
        }

        elsif ( $strand =~ /-/){
            $extendStart = $chromStart ;
            $extendEnd   = $chromEnd + $extendLength;
            print OUT $chrom,"\t", $extendStart,"\t", $extendEnd,"\t",$lineInfojoin,"\n";
        }

    }

    close(IN);
    close(OUT);
    return $outFile;

}

sub downExtend {
    my ( $inputFile, $extendLength ) = @_;
    open( IN, "<$inputFile" ) || die "can't open the $inputFile\n";
    my $outPrefix = basename($inputFile);
    my $outFile = $tmpFileDir."/".$outPrefix.".downExtend" . $extendLength . "nt.txt";
    open( OUT, ">$outFile" ) || die "can't open the $outFile\n";
    while ( my $line = <IN> ) {
        $line =~ s/\s+$//;

        if ( $line =~ /\#/ ) {
            print OUT $line,"\n";
            next;
        }
        my @items  = split /\t/, $line;
        my $chrom      = $items[0];
        my $chromStart = $items[1];
        my $chromEnd   = $items[2];
        my $name       = $items[3];
        my $score      = $items[4];
        my $strand     = $items[5];
        my @lineInfo = @items[3..$#items];
        my $lineInfojoin  = join "\t", @lineInfo;
        my $extendStart = "";
        my $extendEnd   = "";
        if ( $strand =~ /\+/){
            $extendStart     = $chromStart;
            $extendEnd   = $chromEnd + $extendLength;
            print OUT $chrom,"\t", $extendStart,"\t", $extendEnd,"\t",$lineInfojoin,"\n";
        }

        elsif ( $strand =~ /-/){
            $extendStart = $chromStart - $extendLength;
            $extendEnd   = $chromEnd;
            print OUT $chrom,"\t", $extendStart,"\t", $extendEnd,"\t",$lineInfojoin,"\n";
        }

    }

    close(IN);
    close(OUT);
    return $outFile;

}



sub getTrackRange{
	my ( $inputFile,$outDir ) = @_;
	open (INPUT,"<$inputFile");
    my $outPrefix=basename($inputFile);
	open (OUTPUT1,">$tmpFileDir/$final_prefix.trackRange.plus.bed");
	open (OUTPUT2,">$tmpFileDir/$final_prefix.trackRange.minus.bed");
	while ( my $line = <INPUT> ) {
	    my @items = split /\s+/, $line;
	    my $strand = $items[5];

	    if ( $line =~ /\#/ ) {
	        next;
	    }

	    if ( $strand =~ /\+/ ) {
	        print OUTPUT1 $line;
	    }

	    if ( $strand =~ /-/ ) {
	        print OUTPUT2 $line;
	    }
	}

	close(INPUT);
	close(OUTPUT1);
	close(OUTPUT2);	
}



sub trackBrowser{
	my ( $track1,$track2,$outDir ) = @_;
    my $outPrefix=basename($inputFile);
	my $commandLine1 = 
		"./configure/doTrackBrowser.pl " . 
		"--conf trackBrowserConfigure_".$org.".xml " . 
		"--track $track1 " . 
		"--range $tmpFileDir/$final_prefix.trackRange.plus.bed " . 
		"--output $trackSVGdir";
	executeCommand("$commandLine1");

	my $commandLine2 = 
		"./configure/doTrackBrowser.pl " . 
		"--conf trackBrowserConfigure_".$org.".xml " . 
		"--track $track2 " . 
		"--range $tmpFileDir/$final_prefix.trackRange.minus.bed " . 
		"--output $trackSVGdir";
	executeCommand("$commandLine2");
	
}