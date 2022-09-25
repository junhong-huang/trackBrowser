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
use SVG;
use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use Bio::DB::BigWig 'binMean';
use JSON::Parse 'json_file_to_perl';
use JSON::Parse ':all';
use GD;
use GD::Text;
use feature qw( say );
use Bio::DB::HTS::Tabix;

my $help        = 0;
my $rangeFile   = "";
my $confXmlFile = "";
my $trackFile   = "";
my $outputDir   = "./output";
my $prefix      = "prefix";
my $confDir     = "./configure";

if ( scalar(@ARGV) < 2 ) {
    usage();
}

my $optLong = GetOptions(
    "conf=s"   => \$confXmlFile,    # configure file
    "track=s"  => \$trackFile,      # track file
    "range=s"  => \$rangeFile,      # range file
    "output=s" => \$outputDir,      # output dir name
    "help|h"   => \$help
) or usage();

&usage() if $help;

if ( $confXmlFile eq "" ) {
    die "Please set the --conf option\n";
}
if ( $trackFile eq "" ) {
    die "Please set the --track option\n";
}
if ( $rangeFile eq "" ) {
    die "Please set the --range option\n";
}

if ( !( -e $outputDir ) ) {
    warn "ceate output dir\n";
    &executeCommand("mkdir $outputDir");
}

&drawTrackLoop( $confXmlFile, $trackFile, $rangeFile, $outputDir );
##########################################################################
# MAIN FUNCTIONS as follows:
##########################################################################

sub usage
### usage function
{
    print
        "Usage: perl $0 [options] --conf <configure file> --track <track file> --range <range file>\n";
    print "--conf    configure file\n";
    print "--track   track file\n";
    print "--range   range file[chrx\tstart\tend]\n";
    print "--output  output directory\n";
    print "--help    help information\n";
    exit;
}

sub drawTrackLoop {
    my ( $confFile, $trackFile, $rangeFile, $outputDir ) = @_;
    my $extendLen = 0;#20
    my $orgExtendLen = $extendLen;


    open( RANGE, "<$rangeFile" ) || die "can't open the $rangeFile\n";
    while ( my $line = <RANGE> ) {
        $line =~ s/\s+$//;
        next if ( $line =~ /chrom/ );
        my @items      = split /\s+/, $line;
        # $extendLen     = int(($items[2] - $items[1]) * 0.2);
        $extendLen = $orgExtendLen if ($extendLen < $orgExtendLen);
        my $chromStart = $items[1] - $extendLen;
        my $chromEnd   = $items[2] + $extendLen;
        my $strand     = $items[5];
        $chromStart = 0 if ( $chromStart < 0 );
        my $range
            = $items[0] . ":" . $chromStart . "-" . $chromEnd . ":" . $strand;
        my $prefixName = $items[0] . ":" . $items[1] . "-" . $items[2];
        #my @names = split /\|/, $items[3];
        #if (scalar(@names)>=5)
        #{
        #    $prefixName = $prefixName."_".$names[3];
        #} 
        if ($strand eq "+")
        {
            $prefixName .= "_plus";
        }
        else
        {
            $prefixName .= "_minus";
        }
        my $commandLine
            = " perl "
            . $confDir . "/"
            . "trackBrowser.pl --conf "
            . $confDir . "/"
            . $confFile
            . " --track "
            . $trackFile
            . " --range "
            . $range
            . " --output "
            . $outputDir
            . " --prefix "
            . $prefixName;
        &executeCommand($commandLine);
    }
    close(RANGE);
}
