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
use Bio::DB::HTS::Faidx;

my $help        = 0;
my $range       = "";
my $confXmlFile = "";
my $trackFile   = "";
my $outputDir   = "./output";
my $prefix      = "prefix";

if ( scalar(@ARGV) < 2 ) {
    usage();
}

my $optLong = GetOptions(
    "conf=s"   => \$confXmlFile,    # configure file
    "track=s"  => \$trackFile,      # track file
    "range=s"  => \$range,          # range file
    "prefix=s" => \$prefix,         # prefix of output file
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
if ( $range eq "" ) {
    die "Please set the --range option\n";
}

if ( !( -e $outputDir ) ) {
    #warn "ceate output dir\n";
    &executeCommand("mkdir $outputDir");
}

### configure file
my $confHash = XMLin($confXmlFile);
### track file
my $trackHash = &getTracks($trackFile);
our ( $chrom, $chromStart, $chromEnd, $strand ) = &getGenomicPos($range);
our $span          = $chromEnd - $chromStart + 1;
our $labelWidth    = $confHash->{'labelWidth'};
our $trackWidth    = $confHash->{'width'} - $labelWidth;
our $baseWidth     = $trackWidth / $span;
our $trackSpace    = $confHash->{'trackSpace'};
our $tickSize      = $confHash->{'tickSize'};
our $fontFamily    = $confHash->{'fontFamily'};
our $fontSize      = $confHash->{'fontSize'};
our $itemSpace     = 3;
our $lineWidth     = 0.5;
our $tickNum       = 3;
our $minWidth      = 0.1;
our $minHeight     = 0.1;
our $arrowFontSize = $fontSize / 2;
our ( $fontHeight, $fontSpace ) = &getFixFontSize( $fontFamily, $fontSize );
our $defaultFontColor = "black";
our $changeVal     = 0;
our $setMaxVal     = 40000;
our $setMinVal     = 0;

my $genomeFasta = $confHash->{'genomeFile'};
my $seq         = &getSequence($genomeFasta);

#warn $seq, "\t", $seqLen, "\n";

my $svgFile = &drawTracks( $confHash, $trackHash, $outputDir, $prefix );

##########################################################################
# MAIN FUNCTIONS as follows:
##########################################################################

sub usage
### usage function
{
    print
        "Usage: perl $0 [options] --conf <configure file> --track <track file> --range <genomic position>\n";
    print "--conf    configure file\n";
    print "--track   track file\n";
    print "--range   genomic position[chrx:start-end]\n";
    print "--help    help information\n";
    exit;
}

sub getTracks {
    my $jsonFile = shift @_;
    my $jsonHash = json_file_to_perl($jsonFile);

    #print $jsonHash->{'tracks'}->[0]->{'trackType'}, "\n";
    return $jsonHash;
}

sub getSequence {
    my $genomeFasta = shift @_;
    my $seqIndex    = Bio::DB::HTS::Faidx->new($genomeFasta);
    my $location    = $chrom . ":" . $chromStart . "-" . $chromEnd;
    my ( $seq, $seqLen ) = $seqIndex->get_sequence($location);
    $seq = uc($seq);
    if ( $strand eq "-" ) {

        #$seq = reverse($seq);
        $seq =~ tr/ACGTacgt/TGCAtgca/;
    }
    return $seq;
}

sub drawTracks {
    my ( $confHash, $trackHash, $outputDir, $prefix ) = @_;
    my $height    = 20;                      # cummulated height
    my $svgWidth  = $confHash->{'width'};
    my $svgHeight = $confHash->{'height'};

    # create an SVG object
    my $svg = SVG->new( width => $svgWidth, height => $svgHeight );###创建SVG对象
    my $arrowId = 'whiteArrow';
    &drawArrowMarker( $svg, $arrowId, '255,255,255' );

    # &drawTrackLine( $svg, $height, $svgHeight );

    $height = &drawScaleAndPosition( $svg, $height, $tickNum ) + $trackSpace;

    my @tracks   = @{ $trackHash->{'tracks'} };
    my $trackNum = scalar(@tracks);
    for ( my $i = 0; $i < $trackNum; $i++ ) {
        my $track = $tracks[$i];
        if ( $track->{'display'} eq "hide" ) {
            next;
        }
        my $trackType = $track->{'trackType'};
        &drawArrowMarker( $svg, $track->{'trackId'}, $track->{'color'} );
        if ( $trackType eq "bigWig" ) {
            $height = &drawWigSVG( $svg, $track, $height );
            $height += $trackSpace;
        }
        if ( $trackType eq "bed" ) {
            $height = &drawBedSVG( $svg, $track, $height );
            $height += $trackSpace;
        }
    }
    my $svgFile = $outputDir . "/" . $prefix . ".svg";
    open( OUTPUT, ">$svgFile" ) || die "can't open the $svgFile\n";
    print OUTPUT $svg->xmlify, "\n";
    close(OUTPUT);
}

sub drawArrowMarker {
    my ( $svg, $markerId, $markerColor ) = @_;
    my $trackColor = 'rgb(' . $markerColor . ')';
    my $arrowId    = "arrow_" . $markerId;
    my $tag        = $svg->marker(
        id           => $arrowId,
        markerWidth  => 4,
        markerHeight => 4,
        refX         => 2,
        refY         => 2,
        orient       => "auto"
    );
    $tag->tag(
        'path',
        d     => "M0,0 L4,2 0,4",
        id    => $arrowId . "_path",
        style => {
            'fill-opacity' => 0,
            'stroke'       => $trackColor,
            'stroke-width' => $lineWidth,
        }
    );
}

sub drawScaleAndPosition {
    my ( $svg, $height, $expectNum ) = @_;
    my ( $posRef, $gap )
        = &generateAxis( $chromStart, $chromEnd, $expectNum );
    my @positions  = @{$posRef};
    my $viewHeight = $height + $tickSize;
    my $gapPix     = $gap * $baseWidth;
    my $x1         = $labelWidth + ( $trackWidth - $gapPix ) / 2;
    my $x2         = $x1 + $gapPix;

    my $tag = $svg->line(
        x1    => $x1,
        y1    => $height,
        x2    => $x2,
        y2    => $height,
        style => {
            'fill-opacity' => 1,
            'stroke'       => 'rgb(0,0,0)',
            'stroke-width' => $lineWidth
        }
    );
    $tag = $svg->line(
        x1    => $x1,
        y1    => $height - $tickSize,
        x2    => $x1,
        y2    => $height,
        style => {
            'fill-opacity' => 1,
            'stroke'       => 'rgb(0,0,0)',
            'stroke-width' => $lineWidth
        }
    );

    $tag = $svg->line(
        x1    => $x2,
        y1    => $height - $tickSize,
        x2    => $x2,
        y2    => $height,
        style => {
            'fill-opacity' => 1,
            'stroke'       => 'rgb(0,0,0)',
            'stroke-width' => $lineWidth
        }
    );
    my ( $gapWidth, $gapHeight )
        = &getLabelSize( $fontFamily, $fontSize, $gap . " bases" );

    my $labelX = $labelWidth + $trackWidth / 2 + $gapWidth / 2;
    my $labelY = $height - $itemSpace;

    my $text = $svg->text(
        x      => $labelX,
        y      => $labelY,
        -cdata => $gap . " bases",
        style  => {
            'font'              => $fontFamily,
            'font-size'         => $fontSize,
            'fill'              => 'rgb(0,0,0)',
            'text-anchor'       => 'end',
            'dominant-baseline' => "baseline",
        }
    );

    $labelX = $labelWidth + ( $trackWidth - $gapPix ) / 2 - $itemSpace;
    $svg->text(
        x      => $labelX,
        y      => $height,
        -cdata => $chrom,
        style  => {
            'font'              => $fontFamily,
            'font-size'         => $fontSize,
            'fill'              => 'rgb(0,0,0)',
            'text-anchor'       => 'end',
            'dominant-baseline' => "baseline",
        }
    );

    $height     += $gapHeight + $itemSpace;
    $viewHeight += $gapHeight + $itemSpace;

    foreach my $pos (@positions) {
        my $x = $labelWidth + ( $pos - $chromStart ) * $baseWidth;
        my $tag = $svg->line(
            x1    => $x,
            y1    => $height,
            x2    => $x,
            y2    => $viewHeight,
            style => {
                'fill-opacity' => 1,
                'stroke'       => 'rgb(0,0,0)',
                'stroke-width' => $lineWidth
            }
        );
        my ( $nameWidth, $nameHeight )
            = &getLabelSize( $fontFamily, $fontSize, $pos );
        my $labelX = $x + $nameWidth / 2;
        my $labelY = $height;
        my $text   = $svg->text(
            x      => $labelX,
            y      => $labelY,
            -cdata => $pos,
            style  => {
                'font'              => $fontFamily,
                'font-size'         => $fontSize,
                'fill'              => 'rgb(0,0,0)',
                'text-anchor'       => 'end',
                'dominant-baseline' => "baseline",
            }
        );
    }

    my @bases = split //, $seq;
    my ( $nameWidth, $nameHeight )
        = &getLabelSize( $fontFamily, $fontSize, $bases[0] );
    if ( scalar(@bases) * $nameWidth < $trackWidth ) {
        $viewHeight += $gapHeight;
        my %baseColor = (
            'A' => 'rgb(0,204,0)',
            'C' => 'rgb(29,32,136)',
            'G' => 'rgb(243,152,0)',
            'T' => 'rgb(204,0,0)'
        );
        my $idx = 1;
        foreach my $base (@bases) {
            my $x = $labelWidth + $idx * $baseWidth;
            my ( $nameWidth, $nameHeight )
                = &getLabelSize( $fontFamily, $fontSize, $base );
            my $labelX = $x;
            my $labelY = $viewHeight;
            my $color  = $baseColor{$base};
            my $text   = $svg->text(
                x      => $labelX,
                y      => $labelY,
                -cdata => $base,
                style  => {
                    'font'              => $fontFamily,
                    'font-size'         => $fontSize,
                    'fill'              => $color,
                    'text-anchor'       => 'end',
                    'dominant-baseline' => "baseline",
                }
            );
            $idx++;
        }
    }
    return $viewHeight + $itemSpace;
}

sub getOptInterval {
    my $inv = shift @_;
    if ( $inv > 5 ) {
        return $inv;
    }
    elsif ( $inv > 2 ) {
        return $inv * 2;
    }
    else {
        return 10;
    }
}

sub generateAxis {
    my ( $start, $end, $expectNum ) = @_;
    my $dif         = $end - $start + 1;
    my $difGap      = $dif / ( $expectNum - 1 );
    my $exponent    = &log10($difGap) - 1;
    my $exponentInt = int($exponent);

    if ( $exponent < 0 && abs($exponent) > 1e-8 ) {
        $exponentInt--;
    }

    my $step = int( $difGap / ( 10**$exponentInt ) );
    my @fixSteps = ( 10, 20, 25, 50, 100 );
    my $fixStep = 10;
    for ( my $i = scalar(@fixSteps) - 1; $i >= 1; $i-- ) {
        if ( $step > ( $fixSteps[$i] + $fixSteps[ $i - 1 ] ) / 2 ) {
            $fixStep = $fixSteps[$i];
            last;
        }
    }

    my $degreeGap = $fixStep * ( 10**$exponentInt );
    my $start1    = $start / $degreeGap;
    my $start2    = int($start1);
    if ( $start1 < 0 && abs( $start1 - $start2 ) > 1e-8 ) {
        $start2--;
    }
    my $degreeStart = $start2;

    my $end1 = $end / $degreeGap;
    my $end2 = int($end1);
    if ( $end1 >= 0 && abs( $end1 - $end2 ) > 1e-8 ) {
        $end2++;
    }
    my $degreeEnd = $end2;

    my @results = ();
    for ( my $i = 0; $i <= $expectNum; $i++ ) {
        my $value = ( $degreeStart * $degreeGap + $i * $degreeGap );
        if ( $value > $start && $value <= $end ) {
            push @results, $value;
        }
    }

    #warn "@results\n";

    return ( \@results, $degreeGap );
}

sub log10 {
    my $n = shift @_;
    return log($n) / log(10);
}

sub drawBedSVG {
    my ( $svg, $track, $height ) = @_;
    my $bedFile    = $track->{'data'};
    my $trackColor = 'rgb(' . $track->{'color'} . ')';
    my $trackName  = $track->{'name'};
    my $bedItemRef = &getBedItems($bedFile);
    my @bedItems   = @{$bedItemRef};
    my $viewHeight = $height;
    my $geneIdx    = 0;

    @bedItems = sort {
               ( $a->{'start'} <=> $b->{'start'} )
            || ( $a->{'score'} <=> $b->{'score'} )
    } @bedItems;

    for my $bed (@bedItems) {
        $geneIdx++;
        my $drawName   = 1;
        my $chr        = $bed->{'chr'};
        my $start      = $bed->{'start'};
        my $end        = $bed->{'end'};
        my $name       = $bed->{'name'};
        my $strand     = $bed->{'strand'};
        my $thickStart = $bed->{'thickStart'};
        my $thickEnd   = $bed->{'thickEnd'};
        my $blockCount = $bed->{'blockCount'};
        my $arrowId    = "arrow_" . $track->{'trackId'};

        my ( $nameWidth, $nameHeight )
            = &getLabelSize( $fontFamily, $fontSize, $name );

        for ( my $i = 0; $i < $blockCount; $i++ ) {
            my $blockStart = $start + $bed->{'blockStarts'}->[$i] + 1;
            my $blockEnd   = $blockStart + $bed->{'blockSizes'}->[$i] - 1;

            my $blockThickStart = 0;
            my $blockThickEnd   = 0;

            last if ( $blockStart > $chromEnd );
            next if ( $blockEnd < $chromStart );

            $blockStart = $chromStart if ( $blockStart < $chromStart );
            $blockEnd   = $chromEnd   if ( $blockEnd > $chromEnd );

            if (&overlapLength( $blockStart, $blockEnd, $thickStart,
                    $thickEnd ) > 0
                )
            {
                if ( $thickStart < $blockStart ) {
                    $blockThickStart = $blockStart;
                }
                else {
                    $blockThickStart = $thickStart;
                }
                if ( $thickEnd < $blockEnd ) {
                    $blockThickEnd = $thickEnd;
                }
                else {
                    $blockThickEnd = $blockEnd;
                }
            }
            my $startDist = $baseWidth * ( $blockStart - $chromStart );
            my $x         = $labelWidth + $startDist;
            my $w         = $baseWidth * ( $blockEnd - $blockStart + 1 );
            $w = $minWidth if ( $w < $minWidth );
            my $h = $nameHeight / 2;
            my $y = $viewHeight + $nameHeight / 4;
            $svg->rectangle(
                x      => $x,
                y      => $y,
                width  => $w,
                height => $h,
                style  => {
                    'fill'           => $trackColor,
                    'stroke'         => 'black',
                    'stroke-width'   => 0,
                    'stroke-opacity' => 1,
                    'fill-opacity'   => 1
                },
            );

            if ( $blockThickStart != 0 && $blockThickEnd != 0 ) {
                my $thickDist
                    = $baseWidth * ( $blockThickStart - $chromStart );
                my $tx = $labelWidth + $thickDist;
                my $tw
                    = $baseWidth * ( $blockThickEnd - $blockThickStart + 1 );
                $tw = $minWidth if ( $tw < $minWidth );
                my $th = $nameHeight;
                my $ty = $viewHeight;
                $svg->rectangle(
                    x      => $tx,
                    y      => $ty,
                    width  => $tw,
                    height => $th,
                    style  => {
                        'fill'           => $trackColor,
                        'stroke'         => 'black',
                        'stroke-width'   => 0,
                        'stroke-opacity' => 1,
                        'fill-opacity'   => 1
                    },
                );
            }
        }    # for $i blockCount

        if ($drawName) {
            my $labelX = $labelWidth - $itemSpace;
            my $labelY = $viewHeight + $fontSpace / 2;

            #$labelY = $viewHeight + $nameHeight / 4 if ( $blockCount == 1 );
            $start = $start + 1;
            $start = $chromStart if ( $start < $chromStart );
            my $regionStartPix = $baseWidth * ( $start - $chromStart );
            $end = $chromEnd if ( $end > $chromEnd );
            my $regionEndPix = $baseWidth * ( $end - $chromStart + 1 );

            if ( $regionStartPix > $nameWidth ) {
                $labelX = $labelWidth + $regionStartPix - $itemSpace;
            }
            my $text = $svg->text(
                x      => $labelX,
                y      => $labelY,
                -cdata => $name,

                #'transform' => $translateStr,
                style => {
                    'font'              => $fontFamily,
                    'font-size'         => $fontSize,
                    'fill'              => $trackColor,
                    'text-anchor'       => 'end',
                    'dominant-baseline' => "hanging",
                }
            );

            if ( $blockCount > 1 ) {
                my $tag = $svg->line(
                    x1    => $labelWidth + $regionStartPix,
                    y1    => $viewHeight + $nameHeight / 2,
                    x2    => $labelWidth + $regionEndPix,
                    y2    => $viewHeight + $nameHeight / 2,
                    style => {
                        'fill-opacity' => 0,
                        'stroke'       => $trackColor,
                        'stroke-width' => $lineWidth
                    }
                );
            }    # if blockCount
            else {
                $arrowId = "arrow_whiteArrow";
            }
            my @xps = ();
            my @yps = ();
            for (
                my $i = $regionStartPix;
                $i < $regionEndPix;
                $i += $itemSpace * 2
                )
            {
                my $xp = $labelWidth + $i;
                my $yp = $viewHeight + $nameHeight / 2;
                push @xps, $xp;
                push @yps, $yp;
            }
            if ( $strand eq "-" ) {
                @xps = reverse(@xps);
                @yps = reverse(@yps);
            }
            my $points = $svg->get_path(
                x     => \@xps,
                y     => \@yps,
                -type => 'polyline'
            );

            my $tagLine = $svg->polyline(
                %$points,
                'id'         => 'arrowPolyline_' . $name . "_" . $geneIdx,
                'marker-mid' => "url(#" . $arrowId . ")",
                style        => {
                    'fill-opacity'   => 0,
                    'stroke'         => $trackColor,
                    'stroke-opacity' => 0,
                }
            );
            $drawName = 0;
        }    # drawName end

        $viewHeight += $nameHeight + $itemSpace;
    }    # for points
    &drawBedTrackInfo( $svg, $trackName, $height, $viewHeight, $trackColor );
    return $viewHeight;
}

sub minValue {
    my ( $v1, $v2 ) = @_;
    if ( $v1 < $v2 ) {
        return $v1;
    }
    else {
        return $v2;
    }
}

sub maxValue {
    my ( $v1, $v2 ) = @_;
    if ( $v1 > $v2 ) {
        return $v1;
    }
    else {
        return $v2;
    }
}

sub overlapLength {
    my ( $qStart, $qEnd, $tStart, $tEnd ) = @_;
    my $maxStart = 0;
    my $minEnd   = 0;
    $maxStart = maxValue( $qStart, $tStart );
    $minEnd = minValue( $qEnd, $tEnd );
    return ( $minEnd - $maxStart );
}

sub getBedItems {
    my $bedFile  = shift @_;
    my @bedItems = ();
    my $tabix    = Bio::DB::HTS::Tabix->new( filename => $bedFile );

    #say $tabix->header;
    my $region = $chrom . ":" . $chromStart . "-" . $chromEnd;
    my $iter   = $tabix->query($region);
    if ( !defined($iter) ) {
        return \@bedItems;
    }
    while ( my $line = $iter->next ) {

        #say $line;
        my %itemHash = (
            'chr'         => '',
            'start'       => 0,
            'end'         => 0,
            'name'        => '',
            'score'       => 900,
            'strand'      => '+',
            'thickStart'  => 0,
            'thickEnd'    => 0,
            'itemRgb'     => 0,
            'blockCount'  => 0,
            'blockSizes'  => [],
            'blockStarts' => [],
        );

        my @items = split /\s+/, $line;
        my $itemNum = scalar(@items);
        $itemHash{'chr'}    = $items[0];
        $itemHash{'start'}  = $items[1];
        $itemHash{'end'}    = $items[2];
        $itemHash{'name'}   = $items[3] if ( defined( $items[3] ) );
        $itemHash{'score'}  = $items[2] - $items[1];
        $itemHash{'strand'} = $items[5] if ( defined( $items[5] ) );
        if ( $itemNum >= 12 ) {
            $itemHash{'thickStart'} = $items[6];
            $itemHash{'thickEnd'}   = $items[7];
            $itemHash{'itemRgb'}    = $items[8];
            $itemHash{'blockCount'} = $items[9];
            $items[10] =~ s/\,$//;
            $items[11] =~ s/\,$//;
            my @sizes  = split /\,/, $items[10];
            my @starts = split /\,/, $items[11];
            $itemHash{'blockSizes'}  = \@sizes;
            $itemHash{'blockStarts'} = \@starts;
        }
        else {
            $itemHash{'thickStart'}  = $items[1];
            $itemHash{'thickEnd'}    = $items[2];
            $itemHash{'blockCount'}  = 1;
            $itemHash{'blockSizes'}  = [ $items[2] - $items[1] ];
            $itemHash{'blockStarts'} = [0];
        }
        if (overlapLength(
                $itemHash{'start'}, $itemHash{'end'},
                $chromStart - 1,    $chromEnd
            ) > 0
            )
        {
            push @bedItems, \%itemHash;
        }
    }
    $tabix->close;
    return \@bedItems;
}

sub drawBedTrackInfo {
    my ( $svg, $trackName, $height, $viewHeight, $trackColor ) = @_;

    my $tag = $svg->line(
        x1    => $labelWidth - $tickSize,
        y1    => $height,
        x2    => $labelWidth,
        y2    => $height,
        style => {
            'fill-opacity' => 1,
            'stroke'       => $trackColor,
            'stroke-width' => $lineWidth
        }
    );

    $tag = $svg->line(
        x1    => $labelWidth,
        y1    => $height,
        x2    => $labelWidth,
        y2    => $viewHeight,
        style => {
            'fill-opacity' => 0,
            'stroke'       => $trackColor,
            'stroke-width' => $lineWidth
        }
    );

    $tag = $svg->line(
        x1    => $labelWidth - $tickSize,
        y1    => $viewHeight,
        x2    => $labelWidth,
        y2    => $viewHeight,
        style => {
            'fill-opacity' => 0,
            'stroke'       => $trackColor,
            'stroke-width' => $lineWidth
        }
    );

    #track name
    #track
    my ( $nameWidth, $nameHeight )
        = &getLabelSize( $fontFamily, $fontSize, $trackName );

    #my $text = $svg->text(
    #    x      => $labelWidth - $tickSize,
    #    y      => $height + ( $viewHeight - $height ) / 2,
    #    -cdata => $trackName,
    #    style  => {
    #        'font'        => $fontFamily,
    #        'font-size'   => $fontSize,
    #        'fill'        => $trackColor,
    #        'text-anchor' => 'end'
    #    }
    #);
    my $text = $svg->text(
        x      => $labelWidth + $trackWidth / 2 + $nameWidth / 2,
        y      => $height - $itemSpace,
        -cdata => $trackName,
        style  => {
            'font'              => $fontFamily,
            'font-size'         => $fontSize,
            'fill'              => $trackColor,
            'text-anchor'       => 'end',
            'dominant-baseline' => "baseline",
        }
    );
    return $tag;
}

sub drawTrackLine {
    my ( $svg, $height, $svgHeight ) = @_;
    my $tag = $svg->line(
        x1    => $labelWidth,
        y1    => $height,
        x2    => $labelWidth,
        y2    => $svgHeight,
        style => {
            'fill-opacity' => 0,
            'stroke'       => 'rgb(250,123,23)'
        }
    );
    return $tag;
}

sub drawWigTrackInfo {
    my ( $svg, $trackName, $height, $viewHeight, $trackColor, $minVal,
        $maxVal )
        = @_;
    $minVal = sprintf( "%.2f", $minVal );
    $maxVal = sprintf( "%.2f", $maxVal );

    $minVal =~ s/\.0+$//;
    $maxVal =~ s/\.0+$//;

    my $minStr = $minVal;
    my $maxStr = $maxVal;

    my $tag = $svg->line(
        x1    => $labelWidth - $tickSize,
        y1    => $height,
        x2    => $labelWidth,
        y2    => $height,
        style => {
            'fill-opacity' => 0,
            'stroke'       => $trackColor,
            'stroke-width' => $lineWidth
        }
    );

    $tag = $svg->line(
        x1    => $labelWidth,
        y1    => $height,
        x2    => $labelWidth,
        y2    => $viewHeight,
        style => {
            'fill-opacity' => 0,
            'stroke'       => $trackColor,
            'stroke-width' => $lineWidth
        }
    );

    $tag = $svg->line(
        x1    => $labelWidth - $tickSize,
        y1    => $viewHeight,
        x2    => $labelWidth,
        y2    => $viewHeight,
        style => {
            'fill-opacity' => 0,
            'stroke'       => $trackColor,
            'stroke-width' => $lineWidth
        }
    );

    #maxVal
    my ( $maxWidth, $maxHeight )
        = &getLabelSize( $fontFamily, $fontSize, $maxStr );

    #warn "max:" . $maxStr, "\t", $maxWidth, "\t", $maxHeight, "\n";
    my $text = $svg->text(
        x      => $labelWidth - $tickSize,
        y      => $height - $maxHeight / 2 + $fontSpace / 2,
        -cdata => $maxStr,
        style  => {
            'font'              => $fontFamily,
            'font-size'         => $fontSize,
            'fill'              => $trackColor,
            'text-anchor'       => 'end',
            'dominant-baseline' => "hanging",
        }
    );

    #$svg->rectangle(
    #    x      => $labelWidth - $maxWidth - $tickSize,
    #    y      => $height,
    #    width  => $maxWidth,
    #    height => $maxHeight,
    #    style  => {
    #        'fill'           => $trackColor,
    #        'stroke'         => 'black',
    #        'stroke-width'   => 1,
    #        'stroke-opacity' => 1,
    #        'fill-opacity'   => 0,
    #    },
    #);

    #minVal
    my ( $minWidth, $minHeight )
        = &getLabelSize( $fontFamily, $fontSize, $minStr );

    #warn "min:" . $minStr, "\t", $minWidth, "\t", $minHeight, "\n";
    $text = $svg->text(
        x      => $labelWidth - $tickSize,
        y      => $viewHeight - $minHeight / 2 + $fontSpace / 2,
        -cdata => $minStr,
        style  => {
            'font'              => $fontFamily,
            'font-size'         => $fontSize,
            'fill'              => $trackColor,
            'text-anchor'       => 'end',
            'dominant-baseline' => "hanging",
        }
    );

    #track name
    #track
    my ( $nameWidth, $nameHeight )
        = &getLabelSize( $fontFamily, $fontSize, $trackName );
    $text = $svg->text(
        x      => $labelWidth - $tickSize,
        y      => $height + ( $viewHeight - $height ) / 2 - $nameHeight / 2,
        -cdata => $trackName,
        style  => {
            'font'              => $fontFamily,
            'font-size'         => $fontSize,
            'fill'              => $trackColor,
            'text-anchor'       => 'end',
            'dominant-baseline' => "hanging",
        }
    );
    return $tag;
}

sub getLabelWidthHeight {
    my ( $fontFamily, $fontSize, $label ) = @_;
    my $font = lc($fontFamily);
    my $gd_text = GD::Text->new() or die GD::Text::error();
    $gd_text->set_font( $font, $fontSize ) or die $gd_text->error;
    $gd_text->set_text($label);
    my ( $w, $h ) = $gd_text->get( 'width', 'height' );

    return ( $w, $h );
}

sub getWidthAndHeight {
    my $boundRef = shift @_;

    #@bounds[0,1]  Lower left corner (x,y)
    #@bounds[2,3]  Lower right corner (x,y)
    #@bounds[4,5]  Upper right corner (x,y)
    #@bounds[6,7]  Upper left corner (x,y)
    # llx lower left x
    # lly lower left y
    # lrx lower right x
    # lry lower right y
    # etc

    my @bounds = @{$boundRef};
    my ( $llx, $lly, $lrx, $lry, $ulx, $uly, $urx, $ury ) = @bounds;
    my ( $w, $h );
    if ( $lly == $lry ) {

        # the string box is horizontal or vertical
        $w = abs( $lrx - $llx ) - 1;
        $h = abs( $uly - $lly ) - 1;
    }
    else {
        $w = sqrt(
            ( abs( $lrx - $llx ) - 1 )**2 + ( abs( $lry - $lly ) - 1 )**2 );
        $h = sqrt(
            ( abs( $ulx - $llx ) - 1 )**2 + ( abs( $uly - $lly ) - 1 )**2 );
    }
    return ( $w, $h );
}

# Return the width and height of text
# modified from circos, thanks
sub getLabelSize {
    my ( $fontFamily, $fontSize, $label ) = @_;
    my $black = 0;
    my $angle = 0;
    my $x     = 0;
    my $y     = 0;

    my %fontHash = (
        'Arial' =>
            './fonts/arial.ttf',
        'Courier' =>
            './fonts/cour.ttf',
        'CMUBright-Roman' =>
            './fonts/cmunbmr.ttf'
    );

    my $font = $fontHash{$fontFamily};

    #my $font      = "./cour.ttf";
    my @strings    = ('SNRUAC');
    my $meanHeight = 0;
    foreach my $ch (@strings) {
        my @chBounds = GD::Image->stringFT(
            $black, $font,
            $fontSize,
            $angle, $x, $y, $ch,
            {   linespacing => 0.6,
                charmap     => 'Unicode',
            }
        );
        my ( $w, $h ) = &getWidthAndHeight( \@chBounds );
        $meanHeight = $h;
    }

    my @bounds = GD::Image->stringFT(
        $black, $font,
        $fontSize,
        $angle, $x, $y, $label,
        {   linespacing => 0.6,
            charmap     => 'Unicode',
        }
    );

    my ( $width, $height ) = &getWidthAndHeight( \@bounds );

    #warn "@bounds\n", "\n";
    #warn $label . "\twidth:\t" . $width . "\theight:\t", $height,
    #    "\tcharHeight:", $meanHeight, "\n";

    #$h = $fontSize;
    $height = $meanHeight
        if ( $meanHeight > $height / 2 );
    return ( $width, $height );
}

# Return the width and height of text
# modified from circos, thanks
sub getFixFontSize {
    my ( $fontFamily, $fontSize ) = @_;
    my $black = 0;
    my $angle = 0;
    my $x     = 0;
    my $y     = 0;

    my %fontHash = (
        'Arial' =>
            './fonts/arial.ttf',
        'Courier' =>
            './fonts/cour.ttf',
        'CMUBright-Roman' =>
            './fonts/cmunbmr.ttf'
    );

    my $font = $fontHash{$fontFamily};

    #my $font      = "./cour.ttf";
    my @strings    = ( 'SNRUAC', '|||||||||' );
    my @heights    = ();
    my $fontSpace  = 0;
    my $fontHeight = 0;
    foreach my $ch (@strings) {
        my @chBounds = GD::Image->stringFT(
            $black, $font,
            $fontSize,
            $angle, $x, $y, $ch,
            {   linespacing => 0.6,
                charmap     => 'Unicode',
            }
        );
        my ( $w, $h ) = &getWidthAndHeight( \@chBounds );
        push @heights, $h;
    }
    $fontHeight = $heights[0];
    $fontSpace  = $heights[1] - $heights[0];
    $fontSpace  = 0 if ( $fontSpace < 0 );
    return ( $fontHeight, $fontSpace );
}

sub drawWigSVG {
    my ( $svg, $track, $height ) = @_;
    my $wigFile = $track->{'data'};
    my ( $maxPt, $viewPt, $minPt ) = split /\:/, $track->{'maxHeightPixels'};
    my ( $pointRef, $minVal, $maxVal ) = &getWigVals($wigFile);
    $minVal = 0 if ( $minVal > 0 );
    if ( defined( $track->{'viewLimits'} ) && $track->{'autoScale'} eq "off" )
    {
        ( $minVal, $maxVal ) = split /\:/, $track->{'viewLimits'};
    }
    if ($changeVal)
    {
        $minVal = $setMinVal;
        $maxVal = $setMaxVal;
    }
    my $valSpan = $maxVal - $minVal;
    my @points  = @{$pointRef};
    if ( $valSpan == 0 ) {
        $valSpan = 1;
    }
    my $valHeight  = $viewPt / $valSpan;
    my $trackColor = 'rgb(' . $track->{'color'} . ')';
    my $trackName  = $track->{'name'};
    my $viewHeight = $height + $viewPt;

    &drawWigTrackInfo( $svg, $trackName, $height, $viewHeight, $trackColor,
        $minVal, $maxVal );

    for my $p (@points) {
        my $start = $p->start;
        my $end   = $p->end;
        my $val   = abs( $p->score );
        $val = $val - $minVal;
        $val = 0 if ( $val < 0 );
        $val = $valSpan if ( $val > $valSpan );
        my $x = $labelWidth + $baseWidth * ( $start - $chromStart );
        my $w = $baseWidth * ( $end - $start + 1 );
        $w = $minWidth if ( $w < $minWidth );
        my $h = $valHeight * $val;
        $h = $minHeight if ( $h > 0 && $h < $minHeight );
        my $y = $viewHeight - $h;

        $svg->rectangle(
 
            x      => $x,
            y      => $y,
            width  => $w,
            height => $h,
            style  => {
            
                'fill'           => $trackColor,
                'stroke'         => 'black',
                'stroke-width'   => 0,
                'stroke-opacity' => 1,
                'fill-opacity'   => 1,
            },
        );
    }    # for points
    return $viewHeight;
}

sub getGenomicPos {
    my $range = shift @_;
    $range =~ s/\,//g;
    my $chr    = "";
    my $start  = 0;
    my $end    = 0;
    my $strand = '+';
    if ( $range =~ /^(\S+)\:(\d+)\-(\d+)$/ ) {
        $chr   = $1;
        $start = $2;
        $end   = $3;
    }
    else {
        if ( $range =~ /^(\S+)\:(\d+)\-(\d+)\:([+-])$/ ) {
            $chr    = $1;
            $start  = $2;
            $end    = $3;
            $strand = $4;
        }
        else {
            die "Error Format for genomic position: " . $range;
        }
    }
    return ( $chr, $start, $end, $strand );
}

sub getWigVals {
    my $bwFile    = shift @_;
    my @positions = ();
    my $maxVal    = 0 - 1e9;
    my $minVal    = 1e9;
    my $wig       = Bio::DB::BigWig->new( -bigwig => $bwFile );
    my @points    = $wig->features(
        -seq_id => $chrom,
        -start  => $chromStart,
        -end    => $chromEnd
    );

    for my $p (@points) {
        my $start = $p->start;
        my $end   = $p->end;
        my $val   = abs( $p->score );
        $minVal = $val if ( $val < $minVal );
        $maxVal = $val if ( $val > $maxVal );
    }
    return ( \@points, $minVal, $maxVal );
}
