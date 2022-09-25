package annoTools;

BEGIN {
    $VERSION = '0.1';
}

use strict;
use warnings;
use commonTools;

use base 'Exporter';
our @EXPORT = qw(annotateBed6Genes);

sub annotateBed6Genes {
    my ( $program, $options, $exonAnnotation, $targetFile ) = @_;
    if ( $program eq "" ) {
        $program = '/home/jianhua/softwares/bedtools/bin/bedtools';
    }
    my %targetHash = ();
    my $targetBed  = $targetFile . ".anno.bed";
    open( TARGET, "<$targetFile" ) || die "can't open the $targetFile\n";
    open( OUTPUT, ">$targetBed" )
        || die "can't open the $targetBed\n";
    while ( my $line = <TARGET> ) {
        $line =~ s/\s+$//;
        if ($line =~ /\#/)
        {
            next;
        }
        my @items = split /\t/, $line;
        $targetHash{ $items[3] } = $line;
        print OUTPUT $items[0], "\t", $items[1], "\t", $items[2], "\t",
            $items[3], "\t", $items[4], "\t", $items[5], "\n";
    }
    close(TARGET);
    close(OUTPUT);

    my $annotatedOutputFile = $targetFile . ".annotated";
    my $cmdLine
        = $program . " "
        . $options
        . " -a $targetBed -b "
        . $exonAnnotation
        . " >$annotatedOutputFile";
    &executeCommand($cmdLine);

    my %annotatedHash = ();
    open( ANNO, "<$annotatedOutputFile" )
        || die "can't open $annotatedOutputFile\n";
    while ( my $line = <ANNO> ) {
        $line =~ s/\s+$//;
        my @items = split /\t/, $line;
        $annotatedHash{ $items[3] } .= $items[9] . ",";
    }
    close(ANNO);

    my $allAnnoFile = $targetFile . ".anno.allTypes";

    my @suffixType = (
        "knownTranRNAs", "knownSnoRNAs", "knownMiRNAs",  "knownCDS",
        "known3UTR",     "known5UTR",    "proteinExons", "pseudogenes",
        "knownNcRNAs",   "novelSnoRNAs", "novelNcRNAs",  "knownRmsks",
        "intron",        "intergenic"
    );

    my @outfps = ();
    for ( my $i = 0; $i < scalar(@suffixType); $i++ ) {
        my $annoFile
            = $targetFile . ".anno." . $suffixType[$i];
        open( my $fp, ">$annoFile" ) || die "can't open the $annoFile\n";
        push @outfps, $fp;
    }
    open( my $afp, ">$allAnnoFile" ) || die "can't open the $allAnnoFile\n";

    foreach my $name ( sort keys %targetHash ) {
        my $anno
            = "intergenic" . "\t"
            . "intergenic" . "\t"
            . "intergenic" . "\t"
            . "intergenic";
        my $geneTypes   = "intergenic";
        my $geneRegions = "intergenic";
        my $geneId      = "intergenic";
        my $typeTag     = 0;
        my $geneSymbols = "intergenic";

        if ( defined( $annotatedHash{$name} ) ) {
            ( $geneTypes, $geneRegions, $geneId, $geneSymbols, $anno )
                = annotatedResults( $annotatedHash{$name} );
        }
        $typeTag
            = &getGeneType( $geneTypes, $geneRegions, $geneId, $geneSymbols );
        my $outfp = $outfps[$typeTag];
        print $outfp $targetHash{$name}, "\t", $anno, "\n";
        print $afp $targetHash{$name},   "\t", $anno, "\n";
    }
    my $commandLine = "rm -f $targetBed $annotatedOutputFile";
    &executeCommand($commandLine);
    return $annotatedOutputFile;
}

sub annotatedResults {
    my $targetGenes = shift @_;
    my @annos       = split /\,/, $targetGenes;
    my %annoId      = ();
    my %annoType    = ();
    my %annoRegion  = ();
    my $geneSymbols = '';
    my $geneIds     = '';
    my $geneTypes   = '';
    my $geneRegions = '';
    my $refAnno     = "NA" . "\t" . "NA" . "\t" . "NA";
    my $refTag      = 1;
    foreach my $anno (@annos) {
        $anno =~ s/\,$//;
        my ( $tid, $tgene, $gid, $gene, $biotype, $region ) = split /\|/,
            $anno;
        if ( defined( $annoId{ $gene . $region } ) ) {
            $annoId{ $gene . $region } .= "|" . $tid;
        }
        else {
            $annoId{ $gene . $region } .= $tid;
        }
        $annoType{ $gene . $region } = $biotype;
        my $tag = "a";
        if ( $biotype ne "protein_coding" ) {
            $tag = "b";
        }
        my $rtag = "g";
        if ( $region eq "utr3" ) {
            $rtag = "a";
        }
        if ( $region eq "cds" ) {
            $rtag = "b";
        }
        if ( $region eq "utr5" ) {
            $rtag = "c";
        }
        if ( $region eq "exon" ) {
            $rtag = "d";
        }
        $annoRegion{ $tag . "|" . $gene . "|" . $rtag . "|" . $region }++;
    }

    foreach my $key (
        sort { $annoRegion{$b} <=> $annoRegion{$a} || $a cmp $b }
        keys %annoRegion
        )
    {
        my ( $tag, $gene, $rtag, $region ) = split /\|/, $key;
        my $biotype = $annoType{ $gene . $region };
        my $tid     = $annoId{ $gene . $region };
        $geneSymbols .= $gene . ",";
        $geneRegions .= $region . ",";
        $geneIds     .= $tid . ",";
        $geneTypes   .= $biotype . ",";
    }
    $geneSymbols =~ s/\,$//;
    $geneIds =~ s/\,$//;
    $geneTypes =~ s/\,$//;
    $geneRegions =~ s/\,$//;
    return ( $geneTypes, $geneRegions, $geneIds, $geneSymbols,
              $geneSymbols . "\t"
            . $geneTypes . "\t"
            . $geneRegions . "\t"
            . $geneIds );
}

sub getGeneType {
    my ( $geneTypes, $geneRegions, $geneIds, $geneSymbols ) = @_;
    my @transGeneTypes = (
        'rRNA',     'snRNA',          'tRNA',    'scRNA',
        'misc_RNA', 'mascRNA-menRNA', 'Mt_rRNA', 'Mt_tRNA',
        'vaultRNA', 'ribozyme',       'sRNA'
    );
    my @knownSnoTypes
        = ( 'snoRNA', 'CDBox', 'HAcaBox', 'scaRna', 'SNORA', 'SNORD' );
    my @mirnaGeneTypes = ('miRNA');
    my @novelSnoTypes  = ( 'boxCD', 'boxHACA' );
    my @knownNcTypes   = (
        'ncRNA',                    'lncRNA',
        'lincRNA',                  'processed_transcript',
        '3prime_overlapping_ncrna', 'antisense',
        'non_coding',               'sense_intronic',
        'sense_overlapping',        'TEC',
        'known_ncrna',              'macro_lncRNA'
    );
    my @pseudogeneTypes = (
        'transcribed_processed_pseudogene',
        'transcribed_unprocessed_pseudogene',
        'unprocessed_pseudogene',
        'processed_pseudogene',
        'transcribed_unitary_pseudogene',
        'translated_processed_pseudogene',
        'translated_unprocessed_pseudogene',
        'unitary_pseudogene',
        'polymorphic_pseudogene',
        'pseudogene',
        'IG_C_pseudogene',
        'IG_J_pseudogene',
        'TR_V_pseudogene',
        'TR_J_pseudogene'
    );
    my @novelNcTypes        = ('novNcRNA');
    my @utr3                = ('utr3');
    my @utr5                = ('utr5');
    my @cds                 = ('cds');
    my @rmsk                = ( 'rmsk', 'trf', 'dust' );
    my @retainedIntronTypes = ('protein_coding');
    my @intron              = ('intron');
    my @intergenic          = ('intergenic');
    my @types               = (
        \@transGeneTypes,      \@knownSnoTypes,   \@mirnaGeneTypes,
        \@cds,                 \@utr3,            \@utr5,
        \@retainedIntronTypes, \@pseudogeneTypes, \@knownNcTypes,
        \@novelSnoTypes,       \@novelNcTypes,    \@rmsk,
        \@intron,              \@intergenic
    );
    my $typeNum = scalar(@types);
    my $typeTag = $typeNum - 1;

    my @gtypes   = split /\,/, $geneTypes;
    my @gregions = split /\,/, $geneRegions;
    my @gsymbols = split /\,/, $geneSymbols;

    for ( my $i = 0; $i < scalar(@types); $i++ ) {
        my $findTag  = 0;
        my @subTypes = @{ $types[$i] };
        foreach my $geneType (@subTypes) {
            for ( my $j = 0; $j < scalar(@gtypes); $j++ ) {
                my $gtype   = $gtypes[$j];
                my $gregion = $gregions[$j];
                my $gsymbol = $gsymbols[$j];
                if ( $gtype =~ /^$geneType/i && $gregion !~ /^intron/i ) {
                    $typeTag = $i;
                    $findTag = 1;
                    last;
                }
                if ( $gsymbol =~ /^$geneType\d+/i && $gregion !~ /^intron/i )
                {
                    $typeTag = $i;
                    $findTag = 1;
                    last;
                }
                if ( $gregion =~ /^$geneType/i && $gregion !~ /^intron/i ) {
                    $typeTag = $i;
                    $findTag = 1;
                    last;
                }
            }    ### for gtypes
            if ( $findTag == 1 ) {
                last;
            }
        }    ### for subTypes
        if ( $geneIds =~ /$subTypes[0]\d+/i ) {
            $typeTag = $i;
            $findTag = 1;
            last;
        }
        if ( $findTag == 1 ) {
            last;
        }
    }    ### for @types
    ### for intron
    if ( ( $typeTag == $typeNum - 1 ) && $geneRegions =~ /intron/i ) {
        $typeTag -= 1;
    }
    return $typeTag;
}

1;
__END__
