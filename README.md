# trackBrowser

trackBrowser: Visualization of sequencing data within a given genome region.

Overview:
---------

TrackBrowser is a perl script for visualization of sequencing data within a given genome region. Users can input sequencing data and regions of interest, and output visualized SVG images.

Usage:
---------

Usage:  plotTrackBrowser.pl  -org <species> -i <regions of interest> -upExtend <upstream extension> -downExtend <downstream extension> -upOffset <upstream offset> -downOffset <downstream offset> -outDir <output dir> -t1 <trackFile plus strand> -t2 <trackFile minus strand> <BR>
[options]<BR>
-org                    : species<BR>
-i                        : regions of interest<BR>
-upExtend          : upstream extension <BR>
-downExtend      : downstream extension <BR>
-upOffset            : upstream offset<BR>
-downOffset       : downstream offset<BR>
-outDir                : output dir<BR>
-t1                       : trackFile plus strand<BR>
-t2                       : trackFile minus strand<BR>


Installation:<BR>
---------

No installation required.<BR>

System requirements:
---------

(1) Some self-developed perl modules are located in directory ./packages, which could be load by the perl command: <BR>"BEGIN { unshift( @INC, "./packages" ); }"<BR>(2) Some perl modules need to be installed<BR>Users can search for the module names in the website meta::cpan (https://metacpan.org/) and get the installation method in the 'installation instructions'. Generally, it can be installed using the cpanm command<BR>cpanm XML::Simple<BR>
cpanm File::Basename<BR>
cpanm SVG<BR>
cpanm Bio::DB::BigFile<BR>
cpanm Bio::DB::BigFile::Constants<BR>
cpanm Bio::DB::BigWig<BR>
cpanm JSON::Parse<BR>
cpanm GD<BR>
cpanm feature<BR>
cpanm Bio::DB::HTS::Tabix<BR>
cpanm Bio::DB::HTS::Faidx<BR>

Prerequisites:<BR>
---------

(1) Genome:<BR>
You can download genome of interest from UCSC: wget -c 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'<BR>or use the demo genome we provided in demo directory: ./demo/hg38.chrY.fa<BR>

(2)Track data:<BR>
You can use your own sequencing data (bigwig format)<BR>or use the track data we provided in demo directory: ./demo/trackData<BR>

(3) Configuration file:<BR>
If you use the genome of another species, remember to change the contents of the genomeFile in the configuration file<BR>Customize the configuration file trackBrowserConfigure_{$org}.xml (e.g. ./configure/trackBrowserConfigure_hg38.xml) accorrding to your needs<BR>

run trackBrowser:
---------

perl plotTrackBrowser.pl -org hg38 -i ./demo/plot_region.bed -upExtend 20 -downExtend 20 -upOffset 0 -downOffset 0 -outDir . -t1 ./configure/trackFile_human_plusStrand.json -t2 ./configure/trackFile_human_minusStrand.json<BR>

Output:
---------

Visualized SVG images<BR>

<img src=".\results\chrY_14056154-14057020.png" alt="image-20220926095517453" style="zoom: 80%;" />



Acknowledgements:
---------

Thanks a lot to everyone who contributed to the public code used by trackBrowser.<BR>

Contact :
---------

*****************************************************************************************<BR>
 \*	trackBrowser: Visualization of sequencing data within a given genome region.<BR>
 \*<BR>
 \*	Author : Jian-Hua Yang <yangjh7@mail.sysu.edu.cn><BR>
 \* <BR>
 \*	RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
 \*	<BR>
 \*  Create date: 24/09/2022<BR>
 \*  <BR>
 \*  last modified time: 24/09/2022<BR>
 ****************************************************************************************<BR>
