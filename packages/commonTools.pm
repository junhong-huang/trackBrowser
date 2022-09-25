package commonTools;
BEGIN {
    $VERSION = '0.1';
}
use strict;
use warnings;
use base 'Exporter';
our @EXPORT = qw(executeCommand limitProgramNum);

sub executeCommand {
    my $optNum      = scalar(@_);
    my $commandLine = "";
    my $nohup       = "";
    my $logFile     = "";
    if ( $optNum == 1 ) {
        $commandLine = shift @_;
    }
    elsif ( $optNum == 2 ) {
        ( $commandLine, $nohup ) = @_;
    }
    elsif ( $optNum == 3 ) {
        ( $commandLine, $nohup, $logFile ) = @_;
    }
    if ( $logFile ne "" ) {
        $commandLine = $commandLine . " 2>$logFile";
    }
    if ( $nohup eq "nohup" ) {
        $commandLine = "nohup " . $commandLine . " &";
    }
    warn "### Execute command: ", $commandLine, " ###\n";
    system($commandLine) == 0 or die "system # $commandLine # failed: $?\n";
}

sub limitProgramNum {
    my $optNum        = scalar(@_);
    my $program       = "";
    my $user          = "huangjh";
    my $maxProgramNum = 5;
    my $programNum    = 0;
    if ( $optNum == 1 ) {
        $program = shift @_;
    }
    elsif ( $optNum == 2 ) {
        ( $program, $maxProgramNum ) = @_;
    }
    elsif ( $optNum == 3 ) {
        ( $program, $maxProgramNum, $user ) = @_;
    }
    if ( $maxProgramNum > 70 ) {
        $maxProgramNum = 70;
    }
    if ( $maxProgramNum < 1 ) {
        $maxProgramNum = 1;
    }
    my $tmpFile = "program.log";
    while (1) {
        $programNum = 0;
        my $commandLine = "ps -fu " . $user . " >$tmpFile";
        system($commandLine) == 0 or die "system #$commandLine# failed: $?\n";

        open( TMP, "<$tmpFile" ) || die "can't open the $tmpFile\n";
        while ( my $line = <TMP> ) {
            $line =~ s/\s+$//;
            if ( $line =~ /$program/ ) {
                $programNum++;
            }
        }
        close(TMP);
        if ( $programNum <= $maxProgramNum ) {
            last;#跳出循环
        }
        else {
            sleep(60);
        }
        $programNum = 0;
    }
    my $commandLine = "rm -f $tmpFile";
    system($commandLine) == 0 or die "system #$commandLine# failed: $?\n";
}
1;
__END__
