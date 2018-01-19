#!/usr/bin/env perl
# Last changed Time-stamp: <2018-01-19 18:24:17 mtw>
# -*-CPerl-*-
#
# NOTE: Before running this script, it is recommended to set two
# environment variables:
# $KINCONG_KINFOLD_PATH : path of the Kinfold executable
# $KINCONG_RNALILA_PATH : path to the RNAlila executables (ie RNAwalk)
#
# Alternatively, paths to the Kinfold and RNAwalk executables can be
# provided via the command line

use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use version; my $VERSION= qv('0.1');
use Pod::Usage;
use Data::Dumper;
use strict;
use warnings;
use Carp;
use IPC::Cmd qw(can_run run);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($Kinfold_path, $Kinfold_exe, $Kinfold_cmd);
my ($RNAwalk_path, $RNAwalk_exe, $RNAwalk_cmd);
my ($seq, $i,$glen,$time,$grow,$num, $verbose);
my ($logprefix, $logsuffix, $logfile, $kinfold);
my ($want_ctsteps);
my %minima = ();
my %endstrucs = ();
my @ctsteps = ();

if(defined $ENV{'KINCONG_KINFOLD_PATH'}){
    $Kinfold_path = $ENV{'KINCONG_KINFOLD_PATH'};
    $Kinfold_exe  = $Kinfold_path."/Kinfold";
}
else{$Kinfold_exe = "Kinfold";}
if(defined $ENV{'KINCONG_RNALILA_PATH'}){
    $RNAwalk_path = $ENV{'KINCONG_RNALILA_PATH'};
    $RNAwalk_exe  = $RNAwalk_path."/RNAwalk";
}
else{$RNAwalk_exe = "RNAwalk"};

$logprefix    = "foo";
$logsuffix    = ".KinCong.log";
$verbose      = 0;
$want_ctsteps = 0;
$glen         = 10;
$time         = 0;
$grow         = 4000;  # 4000
$num          = 1000;
$kinfold      = undef;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless
  GetOptions("Kinfold=s"          => \$Kinfold_exe,
	     "RNAwalk=s"          => \$RNAwalk_exe,
	     "c|ctsteps"          => sub{$want_ctsteps=1},
	     "g|grow=i"           => \$grow,
	     "l|glen=i"           => \$glen,
	     "n|num=i"            => \$num,
	     "v|verbose"          => sub{$verbose=1},
	     "L|log=s"            => \$logprefix,
	     "man"                => sub{pod2usage(-verbose => 2)},
	     "h|help"             => sub{pod2usage(1)}
	    );

$logfile = $logprefix.$logsuffix;
open (LOG, ">$logfile");
print LOG "# KinCong $VERSION \n";

# test if Kinfold binary accessible
# test if input integers > 0

# read seqence from stdin
$seq = <STDIN>;
chomp($seq);
croak "not an RNA input sequence, exiting\n"
  unless ($seq =~ /^[AUGCaugc]+$/);

# test if sequence length >= $glen
croak "input sequence shorter than $glen, exiting\n"
  unless (length($seq) >= $glen);

print LOG "# input sequence: $seq\n";
print LOG "# Kinfold: $Kinfold_exe \n";
print LOG "# RNAwalk: $RNAwalk_exe \n";
print LOG "# sequence length: ".length($seq)."\n";
print LOG "# seed length: $glen\n";
print LOG "# grow speed: $grow\n";
print LOG "# number of trajectories: $num\n";

# compute time
$time = (length($seq)-$glen)*$grow;
print LOG "# computed time: $time\n";

# build Kinfold command line
$kinfold = can_run($Kinfold_exe) or croak "ERROR: $kinfold not found";
$Kinfold_cmd = "echo '$seq' | $kinfold --met --noShift --time $time --glen=$glen --grow=$grow --num 1";
print LOG "# Kinfold command: ".$Kinfold_cmd."\n";


for($i=0;$i<$num;$i++){
  my ($struct, $gradmin, $minstat, $ekin, $egrad);
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $Kinfold_cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR Call to $kinfold unsuccessful\n";
    print STDERR "ERROR This is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
  my $stdout_buffer = join "", @$stdout_buf;
  my @out = split /\n/, $stdout_buffer;
  push @ctsteps, scalar @out;
  my $kinfendline = $out[-1];

  croak "Cannot parse dot bracket string in Kinfold output"
    unless ($kinfendline =~ /^(\S+)\s+(-?\d+\.\d+)/);
  $struct = $1;
  $ekin   = $2;
  if ($verbose == 1){
    printf LOG ("%s %6.2f KinStop\n",$struct,$ekin);
  }
  $endstrucs{$struct}++;

  # get corresponding gradient minimum
  my $RNAwalk_cmd = "echo \"$seq\n$struct\"  | ".$RNAwalk_exe." -t G -l 10000";
  print LOG "# RNAwalk command: $RNAwalk_cmd\n";
  open(GRAD, $RNAwalk_cmd."| tail -1 |") or die "cannot open RNAwalk pipe";
  while (<GRAD>){
    if ($verbose == 1){
	my $line = $_; chomp $line; printf LOG ("%s\n",$line);
    }
    croak "Cannot parse dot bracket string in RNAwalk output"
      unless ($_ =~ /^([.()]+)\s+(\-?\d+\.\d+)\s([SDT\*])/);
    $gradmin = $1;
    $egrad   = $2;
    $minstat = $3;
    if ($verbose == 1){
	printf LOG ("%s %6.2f gradmin %s\n",$gradmin,$egrad,$minstat);
    }
    $minima{$gradmin}++;
  }
}

if ($want_ctsteps){
  open my $ctsfile, ">", "kincong.cts";
  foreach my $i (@ctsteps){print $ctsfile "$grow\t$i\n"}
  close($ctsfile);
}

if ($verbose == 1) {
    print LOG "# gradient minima and number of hits\n";
    print LOG Dumper(\%minima);
}

# print minima (key from histogram) with probability = value/$num probability
my $sum = 0.;
my $ratio = 0.;
foreach my $k (keys %minima){
  $ratio = $minima{$k}/$num;
  printf ("%s %6.4f\n",$k,$ratio);
  $sum += $ratio;
}

printf LOG ("# Checksum: %10.8f\n",$sum);
carp "WARNING: Population probabilities don't sum up to 1.0"
    if ($sum > 1.0001 || $sum < 0.9999);

close(LOG);
__END__

=head1 NAME

KinCong.pl - Co-transcriptional initialization for kinetic RNA folding

=head1 SYNOPSIS

KinCong.pl [--Kinfold F<EXE>] [--RNAwalk F<EXE>] [-g|--glen I<INT>]
[-n|--num I<INT>] [options]

=head1 DESCRIPTION

This tool computes start populations for (coarse grained) kinetic RNA
folding from co-transcriptional sampling. Summary statistics are
derived from Kinfold trajectories of growing chain fragments. For
this, a gradient walk is performed starting at each fully elongated
Kinfold stop structure.

This approach provides a good estimation for the population density of
different macrostates at the end of transcription. In this line,
KinCong yields biophysically justifiable data and outperforms previous
assumptions about starting conditions (i.e. starting at the unfolded
open chain conformation).

=head1 PREREQUISTES

This script runs two external programs, I<Kinfold> from the ViennaRNA
package and I<RNAwalk> from RNAlila. The script needs to know where to
locate these two third party executables. This can be accomplished via
two environment variables, which need to be set in your environment:

=over

=item I<$KINCONG_KINFOLD_PATH>: path of the Kinfold executable

=item I<$KINCONG_RNALILA_PATH>: path to the RNAlila executables (ie RNAwalk)

=back

Alternatively, the full path to both the Kinfold and RNAwalk
executables can be provided via command line options.

=head1 OPTIONS

=over

=item B<--Kinfold>

Full path of the F<Kinfold> executable

=item B<--RNAwalk>

Full path of the RNAwalk executable

=item B<--num -n>

Number of Kinfold trajectories to sample (default: 1000)

=item B<--glen -l>

Initial size of the growing chain (default: 10)

=item B<--grow -g>

Let the chain grow every I<INT> time units. This differes from the
original Kinfold command line option, which takes a I<FLOAT> here.

=item B<--log>

Log file extension. Default is ".KinCong.log". The log file is created
in the current working directory.

=item B<--debug -d>

Turn on debugging

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
