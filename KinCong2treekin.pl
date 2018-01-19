#!/usr/bin/env perl
# Last changed Time-stamp: <2018-01-19 18:30:33 mtw>
# -*-CPerl-*-

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use RNA;

my $barfile = undef;
my $command;
my $silent=0;
my $t8 = undef;
my %pp = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless
  GetOptions(
	     "b|barfile=s" => \$barfile,
	     "s|silent"    => sub{$silent=1},
	     "t|teight"    => \$t8,
	     "h|help"      => sub{pod2usage(1)},
	     "man"         => sub{pod2usage(-exitstatus => 1, -verbose => 2)}
	    ) and $barfile or pod2usage(-verbose => 1);


while (<>) {
  my @KC = split;
  my $idx = map_to_barfile($KC[0], $barfile, $silent);
  if (defined($pp{$idx})){$pp{$idx}+=$KC[1]}
  else {$pp{$idx}=$KC[1]}
}

foreach my $i (sort {$a <=> $b} keys %pp){
  $command .= " --p0 $i=$pp{$i}";
}

 if (defined $t8){
      $command .= " --t8 $t8 ";
  }

print "treekin --t0 0.001 ".$command." -m I --bin < rates.bin\n";


# ********************************* #
# ~~~~~~~~~~ Subroutines ~~~~~~~~~~ #
# ................................. #

sub map_to_barfile {
  my ($str, $barfile, $silent) = @_;

  my $fh = new IO::File "< $barfile" or die "no file: $barfile found\n";
  my $firstline = <$fh>;
  my $sequence = $1 if ($firstline =~ /(\S+)/);

  my $idx = 0;
  my ($mdx, $min_dist, $mstr);
  while(<$fh>) {
    my @line = split;
    if ($line[1] eq $str) {
      $idx=$line[0];
      last;
    } else {
      my $dist = RNA::bp_distance($line[1], $str);
      ($min_dist, $mdx, $mstr) = ($dist, $line[0], $line[1]) if (!$mdx || $min_dist > $dist);
    }
  } 
  $fh->close;

  if (!$idx) {
    print STDERR "$0 -- No exact match found, minimum bpdist = $min_dist:\n$str\n$mstr ($mdx)\n" unless $silent;
    return $mdx;
  }
  return $idx;
}

__END__

=head1 NAME

KinCong2treekin.pl - construct a treekin call from KinCong output

=head1 SYNOPSIS

KinCong2treekin.pl [-b|--barfile I<FILE>] [-t|--teight = I<INT>]
[-s|--silent] < kinkong.out

=head1 DESCRIPTION

Map a set of local minima (I<lm>) and their occupancy (I<o>)to
corresponding leaves of a barrier tree. This information is printed in
form of a I<treekin> commandline call with the respective --p0
I<lm>=I<o> statements. If no directly corresponding local minimum is
found in the barrier tree, the structure will be mapped to the local
minimum with the least base-pair distance and a warning is issued.

=head1 PREREQUISTES

This script is ment to be used in a pipeline to post-process output
from I<KinCong.pl> and translate it into a system call for
I<treekin>. It requires a .bar file from I<barriers> for mapping
I<KinCong.pl> output to macrostates (i.e. local minima) of the energy
landscape.

=head1 OPTIONS

=over

=item B<-b, --barfile> I<FILE>

Specify a barrier tree file in the format from I<barriers>.

=item B<-s, --silent>

Switch off the warnings if no exact match was found for a local minimum.

=item B<-h, --help>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHORS

=over

=item Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=item Stefan Badelt E<lt>badelt@caltech.eduE<gt>

=back

=cut
