#!/usr/bin/perl
use 5.14.0;
use Carp;
use File::Copy;
use File::Basename;
use lib qw(/home/sco /home/sco/perllib);
use File::Temp qw(tempfile tempdir);
my $tempdir = qw(/home/sco/volatile);
my $template="rcmdXXXXX";
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
    tablistV tablistVE linelistV linelistVE tablistH linelistH
    tablistER tablistVER linelistER linelistVER tabhashER tabhashVER);
use File::Spec;

# {{{ Getopt::Long stuff
use Getopt::Long;
my $testCnt = 0;
my $aperture = 30;      # Even number please.
my $region = 3000;      # Even number please.
my $outfile;
my $quiet = 0;
GetOptions (
"testcnt:i" => \$testCnt,
"regionwidth:i" => \$region,
"outfile=s" => \$outfile,
"aperture:i" => \$aperture,
"quiet" => \$quiet
);

=pod

=begin markdown

# localNormalisation.pl

## Description

Input is a depth file made by the samtools depth command.

Output is a bedgraph file.

## Usage

~~~ {.sh}

perl code/localNormalisation.pl -quiet \
-outfile bedgraph/Z4_ln.bedgraph -- depth/Z4.sorted.depth

~~~

## Options

* -regionwidth: integer. Typically 1000 to 6000
* -aperture: integer. Typically 20 to 50
* -outfile: The name of the output bedgraph file.
* -quiet: Prevent progress messages to STDERR (for parallel).

 
=end markdown

=cut

# }}}

if($region % 2) {
croak("Region width specified ($region) in not an even number");
}


# {{{ constants
my $halfReg = $region/2;
my $halfAp = $aperture/2;
my $regionMinusAperture = $region - $aperture;
# }}}


my @infiles = @ARGV;

# print(@infiles);

# {{{ input and output formats
#
### format of depth files ###
# KY5	1	38
# KY5	2	54
# KY5	3	59
# KY5	4	71
#############################

### format of bedgraph files ###
# track type=bedGraph name=WT3.sorted description=WT3.sorted
# KY5	1	10	96
# KY5	11	20	143
# KY5	21	30	188
# KY5	31	40	236
# KY5	41	50	277
# KY5	51	60	327
# KY5	61	70	378
# KY5	71	80	442
# KY5	81	90	484
################################
# }}}

my $ifn = shift(@infiles);
my ($basename, $directory, $ext)= fileparse($ifn, qr/\..*$/);
my $trackname = $basename . "_ln";
my $lnOut = $outfile; 
open(my $lnfh, ">", $lnOut);
my @cov;
open(INFILE, "<$ifn") or croak("Could not open $ifn");
my $lineCnt = 0;
while(<INFILE>) {
  my $line=$_;
  chomp($line);
  if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
  my @ll=split(/\t/, $line);
  # push(@cov, $ll[2]);
  push(@cov, \@ll);
}
close(INFILE);
unless($quiet) {
linelistER("$ifn read\n");
}

print($lnfh join(" ", "track", "type=bedGraph", "name=$trackname",
"description=$trackname"), "\n");
my $lastdx = $#cov;
my $pos = $halfAp;

my $doneCnt = 0;
while($pos <= ($#cov-$halfAp)) {
  my $regLeft = ($pos - $halfReg) >= 0 ? $pos-$halfReg : 0;
  my $regRight = ($regLeft + $region) <= $#cov ? $regLeft+$region : $#cov;

  my $apLeft = $pos - $halfAp;
  my $apRight = $pos + $halfAp;
  my @apRange = ($apLeft .. $apRight);

  my $apCov = 0;
  foreach my $temp ($apLeft..$apRight) {
    $apCov += $cov[$temp]->[2];
  }
  my $regCov = 0;
  foreach my $temp ($regLeft..($apLeft - 1)) {
    $regCov += $cov[$temp]->[2];
  }
  foreach my $temp (($apRight + 1)..$regRight) {
    $regCov += $cov[$temp]->[2];
  }

  my $apertureDensity = $apCov / $aperture;
  my $regionDensity = $regCov / $regionMinusAperture;

  my $enrich;
  if($regionDensity > 0) {
    $enrich = $apertureDensity / $regionDensity;
  }
  else {
    $enrich = 0;
  }

if($testCnt and (not $quiet)) {
  tablistE($regLeft, $apLeft, $apRight, $regRight, "-", $apCov, $regCov, "-",
  $apertureDensity, $regionDensity, $enrich);
}

#  print($lnfh "$pos\t$enrich\n");
  print($lnfh "$cov[$pos]->[0]\t$apLeft\t$apRight\t$enrich\n");

  $pos += $halfAp;
  unless($pos % 10000 and (not $quiet)) {
    linelistER("$pos           ");
  }
  $doneCnt += 1;
  if($testCnt and $doneCnt >= $testCnt) { last; }
}

close($lnfh);
unless($quiet) {
linelistE();
}


__END__
