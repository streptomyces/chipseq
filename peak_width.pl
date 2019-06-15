#!/usr/bin/perl
use 5.14.0;
use utf8;
use Carp;
use lib qw(/home/sco /home/sco/mnt/smoke/perllib);
use File::Basename;
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
    tablistV tablistVE linelistV linelistVE tablistH linelistH
    tablistER tablistVER linelistER linelistVER tabhashER tabhashVER csvsplit);
use File::Spec;

# {{{ Getopt::Long
use Getopt::Long;
my $outdir;
my $indir;
my $fofn;
my $outex; # extension for the output filename when it is derived on infilename.
my $conffile = qq(local.conf);
my $bedgraphfn;
my $errfile;
my $runfile;
my $outfile;
my $testCnt = 0;
our $verbose;
my $skip = 0;
my $help;
GetOptions (
"outfile:s" => \$outfile,
"dirout:s" => \$outdir,
"indir:s" => \$indir,
"fofn:s" => \$fofn,
"extension:s" => \$outex,
"bedgraphfn:s" => \$bedgraphfn,
"conffile:s" => \$conffile,
"errfile:s" => \$errfile,
"runfile:s" => \$runfile,
"testcnt:i" => \$testCnt,
"skip:i" => \$skip,
"verbose" => \$verbose,
"help" => \$help
);
# }}}

# {{{ POD markdown

=pod

=begin markdown

# script.pl

## Usage

~~~ {.sh}

perl code/peak_width.pl -bedgraph aarout/BldC_10.csv \
-outfile out -verbose -test 0 -- List_for_peak_width_analysis.csv

~~~

## Description

A description of what this script does.

## Options

* -outfile
* -errfile
* -testcnt: Defaults to 0 which means process all input.
* -conffile: Defaults to local.conf
* -help: Display this documentation and exit.

These are the commonly used one. There are a few more.

=end markdown

=cut

# }}}

if($help) {
exec("perldoc $0");
exit;
}

# {{{ open the errfile
if($errfile) {
open(ERRH, ">", $errfile);
print(ERRH "$0", "\n");
close(STDERR);
open(STDERR, ">&ERRH"); 
}
# }}}

# {{{ Populate %conf if a configuration file 
my %conf;
if(-s $conffile ) {
  open(my $cnfh, "<", $conffile);
  my $keyCnt = 0;
  while(my $line = readline($cnfh)) {
    chomp($line);
    if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
    my @ll=split(/\s+/, $line, 2);
    $conf{$ll[0]} = $ll[1];
    $keyCnt += 1;
  }
  close($cnfh);
  linelistE("$keyCnt keys placed in conf.");
}
elsif($conffile ne "local.conf") {
linelistE("Specified configuration file $conffile not found.");
}
# }}}

# {{{ outdir and outfile business.
my $ofh;
my $idofn = 0;    # Flag for input filename derived output filenames. 
if($outfile) {
  my $ofn;
  if($outdir) {
    unless(-d $outdir) {
      unless(mkdir($outdir)) {
        croak("Failed to make $outdir. Exiting.");
      }
    }
    $ofn = File::Spec->catfile($outdir, $outfile);
  }
  else {
    $ofn = $outfile;
  }
  open($ofh, ">", $ofn);
}
elsif($outdir) {
linelistE("Output filenames will be derived from input");
linelistE("filenames and placed in $outdir");
    unless(-d $outdir) {
      unless(mkdir($outdir)) {
        croak("Failed to make $outdir. Exiting.");
      }
    }
$idofn = 1;
}
else {
  open($ofh, ">&STDOUT");
}
select($ofh);
# }}}

# {{{ populate @infiles
my @infiles;
if(-e $fofn and -s $fofn) {
open(FH, "<", $fofn);
while(my $line = readline(FH)) {
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my $fn;
if($indir) {
$fn = File::Spec->catfile($indir, $line);
}
else {
$fn = $line;
}

push(@infiles, $fn);
}
close(FH);
}
else {
@infiles = @ARGV;
}

# }}}

# read the bedgraph file into @bg.
my @bg;
open(my $bgfh, "<$bedgraphfn") or croak("Could not open $bedgraphfn");
my $bghead = readline($bgfh);
while(my $line = readline($bgfh)) {
  chomp($line);
  if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
  my @ll=split(/\t/, $line);
  push(@bg, \@ll);
}
close($bgfh);
my @temp = sort bgsorter @bg;
@bg = @temp; @temp = ();
my $nrowbg = scalar(@bg);
# linelistVE($nrowbg);

# {{{ Cycle through all the infiles.
for my $infile (@infiles) {
my ($noex, $dir, $ext)= fileparse($infile, qr/\.[^.]*/);
my $bn = $noex . $ext;
# tablistE($infile, $bn, $noex, $ext);

if($idofn) {
my $ofn = File::Spec->catfile($outdir, $noex . "_out" . $outex);
open(OFH, ">", $ofn) or croak("Failed to open $ofn");
select(OFH);
}


open(my $ifh, "<$infile") or croak("Could not open $infile");
my $lineCnt = 0;
$skip = 1;
if($skip) {
for (1..$skip) { my $discard = readline($ifh); }
}
# local $/ = ""; # For reading multiline records separated by blank lines.
while(my $line = readline($ifh)) {
chomp($line);
# if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
# my @ll=split(/\t/, $line);
my $position = $line;
# tablist($position);
my $posdx = getdx($position);
unless(defined($posdx)) { 
  tablist($position, split("", "-" x 3));
  next;
}
my $rr = $bg[$posdx];
# tablist($posdx, @{$rr});

my @before = before($posdx, $position);
# for my $rr (@before) {
# tablist(@{$rr});
# }
# linelist("------");
my @after = after($posdx, $position);

my @sorted = sort bgsorter (@before, @after);

my $start = $sorted[0]->[0];
my $end = $sorted[$#sorted]->[0];
tablist($position, $start, $end, $end - $start);


# for my $rr (@sorted) {
# tablist(@{$rr});
# }
# linelist();

$lineCnt += 1;
if($testCnt and $lineCnt >= $testCnt) { last; }
if($runfile and (not -e $runfile)) { last; }
}
close($ifh);
close(OFH);
}
# }}}

exit;

# Multiple END blocks run in reverse order of definition.
END {
close($ofh);
close(STDERR);
close(ERRH);
# $handle->disconnect();
}

# {{{ sub bgsorter
sub bgsorter {
  return($a->[0] <=> $b->[0]);
}
# }}}


sub getdx {
my $pos = shift(@_);
my $dx = 0;
for my $rr (@bg) {
if($rr->[0] == $pos) {
return($dx);
}
$dx += 1;
}
return();
}


# {{{ sub before
sub before {
  my $dx = shift(@_);
  my $pos = shift(@_);
  # $dx += 1;
  my @retlist = $bg[$dx];
  while($dx > 0) {
    my $currpos = $bg[$dx]->[0];
    $dx -= 1;
    if(($currpos - $bg[$dx]->[0]) <= 30) {
      push(@retlist, $bg[$dx]);
    }
    else {
      last;
    }
  }
  return(@retlist);
}
# }}}

# {{{ sub after
sub after {
  my $dx = shift(@_);
  my $pos = shift(@_);
  my @retlist;
  while($dx < ($nrowbg)) {
    my $currpos = $bg[$dx]->[0];
    $dx += 1;
    if($dx >= $nrowbg) { last; }
    if(($bg[$dx]->[0] - $currpos) <= 30) {
      push(@retlist, $bg[$dx]);
    }
    else {
      last;
    }
  }
  return(@retlist);
}
# }}}

__END__

# R code to combine the two csv files output by two runs
# of this script.

cn <- c("pp", "left", "right", "width");
pw10 <- read.table("../BldC_10_peakwidth.csv",
quote = "", stringsAsFactors = FALSE,
header = F, row.names = 1, col.names = cn, sep = "\t");

pw14 <- read.table("../BldC_14_peakwidth.csv",
quote = "", stringsAsFactors = FALSE,
header = F, row.names = 1, col.names = cn, sep = "\t");

pw <- pw10;

colnames(pw) <- c("left10h", "right10h", "width10h");

pw$left14h <- pw14[rownames(pw), "left"]; 
pw$right14h <- pw14[rownames(pw), "right"]; 
pw$width14h <- pw14[rownames(pw), "width"]; 

# Note that all columns are character because some have "-" in
# the csv file.

head(pw);

fh <- file("../peakwidths.csv", open = "w");
h <- paste(c("peakpos", colnames(pw)), collapse = "\t");
writeLines(h, fh);
write.table(pw, file = fh, quote = F, col.names = F,
row.names = T, sep = "\t");
close(fh);
#

