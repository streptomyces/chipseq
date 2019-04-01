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
my $errfile;
my $runfile;
my $outfile;
my $testCnt = 0;
our $verbose;
my $skip = 0;
my $help;
GetOptions (
"outfile:s" => \$outfile,
"outdir:s" => \$outdir,
"indir:s" => \$indir,
"fofn:s" => \$fofn,
"extension:s" => \$outex,
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

# bedgraph_diff.pl

## Usage

~~~ {.sh}

perl code/bedgraph_diff.pl -outfile stuff.bedgraph \
-- bedgraph/J2_ln.bedgraph bedgraph/WT2_ln.bedgraph

~~~

## Description

Subtracts the column 4 values of the second bedgraph file from the
column 4 values of the first bedgraph file and writes a new bedgraph
file with these values. The first three columns are not touched.
No checking is done other than that there should be the same number of
lines in the two files.

## Options

* -outfile
* -errfile
* -testcnt: Defaults to 0 which means process all input.
* -conffile: Defaults to local.conf
* -help: Display this documentation and exit.


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

# {{{ Outdir and outfile business.
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

my @lol;
my @head;
my $head1;
my $head2;

# {{{ Cycle through all the infiles.
my $cycle = 0;
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
my $temp = readline($ifh); chomp($temp);
$head[$cycle] = $temp;



# local $/ = ""; # For reading multiline records separated by blank lines.
while(my $line = readline($ifh)) {
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @ll=split(/\t/, $line);
push(@{$lol[$cycle]}, \@ll);

$lineCnt += 1;
if($testCnt and $lineCnt >= $testCnt) { last; }
if($runfile and (not -e $runfile)) { last; }
}
close($ifh);
close(OFH);

$cycle += 1;
}
# }}}

unless(checklol(@lol)) {
croak("Lists not of the same length\n");
}



$head[0] =~ s/name=(.*?)\s/name=$1_mean /;
$head[0] =~ s/description=(.*?)$/description=$1_mean/;
$head[0] =~ s/A//g;

linelist($head[0]);


my @fl = @{$lol[0]};

for(my $dx = 0; $dx <= $#fl; $dx += 1) {
my $rowtot = 0;
for my $lr (@lol) {
$rowtot += $lr->[$dx]->[3];
}
my $mean = $rowtot / scalar(@lol);
tablist(@{$fl[$dx]}[0..2], $mean);
}



exit;

# Multiple END blocks run in reverse order of definition.
END {
close($ofh);
close(STDERR);
close(ERRH);
# $handle->disconnect();
}


sub checklol {
my @lol = @_;
my @lens;

for my $lr (@lol) {
push(@lens, scalar(@{$lr}));
}

my $first = shift(@lens);

for my $len (@lens) {
  if($len != $first) {
    tablistE($first, @lens);
    return(0);
  }
}

return(1);
}

