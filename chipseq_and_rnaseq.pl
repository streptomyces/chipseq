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
my $chipfn;
my $rnafn;
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
"chipseqfn:s" => \$chipfn,
"rnaseqfn:s" => \$rnafn,
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

=pod

rm BldC_chipseq_rnaseq_combined.csv
perl code/chipseq_and_rnaseq.pl -outfile BldC_chipseq_rnaseq_combined.csv \
-chip BldCRegulon_lndiffFinalTable.csv -rnaseq bldC_rnaseq.txt
cp BldC_chipseq_rnaseq_combined.csv ~/mnt/wstemp/matt/


=cut


# {{{ POD markdown

=pod

=begin markdown

# script.pl

## Usage

~~~ {.sh}

perl script.pl -outfile outfn -- infile1 infile2 infile3 ...

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

my ($noex, $dir, $ext)= fileparse($rnafn, qr/\.[^.]*/);
my $rnabn = $noex . $ext;

my %rnaseq;
$skip = 1;
open(my $rnafh, "<$rnafn") or croak("Could not open $rnafn");
my $lineCnt = 0;
if($skip) {
for (1..$skip) { my $discard = readline($rnafh); }
}
# local $/ = ""; # For reading multiline records separated by blank lines.
while(my $line = readline($rnafh)) {
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @ll=split(/\t/, $line);
my $gene = $ll[0];
my @rlfc = @ll[4,7,10];
$rnaseq{$gene} = \@rlfc;
$lineCnt += 1;
if($testCnt and $lineCnt >= $testCnt) { last; }
if($runfile and (not -e $runfile)) { last; }
}
close($rnafh);



($noex, $dir, $ext)= fileparse($chipfn, qr/\.[^.]*/);
my $chipbn = $noex . $ext;

$skip = 1;
open(my $chipfh, "<$chipfn") or croak("Could not open $chipfn");
$lineCnt = 0;
if($skip) {
  for (1..$skip) {
    my $head = readline($chipfh); 
    chomp($head);
    tablist(split(/\t/, $head),
    "lgene_rnaseq10", "lgene_rnaseq14", "lgene_rnaseq22",
    "igene_rnaseq10", "igene_rnaseq14", "igene_rnaseq22",
    "rgene_rnaseq10", "rgene_rnaseq14", "rgene_rnaseq22"
    );
  }
}
# local $/ = ""; # For reading multiline records separated by blank lines.
while(my $line = readline($chipfh)) {
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @ll=split(/\t/, $line);
my $lgene = $ll[6] eq '-' ? undef : $ll[6];
my $igene = $ll[9] eq '-' ? undef : $ll[9];
my $rgene = $ll[12] eq '-' ? undef : $ll[12];
my @rlfc = split("", ("-" x 9));
if($lgene) {
my $lfcr = $rnaseq{$lgene};
splice(@rlfc, 0, 3, @{$lfcr});
}
if($igene) {
my $lfcr = $rnaseq{$igene};
splice(@rlfc, 3, 3, @{$lfcr});
}
if($rgene) {
my $lfcr = $rnaseq{$rgene};
splice(@rlfc, 6, 3, @{$lfcr});
}

tablist(@ll, @rlfc);


$lineCnt += 1;
if($testCnt and $lineCnt >= $testCnt) { last; }
if($runfile and (not -e $runfile)) { last; }
}
close($chipfh);


exit;

# Multiple END blocks run in reverse order of definition.
END {
close($ofh);
close(STDERR);
close(ERRH);
# $handle->disconnect();
}

