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
my $scofn;
my $vnzfn;
my $namesfn = qq(vnzNames.list);
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
"scofn:s" => \$scofn,
"vnzfn:s" => \$vnzfn,
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

perl code/vnz_sco_bldc_regulon.pl \
-scofn BldCChIP-chipSCONatureCommunications.csv -vnzfn BldCRegulon_lndiffFinalTable.csv

rm BldC_vnz_sco_regulon_overlap.csv
perl code/vnz_sco_bldc_regulon.pl -outfile BldC_vnz_sco_regulon_overlap.csv \
-scofn BldCChIP-chipSCONatureCommunications.csv -vnzfn BldCRegulon_lndiffFinalTable.csv
cp BldC_vnz_sco_regulon_overlap.csv ~/mnt/wstemp/matt/


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

# {{{ Hash of vnz to sco orthologs. %scortho.

my %scortho;
open(my $namesfh, "<$namesfn") or croak("Could not open $namesfn");
while(my $line = readline($namesfh)) {
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @ll=split(/\t/, $line);
if($ll[0] =~ m/^vnz_/ and $ll[3] =~ m/^SCO/) {
$scortho{$ll[0]} = $ll[3];
}
}
close($namesfh);

# }}}

# {{{ Read Sco ChIP-Chip data from $scofn. %left, %right, %in.

my %left;
my %in;
my %right;
my %sco;
$skip = 1;
open(my $scofh, "<$scofn") or croak("Could not open $scofn");
my $lineCnt = 0;
if($skip) {
for (1..$skip) { my $discard = readline($scofh); }
}
# local $/ = ""; # For reading multiline records separated by blank lines.
while(my $line = readline($scofh)) {
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @ll=split(/\t/, $line);
my $lgene = $ll[4];
my $rgene = $ll[8];
my $igene = $ll[12];
if($lgene) {
$left{$lgene} = [@ll[5,6]];
$sco{$lgene} = ["left", @ll[5, 6]];
}
if($rgene) {
$right{$rgene} = [@ll[9,10]];
$sco{$rgene} = ["right", @ll[9, 10]];
}
if($igene) {
$in{$igene} = [@ll[13,14]];
$sco{$igene} = ["in", @ll[13, 14]];
}
$lineCnt += 1;
if($testCnt and $lineCnt >= $testCnt) { last; }
if($runfile and (not -e $runfile)) { last; }
}
close($scofh);

# Below was used for testing that no Sco gene occurs in two
# hashes.
# my %counts;
# for my $key ((keys(%left), keys(%right), keys(%in))) {
# # linelist($key);
# $counts{$key} += 1;
# }
# 
# for my $key (keys %counts) {
# if($counts{$key} == 1) {
# tablist($key, $counts{$key});
# }
# }


# }}}


# {{{ Read vnz ChIP-Seq data from $vnzfn.

$skip = 1;
open(my $vnzfh, "<$vnzfn") or croak("Could not open $vnzfn");
$lineCnt = 0;
if($skip) {
  for (1..$skip) {
    my $head = readline($vnzfh); 
    chomp($head);
    tablist(split(/\t/, $head),
    "lvnzgene", "lscogene", "lscofrompeak", "lscostrand", "lscodist",
    "ivnzgene", "iscogene", "iscofrompeak", "iscostrand", "iscodist",
    "rvnzgene", "rscogene", "rscofrompeak", "rscostrand", "rscodist"
    );
  }
}
# local $/ = ""; # For reading multiline records separated by blank lines.
while(my $line = readline($vnzfh)) {
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @ll=split(/\t/, $line);
my $lgene = $ll[6] eq '-' ? undef : $ll[6];
my $igene = $ll[9] eq '-' ? undef : $ll[9];
my $rgene = $ll[12] eq '-' ? undef : $ll[12];

my @add = split("", ("-" x 15));

my $vgcyc = 0;
for my $vnzgene ($lgene, $igene, $rgene) {
  if($vnzgene) {
      splice(@add, $vgcyc * 5, 1, $vnzgene);
    if(exists($scortho{$vnzgene})) {
      my $scogene = $scortho{$vnzgene};
      splice(@add, $vgcyc * 5 + 1, 1, $scogene);
      if(exists($sco{$scogene})) {
      my $scoref = $sco{$scogene};
      splice(@add, $vgcyc * 5 + 2, 3, @{$scoref});
      }
    }
  }
$vgcyc += 1;
}

# tablist(@ll);
tablist(@ll, @add);
# linelist();

$lineCnt += 1;
if($testCnt and $lineCnt >= $testCnt) { last; }
if($runfile and (not -e $runfile)) { last; }
}
close($vnzfh);

# }}}

exit;

# Multiple END blocks run in reverse order of definition.
END {
close($ofh);
close(STDERR);
close(ERRH);
# $handle->disconnect();
}

