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
my $nrthresh = 10;
my $errfile;
my $progint = 100_000; # Progress indicator interval.
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
"nrthresh:i" => \$nrthresh,
"extension:s" => \$outex,
"conffile:s" => \$conffile,
"progint:i" => \$progint,
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

perl code/pe_sam_redundant_filter.pl -test 100000 -progint 10000 \
-nrthresh 6 \
-dirout samnr \
-- sam/aSepH-WT-22A-IP.sam
#;

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


# {{{ Cycle through all the infiles.
for my $infile (@infiles) {
  my ($noex, $dir, $ext)= fileparse($infile, qr/\.[^.]*/);
  my $bn = $noex . $ext;
# tablistE($infile, $bn, $noex, $ext);

  if($idofn) {
    my $ofn = File::Spec->catfile($outdir, $noex . "_nr" . $ext);
    open(OFH, ">", $ofn) or croak("Failed to open $ofn");
    select(OFH);
  }

my %posCnt;

  open(my $ifh, "<$infile") or croak("Could not open $infile");
  my $lineCnt = 0;
  if($skip) {
    for (1..$skip) { my $discard = readline($ifh); }
  }
# local $/ = ""; # For reading multiline records separated by blank lines.
  while(my $line = readline($ifh)) {
    chomp($line);
    if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
    if($line=~m/^@/) { linelist($line); next;}
    my @ll=split(/\t/, $line);
    if($ll[3] =~ m/\D/ or $ll[7] =~ m/\D/ or $ll[5] =~ m/[^0-9M]/) { next; }
    if(abs($ll[8]) > 0) {
      my $postr = $ll[3] . ":" . $ll[7];
      if($posCnt{$postr} >= $nrthresh) {
        $posCnt{$postr} += 1;
      }
      else {
        linelist($line);
        $posCnt{$postr} += 1;
      }
    }

    $lineCnt += 1;
    if($testCnt and $lineCnt >= $testCnt) { last; }
    if($runfile and (not -e $runfile)) { last; }
    unless($lineCnt % $progint) {
      tablistE($infile, $lineCnt);
    }
  }
  close($ifh);
  close(OFH);
# }}}

linelist(length(keys(%posCnt)));

for my $pos (keys(%posCnt)) {
if($posCnt{$pos} > $nrthresh) {
tablistE($pos, $posCnt{$pos});
}
}

}


exit;

# Multiple END blocks run in reverse order of definition.
END {
close($ofh);
close(STDERR);
close(ERRH);
# $handle->disconnect();
}

__END__

# {{{ For gzipped files.
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use File::Temp qw(tempfile tempdir);
my $tempdir = qq(/mnt/volatile);
my $template = qq(gunzippedXXXXX);
# Now you can use the gunzip function. AutoClose closes the file being
# written to after all the writing has been done.

my($gbfh, $gbfn)=tempfile($template, DIR => $tempdir, SUFFIX => '.gbff');
unless(gunzip $infile => $gbfh, AutoClose => 1) {
  close($gbfh); unlink($gbfn);
  die "gunzip failed: $GunzipError\n";
}
# }}}

