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

=head1 Name

bedgraph_diff.pl

=head2 Examples


perl code/bedgraph_log2ratio.pl -outfile bedgraph/lfc1.bedgraph \
-- bedgraph/tpm_3Flag1.bedgraph bedgraph/tpm_WT1.bedgraph

perl code/bedgraph_log2ratio.pl -outfile bedgraph/lfc.bedgraph \
-- bedgraph/tpm_3Flag.bedgraph bedgraph/tpm_WT.bedgraph

=head2 Description

log2 of column 4 of the first bedgraph file minus the log2 of
column 4 of the second bedgraph file.

The input file could come from bedgraph_mean.pl.

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

my @lr1;
my @lr2;
my $head1;
my $head2;

# {{{ Cycle through all the infiles.
my $cycle = 0;
for my $infile (@infiles) {
  $cycle += 1;
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
if($cycle == 1) {
$head1 = readline($ifh); chomp($head1)
}
if($cycle == 2) {
$head2 = readline($ifh); chomp($head2);
}

# local $/ = ""; # For reading multiline records separated by blank lines.
while(my $line = readline($ifh)) {
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @ll=split(/\t/, $line);
if($cycle == 1) {
push(@lr1, \@ll);
}
elsif($cycle == 2) {
push(@lr2, \@ll);
}

$lineCnt += 1;
if($testCnt and $lineCnt >= $testCnt) { last; }
if($runfile and (not -e $runfile)) { last; }
}
close($ifh);
close(OFH);

}
# }}}

unless($#lr1 == $#lr2) {
croak("Lists not of the same length\n");
}

$head1 =~ s/name=(.*?)\s/name=$1_diff /;
$head1 =~ s/description=(.*?)$/description=$1_lfc/;
$head1 =~ s/_mean//g;

linelist($head1);
for(my $dx = 0; $dx <= $#lr1; $dx += 1) {

my $signal = log2($lr1[$dx]->[3]);
my $baseline = log2($lr2[$dx]->[3]);
my $lfc = $signal - $baseline;

if($signal < 0.5 and $baseline < 0.5) { next; }

my @ol = @{$lr1[$dx]};
$ol[3] = $lfc;
tablist(@ol);
}


exit;

# Multiple END blocks run in reverse order of definition.
END {
close($ofh);
close(STDERR);
close(ERRH);
# $handle->disconnect();
}


sub log2 {
my $num = shift(@_);
$num += 1;
if($num < 0) { return 0; }
else {
return(log($num) / log(2));
}
}
