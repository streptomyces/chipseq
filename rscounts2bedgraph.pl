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
my $outdir = qq(bedgraph);
my $indir;
my $fofn;
my $outpre;
my $conffile = qq(local.conf);
my $errfile;
my $runfile;
my $genolen = 5365318;
my $outfile;
my $testCnt = 0;
my $sectlen = 30;
our $verbose;
my $skip = 0;
my $help;
GetOptions (
"outfile:s" => \$outfile,
"dirout:s" => \$outdir,
"indir:s" => \$indir,
"fofn:s" => \$fofn,
"genolen|length:i" => \$genolen,
"pre|outpre:s" => \$outpre,
"conffile:s" => \$conffile,
"errfile:s" => \$errfile,
"runfile:s" => \$runfile,
"testcnt:i" => \$testCnt,
"sectlen:i" => \$sectlen,
"skip:i" => \$skip,
"verbose" => \$verbose,
"help" => \$help
);
# }}}

# {{{ POD markdown

=head1 Name

rscounts2bedgraph.pl

=head2 Examples

perl code/rscounts2bedgraph.pl -sectlen 30 -genolen 5365318 -outdir bedgraph \
-- rscounts.txt

=head2 Description

Genome length is needed to make sure that the end remains within the
end of the genome.

$sectlen is used to work out the left and right coords for each section
which looks like C<Avin_15>.

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


my $halfsect = int($sectlen / 2);
my @ofh;

# {{{ Cycle through all the infiles.
for my $infile (@infiles) {
my ($noex, $dir, $ext)= fileparse($infile, qr/\.[^.]*/);
my $bn = $noex . $ext;
# tablistE($infile, $bn, $noex, $ext);

open(my $ifh, "<$infile") or croak("Could not open $infile");
my $lineCnt = 0;
my $header = readline($ifh);
chomp($header);
my @hl=split(/\t/, $header);
my @bns = @hl; shift(@bns);

my $cycle = 0;
for my $bn (@bns) {
open($ofh[$cycle], ">", File::Spec->catfile($outdir, $outpre . $bn . ".bedgraph"));
$cycle += 1;
}

$cycle = 0;
for my $ofh (@ofh) {
my $trackline = "track type=bedgraph name=" . $outpre . $bns[$cycle] . " description=" . $outpre . $bns[$cycle];
tablistH($ofh, $trackline);
$cycle += 1;
}


# local $/ = ""; # For reading multiline records separated by blank lines.
while(my $line = readline($ifh)) {
  chomp($line);
  if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
  my @ll=split(/\t/, $line);
  my $section = shift(@ll);
  my ($mol, $start, $end) = section2loc($section);
  my $cycle = 0;
  for my $ofh (@ofh) {
    tablistH($ofh, $mol, $start, $end, $ll[$cycle]);
    $cycle += 1;
  }
  $lineCnt += 1;
  if($testCnt and $lineCnt >= $testCnt) { last; }
  if($runfile and (not -e $runfile)) { last; }
}
close($ifh);
}
# }}}

exit;

# Multiple END blocks run in reverse order of definition.
END {
for my $ofh (@ofh) {
  close($ofh);
}
close(STDERR);
close(ERRH);
# $handle->disconnect();
}

sub section2loc {
my $section = shift(@_);
my ($mol, $mid) = split(/_/, $section);
my $start = $mid - $halfsect;
my $end = $mid + ($halfsect - 1);
if($end > $genolen) { $end = $genolen; }
return($mol, $start, $end);
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

