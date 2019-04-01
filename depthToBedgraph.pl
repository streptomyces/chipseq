#!/usr/bin/perl
use 5.14.0;
use utf8;
use Carp;
use lib qw(/home/sco /home/sco/perllib);
use File::Basename;
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
    tablistV tablistVE linelistV linelistVE tablistH linelistH
    tablistER tablistVER linelistER linelistVER tabhashER tabhashVER csvsplit);
use File::Spec;

# {{{ Getopt::Long
use Getopt::Long;
my $indir;
my $conffile = qq(local.conf);
my $errfile;
my $ofn;
my $aper = 10;
my $testCnt = 0;
our $verbose;
my $skip = 0;
my $help;
GetOptions (
"outfile:s" => \$ofn,
"indir:s" => \$indir,
"aperture:i" => \$aper,
"conffile:s" => \$conffile,
"errfile:s" => \$errfile,
"testcnt:i" => \$testCnt,
"skip:i" => \$skip,
"verbose" => \$verbose,
"help" => \$help
);
# }}}

# {{{ POD Example

=head1 Name

Change me.

=head2 Example

 perl changeme.pl -outfile stuff.out -- inputfile1 inputfile2

Note that input files are always specified as non-option arguments.

=cut

# }}}

# {{{ POD blurb

=head2 Blurb

Some kind of description here.

=cut

# }}}

# {{{ POD Options

=head2 Options

=over 2

=item -help

Displays help and exits. All other arguments are ignored.

=item -outfile

If specified, output is written to this file. Otherwise it
is written to STDOUT. This is affected by the -outdir option
described below.

=item -outdir

The directory in which output files will be placed. If this is
specified without -outfile then the output filenames are derived
from input filenames and placed in this directory.

If this directory does not exist then an attempt is made to make
it. Failure to make this directory is a fatal error (croak is called).

If -outdir is specified with -outfile then the outfile is placed
in this directory.

=item -extension

By default this ($outex) is undefined. This is the extension to use
when output filenames are derived from input filenames. 

=back

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

# {{{ Open the outfile. Defaults to STDOUT.
my $ofh;
if($ofn) {
  open($ofh, ">", $ofn);
}
else {
  open($ofh, ">&STDOUT");
}
select($ofh);
# }}}

# Populate @infiles
my @infiles = @ARGV;


# {{{ Cycle through all the infiles as text files.
for my $infile (@infiles) {
my ($noex, $dir, $ext)= fileparse($infile, qr/\.[^.]*/);
my $bn = $noex . $ext;
# tablistE($infile, $bn, $noex, $ext);
$noex =~ s/\.sorted\.depth$/_GPE/;

linelist("track type=bedGraph name=$noex description=$noex");
open(my $ifh, "<$infile") or croak("Could not open $infile");
my $lineCnt = 0;
if($skip) {
for (1..$skip) { my $discard = readline($ifh); }
}
# local $/ = ""; # For reading multiline records separated by blank lines.

my $apertot = 0;
my $start = 1;

while(my $record = readline($ifh)) {
chomp($record);
if($record=~m/^\s*\#/ or $record=~m/^\s*$/) {next;}
my @ll=split(/\t/, $record);
my $chr = $ll[0];
my $gps = $ll[1];
my $val = $ll[2];

if($gps % $aper) {
$apertot += $val;
}

if($gps % $aper == 0) {
tablist($chr, $start, $gps, int($apertot/$aper));
$start = $gps + 1;
$apertot = 0;
}

$lineCnt += 1;
if($testCnt and $lineCnt >= $testCnt) { last; }
}
close($ifh);
# Uncomment line below for ifn derived ofn.
# close($idoh);
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

# track type=bedGraph name=track_label description=center_label
#     visibility=display_mode color=r,g,b altColor=r,g,b
#     priority=priority autoScale=on|off alwaysZero=on|off gridDefault=on|off
#     maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper
#     yLineMark=real-value yLineOnOff=on|off
#     windowingFunction=maximum|mean|minimum smoothingWindow=off|2-16
