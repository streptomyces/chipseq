#!/usr/bin/perl
use 5.14.0;
use utf8;
use Carp;
use lib qw(/home/sco /home/sco/perllib);
use File::Basename;
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
    tablistV tablistVE linelistV linelistVE tablistH linelistH
    tablistER tablistVER linelistER linelistVER tabhashER tabhashVER);
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
my $cpus = 8;
our $verbose;
my $skip = 0;
my $help;
GetOptions (
"outfile:s" => \$outfile,
"outdir:s" => \$outdir,
"indir:s" => \$indir,
"fofn:s" => \$fofn,
"cpus:i" => \$cpus,
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

# {{{ POD Example

=head1 Name

Change me.

=head2 Example

 perl changeme.pl -outfile out.txt -- inputfile1 inputfile2

Note that input files are always specified as non-option arguments.

# }}}

# {{{ POD Options and blurb

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

=head2 Blurb

Uses F<Sco::Common> for a variety of printing functions.

If neither -outfile nor -outdir are specified then the output
is to STDOUT.

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

opendir(DIR, "bedgraph");
while(my $bn = readdir(DIR)) {
  my $file = File::Spec->catfile("bedgraph", $bn);
  if($file =~ m/mean/) { next; }
  if(-f $file and $file =~ m/BldC-.*\d{2}A.*\.bedgraph/) {
    my $read1 = $file;
    my $read2 = $read1; $read2 =~ s/A/B/;
    my $outfn = $read1; $outfn =~ s/_ln/_ln_mean/;
    $outfn =~ s/A//;
    my $xstr = qq(perl code/bedgraph_mean.pl -outfile $outfn -- $read1 $read2 );
    linelist($xstr);
  }
}

# perl code/bedgraph_mean.pl -outfile stuff -- \
# bedgraph/aBldC-14A-IP_ln.bedgraph bedgraph/aBldC-14B-IP_ln.bedgraph



exit;

# Multiple END blocks run in reverse order of definition.
END {
# close($ofh);
close(STDERR);
close(ERRH);
# $handle->disconnect();
}

