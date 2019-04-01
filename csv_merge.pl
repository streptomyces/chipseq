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
"dirout:s" => \$outdir,
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

my @rowid;

my @row;
my @head;

# {{{ Cycle through all the infiles.
my $filecnt = 0;
for my $infile (@infiles) {
my %done;
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
my $header = readline($ifh);
chomp($header);
my @hl = split(/\t/, $header);
push(@head, "Sample", @hl); 

my @forsort;
while(my $line = readline($ifh)) {
chomp($line);
push(@forsort, $line);
}
my @sortedLines = midgene(@forsort);
tablistE("forsort", scalar(@forsort));
tablistE("sortedLines", scalar(@sortedLines));



# local $/ = ""; # For reading multiline records separated by blank lines.
for my $line (@sortedLines) {
# chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @ll=split(/\t/, $line);
if(exists($done{$ll[1]}) or $ll[1] =~ m/NA.NA.NA/) { next; }
$done{$ll[1]} = 1;
push(@{$row[$filecnt]}, [$noex, @ll]);
unless(grep {$_ eq $ll[1]} @rowid) {
  push(@rowid, $ll[1]);
}
$lineCnt += 1;
if($testCnt and $lineCnt >= $testCnt) { last; }
if($runfile and (not -e $runfile)) { last; }
}
close($ifh);
close(OFH);
  $filecnt += 1;
}
# }}}

tablist("Gene", @head);

my @redid = reducelist(@rowid);


my @merged;
for my $rowid (@redid) {
  # linelistE($rowid);
  my @out;
  for my $lr (@row) {
    my $pushed = 0;
    for my $rr (@{$lr}) {
      if($rr->[2] eq $rowid) {
        push(@out, @{$rr});
        $pushed = 1;
      }
    }
    unless($pushed) {
      for(0..14) {
        push (@out, "-");
      }
    }
  }
  push(@merged, \@out);
}

my @done;
for my $rr (@merged) {
  my $gstr = join(":", $rr->[2], $rr->[17], $rr->[32]);
  my @genes;
  for my $gene (split(/:/, $gstr)) {
    unless(grep {$_ eq $gene} @genes
        or $gene eq 'NA' or $gene eq '-') { push(@genes, $gene); }
  }
  for my $gene (@genes) {
    unless( grep {$_ eq $gene} @done ) {
      tablist($gene, @{$rr});
      push(@done, $gene);
    }
  }
# tablist(@gstr);
}





exit;

# Multiple END blocks run in reverse order of definition.
END {
close($ofh);
close(STDERR);
close(ERRH);
# $handle->disconnect();
}

sub midgene {
  my @lines = @_;
  my (@mg, @sg0, @sg1, @sg2, @sg3);
  for my $line (@lines) {
    my @ll = split(/\t/, $line);
    my @gpl = split(/:/, $ll[1]);
    if($line =~ m/NA:NA/) {
      push(@sg0, $line);
    }
    if($line =~ m/:NA:/) {
      push(@sg1, $line);
    }
    elsif($gpl[0] eq 'NA' and $gpl[2] eq 'NA')  {
      push(@mg, $line);
    }
    elsif($line =~ m/NA/) {
      push(@sg2, $line);
    }
    else {
      push(@sg3, $line);
    }
  }
  my @retlist = (@sg0, @sg1, @sg2, @sg3, @mg);
  return(@retlist);
}

sub reducelist {
  my @rowid = @_;
  my @retlist;
  my @seen;
  for my $rowid (@rowid) {
    my @gpl;
    my @temp = split(/:/, $rowid);
    for my $gene (@temp) {
      unless($gene eq 'NA') { push @gpl, $gene; }
    }
    my $allseen = 1;
    for my $gene (@gpl) {
      if( grep {$_ eq $gene} @seen) { }
      else {
        $allseen = 0;
        push(@seen, $gene);
      }
    }
    if($allseen) { }
    else {
      push(@retlist, $rowid);
    }
  }
  tablistE("Unreduced:", scalar(@rowid));
  tablistE("Reduced:", scalar(@retlist));
  return(@retlist);
}
