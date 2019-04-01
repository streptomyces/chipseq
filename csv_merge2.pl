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
use DBI;

# {{{ Getopt::Long
use Getopt::Long;
my $outdir;
my $indir;
my $fofn;
my $distThresh = 500;
my $sqlitefn = qq(/home/sco/sqlite/vnz.sqlite);
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
"distthresh:i" => \$distThresh,
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

my $handle = DBI->connect("DBI:SQLite:dbname=$sqlitefn", '', '');

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


my @all;
my %alh;
my @noex;

# {{{ Cycle through all the infiles.
my $filecnt = 0;
for my $infile (@infiles) {
my %done;
my ($noex, $dir, $ext)= fileparse($infile, qr/\.[^.]*/);
push(@noex, $noex);
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
# push(@head, "Sample", @hl); 

while(my $line = readline($ifh)) {
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @ll=split(/\t/, $line);
my $pos = $ll[0];
my $lfc = $ll[2];
my $lcpm = $ll[3];
my $fdr = $ll[4];

my $key = $noex . "_" . $pos;
$alh{$key} = [$lfc, $lcpm, $fdr];
push(@all, [$noex, $pos, $lfc, $lcpm, $fdr]);


}
close($ifh);
close(OFH);
  $filecnt += 1;
}
# }}}


my @head = qw(position lfc10 lfc14 lfc22 minFDR closestgene lgene
lproduct ldist igene iproduct idist rgene rproduct rdist);
tablist(@head);

my @donepos;
my @donein;
my @donegene = ("-");

my $outCnt = 0;
### @all is like [$noex, $pos, $lfc, $lcpm, $fdr]
my @sortedall = sort {$b->[2] <=> $a->[2]} @all;
for my $rr (@sortedall) {
my $pos = $rr->[1];

if(pushifnotN(\@donepos, $pos)) {
my($lhr, $ihr, $rhr) = lir($pos);

my $ldist = $pos - $lhr->{end_pos};
my $rdist = $rhr->{start_pos} - $pos;
my $idist;
if($ihr->{id} eq '-') {
$idist = '-';
}
elsif($ihr->{strand} == -1) {
$idist = abs($ihr->{end_pos} - $pos);
}
elsif($ihr->{strand} == 1) {
$idist = abs($pos - $ihr->{start_pos});
}
$lhr->{dist} = $ldist;
$ihr->{dist} = $idist;
$rhr->{dist} = $rdist;
my $mindisthr = mindist($lhr, $ihr, $rhr);

my $gdr = genesDone($lhr, $ihr, $rhr);
# tablistE("gdr:", $gdr);
if($gdr) {
next;
}

my @lfc;
for my $noex (@noex) {
my $key = $noex . "_" . $pos;
if(exists($alh{$key})) {
push(@lfc, $alh{$key}->[0]);
}
else { push(@lfc, "-"); }
}


my @lout;
if($ldist <= $distThresh) {
push(@lout, $lhr->{id}, $lhr->{product}, $ldist);
  pushifnot(\@donegene, $lhr->{id});
}
else { @lout = ("-", "-", "-"); }

my @rout;
if($rdist <= $distThresh) {
push(@rout, $rhr->{id}, $rhr->{product}, $rdist);
  pushifnot(\@donegene, $rhr->{id});
}
else { @rout = ("-", "-", "-"); }

my @iout;
if($idist <= $distThresh) {
push(@iout, $ihr->{id}, $ihr->{product}, $idist);
  pushifnot(\@donein, $ihr->{id});
}
else { @iout = ("-", "-", "-"); }


tablist($pos, @lfc, $rr->[4], $mindisthr->{id},
@lout, @iout, @rout
);

$outCnt += 1;
if($testCnt and $outCnt >= $testCnt) { last; }
}
}


# linelist(@donegene);

exit;

# Multiple END blocks run in reverse order of definition.
END {
close($ofh);
close(STDERR);
close(ERRH);
$handle->disconnect();
}

sub mindist {
  my $retobj;
  my $min = 99999;
  for my $hr (@_) {
    if($hr->{dist} ne '-') {
      if($hr->{dist} < $min) {
        $retobj = $hr;
        $min = $hr->{dist};
      }
    }
  }
  return($retobj);
}


# {{{ sub genesDone
sub genesDone {
  my ($lhr, $ihr, $rhr) = @_;
  my @ids = ($lhr->{id}, $ihr->{id}, $rhr->{id});
  if($lhr->{dist} > $distThresh) { $ids[0] = '-'; }
  if($ihr->{dist} > $distThresh) { $ids[1] = '-'; }
  if($rhr->{dist} > $distThresh) { $ids[2] = '-'; }
  my $retval = 1; # 1 means done.
    if(
        inlist(\@donegene, $ids[0])
        and (inlist(\@donegene, $ids[1]) or inlist(\@donein, $ids[1]))
        and inlist(\@donegene, $ids[2])
      ) { return(1); }
    else {
      return(0);
    }
}
# }}}


# {{{ pushifnot and pushifnotN
sub pushifnot {
my $lr = shift(@_);
my $value = shift(@_);
unless( grep {$_ eq $value} @{$lr} ) {
push(@{$lr}, $value);
return(1);
}
return(0);
}
sub pushifnotN {
my $lr = shift(@_);
my $value = shift(@_);
unless( grep {$_ == $value} @{$lr} ) {
push(@{$lr}, $value);
return(1);
}
return(0);
}
# }}}


# {{{ sub lir
sub lir {
  my $pos = shift(@_);
# left gene
my $qstr = qq/select max(end_pos) from features/;
$qstr .= qq/ where end_pos <= $pos and strand == -1/;
$qstr .= qq/ and (pritag = 'CDS' or pritag like '%rna%')/;
my ($lendpos) = $handle->selectrow_array($qstr);
$qstr = qq/select * from features where end_pos = $lendpos and strand = -1/;
$qstr .= qq/ and (pritag = 'CDS' or pritag like '%rna%')/;
my $stmt = $handle->prepare($qstr);
$stmt->execute();
my $lhr = $stmt->fetchrow_hashref();
$stmt->finish();

# right gene
$qstr = qq/select min(start_pos) from features/;
$qstr .= qq/ where start_pos >= $pos and strand == 1/;
$qstr .= qq/ and (pritag = 'CDS' or pritag like '%rna%')/;
my ($rstartpos) = $handle->selectrow_array($qstr);
$qstr = qq/select * from features where start_pos = $rstartpos and strand = 1/;
$qstr .= qq/ and (pritag = 'CDS' or pritag like '%rna%')/;
$stmt = $handle->prepare($qstr);
$stmt->execute();
my $rhr = $stmt->fetchrow_hashref();
$stmt->finish();

# in gene
$qstr = qq/select start_pos, end_pos, strand from features/;
$qstr .= qq/ where start_pos < $pos and end_pos > $pos/;
$qstr .= qq/ and (pritag = 'CDS' or pritag like '%rna%')/;
my ($istartpos, $iendpos, $istrand) = $handle->selectrow_array($qstr);
my $ihr;
if ($istartpos and $iendpos and $istrand) {
$qstr = qq/select * from features where start_pos = $istartpos and strand = $istrand/;
$qstr .= qq/ and end_pos = $iendpos/;
$qstr .= qq/ and (pritag = 'CDS' or pritag like '%rna%')/;
$stmt = $handle->prepare($qstr);
$stmt->execute();
$ihr = $stmt->fetchrow_hashref();
$stmt->finish();
}
unless(ref($lhr)) {
$lhr->{id} = '-';
$lhr->{product} = '-';
}
unless(ref($ihr)) {
$ihr->{id} = '-';
$ihr->{product} = '-';
}
unless(ref($rhr)) {
$rhr->{id} = '-';
$rhr->{product} = '-';
}

return($lhr, $ihr, $rhr);

}
# }}}

sub inlist {
my $lr = shift(@_);
my $val = shift(@_);
if(grep {$_ eq $val} @{$lr}) {
return(1);
}
else { return(0); }
}
