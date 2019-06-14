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
use List::Util qw(any);

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
my $dbfile = '/home/sco/sqlite/avin.sqlite';
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
"sqlitefile|dbfile:s" => \$dbfile,
"runfile:s" => \$runfile,
"testcnt:i" => \$testCnt,
"skip:i" => \$skip,
"verbose" => \$verbose,
"help" => \$help
);
# }}}

# {{{ POD markdown

=head2 Name

lir.pl

=head2 Examples

perl code/lir.pl -outfile lir.txt -- siglfc.txt


=cut

# }}}

if($help) {
exec("perldoc $0");
exit;
}
my $handle = DBI->connect("DBI:SQLite:dbname=$dbfile", '', '');

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
    my $ofn = File::Spec->catfile($outdir, $noex . "_out" . $outex);
    open(OFH, ">", $ofn) or croak("Failed to open $ofn");
    select(OFH);
  }


  open(my $ifh, "<$infile") or croak("Could not open $infile");
  my $lineCnt = 0;
  if($skip) {
    for (1..$skip) { my $discard = readline($ifh); }
  }
  my $header = readline($ifh); chomp($header);
  my @header = split(/\t/, $header);
  my @genehead = qw/LeftGene LeftDistance LeftProduct/;
  push(@genehead, (qw/InGene InDistance InProduct/));
  push(@genehead, (qw/RightGene RightDistance RightProduct/));
  tablist(@header, @genehead);
  my @donelir;
# local $/ = ""; # For reading multiline records separated by blank lines.
  while(my $line = readline($ifh)) {
    chomp($line);
    if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
    my @ll=split(/\t/, $line);
    my ($pos, $enrich, $apv) = @ll;
    # unless($enrich >= 1.5) { next; }
    my %lir = lir($pos, 0); # The second argument is boolean for direction.

    my (%left, %right, %in);
    if(ref($lir{left})) {
      %left = %{$lir{left}};
    }
    if(ref($lir{right})) {
      %right = %{$lir{right}};
    }
    if(ref($lir{in})) {
      %in = %{$lir{in}};
    }


    my @out;
    my @strlir;
    push(@out, split("", "-" x 9));
    if(%left and $left{strand} == -1) {
      my $dist = dist($pos, \%left);
      splice(@out, 0, 3, $left{id}, $dist, $left{product});
      push(@strlir, $left{id});
    }
    if(%in) {
      my $dist = dist($pos, \%in);
      splice(@out, 3, 3, $in{id}, $dist, $in{product});
      push(@strlir, $in{id});
    }
    if(%right and $right{strand} == 1) {
      my $dist = dist($pos, \%right);
      splice(@out, 6, 3, $right{id}, $dist, $right{product});
      push(@strlir, $right{id});
    }
    my $strlir = join("-", @strlir);
    if(any {$_ eq $strlir} @donelir) {
    }
    else {
      push(@donelir, $strlir);
      tablist(@ll, @out);
    }
    $lineCnt += 1;
    if($testCnt and $lineCnt >= $testCnt) { last; }
    if($runfile and (not -e $runfile)) { last; }
  }
  close($ifh);
  close(OFH);
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


# {{{ sub lir
sub lir {
  my $pos = shift(@_);
  my $direction = shift(@_);
  my %lir;
# {{{ left
  my $qstr = qq/select max(end_pos) from features where end_pos < $pos/;
  $qstr .= qq/ and pritag = 'CDS'/;
  if($direction) {
    $qstr .= qq/ and strand = -1/;
  }
  my ($left_end) = $handle->selectrow_array($qstr);
  if($left_end) {
    $qstr = qq/select id from features where end_pos = $left_end/;
    $qstr .= qq/ and pritag = 'CDS'/;
    my ($left_id) = $handle->selectrow_array($qstr);
    if($left_id) { 
      $qstr = qq/select id, locus_tag, pritag, start_pos, end_pos,/;
      $qstr .= qq/ strand, product, olt/;
      $qstr .= qq/ from features where id = '$left_id'/;
      my $stmt = $handle->prepare($qstr);
      $stmt->execute();
      my $lr = $stmt->fetchrow_hashref();
      $lir{left} = $lr ;
    }
  }
# }}}

# {{{ right
  $qstr = qq/select min(start_pos) from features where start_pos > $pos/;
  $qstr .= qq/ and pritag = 'CDS'/;
  if($direction) {
    $qstr .= qq/ and strand = 1/;
  }
  my ($right_start) = $handle->selectrow_array($qstr);
  if($right_start) {
    $qstr = qq/select id from features where start_pos = $right_start/;
    $qstr .= qq/ and pritag = 'CDS'/;
    my ($right_id) = $handle->selectrow_array($qstr);
    if($right_id) { 
      $qstr = qq/select id, locus_tag, pritag, start_pos, end_pos,/;
      $qstr .= qq/ strand, product, olt/;
      $qstr .= qq/ from features where id = '$right_id'/;
      my $stmt = $handle->prepare($qstr);
      $stmt->execute();
      my $lr = $stmt->fetchrow_hashref();
      $lir{right} = $lr ;
    }
  }
# }}}

# {{{ in
  $qstr = qq/select id, locus_tag, pritag, start_pos, end_pos,/;
  $qstr .= qq/ strand, product, olt/;
  $qstr .= qq/ from features where start_pos <= $pos/;
  $qstr .= qq/ and end_pos >= $pos/;
  $qstr .= qq/ and pritag = 'CDS'/;
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  my $lr = $stmt->fetchrow_hashref();
  $lir{in} = $lr ;
# }}}
  return(%lir);
}
# }}}

sub dist {
my $pos = shift(@_);
my $gr = shift(@_);
my ($codon_start, $dist);
if($gr->{strand} == -1) {
$codon_start = $gr->{end_pos};
$dist = $pos - $codon_start;
}
else {
$codon_start = $gr->{start_pos};
$dist = $codon_start - $pos;
}
return($dist);
}


__END__


my @values = ();
foreach my $key (sort keys(%{$hr})) {
push(@values, $hr->{$key});
}
print(++$serial,"\t",join("\t", @values), "\n");

}
$stmt->finish();

### selectrow_array ###
my @row = $handle->selectrow_array($qstr);

# Temporary tables. Start a transaction.
$handle->begin_work();

 my($tmpfh, $table) = tempfile($template, DIR => undef, SUFFIX => undef);
 close($tmpfh); unlink($table);
 linelistE($table);

# Create the temporary table as usual but with the additional clause
# "on commit drop".

 $handle->do(qq/create temporary table $table (file text, locus_tag text,
  bB1 text, bB2 text, bC1 text, bC2 text)
  on commit drop/);

# Now use the table however you wish and when you are done.
 $handle->commit();

$handle->disconnect();



