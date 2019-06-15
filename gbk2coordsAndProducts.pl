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
use Bio::SeqIO;

# {{{ Getopt::Long
use Getopt::Long;
my $indir;
my $conffile = qq(local.conf);
my $errfile;
my $ofn;
my $testCnt = 0;
our $verbose;
my $skip = 0;
my $help;
GetOptions (
"outfile:s" => \$ofn,
"indir:s" => \$indir,
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

my %lut;

open(VNZ, "<", "vnzNames.list");
while(<VNZ>) {
chomp($_);
my @ll = split(/\t/, $_);
my $vlt = shift(@ll);
if($vlt =~ m/^vnz_/) {
$lut{$vlt} = \@ll;
}
}
close(VNZ);


# Populate @infiles
my @infiles = @ARGV;


my @pritags_todo = qw(CDS tRNA rRNA);
my @tag_as_ids = qw(locus_tag gene);
tablist(qw(pritag start end strand product names));

# {{{ Cycle through @infiles as Bio::SeqIO.
for my $infile (@infiles) {
my ($noex, $dir, $ext)= fileparse($infile, qr/\.[^.]*/);
my $bn = $noex . $ext;
my $seqio=Bio::SeqIO->new(-file => $infile);
my $seqout=Bio::SeqIO->new(-fh => $ofh, -format => 'fasta');

while(my $seqobj=$seqio->next_seq()) {
  foreach my $feature ($seqobj->all_SeqFeatures()) {
    my $start = $feature->start();
    my $end = $feature->end();
    my $strand = $feature->strand();
    my $pritag = $feature->primary_tag();
    if( grep { $_ eq $pritag } @pritags_todo ) {
      my $product;
      my $id;
      my $otherids;
      my @tags = $feature->get_all_tags();
      foreach my $tag (@tags) {
        if($tag=~m/product/i) {
          my @tagstrs = $feature->get_tag_values($tag);
          $product = join("; ", @tagstrs);
        }
        if(grep { $_ eq $tag } @tag_as_ids) {
          my @temp = $feature->get_tag_values($tag);
          $id = shift(@temp);
          if(@temp) { $otherids = join(", ", @temp); }
        }
      }
      if($otherids) { $product .= $otherids; }
      unless($product) { $product = $pritag; }
      my $names = "-";
      if(exists($lut{$id})) {
        my @temp;
        for my $vn (@{$lut{$id}}) {
          unless($vn eq "-") { push(@temp, $vn); }
        }
        $names = join(", ", @temp);
      }
      tablist($id, $pritag, $start, $end, $strand, $product, $names);
    }
  }
# $seqout->write_seq($seqobj);
}
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

