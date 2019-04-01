use 5.14.0;
use Carp;
use utf8;
use lib qw(/home/sco /home/sco/perllib);
use File::Basename;
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
    tablistV tablistVE linelistV linelistVE tablistH linelistH
    tablistER tablistVER linelistER linelistVER tabhashER tabhashVER csvsplit);
use File::Spec;

my $indir = $ARGV[0];
my $outdir = $ARGV[1];

unless(-d $outdir) {
  my $mdret = mkdir($outdir);
  unless($mdret) {
    croak("Failed to make $outdir for output files");
  }
}

opendir(DIR, $indir);
while(my $fn = readdir(DIR)) {
  my $ipath = File::Spec->catfile($indir, $fn);
  if(-f $ipath and $fn =~m/\.depth$/) {
    my ($noex, $dir, $ext)= fileparse($ipath, qr/\.[^.]*/);
    my $bn = $noex . $ext;
    $noex =~ s/\.sorted//;
    my $on = $noex . ".bedgraph";
    my $opath = File::Spec->catfile($outdir, $on);
    linelist("perl code/depthToBedgraph.pl -outfile $opath -- $ipath")
  }
}

close(DIR);
exit;

