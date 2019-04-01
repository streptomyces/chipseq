#!/usr/bin/perl
use 5.14.0;
use Carp;
use File::Copy;
use File::Basename;
use lib qw(/home/sco /home/sco/perllib);
use File::Temp qw(tempfile tempdir);
my $tempdir = qw(/home/sco/volatile);
my $template="rcmdXXXXX";

# {{{ Getopt::Long stuff
use Getopt::Long;
my $ipfile;
my $totfile;
my $outfile;
my $diffout;
my $testCnt = 0;
my $region = 1000;      # Even number please.
my $sdthresh = 3.0;
my $pvthresh = 1e-4;
my $lnoutdir = q(lngr);  # defaults to ./lngr.
GetOptions (
"ipfile=s" => \$ipfile,
"totfile=s" => \$totfile,
"outfile:s" => \$outfile,
"diffout=s" => \$diffout,
"testcnt:i" => \$testCnt,
"sdthresh:f" => \$sdthresh,
"pvthresh:f" => \$pvthresh,
"regionwidth:i" => \$region,
"lnoutdir:s" => \$lnoutdir
);


if($outfile) {
open(OUT, ">$outfile");
}
else {
  open(OUT, ">&STDOUT");
}
# }}}

if($lnoutdir) {
$lnoutdir =~ s/\/$//;   # Remove trailing slash from $lnoutdir.
unless(-d $lnoutdir) {
mkdir $lnoutdir or croak "Failed to make directory $lnoutdir\n";
}
}

if($region % 2) {
croak("Region width specified ($region) in not an even number");
}

my ($basename, $directory, $ext)= fileparse($ipfile, qr/\.[^.]*/);

# {{{ constants
my $halfReg = $region/2;
my $aperture = 50;      # Even number please.
my $halfAp = $aperture/2;
# }}}

my %ln;
my %ln1;

# {{{ read the files and populate %ln and %ln1
FILE: for my $ifn ($ipfile, $totfile) {
my ($basename, $directory, $ext)= fileparse($ifn, qr/\.[^.]*/);
my $suffix;
unless($basename =~ m/_$/) { $suffix = '_'; }
my $lnOut = $lnoutdir . "/" . $basename . $suffix . 'ln.gr';
open(my $lnfh, ">", $lnOut);
my @cov;
open(INFILE, "<$ifn") or croak("Could not open $ifn");
my $lineCnt = 0;
while(<INFILE>) {
my $line=$_;
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @llist=split(/\t/, $line);
push(@cov, $llist[1]);
}
close(INFILE);
print(STDERR "$ifn read\n");

my $lastdx = $#cov;
my $pos = $halfAp;
while($pos <= ($#cov-$halfAp)) {
my $regLeft = ($pos - $halfReg) >= 0 ? $pos-$halfReg : 0;
my $regRight = ($regLeft + $region) <= $#cov ? $regLeft+$region : $#cov;

my $apLeft = $pos - $halfAp;
my $apRight = $pos + $halfAp;

my $apCov = 0;
foreach my $temp ($apLeft..$apRight) {
$apCov += $cov[$temp];
}
my $regCov = 0;
foreach my $temp ($regLeft..$regRight) {
$regCov += $cov[$temp];
}

my $fracInAp;
if($regCov > 0) {
$fracInAp = $apCov / $regCov;
}
else { $fracInAp = 0; }
my $intFrac = $fracInAp * $region;
# the use of $region above is just as a multiplier.
# It is not related to the length of the region in
# any manner. Which means, it can be changed to any
# other value without any other consequences.

# Note that $intFrac is used only to print out to the FH below.
# For further calculations we use $fracInAp.
print($lnfh "$pos\t$intFrac\n");

if($ifn eq $ipfile) {
$ln{$pos} = $fracInAp;
}
elsif($ifn eq $totfile) {
$ln1{$pos} = $fracInAp;
}
$pos += $halfAp;
}
close($lnfh);
}
# }}}


my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.csv'); # input data for R
for my $pos (sort {$a<=>$b} keys(%ln)) {
#my $diff = ($ln{$pos} - $ln1{$pos}) > 0 ? $ln{$pos} - $ln1{$pos} : 0;
my $diff = ($ln{$pos} - $ln1{$pos});
print($fh "$pos\t$diff\n");
}
close($fh);
print(STDERR "Difference data written to $fn\n");
copy($fn, $diffout);

my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.R'); # R commands

my($fh2, $fn2)=tempfile($template, DIR => $tempdir, SUFFIX => '.csv'); # output from R
close($fh2);

# {{{ R code
print $fh1 <<"AAR";
cn <- c("pos", "cov");
temp<-read.table("$fn", header=FALSE, stringsAsFactors=FALSE, sep="\\t", col.names = cn);
df<-temp;
covsd <- sd(df\$cov);
covbar <- mean(df\$cov);
pvals <- pnorm(df\$cov, sd = covsd, mean = covbar, lower.tail = FALSE);
df\$pval <- pvals;

#thresh <- covbar + (covsd * $sdthresh);
#plus <- df[df\$cov >= thresh, ];
#minus <- df[df\$cov < -thresh, ];
#iframe <-rbind(plus, minus);

iframe <- df[df\$pval <= $pvthresh, ];

ordvec <- order(iframe\$pos);
jframe <- iframe[ordvec, ];
outf <- data.frame(pos = jframe\$pos, cov = jframe\$cov);
write.table(df, file="/home/sco/temp/dfpv.csv", sep="\\t", quote=FALSE, col.names=TRUE, row.names=FALSE);
write.table(outf, file="$fn2", sep="\\t", quote=FALSE, col.names=FALSE, row.names=FALSE);
AAR
close($fh1);
print(STDERR "R code written to $fn1\n");
# }}}

### Call to R below.
`/usr/bin/R CMD BATCH --vanilla --silent --slave $fn1 /dev/null`;
copy($fn1, "thresh.R");
unlink($fn, $fn1);

open(ROUT, "<$fn2");
while(<ROUT>) {
# do something with the lines here
  print(OUT $_);
}
close(ROUT);
unlink($fn2);
close(OUT);

__END__

=head1 Name

localNormalPairSig.pl

=head2 Example

    for ipfn in $(ls coverage/*_IP.gr);
      do
    bn=$(basename $ipfn _IP.gr);
    tfn="coverage/"$bn"_T.gr";
    dfn="IP-Tdiff/"$bn"_diff.gr";
    sfn="sigDiff/"$bn"_sigDiff.gr";
        # echo $ipfn $tfn $dfn $sfn
        echo "Starting $bn at " $(date);
        perl code/localNormalPairSig.pl -pvth 1e-4 -region 4000 \
        -ip $ipfn -tot $tfn -out $sfn -lnoutdir lngr \
        -diff $dfn
        echo "Finished $bn at " $(date);
        echo;
        done

=head2 Comments

The names of the local normalised .gr files are decided by the code
based on the names of the input files.

There is a call to R CMD BATCH indicated by a comment

 Call to R below.

The actual call looks like.
 
 /usr/bin/R CMD BATCH --vanilla --silent --slave $fn1 /dev/null

=cut


