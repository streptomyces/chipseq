#!/usr/bin/perl
use 5.14.0;
use Carp;
use lib qw(/home/sco /home/sco/perllib);
use File::Basename;
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals);
use Sco::Genbank;
my $gbk = Sco::Genbank->new();


# {{{ Usage message

my $usage = <<"USAGE";

 -infile:   Input file to read data from.

 -outfile:  File to which output will be written.
            If this is not specified output is to STDOUT.

 -append:   Boolean. If true, the output file is opened
            for appending.

 -testcnt:  May be used to stop before processing
            all the data in the input file. For testing.

 -featlstr: List of features to include. Comma separated, no spaces.

 -molecule: String to be used in the first column of the output file.

USAGE

# }}}

# {{{ Getopt::Long stuff
use Getopt::Long;
my $infile;
my $outfile;
my $featlstr;
my $tagasid = qq(locus_tag);
my $molecule;
my $testCnt = 0;
my $append = 0;
my $verbose;
my $skip = 0;
my $help;
GetOptions (
"infile|filename=s" => \$infile,
"outfile=s" => \$outfile,
"tagasid:s" => \$tagasid,
"features:s" => \$featlstr,
"molecule=s" => \$molecule,
"testcnt:i" => \$testCnt,
"skip:i" => \$skip,
"append" => \$append,
"verbose" => \$verbose,
"help" => \$help
);

if((not $infile) or $help) {
print($usage);
exit;
}
# }}}

# CDS, tRNA, rRNA are done by default. If you need any other
# tags included use the features option.

my @features;
if($featlstr) {
@features = split(/,/, $featlstr);
}
if(@features) {
$gbk->saffile(file => $infile, molecule => $molecule, tagasid => $tagasid, 
              outfilename => $outfile, features => [@features]);
}
else {
$gbk->saffile(file => $infile, molecule => $molecule, tagasid => $tagasid,
              outfilename => $outfile);
}

exit;


# }}}




