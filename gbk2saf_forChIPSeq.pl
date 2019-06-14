#!/usr/bin/perl
use 5.14.0;
use Carp;
use lib qw(/home/sco /home/sco/perllib);
use File::Basename;
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals);
use Sco::Genbank;
my $gbk = Sco::Genbank->new();

# {{{ POD

=head1 Name

gbk2saf_forChIPSeq.pl

=head2 Examples

 perl code/gbk2saf_forChIPSeq.pl -outfile Avin_chipseq.saf -molecule Avin \
 -prefix Avin -- Avin.gbk


 -outfile:  File to which output will be written.
            If this is not specified output is to STDOUT.

 -molecule: String to be used in the first column of the output file.

=cut


# }}}

# {{{ Getopt::Long stuff
use Getopt::Long;
my $outfile;
my $molecule;
my $prefix;
my $by = 30;
my $help;
GetOptions (
"outfile=s" => \$outfile,
"molecule=s" => \$molecule,
"prefix=s" => \$prefix,
"by:i" => \$by,
"help" => \$help
);

my $infile = shift(@ARGV);


if((not $infile) or $help) {
exec("perldoc $0");
exit;
}
# }}}

$gbk->chipseq_saffile(file => $infile, molecule => $molecule, prefix => $prefix,
    by => $by, outfilename => $outfile);

exit;


# }}}




