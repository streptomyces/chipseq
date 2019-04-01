# ChIP-Seq Jul 2018

<!-- vim: set tw=70 nosmartindent: -->

Matt Hutchings emailed on 19 July

    Hi Govind.
    I hope you’re well?
    Rebecca Devine (CC) is a PhD student with me and Barrie (50:50 split) and
    working on the formicamycins project with Qin. She’s done a ChIP
    experiment and the seq data just came in (see below). I was hoping you
    could run it through your pipeline and generate files we can view in IGB?
    The genome is Streptomyces formicae KY5 which is complete and on StrepDB
    and reported here:
    https://www.sciencedirect.com/science/article/pii/S0168165617317443
    If you don’t have time that’s fine, we’ll find another way. I’m sure
    you’re swamped as always.
 
    Thanks.
    Matt
    Begin forwarded message:

     From: "GENEWIZ Next Generation Sequencing Services" <NGS@genewiz.com>
     To: "GENEWIZ Next Generation Sequencing Services" <NGS@genewiz.com>,
     "Matthew Hutchings (BIO - Staff)" <M.Hutchings@uea.ac.uk>
     Cc: "Amy Drummond, GW/UK" <Amy.Drummond@genewiz.com>, "Rebecca Devine
     (BIO - Student)" <R.Devine@uea.ac.uk>
     Subject: RE: Update for Project Quote: KM1711222_R1

     Hello Matt,



     I hope all is well. Your project KM1711222_R1 has completed and was
     uploaded to your account on our sFTP server. Below are your login
     credentials:



     Host:     sftp://ftp3.genewiz.com

     User:     matt.hutchings

     Password:      MjM4ODk4YWFhNTYzZjNiNGIxMzUwMjM2

     Port:      22


#### Fri 20 Jul 2018

Data downloaded into

/mnt/isilon/customers/matt_hutchings/2018_07_20/fromContractor/

~~~ {.sh}
cd fromContractor
md5sum -c md5sum_list.txt
cd ..
~~~

Above, all files checksum OK.

GCA_002556545.1_ASM255654v1_genomic.gbff downloaded from Genbank on
Mon 23 Jul 2018 using a command line FTP client.

~~~ {.sh}
perl ~/code/perl/gbk2fna.pl -infile refseq/ky5.gbk -outfile \
refseq/ky5.fna

cd /mnt/isilon/bowtie2Indexes/ky5
bowtie2-build ky5.fna ky5
cd -

perl code/comgen.pl | parallel --jobs 1

perl ~/code/perl/seq/gbk2bed_forChIPSeq.pl -infile refseq/ky5.gbk \
-prefix KY5 -molecule KY5 -by 30 -outfile ky5_chip.bed


pushd sam
parallel --jobs 4 samtools view -@ 3 -b -o {.}.bam {} \
::: $(ls --color=never ./*.sam)

parallel --jobs 4 samtools sort -@ 3 -o {.}.sorted.bam {} \
::: $(ls --color=never ./*.bam)

parallel --jobs 12 samtools index {} \
::: $(ls --color=never ./*.sorted.bam)
pushd


reportcollect () {
for repo in $(ls --color=never bowtie2reports/*bowtie2); do
  echo $repo
  cat $repo
  echo
done
}

reportcollect > alignment_reports.txt

parallel --dry-run samtools depth -a {} \> depth/{/.}.depth \
::: $(ls --color=never sam/*.sorted.bam)

perl code/bedgraph_comgen.pl depth bedgraph

~~~

#### Wed 01 Aug 2018


~~~ {.sh}

para-ln () {
for depth in $(ls --color=never depth/*.depth); do
bn=$(basename $depth .sorted.depth);
ofn=bedgraph/${bn}_ln.bedgraph;
echo perl code/localNormalisation.pl -quiet -outfile $ofn -- $depth
done
}
perl code/localNormalisation.pl -quiet -outfile stuff -test 100 \
-- depth/Z4.sorted.depth


para-diff () {
for strain in GF J Z; do
for tp in 2 3 4; do
stp=$strain$tp;
strainfn=${stp}_ln.bedgraph;
wtfn=WT${tp}_ln.bedgraph;
outfn=${stp}_lndiff.bedgraph;
echo perl code/bedgraph_diff.pl -outfile bedgraph/$outfn \
-- bedgraph/$strainfn bedgraph/$wtfn
done;
done;
}
#;

~~~

~~~ {.r}

readf <- read.table("../bedgraph/J4_lndiff.bedgraph", sep = "\t",
header = F, skip = 1);

j4 <- readf$V4;

hist(j4, breaks = 30);
x <- rnorm(length(j4), mean = mean(j4), sd = sd(j4))
hist(x, col = "blue",
add = TRUE, breaks = 30)



summary(j4);

pv <- pnorm(j4, mean = mean(j4), sd = sd(j4), lower.tail = F)

apv <- p.adjust(pv, method = "BH")


select <- apv < 1e-5
apvs <- apv[select]

outdf <- cbind(readf[select,], apvs)


# fit <- fitdistr(j4, "cauchy");
# 
# pcv <- pcauchy(j4, location = median(j4),
# scale = sum(quantile(j4, c(0.25, 0.75))), lower.tail = F);
# 
# apcv <- p.adjust(pcv, method = "BH");
# select <- apcv < 1e-5;
# apcvs <- apcv[select];
# 
# outcauchy <- cbind(readf[select,], apcvs)
# nrow(outcauchy);
# #;
# 
# hist(j4);
# hist(rcauchy(30000L, location = fit$estimate["location"],
# scale = fit$estimate["scale"]), breaks = 30);
# 
# 
# fit <- fitdistr(j4, "negative binomial");

~~~
