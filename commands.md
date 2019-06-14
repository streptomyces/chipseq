# Corinne

### Fri 07 Jun 2019

~~~ 

perl code/comgen.pl -cpu 3 -test 200000 | parallel

perl code/comgen.pl -cpu 3 -test 0 | parallel
~~~

### Mon 10 Jun 2019

~~~ 
reportcollect () {
for repo in $(ls --color=never bowtie2reports/*bowtie2); do
  echo $repo
  cat $repo
  echo
done
}
reportcollect > alignment_reports.txt

pushd sam
date
parallel --jobs 4 samtools view -@ 4 -b -o {.}.bam {} \
::: $(ls --color=never ./*.sam)

date
parallel --jobs 4 samtools sort -m 4G -@ 4 -o {.}.sorted.bam {} \
::: $(ls --color=never ./*.bam)

date
parallel --jobs 14 samtools index {} \
::: $(ls --color=never ./*.sorted.bam)
date
pushd
~~~

Depth and local normalisation.

~~~ 

parallel --dry-run samtools depth -a {} \> depth/{/.}.depth \
::: $(ls --color=never sam/*.sorted.bam)

# Test
perl code/localNormalisation.pl -quiet -test 1000 -outfile \
test.bedgraph -- depth/3Flag1.sorted.depth

# Real run
para-ln () {
for depth in $(ls --color=never depth/*.depth); do
bn=$(basename $depth .sorted.depth);
ofn=bedgraph/${bn}_ln.bedgraph;
echo perl code/localNormalisation.pl -aperture 30 -region 3000 \
-quiet -outfile $ofn -- $depth
done
}
para-ln | parallel --dry-run


~~~

~~~ 

perl ~/code/perl/seq/gbk2bed.pl -infile Avin.gbk -outfile Avin.bed \
-molecule Avin
cp Avin.bed /mnt/isilon/bedfiles/
~~~

Below, use Rsubread featureCounts() to get reads for every section.

~~~ {.r}

sam.path <- c("../sam");
saf.file <- c("../Avin_chipseq.saf");
sam.files <- list.files(path = sam.path, pattern = "\\.sorted\\.bam$",
full.names = TRUE);

fc <- featureCounts(annot.ext = saf.file, nthreads = 12,
allowMultiOverlap = TRUE, isPairedEnd = TRUE,
countChimericFragments = TRUE, minMQS = 38, nonSplitOnly = TRUE,
requireBothEndsMapped = TRUE,
files = sam.files);

counts <- fc$counts;

cn <- colnames(counts);
ncn <- sub("^.*sam\\.", "", cn)
ncn <- sub("\\.sorted\\.bam$", "", ncn);

colnames(counts) <- ncn


ofh <- file("../rscounts.txt", open = "wt");
row1 <- c("section", colnames(counts))
ts <- paste(row1, collapse = "\t");
writeLines(ts, con = ofh);
write.table(counts, file = ofh, col.names = F, row.names = T, quote = F, sep = "\t")
close(ofh);

<!-- {{{ R functions scf and tpm -->

~~~ {.r}

# scf. return scaling factors for each number in a vector.
scf <- function(y) {
x <- sort(y);
unit <- x[floor(length(x) / 2)]
scf <- y / unit;
names(scf) <- names(y);
return(scf);
}

# TPM
# count for each gene / length of gene in KB (RPK).
# Sum of all the RPKs / 1e6 (scaling factor).
# RPK of each gene / the scaling factor.
tpm <- function(y) {
y <- counts;
retmat <- y; 
for (cn in colnames(y)) {
cnt <- y[, cn];
rpk <- cnt / (30 / 1000);
sf <- sum(rpk) / 1e6;
retmat[, cn] <- rpk / sf;
}
return(retmat);
}

~~~

<!-- }}} -->

scal <- scf(colSums(counts));
normcnt <- counts;
for (cn in colnames(counts)) {
normcnt[, cn] <- round(counts[, cn] / scal[cn]);
}

ofh <- file("../normcounts.txt", open = "wt");
row1 <- c("section", colnames(normcnt))
ts <- paste(row1, collapse = "\t");
writeLines(ts, con = ofh);
write.table(normcnt, file = ofh, col.names = F, row.names = T, quote = F, sep = "\t")
close(ofh);

tpmat <- tpm(counts);

ofh <- file("../tpm.txt", open = "wt");
row1 <- c("section", colnames(tpmat))
ts <- paste(row1, collapse = "\t");
writeLines(ts, con = ofh);
write.table(tpmat, file = ofh, col.names = F, row.names = T, quote = F, sep = "\t")
close(ofh);


~~~

Below, write bedgraph files for each column in rscounts.txt.
~~~ 
perl code/rscounts2bedgraph.pl -sectlen 30 -genolen 5365318 \
-dirout bedgraph -outpre rs_ -- rscounts.txt

perl code/rscounts2bedgraph.pl -sectlen 30 -genolen 5365318 \
-dirout bedgraph -outpre norm_ -- normcounts.txt

perl code/rscounts2bedgraph.pl -sectlen 30 -genolen 5365318 \
-dirout bedgraph -outpre tpm_ -- tpm.txt
~~~


### Thu 13 Jun 2019

~~~ 
perl code/bedgraph_mean.pl -outfile bedgraph/tpm_3Flag.bedgraph \
-- bedgraph/tpm_3Flag1.bedgraph bedgraph/tpm_3Flag2.bedgraph

perl code/bedgraph_mean.pl -outfile bedgraph/tpm_WT.bedgraph \
-- bedgraph/tpm_WT1.bedgraph bedgraph/tpm_WT2.bedgraph

perl code/bedgraph_log2ratio.pl -outfile bedgraph/lfc.bedgraph \
-- bedgraph/tpm_3Flag.bedgraph bedgraph/tpm_WT.bedgraph

~~~



~~~ 

rm(list=ls());
lfc <- read.table("../bedgraph/lfc.bedgraph", skip = 1, sep = "\t",
stringsAsFactors = F, col.names = c("chr", "start", "end", "lfc"));

pv <- pnorm(lfc$lfc, mean = mean(lfc$lfc), sd = sd(lfc$lfc), lower.tail = F)
apv <- p.adjust(pv, method = "BH");

temp <- lfc[, c(2,3)]
sectmid <- round(rowMeans(temp));

outdf <- data.frame(pos = sectmid, lfc = lfc$lfc, apv = apv);

sigdf <- outdf[outdf$apv <= 0.001,]
ordv <- order(sigdf$apv);
sigdf <- sigdf[ordv,]



write.table(sigdf, file = "../siglfc.txt", quote = F, sep = "\t",
row.names = F, col.names = T);




~~~
