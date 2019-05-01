#!/usr/bin/Rscript --vanilla
library("getopt");
spec = matrix(c(
"infile",          "i",    2,   "character",
"outfile",         "o",    1,   "character"
), byrow=TRUE, ncol=4);
# In column 3 above:
# 0 = no arg,
# 1 = required argument,
# 2 means optional argument.
# Column 4 can be logical, integer, double, complex or character. 
opt <- getopt(spec); # getopt() returns a list.


ifn <- opt$infile;
ofn <- opt$outfile;

ifn <- "wtseph22_scln_mean.bedgraph";
ofn <- "../stuff.csv";


ofh <- file(ofn, open = "w");


ifh <- file(ifn, open = "r");
inhead <- readLines(con = ifh, n = 1);
# outhead <- gsub("_ln", "_scln", inhead);
# writeLines(outhead, con = ofh);


cn <- c("chr", "left", "right", "val")
df <- read.table(ifh, sep = "\t",
col.names = cn, skip = 1, stringsAsFactors = F)
df <- df[, -1]
tmat <- as.matrix(df[, c("left", "right")])
robar <- apply(tmat, 1, mean)
df <- df[, c(-1)]
df$right <- robar
colnames(df)[1] <- c("pos");
cat(mean(df$val), "\n");
cat(sd(df$val), "\n");
pv <- pnorm(df$val, mean = 0, sd = sd(df$val), lower.tail = F)
df$pv <- pv;
df$apv <- p.adjust(pv, method = "fdr");
ordvec <- order(df$pv, decreasing = F);
df <- df[ordvec,]
head(df);


write.table(df, file = ofh, row.names = F, col.names = T,
sep = "\t", quote = F);


close(ifh);
close(ofh);

