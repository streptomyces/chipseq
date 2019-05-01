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

infile <- opt$infile;
outfile <- opt$outfile;
ofh <- file(outfile, open = "w");


ifh <- file(infile, open = "r");
inhead <- readLines(con = ifh, n = 1);
outhead <- gsub("_ln", "_scln", inhead);
writeLines(outhead, con = ofh);


rdf <- read.table(ifh, stringsAsFactors = F, skip = 0);
invec <- rdf[[4]];
minpositive <- min(invec[invec > 0])
cat(minpositive, "\n");
invec[invec < minpositive] <- minpositive;
invec <- log2(invec);
scvec <- scale(invec);


odf <- data.frame(rdf[[1]], rdf[[2]], rdf[[3]], scvec);
write.table(odf, file = ofh, row.names = F, col.names = F,
quote = F, sep = "\t");

close(ifh);
close(ofh);


# plush <- file(description = plusfile, open = "r");
# minush <- file(description = minusfile, open = "r");
# 
# pluswighead <- readLines(con = plush, n = 2);
# minuswighead <- readLines(con = minush, n = 2);
# 
# close(plush);
# close(minush);
# 
# writeLines(pluswighead); # writing to stdout.
# writeLines(minuswighead); # writing to stdout.



