afiles <- list.files(path = "../bedgraph", pattern = ".*A-.*",
full.names = T);

for(afn in afiles) {
bfn <- sub("A", "B", afn);
temp <- read.table(afn, skip = 1, header = F);
a <- temp$V4;
temp <- read.table(bfn, skip = 1, header = F);
b <- temp$V4;
coepred <- 1 - (sum((a - b )^2)/sum((a - mean(a))^2));

corr <- cor(a,b, method = "spearman");

cat(afn, coepred, corr, "\n");

# break;
}

