rm(list = ls());
source("../code/functions.r");

files <- list.files(path = "../bedgraph", pattern = ".*ln_diff\\.bedgraph",
full.names = T);

### Genome features ###
ft <- read.table("../vnz.ft", sep = "\t", stringsAsFactors = F, quote = "");
lut <- ft;
lut[1] <- rownames(ft);
colnames(lut)[1] <- "Gene";

ingene_allow <- c(0);
pg_thresh <- c(500);
fdr_thresh <- c(1e-4);

for(lndfn in files) {
  bn <- basename(lndfn);
  t1 <- sub("_ln_diff\\.bedgraph", "_sig_all.csv", bn);
  outfn <- file.path("../aarout", t1);
  ifh <- file(lndfn, open = "r");
  header <- readLines(ifh, n = 1);
  cat(header, "\n");
  cn <- c("mol", "start", "end", "lndiff");
  indf <- read.table(ifh, header = F, sep = "\t", stringsAsFactors = F,
  col.names = cn);
  close(ifh);
  temp <- paste0("vnz_", indf$start + 15);
  rownames(indf) <- temp;
  lndiff <- indf$lndiff;
  pv <- pnorm(lndiff, mean = mean(lndiff), sd = sd(lndiff), lower.tail = F);
  apv <- p.adjust(pv, method = "BH");
  indf$FDR <- apv;

outfh <- file(outfn, open = "w");
# writeLines("\n", outfh);
write.table(indf, outfh, col.names = TRUE,
row.names = TRUE, sep = "\t", quote = FALSE, na = "-");
close(outfh);

#  break;
}

