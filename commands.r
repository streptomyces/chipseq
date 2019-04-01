# Mon 03 Sep 2018

# See bldc.r and seph.r for actual work. Both of these can be sourced
# at the R prompt to get the output files.

png(file = "/home/sco/mnt/wstemp/matt/bldc_replicates.png", height = 1440,
width = 2560);
par(mfrow = c(2,2));
plot(log2(counts[,1]), log2(counts[,2]), pch = 20, col = "darkblue")
plot(log2(counts[,3]), log2(counts[,4]), pch = 20, col = "darkblue")
plot(log2(counts[,5]), log2(counts[,6]), pch = 20, col = "darkblue")
plot(log2(counts[,7]), log2(counts[,8]), pch = 20, col = "darkblue")
dev.off();
#;

rm(c("out", "outdf"));
out <- list();
cycle <- 1;
for(nam in names(tops_list)) {
tdf <- tops_list[[nam]];
if(cycle == 1) {
rn <- rownames(tdf);
}
out[[nam]] <- tdf[rn, "logFC"];
cycle = cycle + 1;
}
outdf <- as.data.frame(out);
rownames(outdf) <- rn;
rownames(outdf) <- sub("vnz_", "", rownames(outdf));
ford <- order(as.numeric(rownames(outdf)));
lfc <- as.matrix(outdf[ford,]);

head(lfc)
tail(lfc)

#;

source("../code/functions.r");
bc <- read.table("../bldc2.csv", quote = "", stringsAsFactors = FALSE,
header = T, row.names = "position", sep = "\t");

genes <- bc$closestgene;
genes <- unique(genes);
pdfn <- c("../bldc_plots.pdf");

# genes <- sort(genes);
# pdfn <- c("../bldc_plots_by_pos.pdf");

pdf(file = pdfn, paper = "a4", width = 6.5, height = 10);
par(mfrow = c(2,1));
plotCnt <- 0;
for(gene in genes) {
  gdf <- bc[bc$closestgene == gene,]
  gdf <- gdf[1,]
  fdr <- gdf$minFDR;
  cfdr <- sprintf("%0.3E", fdr);
  plotgene(gene, cfdr);
  plotCnt <- plotCnt + 1;
  # if(plotCnt >= 50) { break; }
}
dev.off();
system(paste("cp", pdfn,  "/home/sco/mnt/wstemp/matt/"));
#;


# {{{ making of meanmat

files <- list.files(path = "../bedgraph", pattern = ".*ln_mean\\.bedgraph",
full.names = T);

tl <- list();
for(meanfn in files) {
  t1 <- basename(meanfn);
  regret <- regexpr("BldC(-WT)?-\\d\\d", t1, perl = T);
  t2 <- regmatches(t1, regret);
  # t3 <- sub("-WT-", "-", t2);
  t4 <- gsub("-", "_", t2);
  ifh <- file(meanfn, open = "r");
  header <- readLines(ifh, n = 1);
  cn <- c("mol", "start", "end", "lnmean");
  indf <- read.table(ifh, header = F, sep = "\t", stringsAsFactors = F,
  col.names = cn);
  tl[[t4]] <- indf$lnmean;
  if(exists("gpos")) {
  } else {
    gpos <- indf$start + 15;
  }
  close(ifh);
}

meandf <- data.frame(
control = tl$BldC_14,
b10 = tl$BldC_WT_10,
b14 = tl$BldC_WT_14,
b22 = tl$BldC_WT_22
);
rownames(meandf) <- gpos;

meanmat <- as.matrix(meandf);

saveRDS(meanmat, "meanmat.rds");

# }}}


meanmat <- readRDS("meanmat.rds");

# {{{ calls to plotmean_gene()

source("../code/functions.r");
bc <- read.table("../aarout/BldC_ln.csv", quote = "", stringsAsFactors = FALSE,
header = T, row.names = "position", sep = "\t");

genes <- bc$closestgene;
genes <- unique(genes);
pdfn <- c("../bldc_ln_plots.pdf");

pdf(file = pdfn, paper = "a4", width = 6.5, height = 10);
par(mfrow = c(2,1));
plotCnt <- 0;
for(gene in genes) {
  gdf <- bc[bc$closestgene == gene,]
  gdf <- gdf[1,]
  fdr <- gdf$minFDR;
  cfdr <- sprintf("%0.3E", fdr);
  plotmean_gene(gene, cfdr, meanmat)
  plotCnt <- plotCnt + 1;
  # if(plotCnt >= 50) { break; }
}
dev.off();
system(paste("cp", pdfn,  "/home/sco/mnt/wstemp/matt/"));
#;
# }}}



# {{{ making of diffmat

files <- list.files(path = "../bedgraph", pattern = ".*ln_diff\\.bedgraph",
full.names = T);

tl <- list();
for(difffn in files) {
  t1 <- basename(difffn);
  regret <- regexpr("BldC_\\d\\d", t1, perl = T);
  t2 <- regmatches(t1, regret);
  # t3 <- sub("-WT-", "-", t2);
  t4 <- gsub("-", "_", t2);
  ifh <- file(difffn, open = "r");
  header <- readLines(ifh, n = 1);
  cn <- c("mol", "start", "end", "lndiff");
  indf <- read.table(ifh, header = F, sep = "\t", stringsAsFactors = F,
  col.names = cn);
  close(ifh);
  tl[[t4]] <- indf$lndiff;
  if(exists("gpos")) {
  } else {
    gpos <- indf$start + 15;
  }
}

diffdf <- data.frame(
b10 = tl$BldC_10,
b14 = tl$BldC_14,
b22 = tl$BldC_22
);
rownames(diffdf) <- gpos;

diffmat <- as.matrix(diffdf);

# saveRDS(diffmat, "diffmat.rds");

# }}}


# {{{ calls to plotdiff_gene()

source("../code/functions.r");
bc <- read.table("../aarout/BldC_ln.csv", quote = "", stringsAsFactors = FALSE,
header = T, row.names = "position", sep = "\t");

genes <- bc$closestgene;
genes <- unique(genes);
pdfn <- c("../bldc_ln_diff_plots.pdf");

pdf(file = pdfn, paper = "a4", width = 6.5, height = 10);
par(mfrow = c(2,1));
plotCnt <- 0;
for(gene in genes) {
  gdf <- bc[bc$closestgene == gene,]
  gdf <- gdf[1,]
  fdr <- gdf$minFDR;
  cfdr <- sprintf("%0.3E", fdr);
  plotdiff_gene(gene, cfdr, diffmat)
  plotCnt <- plotCnt + 1;
  # if(plotCnt >= 50) { break; }
}
dev.off();
system(paste("cp", pdfn,  "/home/sco/mnt/wstemp/matt/"));
#;
# }}}


# Making list of genes for the dcw region.
# {{{ dcw region
dcw <- character(0);
for(gn in seq(from = 8495, to = 8570, by = 5)) {
cgn <- sprintf("%05d", gn);
vnz <- paste0("vnz_", cgn);
cat(vnz, cgn, "\n")
dcw <- c(dcw, vnz);
}
#;
# }}}

# {{{ calls to plotmean_gene()

source("../code/functions.r");

pdfn <- c("../dcw_ln_plots.pdf");

pdf(file = pdfn, paper = "a4", width = 6.5, height = 10);
par(mfrow = c(2,1));
plotCnt <- 0;
for(gene in dcw) {
  cfdr <- "";
  plotmean_gene(gene, cfdr, meanmat);
  plotCnt <- plotCnt + 1;
  # if(plotCnt >= 3) { break; }
}
dev.off();
system(paste("cp", pdfn,  "/home/sco/mnt/wstemp/matt/"));
#;


source("../code/functions.r");
pdfn <- c("../dcw_ln_diff_plots.pdf");

pdf(file = pdfn, paper = "a4", width = 6.5, height = 10);
par(mfrow = c(2,1));
plotCnt <- 0;
for(gene in dcw) {
  cfdr <- "";
  plotdiff_gene(gene, cfdr, diffmat);
  plotCnt <- plotCnt + 1;
  # if(plotCnt >= 3) { break; }
}
dev.off();
system(paste("cp", pdfn,  "/home/sco/mnt/wstemp/matt/"));
#;




# }}}

# Thu 11 Oct 2018

diffmat <- readRDS("diffmat.rds");
meanmat <- readRDS("meanmat.rds");
source("../code/functions.r");
lut <- mb.product();
lut <- mb.gene.names();

# {{{ make lut data frame
lutcn <- c("chr", "start", "end", "gene", "ignore1", "strand",
"from", "to", "color");
lut <- read.table("../vnz.bed", sep = "\t", stringsAsFactors = F,
skip = 1, col.names = lutcn);
lut$start = lut$start + 1;
rownames(lut) <- lut$gene;
lut <- lut[, c("start", "end", "strand")];
strand <- numeric(nrow(lut));
strand[lut$strand == '+'] <- 1;
strand[lut$strand == '-'] <- -1;
lut$strand <- strand;
# }}}

# {{{ calls to plotmean_gene_eps()

source("../code/functions.r");
bc <- read.table("../aarout/BldC_ln.csv", quote = "", stringsAsFactors = FALSE,
header = T, row.names = "position", sep = "\t");

genes <- bc$closestgene;
genes <- unique(genes);

plotCnt <- 0;
for(gene in genes) {
  gdf <- bc[bc$closestgene == gene,]
  gdf <- gdf[1,]
  fdr <- gdf$minFDR;
  cfdr <- sprintf("%0.3E", fdr);
  epsfn <- file.path("../eps", paste0(gene, "_lt.eps"));
  plotmean_gene_eps(gene, cfdr, meanmat, epsfn, 1)
  system(paste("cp", epsfn,  "/home/sco/mnt/wstemp/matt/eps/"));
  
  epsfn <- file.path("../eps", paste0(gene, ".eps"));
  plotmean_gene_eps(gene, cfdr, meanmat, epsfn, 0)
  system(paste("cp", epsfn,  "/home/sco/mnt/wstemp/matt/eps/"));
  
  plotCnt <- plotCnt + 1;
  cat(paste(plotCnt, gene, "\n"));
#  if(plotCnt >= 3) { break; }
}

#;
# }}}

# Mon 15 Oct 2018

# {{{ Volcano plot for RNA-Seq data.

bc <- read.table("../aarout/BldC_ln.csv", quote = "", stringsAsFactors = FALSE,
header = T, row.names = "position", sep = "\t");

bc10 <- bc[bc$lndiff10 != "-",]
temp <- bc10$closestgene
bc10genes <- unique(temp);
length(bc10genes);

bc14 <- bc[bc$lndiff14 != "-",]
temp <- bc14$closestgene
bc14genes <- unique(temp);
length(bc14genes);

bc22 <- bc[bc$lndiff22 != "-",]
temp <- bc22$closestgene
bc22genes <- unique(temp);
length(bc22genes);

genes <- list(bc10 = bc10genes, bc14 = bc14genes, bc22 = bc22genes);


red <- read.table("../bldC_rnaseq.txt", header = TRUE, sep = "\t",
stringsAsFactors = FALSE, quote = "", row.names = "gene");


hours <- list(
h10 = c(4,5, 10),
h14 = c(7,8, 14),
h22 = c(10,11, 22)
);

hour <- hours[1];


for(dxv in hours) {
epsfn <- file.path("../volcano", paste0(dxv[3], "_hour_volcano.eps"));
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special");
postscript(file = epsfn);
gln <- paste0("bc", dxv[3]);
glg <- genes[[gln]];

temp1 <- red[glg,]
temp2 <- red[setdiff(rownames(red), glg), ];
nrow(temp1);
nrow(temp2);

plot(temp2[,dxv[1]], -(log10(temp2[,dxv[2]])), pch = 20,
col = "lightgray", xlab = "log fold change",
ylab = expression(paste("-", log[10], " P-value")),
main = paste(dxv[3], "hour")
);
points(temp1[,dxv[1]], -(log10(temp1[,dxv[2]])), pch = 20,
col = "darkred");
abline(v = c(-1, 1), lty = 2);
dev.off();
system(paste("cp", epsfn, "/home/sco/mnt/wstemp/matt/volcano_plots/"));

}
#;
# }}}


# {{{ calls to plotmean_gene_jpg()

source("../code/functions.r");
bc <- read.table("../aarout/BldC_ln.csv", quote = "", stringsAsFactors = FALSE,
header = T, row.names = "position", sep = "\t");

genes <- bc$closestgene;
genes <- unique(genes);

plotCnt <- 0;
for(gene in genes) {
  gdf <- bc[bc$closestgene == gene,]
  gdf <- gdf[1,]
  fdr <- gdf$minFDR;
  cfdr <- sprintf("%0.3E", fdr);
  outfn <- file.path("../jpeg", paste0(gene, "_lt.jpg"));
  plotmean_gene_jpg(gene, cfdr, meanmat, outfn, drawgene = 1, do22 = 1)
  system(paste("cp", outfn,  "/home/sco/mnt/wstemp/matt/jpeg/"));
  
  outfn <- file.path("../jpeg", paste0(gene, ".jpg"));
  plotmean_gene_jpg(gene, cfdr, meanmat, outfn, drawgene = 0, do22 = 1)
  system(paste("cp", outfn,  "/home/sco/mnt/wstemp/matt/jpeg/"));
  
  outfn <- file.path("../jpeg", paste0(gene, "_no22_lt.jpg"));
  plotmean_gene_jpg(gene, cfdr, meanmat, outfn, drawgene = 1, do22 = 0)
  system(paste("cp", outfn,  "/home/sco/mnt/wstemp/matt/jpeg/"));
  
  outfn <- file.path("../jpeg", paste0(gene, "_no22.jpg"));
  plotmean_gene_jpg(gene, cfdr, meanmat, outfn, drawgene = 0, do22 = 0)
  system(paste("cp", outfn,  "/home/sco/mnt/wstemp/matt/jpeg/"));
  
  
  plotCnt <- plotCnt + 1;
  cat(paste(plotCnt, gene, "\n"));
  # if(plotCnt >= 3) { break; }
}

#;
# }}}


# {{{ Volcano plot for RNA-Seq data jpegs.

bc <- read.table("../aarout/BldC_ln.csv", quote = "", stringsAsFactors = FALSE,
header = T, row.names = "position", sep = "\t");

bc10 <- bc[bc$lndiff10 != "-",]
temp <- bc10$closestgene
bc10genes <- unique(temp);
length(bc10genes);

bc14 <- bc[bc$lndiff14 != "-",]
temp <- bc14$closestgene
bc14genes <- unique(temp);
length(bc14genes);

bc22 <- bc[bc$lndiff22 != "-",]
temp <- bc22$closestgene
bc22genes <- unique(temp);
length(bc22genes);

genes <- list(bc10 = bc10genes, bc14 = bc14genes, bc22 = bc22genes);

red <- read.table("../bldC_rnaseq.txt", header = TRUE, sep = "\t",
stringsAsFactors = FALSE, quote = "", row.names = "gene");

hours <- list(
h10 = c(4,5, 10),
h14 = c(7,8, 14),
h22 = c(10,11, 22)
);

hour <- hours[1];


for(dxv in hours) {
outfn <- file.path("../volcano", paste0(dxv[3], "_hour_volcano.jpg"));
jpeg(filename = outfn, width = 600, height = 600);
par(cex.axis = 1.5, cex.lab = 1.5);
par(mar = c(5,4.5,4,2));
gln <- paste0("bc", dxv[3]);
glg <- genes[[gln]];

temp1 <- red[glg,]
temp2 <- red[setdiff(rownames(red), glg), ];
nrow(temp1);
nrow(temp2);

plot(temp2[,dxv[1]], -(log10(temp2[,dxv[2]])), pch = 20,
col = "lightgray", xlab = "log fold change",
ylab = expression(paste("-", log[10], " P-value")),
main = paste(dxv[3], "hour")
);
points(temp1[,dxv[1]], -(log10(temp1[,dxv[2]])), pch = 20,
col = "darkred");
abline(v = c(-1, 1), lty = 2);
dev.off();
system(paste("cp", outfn, "/home/sco/mnt/wstemp/matt/volcano_plots/"));

}
#;
# }}}

# Fri 19 Oct 2018
# For processing the output of peak_width.pl

cn <- c("pp", "left", "right", "width");
pw10 <- read.table("../BldC_10_peakwidth.csv",
quote = "", stringsAsFactors = FALSE,
header = F, row.names = 1, col.names = cn, sep = "\t");

pw14 <- read.table("../BldC_14_peakwidth.csv",
quote = "", stringsAsFactors = FALSE,
header = F, row.names = 1, col.names = cn, sep = "\t");

pw <- pw10;

colnames(pw) <- c("left10h", "right10h", "width10h");

pw$left14h <- pw14[rownames(pw), "left"]; 
pw$right14h <- pw14[rownames(pw), "right"]; 
pw$width14h <- pw14[rownames(pw), "width"]; 

# Note that all columns are character because some have "-" in
# the csv file.

head(pw);

fh <- file("../peakwidths.csv", open = "w");
h <- paste(c("peakpos", colnames(pw)), collapse = "\t");
writeLines(h, fh);
write.table(pw, file = fh, quote = F, col.names = F,
row.names = T, sep = "\t");
close(fh);
#;

# Wed 24 Oct 2018

diffmat <- readRDS("diffmat.rds");
meanmat <- readRDS("meanmat.rds");

# {{{ make lut data frame
lutcn <- c("chr", "start", "end", "gene", "ignore1", "strand",
"from", "to", "color");
lut <- read.table("../vnz.bed", sep = "\t", stringsAsFactors = F,
skip = 1, col.names = lutcn);
lut$start = lut$start + 1;
rownames(lut) <- lut$gene;
lut <- lut[, c("start", "end", "strand")];
strand <- numeric(nrow(lut));
strand[lut$strand == '+'] <- 1;
strand[lut$strand == '-'] <- -1;
lut$strand <- strand;
# }}}

# dcw region is vnz_08495 to vnz_08580.
meanmat <- readRDS("meanmat.rds");

source("../code/functions.r");
leftgene <- c("vnz_08495"); 
rightgene <- c("vnz_08580");
outfn <- c("/home/sco/mnt/wstemp/matt/dcwregion_with22_withgene.jpg");
plotmean_region_jpg(leftgene, rightgene, meanmat, outfn, drawgene = 1, do22 = 1)
outfn <- c("/home/sco/mnt/wstemp/matt/dcwregion_with22.jpg");
plotmean_region_jpg(leftgene, rightgene, meanmat, outfn, drawgene = 0, do22 = 1)
outfn <- c("/home/sco/mnt/wstemp/matt/dcwregion_withgene.jpg"); 
plotmean_region_jpg(leftgene, rightgene, meanmat, outfn, do22 = 0, drawgene = 1)
outfn <- c("/home/sco/mnt/wstemp/matt/dcwregion.jpg"); 
plotmean_region_jpg(leftgene, rightgene, meanmat, outfn, do22 = 0, drawgene = 0)
#;

# Thu 25 Oct 2018


# {{{ make lut data frame
lutcn <- c("chr", "start", "end", "gene", "ignore1", "strand",
"from", "to", "color");
lut <- read.table("../vnz.bed", sep = "\t", stringsAsFactors = F,
skip = 1, col.names = lutcn);
lut$start = lut$start + 1;
rownames(lut) <- lut$gene;
lut <- lut[, c("start", "end", "strand")];
strand <- numeric(nrow(lut));
strand[lut$strand == '+'] <- 1;
strand[lut$strand == '-'] <- -1;
lut$strand <- strand;
# }}}

# dcw region is vnz_08495 to vnz_08580.
meanmat <- readRDS("meanmat.rds");

source("../code/functions.r");
leftgene <- c("vnz_25940"); 
rightgene <- c("vnz_25960");
outfn <- c("/home/sco/mnt/wstemp/matt/hups_with22_withgene.jpg");
plotmean_region_jpg(leftgene, rightgene, meanmat, outfn, drawgene = 1, do22 = 1,
vert = c(5673510, 5673345));
#;


afiles <- list.files(path = "../bedgraph", pattern = ".*ln_diff\\.bedgraph",
full.names = T);

for(afn in afiles) {
bfn <- sub("A", "B", afn);
temp <- read.table(afn, skip = 1, header = F);
a <- temp$V4;

qqnorm(a);
break
}
#;


lne <- read.table("../locallyNormalisedEnrichment.txt", sep = "\t",
header = T, stringsAsFactors = F, quote = "");

up2 <- pmax(lne$WT_10hMinusControl, lne$WT_14hMinusControl);

ford <- order(up2, decreasing = T);
slne <- lne[ford,]
head(slne);

write.table(slne, "../byDiff.txt", col.names = T, row.names = F, quote = F,
sep = "\t");


### SepH related work in 2019 ###

# Fri 29 Mar 2019
# See commands.md for today to see how *_ln_diff.bedgraph files
# were arrived at.


# contrich means enrichment in control has been subtracted from
# the enrichment in the WT.
cn <- c("chr", "left", "right", "contrich")
df22 <- read.table("../bedgraph/wtseph22_ln_diff.bedgraph", sep = "\t",
col.names = cn, skip = 1, stringsAsFactors = F)
df22 <- df22[, -1]
tmat <- as.matrix(df22[, c("left", "right")])
robar <- apply(tmat, 1, mean)
df22 <- df22[, c(-1)]
df22$right <- robar
colnames(df22)[1] <- c("pos");
pv <- pnorm(df22$contrich, mean = 0, sd = sd(df22$contrich), lower.tail = F)
df22$apv <- p.adjust(pv, method = "BH");
ordvec <- order(df22$apv, decreasing = F);
df22 <- df22[ordvec,]
head(df22);
write.table(df22, "../seph22.csv", row.names = F, col.names = T,
sep = "\t", quote = F);


# plot(df22$contrich, -log10(df22$apv))


df14 <- read.table("../bedgraph/wtseph14_ln_diff.bedgraph", sep = "\t",
col.names = cn, skip = 1, stringsAsFactors = F)
df14 <- df14[, -1]
tmat <- as.matrix(df14[, c("left", "right")])
robar <- apply(tmat, 1, mean)
df14 <- df14[, c(-1)]
df14$right <- robar
colnames(df14)[1] <- c("pos");
pv <- pnorm(df14$contrich, mean = 0, sd = sd(df14$contrich), lower.tail = F)
df14$apv <- p.adjust(pv, method = "BH");
ordvec <- order(df14$apv, decreasing = F);
df14 <- df14[ordvec,]
head(df14);
write.table(df14, "../seph14.csv", row.names = F, col.names = T,
sep = "\t", quote = F);



par(mfrow=c(1,2));
hist(df14$contrich)
hist(df22$contrich)


pnorm(1, mean = 0, sd = sd(df22$contrich), lower.tail = F)
sd(df22$contrich)






#;



