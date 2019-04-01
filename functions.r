### Wed 06 Jun 2018
# Copied from /mnt/isilon/customers/matt_bush/2018_02_05/code/

coef2bedgraph <- function(fico, colwant, bgfn, testcnt = 0) {
ofh <- file(bgfn, open = "w");
bn <- basename(bgfn);
bn <- sub("\\.[^.]+$", "", bn, perl = TRUE);
bgline = paste0("track type=bedGraph name=", bn, " description=", bn);
writeLines(bgline, ofh);
rocnt <- 0;
for(rn in rownames(fico)) {
mid <- gps(rn);
left <- mid - 15;
right <- left + 30;
coef <- fico[rn, colwant];
cat("vnz", sprintf("%d", left), sprintf("%d", right), coef, file = ofh, sep = "\t");
cat("\n", file = ofh);
rocnt <- rocnt + 1;
if(testcnt > 0 & rocnt >= testcnt) { break; }
}
close(ofh);
}


# {{{ gps() to convert ChIP section identifier to genomic position.
gps <- function(id) {
f_pos <- sub("^vnz_", "", id);
f_int <- as.integer(f_pos);
return(f_int);
}
# }}}

# {{{ mb.sam2cov()
mb.sam2cov <- function(saf.file, sam.path) {
sam.files <- list.files(path = sam.path, pattern = "\\.sorted\\.bam$",
full.names = TRUE);

fc <- featureCounts(annot.ext = saf.file, nthreads = 12,
allowMultiOverlap = TRUE, isPairedEnd = TRUE,
countChimericFragments = TRUE, minMQS = 38, nonSplitOnly = TRUE,
requireBothEndsMapped = FALSE,
files = sam.files);

counts <- fc$counts;

cn <- colnames(counts);

ncn <- sub("^\\.\\.\\..*sam\\.", "", cn)

ncn <- sub("\\.sorted\\.bam$", "", ncn);
# ncn <- sub("8", "08", ncn);

# cn
# ncn

colnames(counts) <- ncn

# for(cn in colnames(counts)) {
# ofn = file.path("../coverage", paste(cn, ".tdf", sep = ""));
# odf = as.data.frame(counts[, cn]);
# colnames(odf) <- c(cn);
# write.table(odf, file = ofn, quote = F, sep = "\t", col.names = F,
# row.names = T)
# }

return(counts);

}

# }}}

# {{{ mb.sort
mb.sort <- function(indf) {
temp = colnames(indf);
cn <- grep("\\.apv$", temp, value = TRUE)
cat(cn);

df1 <- indf[, cn]

minvec <- apply(df1, 1, min);
ordvec <- order(minvec)

outdf <- indf[ordvec, ];

return(outdf);
}
# }}}

# {{{ mb.product()
mb.product <- function(infile = "../vnz_products.txt") {
rdf <- read.table(file = infile, header = F, stringsAsFactors = F,
sep = "\t", quote = "", col.names = c("gene", "product"), row.names = "gene")
return(rdf);
}
# }}}

# {{{ mb.gene.names()
mb.gene.names <- function(infile = "../vnzNames.list") {
  rdf <- read.table(file = infile, header = F, stringsAsFactors = F,
      sep = "\t", quote = "");
  vnz <- c();
  gene <- c();
  for(rn in rownames(rdf)) {
    rodf <- rdf[rn,]
      romat <- as.matrix(rodf);
    if(grepl("^vnz", romat[1]) & grepl("[a-zA-Z]", romat[5])) {
      vnz <- c(vnz, romat[1])
      gene <- c(gene, romat[5])
    }
  }
  retdf <- data.frame(canonical = gene);
  rownames(retdf) <- vnz;
  return(retdf);
}
# }}}

# {{{ plotpos
plotpos <- function(gpos) {
left <- gpos - 1500;
right <- gpos + 1500;
tmat <- fico[
as.numeric(rownames(fico)) >= left
& as.numeric(rownames(fico)) <= right,
]

rut <- lut[lut$end > left & lut$start < right, ]
rut <- rut[, c("start", "end", "strand")]
rut <- as.matrix(rut);

tmat <- tmat[, -1]
ymin <- min(tmat);
ymax <- max(tmat);

plot(as.numeric(rownames(tmat)), tmat[,1], ylim = c(ymin, ymax),
type = "l", col = "darkred"
);
lines(as.numeric(rownames(tmat)), tmat[,2], ylim = c(ymin, ymax),
type = "l", col = "darkgreen"
);
lines(as.numeric(rownames(tmat)), tmat[,3], ylim = c(ymin, ymax),
type = "l", col = "darkblue"
);

gcol <- c("red", "blue", "green");
for(gene in rownames(rut)) {
lines(x = c(rut[gene, "start"], rut[gene, "end"]),
y = c(0,0), col = gcol[rut[gene, "strand"] + 2]
);
}

text(rut[, "start"], y = 0, labels = rownames(rut));
}
# }}}

# {{{ plotgene
plotgene <- function(gene, fdr) {
temp <- lut[gene, c("start", "end", "strand")];
if(temp$strand[1] == -1) {
gpos <- temp$end;
} else {
  gpos <- temp$start;
}


left <- gpos - 2500;
right <- gpos + 2500;
tmat <- lfc[
as.numeric(rownames(lfc)) >= left
& as.numeric(rownames(lfc)) <= right,
]

rut <- lut[lut$end > left & lut$start < right, ]
rut <- rut[, c("start", "end", "strand")]
rut <- as.matrix(rut);

# tmat <- tmat[, -1]
ymin <- min(tmat);
ymax <- max(tmat);
if(ymin > -0.5) { ymin <- -0.5; }

plot(as.numeric(rownames(tmat)), tmat[,1], ylim = c(ymin, ymax),
type = "l", col = "darkred", main = paste(gene, fdr), ylab = "logFC",
xlab = "Genomic position", lwd = 2
);
lines(as.numeric(rownames(tmat)), tmat[,2], ylim = c(ymin, ymax),
type = "l", col = "darkgreen", lwd = 2
);
lines(as.numeric(rownames(tmat)), tmat[,3], ylim = c(ymin, ymax),
type = "l", col = "darkblue", lwd = 2
);

gcol <- c("red", "blue", "green");
for(gene in rownames(rut)) {
lines(x = c(rut[gene, "start"], rut[gene, "end"]),
y = c(0,0), col = gcol[rut[gene, "strand"] + 2]
);
}

text(rut[, "start"], y = 0, labels = sub("vnz_", "", rownames(rut)),
pos = 4);
}
# }}}

# {{{ plotmean_gene
plotmean_gene <- function(gene, fdr, vmat) {
temp <- lut[gene, c("start", "end", "strand")];
if(temp$strand[1] == -1) {
gpos <- temp$end;
} else {
  gpos <- temp$start;
}


left <- gpos - 2500;
right <- gpos + 2500;
tmat <- vmat[
as.numeric(rownames(vmat)) >= left
& as.numeric(rownames(vmat)) <= right,
]

rut <- lut[lut$end > left & lut$start < right, ]
rut <- rut[, c("start", "end", "strand")]
rut <- as.matrix(rut);

# tmat <- tmat[, -1]
ymin <- min(tmat);
ymax <- max(tmat);
if(ymin > -0.5) { ymin <- -0.5; }

plot(as.numeric(rownames(tmat)), tmat[,1], ylim = c(ymin, ymax),
type = "n", col = "darkred", main = paste(gene, fdr), ylab = "Normalised enrichment",
xlab = "Genomic position", lwd = 2
);

plotcols <- c("black", "darkred", "darkgreen", "darkblue")
names(plotcols) <- colnames(tmat);

for(tp in colnames(tmat)) {
lines(as.numeric(rownames(tmat)), tmat[,tp], ylim = c(ymin, ymax),
type = "l", col = plotcols[tp], lwd = 2
);
}

gcol <- c("red", "blue", "green");
for(gene in rownames(rut)) {
lines(x = c(rut[gene, "start"], rut[gene, "end"]),
y = c(0,0), col = gcol[rut[gene, "strand"] + 2]
);
}

text(rut[, "start"], y = 0, labels = sub("vnz_", "", rownames(rut)),
pos = 4);
}
# }}}

# {{{ plotdiff_gene
plotdiff_gene <- function(gene, fdr, vmat) {
temp <- lut[gene, c("start", "end", "strand")];
if(temp$strand[1] == -1) {
gpos <- temp$end;
} else {
  gpos <- temp$start;
}


left <- gpos - 2500;
right <- gpos + 2500;
tmat <- vmat[
as.numeric(rownames(vmat)) >= left
& as.numeric(rownames(vmat)) <= right,
]

rut <- lut[lut$end > left & lut$start < right, ]
rut <- rut[, c("start", "end", "strand")]
rut <- as.matrix(rut);

# tmat <- tmat[, -1]
ymin <- min(tmat);
ymax <- max(tmat);
if(ymin > -0.5) { ymin <- -0.5; }

plot(as.numeric(rownames(tmat)), tmat[,1], ylim = c(ymin, ymax),
type = "n", col = "darkred", main = paste(gene, fdr),
ylab = expression(paste("Normalised WT - ", Delta, "BldC enrichment")),
xlab = "Genomic position", lwd = 2
);

plotcols <- c("darkred", "darkgreen", "darkblue")
names(plotcols) <- colnames(tmat);

for(tp in colnames(tmat)) {
lines(as.numeric(rownames(tmat)), tmat[,tp], ylim = c(ymin, ymax),
type = "l", col = plotcols[tp], lwd = 2
);
}

gcol <- c("red", "blue", "green");
for(gene in rownames(rut)) {
lines(x = c(rut[gene, "start"], rut[gene, "end"]),
y = c(0,0), col = gcol[rut[gene, "strand"] + 2]
);
}

text(rut[, "start"], y = 0, labels = sub("vnz_", "", rownames(rut)),
pos = 4);
}
# }}}


# {{{ plotmean_gene_eps
plotmean_gene_eps <- function(gene, fdr, vmat, epsfn, drawgene = 1) {
temp <- lut[gene, c("start", "end", "strand")];
if(temp$strand[1] == -1) {
gpos <- temp$end;
} else {
  gpos <- temp$start;
}


left <- gpos - 2500;
right <- gpos + 2500;
tmat <- vmat[
as.numeric(rownames(vmat)) >= left
& as.numeric(rownames(vmat)) <= right,
]

rut <- lut[lut$end > left & lut$start < right, ]
rut <- rut[, c("start", "end", "strand")]
rut <- as.matrix(rut);

# tmat <- tmat[, -1]
ymin <- min(tmat);
ymax <- max(tmat);
if(ymin > -0.5) { ymin <- -0.5; }

setEPS(horizontal = FALSE, onefile = FALSE, paper = "special");
postscript(file = epsfn);

plot(as.numeric(rownames(tmat)), tmat[,1], ylim = c(ymin, ymax),
type = "n", col = "darkred", main = paste(gene, fdr), ylab = "Normalised enrichment",
xlab = "Genomic position (Mbp)", lwd = 2, xaxt = "n"
);

plotcols <- c("black", "darkred", "darkgreen", "darkblue")
names(plotcols) <- colnames(tmat);

for(tp in colnames(tmat)) {
lines(as.numeric(rownames(tmat)), tmat[,tp], ylim = c(ymin, ymax),
type = "l", col = plotcols[tp], lwd = 2
);
}

if(drawgene) {
pup <- par("usr");
interval <- pup[4] - pup[3];
geh <- interval/90; # gene_edge_height
# cat(geh, "\n");
# cat(pup, "\n");
gcol <- c("red", "blue", "green");
for(gene in rownames(rut)) {
lines(x = c(rut[gene, "start"], rut[gene, "end"]),
y = c(0,0), col = gcol[rut[gene, "strand"] + 2]
);
lines(x = c(rut[gene, "start"], rut[gene, "start"]),
y = c(0,geh), col = gcol[rut[gene, "strand"] + 2]
);
}
text(rut[, "start"], y = 0, labels = sub("vnz_", "", rownames(rut)),
pos = 4, cex = 0.6);
}

xtip <- axTicks(1);
axis(1, at = xtip, labels = format(xtip/1e6, digits = 4),
cex.axis = 1.2);

dev.off();

}
# }}}


# {{{ plotmean_gene_jpg
plotmean_gene_jpg <- function(gene, fdr, vmat, outfn,
drawgene = 1, do22 = 1) {
temp <- lut[gene, c("start", "end", "strand")];
if(temp$strand[1] == -1) {
gpos <- temp$end;
} else {
  gpos <- temp$start;
}


left <- gpos - 2500;
right <- gpos + 2500;
tmat <- vmat[
as.numeric(rownames(vmat)) >= left
& as.numeric(rownames(vmat)) <= right,
]

rut <- lut[lut$end > left & lut$start < right, ]
rut <- rut[, c("start", "end", "strand")]
rut <- as.matrix(rut);

# tmat <- tmat[, -1]
ymin <- min(tmat);
ymax <- max(tmat);
if(ymin > -0.5) { ymin <- -0.5; }

jpeg(filename = outfn, width = 600, height = 600);
par(cex.axis = 1.5, cex.lab = 1.5);
par(mar = c(5,4.5,4,2));
plot(as.numeric(rownames(tmat)), tmat[,1], ylim = c(ymin, ymax),
type = "n", col = "darkred", main = paste(gene, fdr), ylab = "Normalised enrichment",
xlab = "Genomic position (Mbp)", lwd = 2, xaxt = "n"
);

plotcols <- c("black", "darkred", "darkgreen", "darkblue")
names(plotcols) <- colnames(tmat);

if(do22) {
tobedone <- colnames(tmat);
} else {
temp <- colnames(tmat);
tobedone <- temp[! grepl("22$", temp)];
}

for(tp in tobedone) {
lines(as.numeric(rownames(tmat)), tmat[,tp], ylim = c(ymin, ymax),
type = "l", col = plotcols[tp], lwd = 2
);
}

if(drawgene) {
pup <- par("usr");
interval <- pup[4] - pup[3];
geh <- interval/100; # gene_edge_height
# cat(geh, "\n");
# cat(pup, "\n");
gcol <- c("red", "blue", "green");
for(gene in rownames(rut)) {
lines(x = c(rut[gene, "start"], rut[gene, "end"]),
y = c(0,0), col = gcol[rut[gene, "strand"] + 2]
);
lines(x = c(rut[gene, "start"], rut[gene, "start"]),
y = c(0,geh), col = gcol[rut[gene, "strand"] + 2]
);
}
text(rut[, "start"], y = 0, labels = sub("vnz_", "", rownames(rut)),
pos = 4, cex = 0.6);
}

xtip <- axTicks(1);
axis(1, at = xtip, labels = format(xtip/1e6, digits = 4),
cex.axis = 1.5);

dev.off();

}
# }}}


# plotmean_region_jpg(leftgene, rightgene, meanmat, outfn, drawgene = 1, do22 = 1)
# {{{ plotmean_region_jpg
plotmean_region_jpg <- function(leftgene, rightgene, vmat, outfn,
drawgene = 1, do22 = 1, vert) {

temp <- lut[leftgene, c("start", "end", "strand")];
left <- temp$start - 200;

temp <- lut[rightgene, c("start", "end", "strand")];
right <- temp$end + 200;


tmat <- vmat[
as.numeric(rownames(vmat)) >= left
& as.numeric(rownames(vmat)) <= right,
]

rut <- lut[lut$end > left & lut$start < right, ]
rut <- rut[, c("start", "end", "strand")]
rut <- as.matrix(rut);

# tmat <- tmat[, -1]
ymin <- min(tmat);
ymax <- max(tmat);
if(ymin > -0.5) { ymin <- -0.5; }

jpeg(filename = outfn, width = 1800, height = 600);
par(cex.axis = 1.5, cex.lab = 1.5);
par(mar = c(5,4.5,4,2));
plot(as.numeric(rownames(tmat)), tmat[,1], ylim = c(ymin, ymax),
type = "n", col = "darkred", main = paste(leftgene, "to", rightgene),
ylab = "Normalised enrichment",
xlab = "Genomic position (Mbp)", lwd = 2, xaxt = "n"
);

plotcols <- c("black", "darkred", "darkgreen", "darkblue")
names(plotcols) <- colnames(tmat);

if(do22) {
tobedone <- colnames(tmat);
} else {
temp <- colnames(tmat);
tobedone <- temp[! grepl("22$", temp)];
}

for(tp in tobedone) {
lines(as.numeric(rownames(tmat)), tmat[,tp], ylim = c(ymin, ymax),
type = "l", col = plotcols[tp], lwd = 2
);
}

if(drawgene) {
pup <- par("usr");
interval <- pup[4] - pup[3];
geh <- interval/100; # gene_edge_height
# cat(geh, "\n");
# cat(pup, "\n");
gcol <- c("red", "blue", "green");
for(gene in rownames(rut)) {
lines(x = c(rut[gene, "start"], rut[gene, "end"]),
y = c(0,0), col = gcol[rut[gene, "strand"] + 2]
);
lines(x = c(rut[gene, "start"], rut[gene, "start"]),
y = c(0,geh), col = gcol[rut[gene, "strand"] + 2]
);
}
text(rut[, "start"], y = 0, labels = sub("vnz_", "", rownames(rut)),
pos = 4, cex = 0.6);
}

xtip <- axTicks(1);
axis(1, at = xtip, labels = format(xtip/1e6, digits = 4),
cex.axis = 1.5);

if(vert) {
  abline(v = vert);
}

dev.off();

}
# }}}
