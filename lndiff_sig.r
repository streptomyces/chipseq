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
# fdr_thresh <- c(1e-4);
fdr_thresh <- c(5e-3);
fdra <- sprintf("%1.0e", 5e-3);
fdrc <- gsub("0|-", "", fdra, perl = T);
outapp <- paste0("_", fdrc, ".csv");
outapp;

for(lndfn in files) {
  bn <- basename(lndfn);
  t1 <- sub("_ln_diff\\.bedgraph", outapp, bn);
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

  select <- indf[indf$FDR <= fdr_thresh,];
  select <- select[order(select$FDR),]

    outflag = FALSE;
  for(iddo in rownames(select)) {
    id_pos <- gps(iddo);
### left
  temp <- lut[((lut$end - ingene_allow) < id_pos & lut$strand == -1) , ]
    left_gene <- temp[nrow(temp),]
    left_gene$distance <- id_pos - left_gene$end;
    if(left_gene$distance > pg_thresh) {
    left_gene <- data.frame(Gene = NA, start = NA, end = NA, strand = NA,
        product = NA, names = NA, distance = NA); 
    }
    colnames(left_gene) <- paste0("l_", colnames(left_gene));
  
### right
  temp <- lut[((lut$start + ingene_allow) > id_pos & lut$strand == 1) , ]
    right_gene <- temp[1,]
    right_gene$distance <- right_gene$start - id_pos;
    if(right_gene$distance > pg_thresh) {
    right_gene <- data.frame(Gene = NA, start = NA, end = NA, strand = NA,
        product = NA, names = NA, distance = NA); 
    }
  colnames(right_gene) <- paste0("r_", colnames(right_gene));


  if(exists("in_gene")) { rm(in_gene); }
  in_gene <- data.frame(Gene = NA, start = NA, end = NA, strand = NA,
      product = NA, names = NA, distance = NA); 
  temp <- lut[(lut$start <= id_pos & lut$end >= id_pos) , ]
  if(nrow(temp) >= 1) {
    in_gene <- temp[1,]
    if(in_gene$strand == 1) {
    in_gene$distance <- in_gene$start - id_pos;
    } else {
    in_gene$distance <- id_pos - in_gene$end;
    }
  }

  if(is.numeric(in_gene$distance)) {
  if(abs(in_gene$distance) > pg_thresh) {
    in_gene <- data.frame(Gene = NA, start = NA, end = NA, strand = NA,
    product = NA, names = NA, distance = NA); 
  }
  }
    colnames(in_gene) <- paste0("i_", colnames(in_gene));
  gene_peak <- paste(left_gene[1, "l_Gene"], in_gene[1, "i_Gene"],
                     right_gene[1, "r_Gene"], sep = ":");
  onedf <- cbind(gene_peak, select[iddo,], id_pos, left_gene, in_gene, right_gene);

  if(outflag) {
    outdf <- rbind(outdf,onedf);
  } else {
    outdf <- onedf;
    outflag <- TRUE;
  }

}

collog <- grepl("mol|id_pos|Gene|start|end|strand", colnames(outdf), perl = T,
ignore.case = F);

collog <- !collog;

outdf <- outdf[, collog];

rownames(outdf) <- sub("^vnz_", "", rownames(outdf));
temp <- c("position", colnames(outdf));
headline <- paste(temp, collapse = "\t");

outfh <- file(outfn, open = "w");
writeLines(headline, outfh);
# writeLines("\n", outfh);
write.table(outdf, outfh, col.names = FALSE,
row.names = TRUE, sep = "\t", quote = FALSE, na = "-");
close(outfh);

#  break;
}

