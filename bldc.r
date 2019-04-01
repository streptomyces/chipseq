# Mon 03 Sep 2018

library("Rsubread");
library("edgeR");
rm(list=ls());

source("../code/functions.r");

# saf.file <- c("../vnz_chipseq.saf");
# counts <- mb.sam2cov(saf.file, "../sam");
# saveRDS(counts, "counts.rds");

all.counts <- readRDS("counts.rds");
bldc.counts <- all.counts[, grepl("BldC", colnames(all.counts))];

### edgeR ###

counts <- bldc.counts;

group <- factor(rep(c("B14", "BW10", "BW14", "BW22"), each = 2))
# group <- factor(rep(c("dS22", "S14", "S22", each = 2))
design <- model.matrix(~group)

y <- DGEList(counts = counts, group = group)

y <- calcNormFactors(y)
y <- estimateDisp(y,design)

# To perform quasi-likelihood F-tests:
gr2coef <- c(2, 3, 4);
names(gr2coef) <- c("BW10", "BW14", "BW22");
fit <- glmQLFit(y, design)


# {{{ list of data frames of top tags.
tops_list <- list();
for(gr in c("BW10", "BW14", "BW22")) {
qlf <- glmQLFTest(fit, coef = gr2coef[gr])
tt.qlf <- topTags(qlf, n = nrow(y));
tops_list[[gr]] <- tt.qlf$table;
}
# }}}
#;


### Genome features ###
ft <- read.table("../vnz.ft", sep = "\t", stringsAsFactors = F, quote = "");
lut <- ft;
lut[1] <- rownames(ft);
colnames(lut)[1] <- "Gene";




sig_list <- list();

for(strainTime in names(tops_list)) {
outfn <- file.path("../aarout", paste0(strainTime, ".csv"));

topf <- tops_list[[strainTime]];

mb_rn <- rownames(topf);
mb_pos <- sub("^vnz_", "", mb_rn);
mb_pos <- as.integer(mb_pos);
mb_ord <- order(mb_pos);

posf <- topf[mb_ord, ]
head(posf);


# {{{ Find peaks. Writes sigdf.
outflag <- FALSE;
lfcThresh <- log2(3);
lcpmThresh <- 1;

# for(rownum in 3:(nrow(posf)-2)) {
#   st <- rownum - 2;
#   en <- rownum + 2;
#   lfc <- posf[st:en, "logFC"];
#   lcpm <- posf[rownum, "logCPM"];
# 
#   if((max(lfc) == lfc[3]) & (lfc[3] >= lfcThresh) & (lcpm >= lcpmThresh)) {
#     if((lfc[2] > lfc[1]) & (lfc[4] > lfc[5])) {
#       subdf <- posf[st:en,];
#       cat(lfc, "\n");
# 
#       if(outflag) {
#         sigdf <- rbind(sigdf, posf[rownum,])
#       } else {
#         sigdf <- posf[rownum,];
#         outflag <- TRUE;
#       }
#     }
#   }
# }
sigdf <- posf[(posf$logFC >= lfcThresh) & (posf$logCPM >= lcpmThresh),]
sig_ord <- sigdf[order(sigdf$FDR),]
# }}}

pg_thresh <- 500;

# {{{ write output file
outflag = FALSE;
for(iddo in rownames(sig_ord)) {
  id_pos <- gps(iddo);
  ingene_allow <- 0;

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
  onedf <- cbind(gene_peak, sig_ord[iddo,], id_pos, left_gene, in_gene, right_gene);

  if(outflag) {
    outdf <- rbind(outdf,onedf);
  } else {
    outdf <- onedf;
    outflag <- TRUE;
  }
}

collog <- grepl("PValue|id_pos|^F$|Gene|start|end|strand", colnames(outdf), perl = T,
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
sig_list[[strainTime]] <- outdf;
#;
# }}}

}

#;

postprocess <- c(
"cd aarout",
"perl ../code/csv_merge.pl -outfile bldc.csv -- BW10.csv BW14.csv BW22.csv",
"cp bldc.csv ~/mnt/wstemp/",
""
);

writeLines(postprocess);

# fico <- fit$coefficients;
# coef2bedgraph(fico, 2, "../bedgraph/bw10_fico.bedgraph");
# coef2bedgraph(fico, 3, "../bedgraph/bw14_fico.bedgraph");
# coef2bedgraph(fico, 4, "../bedgraph/bw22_fico.bedgraph");
