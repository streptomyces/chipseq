# Matt Bush BldC and SepH ChIP-Seq

#### Fri 31 Aug 2018

Used index /mnt/isilon/bowtie2Indexes/vnz/chr


~~~ {.sh}

perl code/comgen.pl -test 12000 -cpu 14 | head -n 1 | parallel

perl code/comgen.pl -cpu 14 | parallel --jobs 1

~~~

#### Sat 01 Sep 2018

~~~ {.sh}

reportcollect () {
for repo in $(ls --color=never bowtie2reports/*bowtie2); do
  echo $repo
  cat $repo
  echo
done
}
reportcollect > alignment_reports.txt

pushd sam
parallel --jobs 4 samtools view -@ 4 -b -o {.}.bam {} \
::: $(ls --color=never ./*.sam)

parallel --jobs 4 samtools sort -m 4G -@ 4 -o {.}.sorted.bam {} \
::: $(ls --color=never ./*.bam)

parallel --jobs 14 samtools index {} \
::: $(ls --color=never ./*.sorted.bam)
pushd



parallel --dry-run samtools depth -a {} \> depth/{/.}.depth \
::: $(ls --color=never sam/*.sorted.bam)


para-ln () {
for depth in $(ls --color=never depth/*.depth); do
bn=$(basename $depth .sorted.depth);
ofn=bedgraph/${bn}_ln.bedgraph;
echo perl code/localNormalisation.pl -quiet -outfile $ofn -- $depth
done
}

perl code/localNormalisation.pl -quiet -test 1000 -outfile \
test.bedgraph -- depth/aBldC-14A-IP.sorted.depth

~~~

#### Sun 02 Sep 2018

bldC is vnz_18945.

~~~ {.sh}

perl code/gbk2coordsAndProducts.pl -outfile vnz.ft -- vnz.gbk

~~~

#### Wed 19 Sep 2018

Matt's email of 4 Sep 2018.

~~~ 
   Hi Govind,

    

   I know you are currently preparing for the Summer School so
   please feel free to just file this e-mail away until your return
   from Croatia. I have plenty to get on with in the meantime.

    

   I have looked at the data and tables and looking ahead to
   publication (aim to submit end of October), I have the following
   queries/requests as the next step:

    

    1. Tables:

         a. Current Table – questions and changes needed

                                                                 
   i.      What does each column represent? logFCv CPM, p-value v
   FDR, F? Are some of these values superfluous?

                                                               
   ii.      All time points in one table to give a single “regulon”
   – for each time point could indicate where peak is significant
   by True/False as previously e.g.WhiA/B.

                                                             
   iii.      In previous tables e.g. WhiA/B, peak distance to gene
   has also been recorded and then only genes which have at least
   one significant position a certain distance (e.g. <300bp
   upstream) of the peak are listed. Currently Left and Right genes
   are always listed (irrespective of distance).

                                                             
   iv.      The strand (+1/-1) is listed. If a peak sits between
   two genes on the same strand then only the downstream gene
   should be listed. In some cases, this also means that genes
   e.g.vnz_27205 are listed multiple times

                                                               
   v.      Sometimes the genes recorded on the left or right of a
   peak don’t make sense e.g. vnz_18940 – peak is clearly upstream
   of vnz_18945 but the “right” gene is instead listed as
   vnz_189500 and 18945 is not listed in the table. I think there
   are several examples of this happening.

                                                             
   vi.      Was significance analysed in 25bp segments as
   previously? If this is true, is the position listed, the most
   significant?

                                                            
   vii.      Some peaks are recorded more than once e.g. vnz_13645
   – only one position should be listed for each gene (the most
   significant?)

                                                          
   viii.      Some peaks are not recorded at some time points where
   (at least visually) in IGB they are clearly present e.g.
   upstream of bldC (vnz_18945). How is significance calculated and
   what is the cut-off?

     b. Additional Tables

                                                                 
   i.      Once the above changes have been made to the ChIP-seq
   table, I’ll need an additional table with BldC RNAseq analysis
   integrated (10 14 and 22 hr time points only) as done previously
   for Affy data with other regulons i.e. each gene listed also has
   the expression information

                                                               
   ii.      Separate table detailing significant peaks in the BldC
   mutant – indicating those that are also present in the WT and
   those that are not present in the WT (from IGB it appears that
   the majority of peaks in the mutant are not present in the WT!)

                                                             
   iii.      Separate table detailing the overlap between the SVEN
   regulon (this work) and SCO regulon (BldC Nature Communications)
   – which genes in this analysis are also identified in the SCO
   BldC work?

    

    2. ChIP-seq plots –plots for 10, 14 and 22 hr time points v
       bldC mutant for each of the significant peaks identified.
       The key should not be integrated into the images. The text
       size/axes should be large for Figures
    3. RNAseq plots:

         a. separate volcano plots for the BldC 10, 14 and 22 hr
            time points (with and without genes identified in
            ChIP-seq highlighted e.g. red dots) – logFC v apv.
            lines crossing the horizontal axes should indicate when
            change in expression e.g. logFC>+1 or <-1
         b. Individual line-plots (showing the absolute expression
            values of WT v mutant at 10hr 14hr and 22hr) for each
            of the genes identified in the ChIP-seq analysis

    4. Whole Genome Plots for each time point v Mutant – should be
       on the same scale

    

   Happy to chat about this when you are available again 😊

    

   As always, offering my thanks doesn’t seem to go far enough – It
   would be impossible for me to do any of this well without your
   help!

    

   Matt

~~~

Script seph.r edited. Then it was sourced in an R session. It writes
SW14.csv and SW22.csv in aarout/.

~~~ 

perl ../code/csv_merge.pl -outfile stuff.csv -- SW14.csv SW22.csv

perl ../code/csv_merge.pl -outfile stuff -- BW10.csv BW14.csv BW22.csv

rm -rf bldc_chipseq bldc_chipseq.zip
mkdir bldc_chipseq
cp aarout/bldc2.csv bldc_chipseq/
cp bedgraph/bw*_fico.bedgraph bldc_chipseq/
cp bedgraph/vnz.{bed,fna} bldc_chipseq/
zip --recurse-paths --to-crlf bldc_chipseq bldc_chipseq
unzip -l bldc_chipseq.zip
cp bldc_chipseq.zip ~/mnt/wstemp/matt/
rm -rf bldc_chipseq
#;


perl code/gbk2sqlite.pl
# Above writes /home/sco/sqlite/

~~~

#### Mon 01 Oct 2018

~~~ 

# For bedgraph_mean.pl
perl code/comgen_mean.pl | parallel

~~~

The above writes \*mean\*.bedgraph files in bedgraph/ for all of the
BldC files. The work is done by bedgraph_mean.pl.

~~~ 

perl code/bedgraph_diff.pl -outfile BldC_10_ln_diff.bedgraph -- \
bedgraph/aBldC-WT-10-IP_ln_mean.bedgraph \
bedgraph/aBldC-14-IP_ln_mean.bedgraph

perl code/bedgraph_diff.pl -outfile BldC_14_ln_diff.bedgraph -- \
bedgraph/aBldC-WT-14-IP_ln_mean.bedgraph \
bedgraph/aBldC-14-IP_ln_mean.bedgraph

perl code/bedgraph_diff.pl -outfile BldC_22_ln_diff.bedgraph -- \
bedgraph/aBldC-WT-22-IP_ln_mean.bedgraph \
bedgraph/aBldC-14-IP_ln_mean.bedgraph

pushd aarout
perl ../code/csv_merge3.pl -outfile BldC_ln.csv -- BldC_10.csv BldC_14.csv BldC_22.csv

cp BldC_ln.csv ~/mnt/wstemp/matt/
pushd

~~~

Script bedgraph_mean.pl was used to get the means of the replicates of
the ln files. This leads to the following 4 files.

    aBldC-14-IP_ln_mean.bedgraph
    aBldC-WT-14-IP_ln_mean.bedgraph
    aBldC-WT-10-IP_ln_mean.bedgraph
    aBldC-WT-22-IP_ln_mean.bedgraph

Then bedgraph_diff.pl gives the following three files from the above 4
files.

    BldC_10_ln_diff.bedgraph
    BldC_14_ln_diff.bedgraph
    BldC_22_ln_diff.bedgraph

The above were processed in R using lndiff_sig.r to produce

    BldC_10.csv
    BldC_14.csv
    BldC_22.csv

Then csv_merge3.pl was used to merge the above 3 files to the following
single file which was emailed to Matt.

    BldC_ln.csv



#### Tue 02 Oct 2018

Matt emailed.

    Hi Govind - thanks so much for this. I think this really clears
    things up. As the next step can I ask for....?
    
    1. 2 sets of ChIP-seq plots – 1) ln only and 2) ln_Diff (this will
    enable me to identify false positives)
    
    2. Can we also look specifically at the dcw cluster
    (vnz_08495-08570).  In the original locally normalised only
    analysis, ftsZ was analysed as a target and other genes in this
    cluster also seem to be bound by BldC – can we extract the
    information for positions in these regions. Do they just fall
    short of the threshold?

In response to the above the following 4 files were sent to Matt.
They were made by calling plotdiff_gene() and plotmean_gene()
functions in functions.r. The actual commands are mentioned in
commands.r.

    2.0M Oct  2 12:41 bldc_ln_plots.pdf
    1.6M Oct  2 13:06 bldc_ln_diff_plots.pdf
    123K Oct  2 17:43 dcw_ln_plots.pdf
    101K Oct  2 18:17 dcw_ln_diff_plots.pdf



#### Wed 03 Oct 2018

    Thanks so much for your work on this. The good news is that I am
    happy that we have settled on a final table that is an accurate
    and objective analysis of the BldC regulon. I have attached this
    table - for all future analysis, it might be better if you use
    this list of genes as reference (I have deleted a few duplicates
    and false positives). I can now start to think about the
    manuscript.
    
    I know you are incredibly busy - so thanks for taking the time
    with this.
    
    Below I have summarised what the next steps are for the analysis.
    We can talk about each step as/when you are available - no rush
    -no pressure!
    
    1. Additional analysis
    
    • Final ChIP-seq table with BldC RNAseq analysis integrated (10 14
    and 22 hr time points only) i.e. each gene listed also has the
    expression information
    
    • Separate table detailing the overlap between the SVEN regulon
    (this work) and SCO regulon (BldC Nature Communications, attached)
    – which genes in this analysis are also identified in the SCO BldC
    work?
    
    • How many “peaks” in the BldC mutant are there - how many of
    these are also in the WT samples (based on local normalisation
    only)?
    
    • To me based on the IGB tracks, the binding of BldC appears
    primarily centred on intergenic regions or at least close to the
    start of a gene. Is this fair?
    
    
    
    2. Final ChIP-seq plots –plots for 10, 14 and 22 hr time points v
    bldC mutant for each of the significant peaks identified. Can the
    dimensions be changed to give “square” panels as in previous
    publications? (easier for building composite figures). The key
    should not be integrated into the images. Produce 2 sets (with and
    without gene identifiers). Note: sometimes it is difficult to see
    where one gene starts and another ends, especially if they are in
    an operon. Can we fix this somehow?
    
    3. RNAseq plots:
    
    a. separate volcano plots for the BldC 10, 14 and 22 hr time
    points (with and without genes identified in ChIP-seq [final
    table] highlighted e.g. red dots) – logFC v apv. lines crossing
    the horizontal axes should indicate when change in expression e.g.
    logFC>+1 or <-1
    
    b. Individual line-plots (showing the absolute expression
    values of WT v mutant at 10hr 14hr and 22hr) for selected genes
    identified in the ChIP-seq analysis (final table). Ideally these
    figures would have the same dimensions as the ChIP-seq figures
    (will help figure construction)
    
    4. Whole Genome Plots  Should be on the same (appropriate scale).
    Produce 2 sets of plots 1) local normalisation only (WT samples
    and mutant) 2) ln_diff
    
    5. Genome plot to include entire dcw cluster region (08495-08580).
    *Does the TSS data correlate with the ChIP-seq data? (Matt will
    look at this)
    
    6. Peak “width” Analysis  - extract sequence beneath “peaks”/most
    Left+Right significant positions and list width in bp. Is there
    any correlation between size of enriched region and the level of
    FC in RNAseq or the number of consensus sites? Predict consensus
    sites…..

#### Thu 04 Oct 2018

BldCRegulon_lndiffFinalTable.csv is the hand edited (by Matt)
version of BldC_ln.csv.

The symlink below was created.

    bldC_rnaseq.txt -> /mnt/isilon/customers/matt_bush/2018_02_05/bldC_rnaseq.txt

Then

~~~ 

rm BldC_chipseq_rnaseq_combined.csv
perl code/chipseq_and_rnaseq.pl -outfile BldC_chipseq_rnaseq_combined.csv \
-chip BldCRegulon_lndiffFinalTable.csv -rnaseq bldC_rnaseq.txt
cp BldC_chipseq_rnaseq_combined.csv ~/mnt/wstemp/matt/

rm BldC_vnz_sco_regulon_overlap.csv
perl code/vnz_sco_bldc_regulon.pl -outfile BldC_vnz_sco_regulon_overlap.csv \
-scofn BldCChIP-chipSCONatureCommunications.csv -vnzfn BldCRegulon_lndiffFinalTable.csv
cp BldC_vnz_sco_regulon_overlap.csv ~/mnt/wstemp/matt/

~~~

#### Mon 08 Oct 2018

~~~ 

rm BldC_vnz_sco_regulon_overlap.csv
perl code/vnz_sco_bldc_regulon2.pl -outfile BldC_vnz_sco_regulon_overlap.csv \
-scofn BldCChIP-chipSCONatureCommunications.csv -vnzfn BldCRegulon_lndiffFinalTable.csv
cp BldC_vnz_sco_regulon_overlap.csv ~/mnt/wstemp/matt/

rm BldC_vnz_sco_regulon_overlap_log2.csv
perl code/vnz_sco_bldc_regulon2.pl -outfile BldC_vnz_sco_regulon_overlap_log2.csv \
-log -scofn BldCChIP-chipSCONatureCommunications.csv \
-vnzfn BldCRegulon_lndiffFinalTable.csv
cp BldC_vnz_sco_regulon_overlap_log2.csv ~/mnt/wstemp/matt/

~~~

#### Wed 10 Oct 2018

Matt emailed yesterday

    Hi Govind,
 
    Thanks for you input over the last few weeks – I’m nearing
    completing the draft of the manuscript. Based on this, I have
    modified the list of priorities for the subsequent analysis.
    I’ve summarised this below. For me now the priority is the
    ChIP-seq panels and RNA-seq volcano plots, then the genome
    plots (whole genome and dcw) and methods. I also want to look
    at the peak “widths” and extract the nucleotide sequences to
    see if I can make any conclusions about the nature of BldC
    binding at individual gene targets.
 
    If you want to chat about any of this at any stage, just let me
    know.
 
    Thanks again for all your help,
 
    Matt


    Additional analysis

     o How many “peaks” in the BldC mutant are not found in the WT
     samples (based on local normalisation only)?

     o To me based on the IGB tracks, the binding of BldC appears
     primarily centred on intergenic regions or at least close to
     the start of a gene. Is this fair?


    Final ChIP-seq plots  –plots for 10, 14 and 22 hr time points v
    bldC mutant for each of the significant peaks identified. Can
    the dimensions be changed to give “square” panels as in
    previous publications? (easier for building composite figures).
    The key should not be integrated into the images. Produce 2
    sets (with and without gene identifiers). Note: sometimes it is
    difficult to see where one gene starts and another ends,
    especially if they are in an operon. Can we fix this somehow?
 
 
    RNAseq plots:
 
    Separate volcano plots for the BldC 10, 14 and 22 hr time
    points (with and without genes identified in ChIP-seq [final
    table] highlighted e.g. red dots) – logFC v apv. lines
    crossing the horizontal axes should indicate when change in
    expression e.g.  logFC>+1 or <-1
 
    Whole Genome Plots  Should be on the same (appropriate scale).
    Produce 2 sets of plots 1) local normalisation only (WT samples
    and mutant) 2) ln_diff
 
    dcw plot  to include ChIP-seq data on entire dcw cluster region
    (08495-08580). *Does the TSS data correlate with the ChIP-seq
    data? (to answer this question it may be necessary to realign
    the TSS data to the vnz chromosome – currently based on Sven15)
 
    Peak “width” Analysis  - extract sequence beneath “peaks”/most
    Left+Right significant positions and list width in bp.
 
    Methods – ChIP-seq and RNA-seq analysis. How was the table
    finally arrived at?
 
#### Fri 19 Oct 2018

Matt emailed today.

> Hi Govind,
> 
> If it's helpful for you, my list of priorities has shifted a
> little (e.g. Whole genome plots whilst useful are not essential
> for submission). Before I can submit this is what I think needs
> to be done (in order of priority):
> 
> 1. Peak “width” Analysis  - extract sequence beneath
> “peaks”/most Left+Right significant positions (only for
> positions identified in 10hr and 14hr timepoints, see attached
> for list of most significant positions)- list: 1) most
> significant position 2) length of sequence below the "peak" in
> bp and 3) the nucleotide sequence itself.
> 
> 2. dcw plot to include ChIP-seq data on entire dcw cluster
> region (08495-08580). This panel can be larger than the standard
> 600 x 600 pixel. Perhaps 600 x 1800 (landscape)?
> 
> 3. hups/ vnz25950 -  can you calculate the minFDR and lndiff
> values based on the 10 and 14 hr timepoints only (falls below
> the cut off at these timepoints)?
> 
> 4. Methods – ChIP-seq and RNA-seq data analysis analysis. As
> similarly generated previously e.g. WhiA/WhiB paper - the
> methods have significantly changed since then
> 
> 5. Upload the RNA-seq and ChIP-seq data to Array Express - 10hr
> and 14hr timepoints only
> 
> 
> Thanks again!
> 
> Matt

The trouble with vnz_25950 is that in 10 and 14 hours there is a
plateau instead of a peak. In 22h there is a peak and at the position
if this peak 10 and 14 hours are below the FDR cut off of 1e-4. See
the figure attached and the table of significances in the csv file in
this region. BldC_10_14_ln.csv is also attached. It excludes 22h and
you can find vnz_25950 near the bottom.


~~~ 

perl code/peak_width.pl -bedgraph aarout/BldC_10.csv -verbose \
-outfile BldC_10_peakwidth.csv -errfile err -- List_for_peak_width_analysis.csv 

less BldC_10_peakwidth.csv

perl code/peak_width.pl -bedgraph aarout/BldC_14.csv -verbose \
-outfile BldC_14_peakwidth.csv -errfile err -- List_for_peak_width_analysis.csv 

less BldC_14_peakwidth.csv

perl code/peak_width.pl -bedgraph aarout/BldC_22.csv -verbose \
-outfile BldC_22_peakwidth.csv -errfile err -- List_for_peak_width_analysis.csv 

less BldC_22_peakwidth.csv

~~~

Used R to combine the above to peakwidths.csv.

~~~ {.bash}

perl code/peak_seq.pl -seqfile vnz.gbk -out peakwidth_seq.csv \
-test 0 \
-- peakwidths.csv


grep '^>' BldC_peaks_max.fna | wc -l
grep '^>' BldC_10and14_peaks.fna | wc -l
grep '^>' BldC_10_peaks.fna | wc -l
grep '^>' BldC_14_peaks.fna | wc -l

wc -l peakwidth_seq.csv

for fn in $(ls --color=never *fna); do
  grep '^>' $fn | wc -l;
done

perl code/peak_seq_22.pl -seqfile vnz.gbk -out peakwidth_seq_22h.csv \
-test 0 \
-- BldC_22_peakwidth.csv

less BldC_22_peakwidth.csv

peakwidth_seq_22h.fna


rm peak_width_analysis.zip $fn
for fn in \
BldC_peaks_max.fna \
BldC_10and14_peaks.fna \
BldC_14_peaks.fna \
BldC_10_peaks.fna \
peakwidth_seq.csv \
BldC_22_peaks.fna \
peakwidth_seq_22h.csv
do
  zip --to-crlf peak_width_analysis.zip $fn
done
  unzip -l peak_width_analysis.zip
#;

re='fna$'
for fn in \
BldC_peaks_max.fna \
BldC_10and14_peaks.fna \
BldC_14_peaks.fna \
BldC_10_peaks.fna \
peakwidth_seq.csv \
BldC_22_peaks.fna \
peakwidth_seq_22h.csv
do


if [[ $fn =~ $re ]]
then
ln=$(grep '^>' $fn | wc -l)
else
ln=$(wc -l $fn)
fi

echo "$fn --- $ln";

done

~~~

#### Thu 25 Oct 2018

~~~ 

pushd aarout
perl ../code/csv_merge3.pl -outfile BldC_10_14_ln.csv -- BldC_10.csv BldC_14.csv
pushd

pushd aarout
perl ../code/csv_merge3.pl -outfile BldC_ln_10_14_5e3.csv -- BldC_{10,14}_5e3.csv
pushd

wc -l BldC_ln_10_14_5e3.csv BldC_10_14_ln.csv \
BldC_ln_5e3.csv BldC_ln.csv


grep '^5673345' BldC_ln_10_14_5e3.csv BldC_10_14_ln.csv \
BldC_ln_5e3.csv BldC_ln.csv

~~~

#### Wed 07 Nov 2018

#### Submission to ArrayExpress

All work related to submission is being done in submission/. It has
its own code/.

~~~ 

mkdir submission

ls fqgz/*{10,14}*


pushd submission
cp -d ../fqgz/*BldC*{10,14}* ./
cp -s ../bedgraph/*{10,14}*_mean.bedgraph ./
popd

~~~

#### Tue 13 Nov 2018

bedgraph_merge.r sourced in an R session. Writes
locallyNormalisedEnrichment.txt.


### SepH ChIP-Seq analysis in 2019.

#### Fri 29 Mar 2019

~~~ 
pushd bedgraph

# This is the deletion strain. So it is the control since polyclonal
# antibody was used.
perl ../code/bedgraph_mean.pl -outfile dseph22_ln_mean.bedgraph -- \
aSepH-22A-IP_ln.bedgraph aSepH-22B-IP_ln.bedgraph

perl ../code/bedgraph_mean.pl -outfile wtseph14_ln_mean.bedgraph -- \
aSepH-WT-14A-IP_ln.bedgraph aSepH-WT-14B-IP_ln.bedgraph

perl ../code/bedgraph_mean.pl -outfile wtseph22_ln_mean.bedgraph -- \
aSepH-WT-22A-IP_ln.bedgraph aSepH-WT-22B-IP_ln.bedgraph



perl ../code/bedgraph_diff.pl -outfile wtseph22_ln_diff.bedgraph -- \
wtseph22_ln_mean.bedgraph  dseph22_ln_mean.bedgraph


perl ../code/bedgraph_diff.pl -outfile wtseph14_ln_diff.bedgraph -- \
wtseph14_ln_mean.bedgraph  dseph22_ln_mean.bedgraph

zip ../seph_IGB *seph*{mean,diff}*bedgraph vnz.bed vnz.fna


pushd
~~~

The zip file made above was sent to Susan.

#### Mon 01 Apr 2019

After doing the stuff in commands.r we get seph22.csv and seph14.csv.

~~~ 
perl code/lir.pl -test 0 -outfile seph22_lir.csv -- seph22.csv
perl code/lir.pl -test 0 -outfile seph14_lir.csv -- seph14.csv
~~~

seph22_lir.csv and seph14_lir.csv were emailed to Susan.


