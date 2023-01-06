#-----------------------------
# 1/6/2023
# Sohyoung Kim
# kimsohy@mail.nih.gov
# generate_bubble_plot.R v1.0.0
# tested on R/4.1.2 and R/4.2.2
# This script generates bubble plot based on the inputs - provided data on the script and input files.
#
# Library required:
#    library(rtracklayer);
#    library(GenomicRanges)
#    library(scales);
#    library(wesanderson);
#    library(viridis);
#    library(ggplot2)
#-----------------------------
ALLOW_LIBRARY_INSTALL=FALSE;
if(ALLOW_LIBRARY_INSTALL){
    list.of.packages<-c('rtracklayer','GenomicRanges','scales','wesanderson','viridis','ggplot2'); 
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)

    ## if library install fails for rtracklayer and GenomicRanges, try:
    # install.packages("BiocManager")
    # BiocManager::install("rtracklayer")
}
#~~~~~~~~~~~~~~~~~~~~~~~~~
# required inputs:
#~~~~~~~~~~~~~~~~~~~~~~~~~
plot.gene<-"ENSMMUG00000022601";
creb1.chr<-'chr12';
# animal names to analyze and their TOA data
a.names<-c("06BA","T63","T96","TBL","Ti7","TiD","TV9","TK1","TK3","TT7","TiJ");
toa<-c(7, 12,  1,  3,  5, 12, 10, 12, 11,  3,  7);
names(toa)<-a.names;

#~~~~~~~~~~~~~~~~~~~~~~~~~
# FILE required: gene annotation file used
#   gtfname<-"Macaca_mulatta.Mmul_10.103.chr.gtf"
#~~~~~~~~~~~~~~~~~~~~~~~~~
gtfname<-"example_subset_CREB1_only_Macaca_mulatta.Mmul_10.103.chr.gtf";
if(plot.gene != "ENSMMUG00000022601"){
   stop('Please use full gtf data, Macaca_mulatta.Mmul_10.103.chr.gtf');
   stop('Please use full ATAC-Seq data, if the gene of interest does not on chr12.');
}
#~~~~~~~~~~~~~~~~~~~~~~~~~
# FILE required: atac peak file (coordinates) and normalized peak files
#~~~~~~~~~~~~~~~~~~~~~~~~~
fname.atac.coord.only<-'atac_110098_peak_list.csv';
fname.basal<-'atac.norm.counts.basal.chr12.csv';
fname.post<-'atac.norm.counts.post.chr12.csv';

#===============================
# Get CREB1 gene annotation data 
#===============================
library(rtracklayer);
gtfd<-readGFF(gtfname);
#------------
# Define TSS for CREB1, ENSMMUG00000022601
# consider 5' boundary of all genomic features as potential TSS for this gene
# Define -250kb/+250kb range from each TSS
#------------
gtfd.creb1<-gtfd[gtfd$gene_id == plot.gene,c('seqid','start','end','strand','gene_id','source')];
colnames(gtfd.creb1)<-c('chr','st','ed','strand','gene.name','source');
gtfd.creb1$chr<-paste0('chr',gtfd.creb1$chr);
gtfd.creb1[gtfd.creb1$strand == '+','ed']<-gtfd.creb1[gtfd.creb1$strand == '+','st'];
gtfd.creb1[gtfd.creb1$strand == '-','st']<-gtfd.creb1[gtfd.creb1$strand == '-','ed'];
gtfd.creb1u.500kb<-unique(gtfd.creb1);
gtfd.creb1u.500kb$tss<-gtfd.creb1u.500kb$st;
gtfd.creb1u.500kb$st<-gtfd.creb1u.500kb$st-250000;
gtfd.creb1u.500kb$ed<-gtfd.creb1u.500kb$ed+250000;

#===============================
# Load ATAC peak coordinates
#===============================
cdat.m<-read.csv(file=fname.atac.coord.only, check.names=FALSE, row.names=1);
cdat.m$id<-rownames(cdat.m);
cdat.m.chr12<-cdat.m[cdat.m$chr %in% creb1.chr, c('id','chr','st','ed')];


#===============================
# Intersect ATAC peaks and potential TSS -/+250kb ranges.
#===============================
library(GenomicRanges)
gr1 <- GRanges(seqnames=cdat.m.chr12$chr,
                ranges=IRanges(start=cdat.m.chr12$st, end=cdat.m.chr12$ed, names=cdat.m.chr12$id))
gr2<- GRanges(seqnames=gtfd.creb1u.500kb$chr,
                ranges=IRanges(start=gtfd.creb1u.500kb$st , end=gtfd.creb1u.500kb$ed, names=rownames(gtfd.creb1u.500kb)) );

fo <- GenomicRanges::findOverlaps(gr1,gr2)
atac.ind<-data.frame(fo)$queryHits;
g.ind<-data.frame(fo)$subjectHits;

#-----------------
# organize the data into matched.500kb object.
# for multiple TSS mapped for each ATAC peak, pick the most upsrteam TSS site.
#-----------------
matched.500kb<-data.frame(dat1.ind=atac.ind, dat2.ind=g.ind);
matched.500kb$atac.id<-cdat.m.chr12[matched.500kb$dat1.ind,'id'];
matched.500kb$atac.st<-cdat.m.chr12[matched.500kb$dat1.ind,'st'];
matched.500kb$atac.ed<-cdat.m.chr12[matched.500kb$dat1.ind,'ed'];
matched.500kb$gene.name<-gtfd.creb1u.500kb[matched.500kb$dat2.ind, 'gene.name'];
matched.500kb$chr<-cdat.m.chr12[matched.500kb$dat1.ind,'chr'];
matched.500kb$tss<-gtfd.creb1u.500kb[matched.500kb$dat2.ind, 'tss'];
matched.500kb$strand<-gtfd.creb1u.500kb[matched.500kb$dat2.ind, 'strand'];

matched.500kb.s<-unique(matched.500kb[,c('gene.name','atac.id','chr','atac.st','atac.ed','tss','strand')]);
tok<-paste0(matched.500kb.s$gene.name,'.',matched.500kb.s$atac.id);

matched.500kb.input<-c();
for(tii in unique(sort(tok) )){
   wii<-which(tok %in% tii);
   matched.500kb.input<-rbind(matched.500kb.input, matched.500kb.s[wii[1],]);
   if(matched.500kb.s[wii[1],'strand']=='+'){
        matched.500kb.input[nrow(matched.500kb.input),'tss']<-min(matched.500kb.s[wii,'tss']);
   } else if(matched.500kb.s[wii[1],'strand']=='-'){
        matched.500kb.input[nrow(matched.500kb.input),'tss']<-max(matched.500kb.s[wii,'tss']);
   } else {
        stop('Update code to handle strand designation other than + or -.');
   }
}

#===============================
# Load normalized ATAC-Seq accessibility data
# and correlated with the accessibity change between post vs. basal and TOA
#===============================
atac.id<-matched.500kb.input$atac.id;
#------------
# Load normalized ATAC count data
#------------
atac.cnts.n.basal<-read.csv(file=fname.basal, check.names=FALSE, row.names=1);
atac.cnts.n.post<-read.csv(file=fname.post, check.names=FALSE, row.names=1);

atac.cnts.n.basal.creb1<-atac.cnts.n.basal[atac.id,a.names];
atac.cnts.n.post.creb1<-atac.cnts.n.post[atac.id,a.names];
atac.cnts.n.input<-atac.cnts.n.basal.creb1+atac.cnts.n.post.creb1;
# atac.cnts.n.input to be used for x-axis of bubble plots to reflect signal stregth in basal and post

#-----------------
# accessibility change defined as log2 fold change at each site for each animal to correlate with TOA
#-----------------
atac.change<-log2(atac.cnts.n.post.creb1+1)-log2(atac.cnts.n.basal.creb1+1);
cval<-cor(t(rbind(toa[a.names],atac.change[atac.id,a.names])), method="spearman")[(1:nrow(atac.change))+1,1];
if(sum(duplicated(atac.id))!=0){
   stop('atac.id should be unique. Cannot be used for the analysis.');
}
if( sum(matched.500kb.input$atac.id!=atac.id)==0){
   matched.500kb.input$cor<-cval[atac.id];
} else {
   stop('matched.500kb.input$atac.id and atac.id must be in the same order.');
}

#======================
# function definition for bubble plot
#======================
bubble_plot_for_one_gene<-function(tmatch.res=matched.500kb.input, multi_tss=TRUE, tag.thr=100, gtitl='CREB1'){
    if(length(unique(sort(tmatch.res$gene.name)))!=1){
        stop('Please provide data with one single gene.');
    }
    if( sum(c("cor","gene.name","tss","strand","atac.id") %in% colnames(matched.500kb.input)) != 5){
	stop('Please check the format of input data. See matched.500kb.input obj for example input data.');
    } 

    # use only if cor is not na.
    tmatch.res.s<-tmatch.res[!is.na(tmatch.res$cor),];
    if(nrow(tmatch.res.s)==0){
        stop('No non-NA cor values were supplied in the input data.');
    }
    if(!('st' %in% colnames(tmatch.res.s))){
       colnames(tmatch.res.s)<-gsub('atac\\.st','st',colnames(tmatch.res.s));
    }
    if(!('ed' %in% colnames(tmatch.res.s))){
       colnames(tmatch.res.s)<-gsub('atac\\.ed','ed',colnames(tmatch.res.s));
    }
    if(!('st' %in% colnames(tmatch.res.s))){
        stop('Please provide ATAC peak start position either as st or atac.st column head.');
    }
    if(!('ed' %in% colnames(tmatch.res.s))){
        stop('Please provide ATAC peak edart position either as ed or atac.ed column head.');
    }

    if(unique(tmatch.res.s$strand) == 0 | unique(tmatch.res.s$strand)=='+'){ # forward strand
        dis= -(tmatch.res.s$tss - ceiling((tmatch.res.s$ed - tmatch.res.s$st)/2 + tmatch.res.s$st))/1000; # distance toward upstream as -.
    } else if (unique(tmatch.res.s$strand) == 1 | unique(tmatch.res.s$strand)=='-') {
        dis= (tmatch.res.s$tss - ceiling((tmatch.res.s$ed - tmatch.res.s$st)/2 + tmatch.res.s$st))/1000; # distance toward upstream as -.
    } else {
        stop('Error in strand information: see match.res');
    }

    if(unique(tmatch.res.s$strand) == 0 | unique(tmatch.res.s$strand)=='+'){ # forward strand
        dis2= -(min(tmatch.res.s$tss) - ceiling((tmatch.res.s$ed - tmatch.res.s$st)/2 + tmatch.res.s$st))/1000; # distance toward upstream as -.
    } else if (unique(tmatch.res.s$strand) == 1 | unique(tmatch.res.s$strand)=='-') {
        dis2= (max(tmatch.res.s$tss) - ceiling((tmatch.res.s$ed - tmatch.res.s$st)/2 + tmatch.res.s$st))/1000; # distance toward upstream as -.
    } else {
        stop('Error in strand information: see match.res');
    }

    if(multi_tss==TRUE){
       tmatch.res.s2<-data.frame(tmatch.res.s, dist=dis);
    } else {
       tmatch.res.s2<-data.frame(tmatch.res.s, dist=dis2);
    }

    #---------------
    # bin correlation and color code
    #---------------
    library(scales);
    library(wesanderson);
    library(viridis);

    bins<-c(-1, seq(-.5,.5,by=1/10), 1);
    col.toa12.p = c('darkgreen', 'blue', 'royalblue',  as.character( wes_palette("Zissou1", 7, type = "continuous") ), 'red3','purple');
    col.toa12 = alpha(col.toa12.p, abs(bins[ -which(bins==0)]) + 0.3 );

    bins[1]= -1.1;
    bins[length(bins)]<-1.1;
    bin.span<-paste0(format(bins, digit=1)[1:12],'<cor<=',format(bins, digit=1)[2:13]);
    bin.span[1]<-"-1<=cor<=-0.5";
    bin.span[12]<-"0.5<cor<=1";

    oi3<-rev(order(tmatch.res.s$cor));
    tmatch.res.s4<-data.frame(tmatch.res.s2[oi3,], tags=rowSums(atac.cnts.n.input[ tmatch.res.s[oi3,'atac.id'], a.names]));

    tmatch.res.s4$col<-rep(NA,nrow(tmatch.res.s4));
    tmatch.res.s4$range<-rep(NA,nrow(tmatch.res.s4));

    for(ii in 1:nrow(tmatch.res.s4))
    {
        wii<-max(which(tmatch.res.s4$cor[ii] - bins > 0));
        tmatch.res.s4$col[ii]<-col.toa12[wii];
        tmatch.res.s4$range[ii]<-bin.span[wii];
    }
    tmatch.res.s4$range<-factor(tmatch.res.s4$range, levels=bin.span[ bin.span %in% tmatch.res.s4$range ]);
    tmatch.res.s4$col<-factor(tmatch.res.s4$col, levels=col.toa12[ col.toa12 %in% tmatch.res.s4$col ]);

    #----------
    # use the thresholded data to identify best atac peak based on the strongest correlation value.
    #----------
    tmatch.res.s4.s<-tmatch.res.s4[tmatch.res.s4$tags>tag.thr,];
    rn.best<-tmatch.res.s4.s$atac.id[which.max(abs(tmatch.res.s4.s$cor))];

    #----------
    # plot bubbles
    #----------
    library(ggplot2)
    p<-ggplot(tmatch.res.s4, aes(x=tags, y=dist, color=range)) + geom_point( aes(size = abs(cor)) );

    bquote.xlab.exp<-bquote('' %<-% Upsteam - Distance~~from~~TSS~~(kb) - Downstream %->% '' );
    p2<-p +  scale_x_continuous(name='Sum of normalized tag counts in basal and post') +
        scale_y_continuous(name=bquote.xlab.exp) +
      scale_color_manual(name='Correlation', labels=levels(tmatch.res.s4$range), values=levels(tmatch.res.s4$col) ) +
      theme_bw(base_size = 14) + ggtitle( gtitl );

    wii<-which(tmatch.res.s4$atac.id %in% rn.best);
    p3<-p2 + geom_point(aes(x=tmatch.res.s4[wii, 'tags'], y=tmatch.res.s4[wii, 'dist']), shape=0, color="magenta", size=9 ) +
                geom_vline(xintercept = tag.thr, linetype="dashed",
                color = "blue", size=1)

    re<-{};
    re$data<-tmatch.res.s4;
    re$plot<-p3;
    re$atac.best<-tmatch.res.s4[tmatch.res.s4$atac.id %in% rn.best,];
    re;
}


re<-bubble_plot_for_one_gene(tmatch.res=matched.500kb.input, multi_tss=TRUE, tag.thr=100, gtitl='CREB1')

#---------------
# plot figure using log scale
#---------------
p<-re$plot + scale_x_continuous(trans = 'log10', name='Sum of normalized tag counts in basal and post');

dev.new()
print(p);

head(re$data)
message('ATAC site with the strongest correlation to TOA: ', re$atac.best$chr,':',
		re$atac.best$st,'-',  re$atac.best$ed);

ggsave(file='bubble-test.pdf', p);




