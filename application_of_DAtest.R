###
# running DAtest and plotting results
###

# load packages
require(DAtest)
require(phyloseq)
require(ggplot2)

# load data (phyloseq)
load("phylo_data.RData")

# select well 1,2,3
phylo_data_1<-phylo_data
phylo_data_1<-prune_samples(grepl("1|2|3", sample_data(phylo_data_1)$sample_name), phylo_data_1)
phylo_data_1<-filter_taxa(phylo_data_1, function(x) sum(x)>0, TRUE)

# agglomerate taxa: rank=phylum
phylo_data_1<-tax_glom(phylo_data_1, taxrank=rank_names(phylo_data_1)[2], NArm=F)

# calculate relative abundances
phylo_data_2<-phylo_data_1
phylo_data_2<-transform_sample_counts(phylo_data_1, function(x) x/sum(x))

# create EAA using flow cytometric cell counts (some methods require integers -> round)
otabEAA<-t(diag(data.frame(sample_data(phylo_data_2))$fc_cells_per_ml) %*% t(otu_table(phylo_data_2)))
colnames(otabEAA)<-rownames(sample_data(phylo_data_2))
phylo_data_3<-phylo_data_2
otu_table(phylo_data_3)<-otu_table(round(otabEAA), taxa_are_rows=T)

# filter low abundance taxa (added to group others)
tot_sum<-sum(rowSums(otu_table(phylo_data_3)))
phylo_data_3<-preDA(phylo_data_3, min.samples=0, min.reads=1, min.abundance=0.0001)

# print summary of resulting data
print(phylo_data_3)
print(paste("Added to others: ", sum(otu_table(phylo_data_3)[rownames(otu_table(phylo_data_3))=="Others"])/tot_sum))

#
# DAtest
#

# spiking pattern
k<-list(c(1,0,0), c(0,1,0), c(0,0,1), c(1,1,1),
       c(2,2,2),
       c(3,0,0), c(0,3,0), c(0,0,3), c(3,3,3),
       c(4,4,4),
       c(5,0,0), c(0,5,0), c(0,0,5), c(5,5,5),
       c(6,6,6))

# effect size
eS<-c(2,4,8,16,32,64,128)

# number of runs for each single testing
nrun<-100

# select which tests to perform
tests<-c("neb", "rai", "per", "bay", "adx", "sam", "qua", "fri", "zpo",
         "znb", "vli", "qpo", "poi", "pea", "wil", "ttt", "ltt", "ltt2",
         "erq", "erq2","ere", "ere2", "msf", "zig", "ds2", "ds2x", "lim",
         "lli", "lli2", "aov", "lao", "lao2", "kru", "lrm", "llm", "llm2",
         "spe","anc", "mva")

# run all tests for all spiking patterns and effect sizes
datest_1<-list()
for(i in eS) {
  eSruns<-lapply(k, function(x) {
    print(paste(i,paste(x, collapse="|")))
    tda=testDA(phylo_data_3, predictor="date_sampling",
               effectSize=i, k=x, R=nrun, core.check=F, relative=F, # relative is set to false
               tests=tests)
    tda=summary(tda)
    tda$spiking_pattern=paste(x, collapse="_")
    tda$effect_size=i
    tda$relative="false"
    tda$paired="false"
    tda$covars="false"
    return(tda)
  })
  datest_1<-c(datest_1, list(eSruns))
}

datest_2=list()
for(i in eS) {
  eSruns=lapply(k, function(x) {
    print(paste(i,paste(x, collapse="|")))
    tda=testDA(phylo_data_3, predictor="date_sampling", paired="sample_name", # using paired
               effectSize=i, k=x, R=nrun, core.check=F, relative=F, # relative is set to false
               tests=tests)
    tda=summary(tda)
    tda$spiking_pattern=paste(x, collapse="_")
    tda$effect_size=i
    tda$relative="false"
    tda$paired="true"
    tda$covars="false"
    return(tda)
  })
  datest_2=c(datest_2, list(eSruns))
}

#
# make plot
#

pldf<-c(datest_1, datest_2)
pldf<-do.call(rbind, do.call(rbind, pldf))
pldf$spiking_pattern<-factor(pldf$spiking_pattern, levels=unique(pldf$spiking_pattern))
#
p1<-ggplot(pldf, aes(effect_size, Score, color=spiking_pattern)) +
  facet_wrap(~paste0(Method, "\n", "relative: ", relative, ", paired: ", paired, ", covars:", covars), ncol = 2) +
  geom_line(alpha=0.5) +
  geom_point(size=2, alpha=0.7) +
  ggtitle("Score = (AUC-0.5) * Spike.detect.rate - FDR") +
  scale_x_continuous(trans="log2") +
  coord_cartesian(ylim=c(0,0.5))
print(p1)


