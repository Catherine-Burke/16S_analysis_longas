#Asaia Beta diversity
library("phyloseq")
library("ggplot2")
#import data to phyloseq
setwd("~/Documents/UTS/skin_microbiome/16S_project/DPI/seq_data/run2_nectar/Asaia_betadiv")
#Define file paths
map<-import_qiime_sample_data("mapping3.txt")
biom_otu_table<- import_biom("otu_table_mc1_w_tax_no_pynast_failures_json.biom","rep_set.tre")
biom_otu_table_v4<-import_biom("otu_table_asaia_v4_no_pynast_failures_json.biom","rep_set_v4.tre")
asaia_data_long <- merge_phyloseq(biom_otu_table,map)
colnames(tax_table(asaia_data_long)) = c("Domain", "Phylum", "Class", "Order", 
                                    "Family", "Genus", "Species")
asaia_data_long
asaia_data_v4 <- merge_phyloseq(biom_otu_table_v4,map)
colnames(tax_table(asaia_data_v4)) = c("Domain", "Phylum", "Class", "Order", 
                                         "Family", "Genus", "Species")
asaia_data_v4
sample_sums(asaia_data_long)
sample_sums(asaia_data_v4)
#Analysis of long seqs rarefied to 100 seqs per sample
asaia_rare100_long<-rarefy_even_depth(asaia_data_long, rngseed=511, replace=FALSE, sample.size=100)
plot_richness(asaia_rare100_long, measures="Observed", x="location", title="observed richness long 16S, 100 seqs per sample")
asaia_wu_ord<-ordinate(asaia_rare100_long, method="PCoA", distance="wunifrac")
asaia_uu_ord<-ordinate(asaia_rare100_long, method="PCoA", distance="unifrac")
plot_ordination(asaia_rare100_long, asaia_wu_ord, color="Source", title="PCoA of weighted unifrac all")+
  scale_color_brewer(palette="Paired")
plot_ordination(asaia_rare100_long, asaia_uu_ord, color="Source", title="PCoA of unweighted unifrac all")+
  scale_color_brewer(palette="Paired")
asaia_rare_sub<-subset_samples(asaia_rare100_long, Source!="Buxton_Peach_A")
sample_names(asaia_rare_sub)
asaia_sub_wu_ord<-ordinate(asaia_rare_sub, method="PCoA", distance="wunifrac")
plot_ordination(asaia_rare_sub, asaia_sub_wu_ord, color="Source", title="PCoA of weighted unifrac without Buxton Peach A")+
  scale_color_brewer(palette="Paired")
asaia_sub_uu_ord<-ordinate(asaia_rare_sub, method="PCoA", distance="unifrac")
plot_ordination(asaia_rare_sub, asaia_sub_uu_ord, color="Source", title="PCoA of unweighted unifrac without Buxton Peach A")+
  scale_color_brewer(palette="Paired")
asaia_rare_sub2<-subset_samples(asaia_rare_sub, !(Source%in%c("Buxton_Peach_A","Buxton_Peach_E")))
sample_names(asaia_rare_sub2)
asaia_sub2_wunifrac<-distance(asaia_rare_sub2, "wunifrac")
asaia_sub2_wu_ord<-ordinate(asaia_rare_sub2, method="PCoA", distance=asaia_sub2_wunifrac)
plot_ordination(asaia_rare_sub2, asaia_sub2_wu_ord, color="Source", title="PCoA of weighted unifrac without Buxton Peach A or E")+
  scale_color_brewer(palette="Paired")
#Analysis of long seqs rarefied to 50 seqs per sample
asaia_rare50<-rarefy_even_depth(asaia_data, rngseed=511, replace=FALSE, sample.size=50)
asaia_wunifrac<-distance(asaia_rare50, "wunifrac")
asaia_wu_ord<-ordinate(asaia_rare50, method="PCoA", distance=asaia_wunifrac)
asaia_uu_ord<-ordinate(asaia_rare50, method="PCoA", distance="unifrac")
plot_ordination(asaia_rare50, asaia_wu_ord, color="Source", title="PCoA of weighted unifrac all\n50 seqs per sample")+
  scale_color_brewer(palette="Paired")
plot_ordination(asaia_rare50, asaia_uu_ord, color="Source", title="PCoA of unweighted unifrac all\n50 seqs per sample")+
  scale_color_brewer(palette="Paired")
asaia_rare_sub<-subset_samples(asaia_rare50, Source!="Buxton_Peach_A")
sample_names(asaia_rare_sub)
asaia_sub_wunifrac<-distance(asaia_rare_sub, "wunifrac")
asaia_sub_wu_ord<-ordinate(asaia_rare_sub, method="PCoA", distance=asaia_sub_wunifrac)
plot_ordination(asaia_rare_sub, asaia_sub_wu_ord, color="Source", title="PCoA of weighted unifrac without Buxton Peach A\n50 seqs per sample")+
  scale_color_brewer(palette="Paired")
asaia_sub_uu_ord<-ordinate(asaia_rare_sub, method="PCoA", distance="unifrac")
plot_ordination(asaia_rare_sub, asaia_sub_uu_ord, color="Source", title="PCoA of unweighted unifrac without Buxton Peach A\n50 seqs per sample")+
  scale_color_brewer(palette="Paired")
asaia_rare_sub2<-subset_samples(asaia_rare_sub, !(Source%in%c("Buxton_Peach_A","Buxton_Peach_E")))
sample_names(asaia_rare_sub2)
asaia_sub2_wunifrac<-distance(asaia_rare_sub2, "wunifrac")
asaia_sub2_wu_ord<-ordinate(asaia_rare_sub2, method="PCoA", distance=asaia_sub2_wunifrac)
plot_ordination(asaia_rare_sub2, asaia_sub2_wu_ord, color="Source", title="PCoA of weighted unifrac without Buxton Peach A or E\n50 seqs per sample")+
  scale_color_brewer(palette="Paired")
asaia_sub2_uunifrac<-distance(asaia_rare_sub2, "unifrac")
asaia_sub2_uu_ord<-ordinate(asaia_rare_sub2, method="PCoA", distance=asaia_sub2_uunifrac)
plot_ordination(asaia_rare_sub2, asaia_sub2_uu_ord, color="Source", title="PCoA of unweighted unifrac without Buxton Peach A or E\n50 seqs per sample")+
  scale_color_brewer(palette="Paired")
#Analysis of long seqs rarefied to 100 seqs per sample
asaia_rare100_v4<-rarefy_even_depth(asaia_data_v4, rngseed=511, replace=FALSE, sample.size=100)
plot_richness(asaia_rare100_v4, measures="Observed", x="location", title="Observed unique sequences of Asaia sp.\nV4 region 100 seqs")
asaia_v4_wu_ord<-ordinate(asaia_rare100_v4, method="PCoA", distance="wunifrac")
asaia_v4_uu_ord<-ordinate(asaia_rare100_v4, method="PCoA", distance="unifrac")
plot_ordination(asaia_rare100_v4, asaia_v4_wu_ord, color="Source", title="PCoA of weighted unifrac all\nV4 region 100 seqs")+
  scale_color_brewer(palette="Paired")
plot_ordination(asaia_rare100_v4, asaia_v4_uu_ord, color="Source", title="PCoA of unweighted unifrac all\nV4 region 100 seqs")+
  scale_color_brewer(palette="Paired")
asaia_v4_rare_sub<-subset_samples(asaia_rare100_v4, Source!="Buxton_Peach_A")
sample_names(asaia_v4_rare_sub)
asaia_v4_sub_wu_ord<-ordinate(asaia_v4_rare_sub, method="PCoA", distance="wunifrac")
plot_ordination(asaia_v4_rare_sub, asaia_v4_sub_wu_ord, color="Source", title="PCoA of weighted unifrac without Buxton Peach A\nV4 region 100 seqs")+
  scale_color_brewer(palette="Paired")
asaia_v4_sub_uu_ord<-ordinate(asaia_v4_rare_sub, method="PCoA", distance="unifrac")
plot_ordination(asaia_v4_rare_sub, asaia_v4_sub_uu_ord, color="Source", title="PCoA of unweighted unifrac without Buxton Peach A\nV4 region 100 seqs")+
  scale_color_brewer(palette="Paired")
asaia_v4_rare_sub2<-subset_samples(asaia_v4_rare_sub, !(Source%in%c("Buxton_Peach_A","Buxton_Peach_E")))
sample_names(asaia_v4_rare_sub2)
asaia_v4_sub2_wu_ord<-ordinate(asaia_rare_sub2, method="PCoA", distance="wunifrac")
plot_ordination(asaia_v4_rare_sub2, asaia_sub2_wu_ord, color="Source", title="PCoA of weighted unifrac without Buxton Peach A or E\nV4 region 100 seqs")+
  scale_color_brewer(palette="Paired")
