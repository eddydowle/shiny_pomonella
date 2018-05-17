

library(gplots)
library(ggplot2)
library(reshape)
library(heatmap3)
library(RColorBrewer)
library(WGCNA)
library(shiny)
library(pool)
library(dplyr)
library(RMySQL)
library(reshape2)
library(shinyWidgets)

#library(RMySQL)
con <- dbConnect(MySQL(),user="raglandlab", password="pomonella",dbname="PomUrbanaGrant", host="localhost")


#to start with lets just work on the RNAtable
apple<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/prelimResults/Table.apple.JuneEJD",header=T,row.names=1,sep="\t",stringsAsFactors = F)
haw<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/prelimResults/Table.haw.JuneEJD",header=T,row.names=1,sep="\t",stringsAsFactors = F)
#my connection through uni and the vpn is super slow so Im going to start with running on my machine so I can work out shiny 
#apple<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/Table.apple.JuneEJD",header=T,row.names=1,sep="\t",stringsAsFactors = F)
colnames(apple)
WithinMonthsDE.apple<-data.frame(apple[c(1:5,44)])
#haw<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/Table.haw.JuneEJD",header=T,row.names=1,sep="\t",stringsAsFactors = F)
colnames(haw)
WithinMonthsDE.haw<-data.frame(haw[c(33,1:4,32)])
WithinMonthsDE.apple<-tibble::rownames_to_column(WithinMonthsDE.apple, "gene_id")
WithinMonthsDE.haw<-tibble::rownames_to_column(WithinMonthsDE.haw, "gene_id")
WithinMonthsDE.haw$month2vs2_haw_logFC.x<-0
WithinMonthsDE.apple$month2vs2_apple_logFC.x<-0
WithinMonthsDE.haw<-WithinMonthsDE.haw[,c(1,2,8,3,4,5,6,7)]
WithinMonthsDE.apple<-WithinMonthsDE.apple[,c(1,2,8,3,4,5,6,7)]
colnames(WithinMonthsDE.haw) <- c("gene_id","flybase", "2M","3M","4M","5M","6M","FDR")
colnames(WithinMonthsDE.apple) <- c("gene_id", "flybase", "2M","3M","4M","5M","6M","FDR")
WithinMonthsDE.haw.melt <- melt(WithinMonthsDE.haw,  c("gene_id","flybase","FDR"))
WithinMonthsDE.apple.melt <- melt(WithinMonthsDE.apple, c("gene_id","flybase","FDR"))
WithinMonthsDE.haw.melt$pop<-"haw"
WithinMonthsDE.apple.melt$pop<-"apple"
#together<-rbind(WithinMonthsDE.haw.melt,WithinMonthsDE.apple.melt)
together<-rbind(WithinMonthsDE.apple.melt,WithinMonthsDE.haw.melt)


AcrossMonthsDE.apple<-data.frame(apple[c(1:5,10)])
AcrossMonthsDE.haw<-data.frame(haw[c(33,1:4,9)])
AcrossMonthsDE.apple<-tibble::rownames_to_column(AcrossMonthsDE.apple, "gene_id")
AcrossMonthsDE.haw<-tibble::rownames_to_column(AcrossMonthsDE.haw, "gene_id")
AcrossMonthsDE.haw$month2vs2_haw_logFC.x<-0
AcrossMonthsDE.apple$month2vs2_apple_logFC.x<-0
AcrossMonthsDE.haw<-AcrossMonthsDE.haw[,c(1,2,8,3,4,5,6,7)]
AcrossMonthsDE.apple<-AcrossMonthsDE.apple[,c(1,2,8,3,4,5,6,7)]
colnames(AcrossMonthsDE.haw) <- c("gene_id","flybase", "2M","3M","4M","5M","6M","FDR")
colnames(AcrossMonthsDE.apple) <- c("gene_id", "flybase", "2M","3M","4M","5M","6M","FDR")
AcrossMonthsDE.haw.melt <- melt(AcrossMonthsDE.haw,  c("gene_id","flybase","FDR"))
AcrossMonthsDE.apple.melt <- melt(AcrossMonthsDE.apple, c("gene_id","flybase","FDR"))
AcrossMonthsDE.haw.melt$pop<-"haw"
AcrossMonthsDE.apple.melt$pop<-"apple"
#together<-rbind(WithinMonthsDE.haw.melt,WithinMonthsDE.apple.melt)
together.across<-rbind(AcrossMonthsDE.apple.melt,AcrossMonthsDE.haw.melt)

#what do I want

#table of indels

#table of snps
#table will have total snps in significant and non-significant genes
#total average snp per genes
#number of snp/length of gene?
#number of snp in each category


#then a heatmap of the snps from sig and snps from non-sig

snpeff_effects<-c('3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant', '5_prime_UTR_variant', 'downstream_gene_variant', 'initiator_codon_variant', 'initiator_codon_variant&non_canonical_start_codon', 'initiator_codon_variant&splice_region_variant', 'intergenic_region', 'intragenic_variant', 'intron_variant', 'missense_variant', 'missense_variant&splice_region_variant', 'non_coding_transcript_exon_variant', 'non_coding_transcript_variant', 'splice_acceptor_variant&intron_variant', 'splice_acceptor_variant&splice_donor_variant&intron_variant', 'splice_donor_variant&intron_variant', 'splice_region_variant', 'splice_region_variant&initiator_codon_variant&non_canonical_start_codon', 'splice_region_variant&intron_variant', 'splice_region_variant&non_coding_transcript_exon_variant', 'splice_region_variant&stop_retained_variant', 'splice_region_variant&synonymous_variant', 'start_lost', 'start_lost&splice_region_variant', 'stop_gained', 'stop_gained&splice_region_variant', 'stop_lost', 'stop_lost&splice_region_variant', 'stop_retained_variant', 'synonymous_variant', 'upstream_gene_variant')


#all posible gene ids in mysql databases
con <- dbConnect(MySQL(),user="raglandlab", password="pomonella",dbname="PomUrbanaGrant", host="localhost")

genelist <- dbGetQuery(con,"SELECT feature_alias.gene_id FROM feature_alias;") %>% arrange(.) %>% distinct(.)

together.test<-together %>% filter(.,FDR < 0.05) %>% distinct(.,gene_id,.keep_all = TRUE)  %>% select(.,gene_id)
together.test
avector <- together.test[['gene_id']]
avector.str<-gsub('^',"gene_id='",avector) %>% gsub('$',"' OR ",.) %>% paste(.,collapse = '') %>%  gsub(' OR $',"",.)

sql2<-"SELECT feature_alias.gene_id, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect,snpFisherurbana.urbana_appleave_hawave_fisher_pvalue FROM feature_alias, annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND (%s) AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue < '%f' ;"
sql2 <- sprintf(sql2,avector.str,0.05)

#this takes some time to run
query2b <- dbGetQuery(con, sql2)

#query for genes not in list


#building table
df.results.snp <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("Genes significant FDR RNA", "Genes non-significant FDR RNA", "Difference")
colnames(df.results.snp) <- x

df.results.snp <- rbind(df.results.snp, 'Total_SNPS' = c(nrow(query2b), nrow(query2b.nonsig),"NA"))

df.results.snp <- rbind(df.results.snp, 'SNPS_per_gene' = c(mean(count(query2b,gene_id)[[2]]), mean(count(query2b.nonsig,gene_id)[[2]]),"NA"))

df.results.snp <- rbind(df.results.snp, 'SNPS_per_gene' = c(mean(count(query2b,gene_id)[[2]]), mean(count(query2b.nonsig,gene_id)[[2]]),"NA"))



sql2<-"SELECT feature_alias.gene_id, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect,snpFisherurbana.urbana_appleave_hawave_fisher_pvalue FROM feature_alias, annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) > 1 AND char_length(snpFisherurbana.alt) > 1) AND (%s) AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue < '%f' ;"
sql2 <- sprintf(sql2,avector.str,0.05)

#this takes some time to run
query2_indel <- dbGetQuery(con, sql2)




notsiggenes<-genelist$gene_id[!(genelist$gene_id %in% together.test$gene_id)]
bvector.str<-gsub('^',"gene_id='",notsiggenes) %>% gsub('$',"' OR ",.) %>% paste(.,collapse = '') %>%  gsub(' OR $',"",.)
sql2<-"SELECT feature_alias.gene_id, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect,snpFisherurbana.urbana_appleave_hawave_fisher_pvalue FROM feature_alias, annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND (%s) AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue < '%f' ;"
sql2 <- sprintf(sql2,bvector.str,0.05)

#this takes some time to run
query2b.nonsig <- dbGetQuery(con, sql2)

#indels

sql2<-"SELECT feature_alias.gene_id, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect,snpFisherurbana.urbana_appleave_hawave_fisher_pvalue FROM feature_alias, annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) > 1 AND char_length(snpFisherurbana.alt) > 1) AND (%s) AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue < '%f' ;"
sql2 <- sprintf(sql2,bvector.str,0.05)

#this takes some time to run
query2_indel.nonsig <- dbGetQuery(con, sql2)





test<- dbGetQuery(con,"SELECT feature_alias.gene_id, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect,snpFisherurbana.urbana_appleave_hawave_fisher_pvalue FROM feature_alias, annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND gene_id='gene10001' AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue < '0.05' ;")


test2<- dbGetQuery(con,"SELECT feature_alias.gene_id, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect,snpFisherurbana.urbana_appleave_hawave_fisher_pvalue FROM feature_alias, annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND gene_id='gene17606' AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue < '0.05' ;")








