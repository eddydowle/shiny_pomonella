################snp estimates compared to across genome#############

#across genome estiamtes from snpeff run of all sites

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


#sql2<-"SELECT feature_alias.gene_id, snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM feature_alias, annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"

#drop the feature_alias as intergenic regions have no loc id
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"

#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 ;"

#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 ;"

#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05;"

sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05;"


sql2 <- sprintf(sql2)

query <- dbGetQuery(con, sql2)

#appleave_hawave<-query
#appleearly_applelate<-query
#hawearly_hawlate<-query
#hawearly_hawlate_appleearly_applelate<-query
hawearly_hawlate_appleave_hawave<-query



#getting all our effects out and merging for splice/exon/5prime
query.splice_acceptor_variant_intron_variant <-query %>% filter(., effect=="splice_acceptor_variant&intron_variant")
query.splice_acceptor_variant_splice_donor_variant_intron_variant <-query %>% filter(., effect=="splice_acceptor_variant&splice_donor_variant&intron_variant")
query.splice_acceptor_variant_splice_region_variant_intron_variant <-query %>% filter(., effect=="splice_acceptor_variant&splice_region_variant&intron_variant")
query.splice_donor_variant_intron_variant <-query %>% filter(., effect=="splice_donor_variant&intron_variant")
query.splice_donor_variant_splice_region_variant_intron_variant <-query %>% filter(., effect=="splice_donor_variant&splice_region_variant&intron_variant")
query.splice_region_variant_stop_retained_variant <-query %>% filter(., effect=="splice_region_variant&stop_retained_variant")
query.start_lost_splice_region_variant <-query %>% filter(., effect=="start_lost&splice_region_variant")
query.stop_gained_splice_region_variant <-query %>% filter(., effect=="stop_gained&splice_region_variant")
query.stop_lost_splice_region_variant<-query %>% filter(., effect=="stop_lost&splice_region_variant")
query.initiator_codon_variant_splice_region_variant<-query %>% filter(., effect=="initiator_codon_variant&splice_region_variant")
query.splice_region_variant_intron_variant <-query %>% filter(., effect=="splice_region_variant&intron_variant")
query.missense_variant_splice_region_variant <-query %>% filter(., effect=="missense_variant&splice_region_variant")
query.splice_region_variant_synonymous_variant <-query %>% filter(., effect=="splice_region_variant&synonymous_variant")
query.splice_region_variant_non_coding_transcript_exon_variant <-query %>% filter(., effect=="splice_region_variant&non_coding_transcript_exon_variant")
query.splice_region_variant <-query %>% filter(., effect=="splice_region_variant")

splice<-rbind(query.splice_acceptor_variant_intron_variant,query.splice_acceptor_variant_splice_donor_variant_intron_variant,query.splice_acceptor_variant_splice_region_variant_intron_variant,query.splice_donor_variant_intron_variant,query.splice_donor_variant_splice_region_variant_intron_variant,query.splice_region_variant_stop_retained_variant,query.start_lost_splice_region_variant,query.stop_gained_splice_region_variant,query.stop_lost_splice_region_variant,query.initiator_codon_variant_splice_region_variant,query.splice_region_variant_intron_variant,query.missense_variant_splice_region_variant,query.splice_region_variant_synonymous_variant,query.splice_region_variant_non_coding_transcript_exon_variant,query.splice_region_variant)

splice<-splice %>% distinct()

query.missense_variant <-query %>% filter(., effect=="missense_variant")
query.synonymous_variant <-query %>% filter(., effect=="synonymous_variant")
query.stop_gained <-query %>% filter(., effect=="stop_gained")
query.stop_retained_variant <-query %>% filter(., effect=="stop_retained_variant")
query.rare_amino_acid_variant <-query %>% filter(., effect=="rare_amino_acid_variant")
query.initiator_codon_variant <-query %>% filter(., effect=="initiator_codon_variant")
query.initiator_codon_variant_non_canonical_start_codon <-query %>% filter(., effect=="initiator_codon_variant&non_canonical_start_codon")
query.stop_lost <-query %>% filter(., effect=="stop_lost")
query.start_lost <-query %>% filter(., effect=="start_lost")
query.non_coding_transcript_exon_variant <-query %>% filter(., effect=="non_coding_transcript_exon_variant")

exon<-rbind(query.missense_variant,query.stop_lost,query.start_lost,query.stop_gained,query.stop_retained_variant,query.synonymous_variant,query.rare_amino_acid_variant,query.initiator_codon_variant,query.initiator_codon_variant_non_canonical_start_codon,query.non_coding_transcript_exon_variant)
exon<-exon %>% distinct()

query.non_coding_transcript_variant <-query %>% filter(., effect=="non_coding_transcript_variant")

noncodingtranscript<-query.non_coding_transcript_variant %>% distinct()

query.5_prime_UTR_variant <-query %>% filter(., effect=="5_prime_UTR_variant")
query.5_prime_UTR_premature_start_codon_gain_variant <-query %>% filter(., effect=="5_prime_UTR_premature_start_codon_gain_variant")
prime5<-rbind(query.5_prime_UTR_variant,query.5_prime_UTR_premature_start_codon_gain_variant)
prime5<-prime5 %>% distinct()


query.3_prime_UTR_variant <-query %>% filter(., effect=="3_prime_UTR_variant")
prime3<-query.3_prime_UTR_variant %>% distinct()

query.intron_variant <-query %>% filter(., effect=="intron_variant")
intron<-query.intron_variant  %>% distinct()

query.intragenic_variant <-query %>% filter(., effect=="intragenic_variant")
intragenic<-query.intragenic_variant %>% distinct()

query.upstream_gene_variant <-query %>% filter(., effect=="upstream_gene_variant")
upstream<-query.upstream_gene_variant %>% distinct()

query.downstream_gene_variant <-query %>% filter(., effect=="downstream_gene_variant")
downstream<-query.downstream_gene_variant %>% distinct()

query.intergenic <-query %>% filter(., effect=="intergenic_region")
intergenic<-query.intergenic %>% distinct()

#total amount should be:
query[1:4] %>% unique() %>% nrow #12614

#across total snpeff valid genome sites:
#Splice:              2005240
#Exon:                37736339
#Noncodingtranscript: 219131
#5prime:              3968953
#3prime:              4951758
#intron:              222613089
#intragenic:          1599585
#upstream:            76802582
#downstream:          58072937
#intergenic:          705992930

#total:               1113962544

#global variables across genome:
#tests:
totalsig<-query[1:4] %>% unique() %>% nrow
total<-1113962544

#splice:
testeffect<-2005240
splice[1:4] %>% unique() %>% nrow() #1200
testsig<-splice[1:4] %>% unique() %>% nrow() 
chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
chitest<-chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-as.data.frame(cbind('splice',chitest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#splice: 


#exon
testeffect<-37736339
exon[1:4] %>% setdiff (.,splice[1:4]) %>% unique() %>% nrow() #23627
testsig<-exon[1:4] %>% setdiff (.,splice[1:4]) %>% unique() %>% nrow() 
chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
chitest<-chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('exon',chitest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#exon: 

#noncodingtranscript
testeffect<-219131
test<-rbind(splice[1:4],exon[1:4])
noncodingtranscript[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #219
testsig<-noncodingtranscript[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() 
chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
chitest<-chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('noncodingtranscript',chitest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#noncodingtranscript: p-value = 1

#prime5
testeffect<-3968953
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4])
prime5[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #4084
testsig<-prime5[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() 
chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
chitest<-chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('prime5',chitest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#prime5 p-value = 7.579e-11

#prime3
testeffect<-4951758
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4])
prime3[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #4798
testsig<-prime3[1:4] %>% setdiff (.,test) %>% unique() %>% nrow()
chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
chitest<-chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('prime3',chitest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#prime3 pvalue: p-value = 8.193e-05

#intron
testeffect<-222613089
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4],prime3[1:4])
intron[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #155481
testsig<-intron[1:4] %>% setdiff (.,test) %>% unique() %>% nrow()
chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
chitest<-chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('intron',chitest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#intron: 

#intragenic
testeffect<-1599585
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4],prime3[1:4],intron[1:4])
intragenic[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #614
testsig<-intragenic[1:4] %>% setdiff (.,test) %>% unique() %>% nrow()
chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
chitest<-chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('intragenic',chitest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#intragenic: p-value = 0.08241


#upstream
testeffect<-76802582  
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4],prime3[1:4],intron[1:4],intragenic[1:4])
upstream[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #41157
testsig<-upstream[1:4] %>% setdiff (.,test) %>% unique() %>% nrow()
chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
chitest<-chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('upstream',chitest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#upstream: 



#downstream
testeffect<-58072937    
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4],prime3[1:4],intron[1:4],intragenic[1:4],upstream[1:4])
downstream[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #32781
testsig<-downstream[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() 
chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
chitest<-chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('downstream',chitest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#downstream:

#intergenic
testeffect<-705992930
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4],prime3[1:4],intron[1:4],intragenic[1:4],upstream[1:4],downstream[1:4])
intergenic[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #338476
testsig<-intergenic[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() 
chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
chitest<-chisq.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('intergenic',chitest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
names(output) <- c("effect", "pvalue","percent_genome","percent_sig_snps")
#intergenic: 

write.table(output,'/media/raglandlab/ExtraDrive3/poolseqPom/comparecountssnpsgenome/hawearly_hawlate_appleave_hawave',sep='\t',quote=F,row.names=F)

#appleave_hawave
#appleearly_applelate
#hawearly_hawlate
#hawearly_hawlate_appleearly_applelate
#hawearly_hawlate_appleave_hawave


#across total snpeff valid genome sites:
#Splice:              2005240
#Exon:                37736339
#Noncodingtranscript: 219131
#5prime:              3968953
#3prime:              4951758
#intron:              222613089
#intragenic:          1599585
#upstream:            76802582
#downstream:          58072937
#intergenic:          705992930

#total:               1113962544

#test:
#sig effect       effect-sig effect
#sig-sig effect   total-total sig-(effect-sig effect)
