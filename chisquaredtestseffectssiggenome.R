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

#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05;"

#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05;"

#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 ;"


sql2 <- sprintf(sql2)

query <- dbGetQuery(con, sql2)

#appleave_hawave<-query
#appleearly_applelate<-query
#hawearly_hawlate<-query
#hawearly_hawlate_appleearly_applelate<-query
#hawearly_hawlate_appleave_hawave<-query
appleearly_appleave_hawave<-query
#hawearly_hawlate_appleave_hawave_appleearly_applelate<-query



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
query[1:4] %>% unique() %>% nrow #

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
splice[1:4] %>% unique() %>% nrow() #
testsig<-splice[1:4] %>% unique() %>% nrow() 
fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
fishtest<-fisher.test(matrix(c(testsig,(testeffect-testsig),(totalsig-testsig),(total-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-as.data.frame(cbind('splice',fishtest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#splice: 

#matrix(c(testsig,(testeffect-testsig),(totalsig-testsig),(total-(testsig-testeffect)-(totalsig-testsig))))
#exon
testeffect<-37736339
exon[1:4] %>% setdiff (.,splice[1:4]) %>% unique() %>% nrow() #
testsig<-exon[1:4] %>% setdiff (.,splice[1:4]) %>% unique() %>% nrow() 
fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
fishtest<-fisher.test(matrix(c(testsig,(testeffect-testsig),(totalsig-testsig),(total-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('exon',fishtest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#exon: 

#noncodingtranscript
testeffect<-219131
test<-rbind(splice[1:4],exon[1:4])
noncodingtranscript[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #
testsig<-noncodingtranscript[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() 
fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
fishtest<-fisher.test(matrix(c(testsig,(testeffect-testsig),(totalsig-testsig),(total-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('noncodingtranscript',fishtest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#noncodingtranscript: 

#prime5
testeffect<-3968953
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4])
prime5[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #4084
testsig<-prime5[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() 
fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
fishtest<-fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('prime5',fishtest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#prime5 

#prime3
testeffect<-4951758
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4])
prime3[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #4798
testsig<-prime3[1:4] %>% setdiff (.,test) %>% unique() %>% nrow()
fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
fishtest<-fisher.test(matrix(c(testsig,(testeffect-testsig),(totalsig-testsig),(total-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('prime3',fishtest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#prime3 pvalue: 

#intron
testeffect<-222613089
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4],prime3[1:4])
intron[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #155481
testsig<-intron[1:4] %>% setdiff (.,test) %>% unique() %>% nrow()
fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
fishtest<-fisher.test(matrix(c(testsig,(testeffect-testsig),(totalsig-testsig),(total-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('intron',fishtest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#intron: 

#intragenic
testeffect<-1599585
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4],prime3[1:4],intron[1:4])
intragenic[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #614
testsig<-intragenic[1:4] %>% setdiff (.,test) %>% unique() %>% nrow()
fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
fishtest<-fisher.test(matrix(c(testsig,(testeffect-testsig),(totalsig-testsig),(total-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('intragenic',fishtest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#intragenic: 


#upstream
testeffect<-76802582  
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4],prime3[1:4],intron[1:4],intragenic[1:4])
upstream[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #41157
testsig<-upstream[1:4] %>% setdiff (.,test) %>% unique() %>% nrow()
fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
fishtest<-fisher.test(matrix(c(testsig,(testeffect-testsig),(totalsig-testsig),(total-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('upstream',fishtest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#upstream: 



#downstream
testeffect<-58072937    
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4],prime3[1:4],intron[1:4],intragenic[1:4],upstream[1:4])
downstream[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #32781
testsig<-downstream[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() 
fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
fishtest<-fisher.test(matrix(c(testsig,(testeffect-testsig),(totalsig-testsig),(total-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('downstream',fishtest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
#downstream:

#intergenic
testeffect<-705992930
test<-rbind(splice[1:4],exon[1:4],noncodingtranscript[1:4],prime5[1:4],prime3[1:4],intron[1:4],intragenic[1:4],upstream[1:4],downstream[1:4])
intergenic[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() #338476
testsig<-intergenic[1:4] %>% setdiff (.,test) %>% unique() %>% nrow() 
fisher.test(matrix(c(testsig,testeffect,(totalsig-testsig),(total-testsig-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
fishtest<-fisher.test(matrix(c(testsig,(testeffect-testsig),(totalsig-testsig),(total-testeffect-(totalsig-testsig))),nrow=2,byrow=T))
output<-rbind(output,cbind('intergenic',fishtest$p.value,(testeffect/total)*100,(testsig/totalsig)*100))
names(output) <- c("effect", "pvalue","percent_genome","percent_sig_snps")
#intergenic: 

#write.table(output,'/media/raglandlab/ExtraDrive3/poolseqPom/comparecountssnpsgenome/appleearly_applelate_appleave_hawave',sep='\t',quote=F,row.names=F)

#appleave_hawave
#appleearly_applelate
#hawearly_hawlate
#hawearly_hawlate_appleearly_applelate
#hawearly_hawlate_appleave_hawave
#appleearly_applelate_appleave_hawave
#hawearly_hawlate_appleave_hawave_appleearly_applelate

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


#is there an elevated tajD in snp regions

#so for apple/haw ave what is the average taj D across all genes

#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"


#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) = 1 AND char_length(snpFisherurbana.alt) = 1) AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"
#LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id

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


sql <- c( "SELECT feature_alias.gene_id,feature_alias.loc,feature_alias.Flybase_tophit,feature_alias.Flybase_FBgn,feature_alias.Flybase_gene_symbol,urbanapoolTajDgene.AppleAve_TajD,urbanapoolTajDgene.AppleEarly_TajD,urbanapoolTajDgene.AppleLate_TajD,urbanapoolTajDgene.HawAve_TajD,urbanapoolTajDgene.HawEarly_TajD,urbanapoolTajDgene.HawLate_TajD FROM feature_alias LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id;")
query <- dbGetQuery(con, sql)

#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 ;"
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 ;"
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 ;"
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"
sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05;"

query2 <- dbGetQuery(con, sql2)

#should it just be if its in a exon? or just associated?




#getting all our effects out and merging for splice/exon/5prime
query.splice_acceptor_variant_intron_variant <-query2 %>% filter(., effect=="splice_acceptor_variant&intron_variant")
query.splice_acceptor_variant_splice_donor_variant_intron_variant <-query2 %>% filter(., effect=="splice_acceptor_variant&splice_donor_variant&intron_variant")
query.splice_acceptor_variant_splice_region_variant_intron_variant <-query2 %>% filter(., effect=="splice_acceptor_variant&splice_region_variant&intron_variant")
query.splice_donor_variant_intron_variant <-query2 %>% filter(., effect=="splice_donor_variant&intron_variant")
query.splice_donor_variant_splice_region_variant_intron_variant <-query2 %>% filter(., effect=="splice_donor_variant&splice_region_variant&intron_variant")
query.splice_region_variant_stop_retained_variant <-query2 %>% filter(., effect=="splice_region_variant&stop_retained_variant")
query.start_lost_splice_region_variant <-query2 %>% filter(., effect=="start_lost&splice_region_variant")
query.stop_gained_splice_region_variant <-query2 %>% filter(., effect=="stop_gained&splice_region_variant")
query.stop_lost_splice_region_variant<-query2 %>% filter(., effect=="stop_lost&splice_region_variant")
query.initiator_codon_variant_splice_region_variant<-query2 %>% filter(., effect=="initiator_codon_variant&splice_region_variant")
query.splice_region_variant_intron_variant <-query2 %>% filter(., effect=="splice_region_variant&intron_variant")
query.missense_variant_splice_region_variant <-query2 %>% filter(., effect=="missense_variant&splice_region_variant")
query.splice_region_variant_synonymous_variant <-query2 %>% filter(., effect=="splice_region_variant&synonymous_variant")
query.splice_region_variant_non_coding_transcript_exon_variant <-query2 %>% filter(., effect=="splice_region_variant&non_coding_transcript_exon_variant")
query.splice_region_variant <-query2 %>% filter(., effect=="splice_region_variant")

splice<-rbind(query.splice_acceptor_variant_intron_variant,query.splice_acceptor_variant_splice_donor_variant_intron_variant,query.splice_acceptor_variant_splice_region_variant_intron_variant,query.splice_donor_variant_intron_variant,query.splice_donor_variant_splice_region_variant_intron_variant,query.splice_region_variant_stop_retained_variant,query.start_lost_splice_region_variant,query.stop_gained_splice_region_variant,query.stop_lost_splice_region_variant,query.initiator_codon_variant_splice_region_variant,query.splice_region_variant_intron_variant,query.missense_variant_splice_region_variant,query.splice_region_variant_synonymous_variant,query.splice_region_variant_non_coding_transcript_exon_variant,query.splice_region_variant)

splice<-splice %>% distinct()

query.missense_variant <-query2 %>% filter(., effect=="missense_variant")
query.synonymous_variant <-query2 %>% filter(., effect=="synonymous_variant")
query.stop_gained <-query2 %>% filter(., effect=="stop_gained")
query.stop_retained_variant <-query2 %>% filter(., effect=="stop_retained_variant")
query.rare_amino_acid_variant <-query2 %>% filter(., effect=="rare_amino_acid_variant")
query.initiator_codon_variant <-query2 %>% filter(., effect=="initiator_codon_variant")
query.initiator_codon_variant_non_canonical_start_codon <-query2 %>% filter(., effect=="initiator_codon_variant&non_canonical_start_codon")
query.stop_lost <-query2 %>% filter(., effect=="stop_lost")
query.start_lost <-query2 %>% filter(., effect=="start_lost")
query.non_coding_transcript_exon_variant <-query2 %>% filter(., effect=="non_coding_transcript_exon_variant")

exon<-rbind(query.missense_variant,query.stop_lost,query.start_lost,query.stop_gained,query.stop_retained_variant,query.synonymous_variant,query.rare_amino_acid_variant,query.initiator_codon_variant,query.initiator_codon_variant_non_canonical_start_codon,query.non_coding_transcript_exon_variant)
exon<-exon %>% distinct()

query.non_coding_transcript_variant <-query2 %>% filter(., effect=="non_coding_transcript_variant")

noncodingtranscript<-query.non_coding_transcript_variant %>% distinct()

query.5_prime_UTR_variant <-query2 %>% filter(., effect=="5_prime_UTR_variant")
query.5_prime_UTR_premature_start_codon_gain_variant <-query2 %>% filter(., effect=="5_prime_UTR_premature_start_codon_gain_variant")
prime5<-rbind(query.5_prime_UTR_variant,query.5_prime_UTR_premature_start_codon_gain_variant)
prime5<-prime5 %>% distinct()


query.3_prime_UTR_variant <-query2 %>% filter(., effect=="3_prime_UTR_variant")
prime3<-query.3_prime_UTR_variant %>% distinct()

query.intron_variant <-query2 %>% filter(., effect=="intron_variant")
intron<-query.intron_variant  %>% distinct()

query.intragenic_variant <-query2 %>% filter(., effect=="intragenic_variant")
intragenic<-query.intragenic_variant %>% distinct()

query.upstream_gene_variant <-query2 %>% filter(., effect=="upstream_gene_variant")
upstream<-query.upstream_gene_variant %>% distinct()

query.downstream_gene_variant <-query2 %>% filter(., effect=="downstream_gene_variant")
downstream<-query.downstream_gene_variant %>% distinct()

query.intergenic <-query2 %>% filter(., effect=="intergenic_region")
intergenic<-query.intergenic %>% distinct()

#maybe just ignore intergenic snps as there loc dont mean much
test<-rbind(splice,exon,noncodingtranscript,prime5,prime3,intron,intragenic,upstream,downstream)

#apple ave & haw ave
#Across everything:
query %>% dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(abs(AppleAve_TajD),na.rm = TRUE),HawAve_TajD_mean=mean(abs(HawAve_TajD),na.rm = TRUE),AppleEarly_TajD_mean=mean(abs(AppleEarly_TajD),na.rm = TRUE), AppleLate_TajD_mean=mean(abs(AppleLate_TajD),na.rm = TRUE), HawEarly_TajD_mean=mean(abs(HawEarly_TajD),na.rm = TRUE),HawLate_TajD_mean=mean(abs(HawLate_TajD),na.rm = TRUE))
#mean abs tajima D value
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.5265845        0.5454203            0.6013117           0.5603815
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5373945         0.5180558
#negative
query %>% filter (.,AppleAve_TajD < 0 & HawAve_TajD < 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          -0.612734       -0.6495739            -0.699358          -0.6060358
#HawEarly_TajD_mean HawLate_TajD_mean
#1         -0.6073344        -0.5523958
#positive
query %>% filter (.,AppleAve_TajD > 0 & HawAve_TajD > 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.5516051        0.5282399            0.3336589           0.3155432
#HawEarly_TajD_mean HawLate_TajD_mean
#1           0.279235         0.2775189


#baseline abs
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5265845        0.5454203            0.6013117           0.5603815  0.5373945         0.5180558
#baseline neg
#1          -0.612734       -0.6495739            -0.699358          -0.6060358 -0.6073344        -0.5523958
#baseline pos
#1          0.5516051        0.5282399            0.3336589           0.3155432 0.279235         0.2775189


#across loc associated sig snps
#abs goes down
query %>% filter (.,loc %in% test$loc) %>% dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(abs(AppleAve_TajD),na.rm = TRUE),HawAve_TajD_mean=mean(abs(HawAve_TajD),na.rm = TRUE),AppleEarly_TajD_mean=mean(abs(AppleEarly_TajD),na.rm = TRUE), AppleLate_TajD_mean=mean(abs(AppleLate_TajD),na.rm = TRUE), HawEarly_TajD_mean=mean(abs(HawEarly_TajD),na.rm = TRUE),HawLate_TajD_mean=mean(abs(HawLate_TajD),na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.4647875        0.4924847            0.5434646           0.5034609
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.4812118         0.4670734
#negative goes up
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD < 0 & HawAve_TajD < 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1         -0.5356646       -0.5793964           -0.6339443          -0.5604866
#HawEarly_TajD_mean HawLate_TajD_mean
#1         -0.5525746        -0.5040419
#positive goes down
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD > 0 & HawAve_TajD > 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.5302562        0.5001141            0.3443894           0.2756307
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.2567729         0.2501515

#baseline abs
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5265845        0.5454203            0.6013117           0.5603815  0.5373945         0.5180558
#baseline neg
#1          -0.612734       -0.6495739            -0.699358          -0.6060358 -0.6073344        -0.5523958
#baseline pos
#1          0.5516051        0.5282399            0.3336589           0.3155432 0.279235         0.2775189

#Apple Early/Late
#abs sig snps
query %>% filter (.,loc %in% test$loc) %>% dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(abs(AppleAve_TajD),na.rm = TRUE),HawAve_TajD_mean=mean(abs(HawAve_TajD),na.rm = TRUE),AppleEarly_TajD_mean=mean(abs(AppleEarly_TajD),na.rm = TRUE), AppleLate_TajD_mean=mean(abs(AppleLate_TajD),na.rm = TRUE), HawEarly_TajD_mean=mean(abs(HawEarly_TajD),na.rm = TRUE),HawLate_TajD_mean=mean(abs(HawLate_TajD),na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.5135464        0.5333521            0.5892989           0.5463974
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5256378         0.5058517
#negative sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD < 0 & HawAve_TajD < 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1         -0.5925766        -0.632095           -0.6844495          -0.5904917
#HawEarly_TajD_mean HawLate_TajD_mean
#1         -0.5946681        -0.5413874
#positive sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD > 0 & HawAve_TajD > 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#  AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.5461312        0.5252799            0.3341139           0.3107599
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.2753513         0.2729742

#baseline abs
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5265845        0.5454203            0.6013117           0.5603815  0.5373945         0.5180558
#baseline neg
#1          -0.612734       -0.6495739            -0.699358          -0.6060358 -0.6073344        -0.5523958
#baseline pos
#1          0.5516051        0.5282399            0.3336589           0.3155432 0.279235         0.2775189

#Haw Early/Late
#abs sig snps
query %>% filter (.,loc %in% test$loc) %>% dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(abs(AppleAve_TajD),na.rm = TRUE),HawAve_TajD_mean=mean(abs(HawAve_TajD),na.rm = TRUE),AppleEarly_TajD_mean=mean(abs(AppleEarly_TajD),na.rm = TRUE), AppleLate_TajD_mean=mean(abs(AppleLate_TajD),na.rm = TRUE), HawEarly_TajD_mean=mean(abs(HawEarly_TajD),na.rm = TRUE),HawLate_TajD_mean=mean(abs(HawLate_TajD),na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.5109921        0.5305701            0.5860236           0.5436952
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5216462         0.5019405
#negative sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD < 0 & HawAve_TajD < 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1         -0.5876549       -0.6266445           -0.6799112          -0.5913814
#HawEarly_TajD_mean HawLate_TajD_mean
#1         -0.5901502        -0.5382934
#positive sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD > 0 & HawAve_TajD > 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.5519408        0.5279052            0.3464879           0.3126028
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.2780832         0.2692894


#baseline abs
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5265845        0.5454203            0.6013117           0.5603815  0.5373945         0.5180558
#baseline neg
#1          -0.612734       -0.6495739            -0.699358          -0.6060358 -0.6073344        -0.5523958
#baseline pos
#1          0.5516051        0.5282399            0.3336589           0.3155432 0.279235         0.2775189

#sig AppleEarlyLate HawEarlyLate
query %>% filter (.,loc %in% test$loc) %>% dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(abs(AppleAve_TajD),na.rm = TRUE),HawAve_TajD_mean=mean(abs(HawAve_TajD),na.rm = TRUE),AppleEarly_TajD_mean=mean(abs(AppleEarly_TajD),na.rm = TRUE), AppleLate_TajD_mean=mean(abs(AppleLate_TajD),na.rm = TRUE), HawEarly_TajD_mean=mean(abs(HawEarly_TajD),na.rm = TRUE),HawLate_TajD_mean=mean(abs(HawLate_TajD),na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.4757942         0.482831            0.5140443           0.4651442
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.4579793         0.4375402
#negative sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD < 0 & HawAve_TajD < 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1         -0.5162645       -0.5578691           -0.6005428          -0.5168822
#HawEarly_TajD_mean HawLate_TajD_mean
#1          -0.539155         -0.485899
#positive sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD > 0 & HawAve_TajD > 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.5852416        0.5572189            0.4081064           0.2514143
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.2931511          0.213704

#baseline abs
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5265845        0.5454203            0.6013117           0.5603815  0.5373945         0.5180558
#baseline neg
#1          -0.612734       -0.6495739            -0.699358          -0.6060358 -0.6073344        -0.5523958
#baseline pos
#1          0.5516051        0.5282399            0.3336589           0.3155432 0.279235         0.2775189

#sig HawEarlyLate AppleAveHawAve
query %>% filter (.,loc %in% test$loc) %>% dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(abs(AppleAve_TajD),na.rm = TRUE),HawAve_TajD_mean=mean(abs(HawAve_TajD),na.rm = TRUE),AppleEarly_TajD_mean=mean(abs(AppleEarly_TajD),na.rm = TRUE), AppleLate_TajD_mean=mean(abs(AppleLate_TajD),na.rm = TRUE), HawEarly_TajD_mean=mean(abs(HawEarly_TajD),na.rm = TRUE),HawLate_TajD_mean=mean(abs(HawLate_TajD),na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.4662917        0.4809814            0.5163435           0.4576273
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.4400501         0.4311533
#negative sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD < 0 & HawAve_TajD < 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1         -0.5219608       -0.5752264           -0.5950885          -0.5361535
#HawEarly_TajD_mean HawLate_TajD_mean
#1         -0.5400254        -0.5057884
#positive sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD > 0 & HawAve_TajD > 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.5427907        0.5105214            0.5142971           0.1836588
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.2511309         0.1357011

#baseline abs
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5265845        0.5454203            0.6013117           0.5603815  0.5373945         0.5180558
#baseline neg
#1          -0.612734       -0.6495739            -0.699358          -0.6060358 -0.6073344        -0.5523958
#baseline pos
#1          0.5516051        0.5282399            0.3336589           0.3155432 0.279235         0.2775189


#sig AppleEarlyLate AppleAveHawAve
query %>% filter (.,loc %in% test$loc) %>% dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(abs(AppleAve_TajD),na.rm = TRUE),HawAve_TajD_mean=mean(abs(HawAve_TajD),na.rm = TRUE),AppleEarly_TajD_mean=mean(abs(AppleEarly_TajD),na.rm = TRUE), AppleLate_TajD_mean=mean(abs(AppleLate_TajD),na.rm = TRUE), HawEarly_TajD_mean=mean(abs(HawEarly_TajD),na.rm = TRUE),HawLate_TajD_mean=mean(abs(HawLate_TajD),na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.4629547        0.5034134            0.5250669           0.4458792
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.4830394         0.4596034
#negative sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD < 0 & HawAve_TajD < 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1         -0.4884025       -0.5536116           -0.5864584          -0.5326133
#HawEarly_TajD_mean HawLate_TajD_mean
#1         -0.5554913        -0.5309627
#positive sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD > 0 & HawAve_TajD > 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.6197486        0.6281556            0.5183894           0.2621068
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.4564026         0.2837647

#baseline abs
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5265845        0.5454203            0.6013117           0.5603815  0.5373945         0.5180558
#baseline neg
#1          -0.612734       -0.6495739            -0.699358          -0.6060358 -0.6073344        -0.5523958
#baseline pos
#1          0.5516051        0.5282399            0.3336589           0.3155432 0.279235         0.2775189


#sig all three
#abs sig snps
query %>% filter (.,loc %in% test$loc) %>% dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(abs(AppleAve_TajD),na.rm = TRUE),HawAve_TajD_mean=mean(abs(HawAve_TajD),na.rm = TRUE),AppleEarly_TajD_mean=mean(abs(AppleEarly_TajD),na.rm = TRUE), AppleLate_TajD_mean=mean(abs(AppleLate_TajD),na.rm = TRUE), HawEarly_TajD_mean=mean(abs(HawEarly_TajD),na.rm = TRUE),HawLate_TajD_mean=mean(abs(HawLate_TajD),na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.5186061        0.5362727            0.6104837           0.4644127
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.5334263         0.4835692
#negative sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD < 0 & HawAve_TajD < 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1         -0.5290029       -0.5157221           -0.5803385          -0.4304185
#HawEarly_TajD_mean HawLate_TajD_mean
#1         -0.6741393         -0.449426
#positive sig snps
query %>% filter (.,loc %in% test$loc) %>% filter (.,AppleAve_TajD > 0 & HawAve_TajD > 0) %>%dplyr::select(.,loc,Flybase_gene_symbol,Flybase_FBgn,AppleAve_TajD,HawAve_TajD,AppleEarly_TajD,AppleLate_TajD,HawEarly_TajD,HawLate_TajD) %>% dplyr::summarise( AppleAve_TajD_mean=mean(AppleAve_TajD,na.rm = TRUE),HawAve_TajD_mean=mean(HawAve_TajD,na.rm = TRUE),AppleEarly_TajD_mean=mean(AppleEarly_TajD,na.rm = TRUE), AppleLate_TajD_mean=mean(AppleLate_TajD,na.rm = TRUE), HawEarly_TajD_mean=mean(HawEarly_TajD,na.rm = TRUE),HawLate_TajD_mean=mean(HawLate_TajD,na.rm = TRUE))
#AppleAve_TajD_mean HawAve_TajD_mean AppleEarly_TajD_mean AppleLate_TajD_mean
#1          0.6083999        0.5984935            0.6308199           0.1975721
#HawEarly_TajD_mean HawLate_TajD_mean
#1          0.4109765         0.1894574



#scanning sig snps through david


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

#mysql -u raglandlab -p
#pomonella

#query <- dbGetQuery(con, sql)

sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, feature_alias.gene_id,feature_alias.Flybase_tophit,feature_alias.Flybase_FBgn,feature_alias.Flybase_gene_symbol, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, feature_alias, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND annotation.loc = feature_alias.loc AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"

#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 ;"
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 ;"
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 ;"
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"
#sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05;"

#biologically signficant
sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, feature_alias.gene_id,feature_alias.Flybase_tophit,feature_alias.Flybase_FBgn,feature_alias.Flybase_gene_symbol, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust, urbanapoolmaf.urbana_appleave_Maj, urbanapoolmaf.urbana_hawave_Maj FROM annotation, feature_alias, snpFisherurbana, urbanapoolmaf WHERE annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.snpId = urbanapoolmaf.snpId AND annotation.loc = feature_alias.loc AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05  AND ABS(urbanapoolmaf.urbana_appleave_Maj - urbanapoolmaf.urbana_hawave_Maj) > 0.5;"


#same direction
sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, feature_alias.gene_id,feature_alias.Flybase_tophit,feature_alias.Flybase_FBgn,feature_alias.Flybase_gene_symbol, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust, urbanapoolmaf.urbana_appleave_Maj, urbanapoolmaf.urbana_hawave_Maj, urbanapoolmaf.urbana_appleearly_Maj, urbanapoolmaf.urbana_applelate_Maj, urbanapoolmaf.urbana_hawearly_Maj, urbanapoolmaf.urbana_hawlate_Maj FROM annotation, feature_alias, snpFisherurbana,urbanapoolmaf WHERE annotation.snpId = snpFisherurbana.snpId AND annotation.loc = feature_alias.loc AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05;"
sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, feature_alias.gene_id,feature_alias.Flybase_tophit,feature_alias.Flybase_FBgn,feature_alias.Flybase_gene_symbol, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust, snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust,snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust,urbanapoolmaf.urbana_appleave_Maj, urbanapoolmaf.urbana_hawave_Maj, urbanapoolmaf.urbana_appleearly_Maj, urbanapoolmaf.urbana_applelate_Maj, urbanapoolmaf.urbana_hawearly_Maj, urbanapoolmaf.urbana_hawlate_Maj FROM annotation, feature_alias, snpFisherurbana,urbanapoolmaf WHERE snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 AND urbanapoolmaf.snpId = snpFisherurbana.snpId AND annotation.snpId = snpFisherurbana.snpId AND annotation.loc = feature_alias.loc ;"

#DE and sig genes
sql2<-"SELECT snpFisherurbana.scaffold, snpFisherurbana.position, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect, feature_alias.gene_id,feature_alias.Flybase_tophit,feature_alias.Flybase_FBgn,feature_alias.Flybase_gene_symbol, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust FROM annotation, feature_alias, snpFisherurbana WHERE annotation.snpId = snpFisherurbana.snpId AND annotation.loc = feature_alias.loc AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 ;"


query2 <- dbGetQuery(con, sql2)
query2 <- query2 %>% unique()


wnt<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/wntSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
tor<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/torSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
insulin<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/insulinSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)



wnttest<-query2 %>% filter(., Flybase_FBgn %in% wnt$flyid) #%>% distinct(.,loc,Flybase_gene_symbol,Flybase_FBgn,.keep_all = TRUE)

wnttest %>% select(.,Flybase_FBgn) %>% unique()
  
query2 %>% select(.,Flybase_FBgn) %>% unique()#3813
#to many to run through david's clustering algorythim but for GO terms BP:
#axon giudance, transcription, imaginal disc derived wing morphogenesis, dendrite morphogenesis, tracheal system development
#cheimcal synaptic transmission, central cord development, regulation of transcription,
#picked out of significants: motor neuron axon quidance, regulation of glucose metabolic process, olfactory learning, neuron/muscle devlopment

#what if we look for 'exons' including splicey exons here

query.splice_region_variant_stop_retained_variant <-query2 %>% filter(., effect=="splice_region_variant&stop_retained_variant")
query.start_lost_splice_region_variant <-query2 %>% filter(., effect=="start_lost&splice_region_variant")
query.stop_gained_splice_region_variant <-query2 %>% filter(., effect=="stop_gained&splice_region_variant")
query.stop_lost_splice_region_variant<-query2 %>% filter(., effect=="stop_lost&splice_region_variant")
query.initiator_codon_variant_splice_region_variant<-query2 %>% filter(., effect=="initiator_codon_variant&splice_region_variant")
query.missense_variant_splice_region_variant <-query2 %>% filter(., effect=="missense_variant&splice_region_variant")
query.splice_region_variant_synonymous_variant <-query2 %>% filter(., effect=="splice_region_variant&synonymous_variant")
query.missense_variant <-query2 %>% filter(., effect=="missense_variant")
query.synonymous_variant <-query2 %>% filter(., effect=="synonymous_variant")
query.stop_gained <-query2 %>% filter(., effect=="stop_gained")
query.stop_retained_variant <-query2 %>% filter(., effect=="stop_retained_variant")
query.rare_amino_acid_variant <-query2 %>% filter(., effect=="rare_amino_acid_variant")
query.initiator_codon_variant <-query2 %>% filter(., effect=="initiator_codon_variant")
query.initiator_codon_variant_non_canonical_start_codon <-query2 %>% filter(., effect=="initiator_codon_variant&non_canonical_start_codon")
query.stop_lost <-query2 %>% filter(., effect=="stop_lost")
query.start_lost <-query2 %>% filter(., effect=="start_lost")
query.non_coding_transcript_exon_variant <-query2 %>% filter(., effect=="non_coding_transcript_exon_variant")

test<-rbind(query.splice_region_variant_stop_retained_variant,query.start_lost_splice_region_variant,query.stop_gained_splice_region_variant,query.stop_lost_splice_region_variant,query.initiator_codon_variant_splice_region_variant,query.missense_variant_splice_region_variant,query.splice_region_variant_synonymous_variant,query.missense_variant,query.synonymous_variant,query.stop_gained,query.stop_retained_variant,query.rare_amino_acid_variant,query.rare_amino_acid_variant,query.initiator_codon_variant,query.initiator_codon_variant_non_canonical_start_codon,query.stop_lost,query.start_lost,query.non_coding_transcript_exon_variant )


test %>% select(.,Flybase_FBgn) %>% unique()
#nothing sig in clusters

#what if we add in splicing 
test <-query2 %>% filter(., effect!="start_lost" & effect!="intron_variant" & effect!="intragenic_variant" & effect!="upstream_gene_variant" & effect!="downstream_gene_variant" & effect!="intergenic_region")
#960 597 flybase unique
#t
test %>% select(.,Flybase_FBgn) %>% unique()

query2 %>% filter(., effect=="missense_variant") 
#162 sites
query2 %>% filter(., effect=="missense_variant") %>%  select(.,Flybase_FBgn) %>% unique() %>% nrow()
#108 genes
query2 %>% filter(., effect=="missense_variant") %>% group_by(loc) %>% dplyr::summarise( length_loc=length(loc)) %>% arrange(length_loc)  %>% print(n = Inf) 


#with 3 sig missense snps
query2 %>% filter(., loc=="LOC108367089")#eyg
query2 %>% filter(., loc=="LOC108362929")#NA hit
query2 %>% filter(., loc=="LOC108361127")#CG30428
query2 %>% filter(., loc=="LOC108358604")#NA hit
#2 sig snps missense
query2 %>% filter(., loc=="LOC108382942")#Uro
query2 %>% filter(., loc=="LOC108372967")#dsx-c73A
query2 %>% filter(., loc=="LOC108368022")#NA
query2 %>% filter(., loc=="LOC108367132")#msd1
query2 %>% filter(., loc=="LOC108365438")#CG13894#flybase went down havent looked at these ones
query2 %>% filter(., loc=="LOC108364980")#Mat89Ba
query2 %>% filter(., loc=="LOC108363436")#nopo
query2 %>% filter(., loc=="LOC108363315")#AdenoK
query2 %>% filter(., loc=="LOC108362933")#Or69a #odarant receptor
query2 %>% filter(., loc=="LOC108362931")#NA
query2 %>% filter(., loc=="LOC108361312")#Nlg1

query2 %>% filter(., effect=="missense_variant") %>% select(.,scaffold, position,ref,alt,loc,Flybase_FBgn,Flybase_gene_symbol)

#scaffold	position	ref	alt	loc	Flybase_FBgn	Flybase_gene_symbol
#NW_016159616.1	14323	A	C	LOC108367089	FBgn0000625	eyg
#NW_016159616.1	14337	C	A	LOC108367089	FBgn0000625	eyg
#NW_016159616.1	14338	A	T	LOC108367089	FBgn0000625	eyg
#NW_016160900.1	10513	C	A	LOC108369718	FBgn0001253	ImpE1
#NW_016169357.1	17312	A	G	LOC108377548	FBgn0002354	l(3)87Df
#NW_016157114.1	401845	T	C	LOC108364061	FBgn0002431	hyd
#NW_016158600.1	14291	G	A	LOC108364274	FBgn0002441	l(3)mbt
#NW_016175384.1	9645	G	C	LOC108380127	FBgn0002931	net
#NW_016186674.1	3258	C	G	LOC108382942	FBgn0003961	Uro
#NW_016186674.1	3259	A	C	LOC108382942	FBgn0003961	Uro
#NW_016157205.1	369200	A	C	LOC108378908	FBgn0004198	ct
#NW_016157390.1	90320	A	G	LOC108358862	FBgn0004379	Klp67A
#NW_016162175.1	27845	C	T	LOC108371588	FBgn0004396	CrebA
#NW_016166953.1	12199	G	A	LOC108376092	FBgn0005659	Ets98B
#NW_016158536.1	101885	T	G	LOC108364058	FBgn0011576	Cyp4d2
#NW_016157132.1	361572	A	T	LOC108367859	FBgn0013813	Dhc98D
#NW_016185194.1	97	C	T	LOC108382676	FBgn0014396	tim
#NW_016157364.1	95438	C	A	LOC108358647	FBgn0014903	CG14630
#NW_016162834.1	24114	A	T	LOC108372437	FBgn0016756	Usp47
#NW_016161892.1	19766	G	C	LOC108371232	FBgn0025697	santa-maria
#NW_016161910.1	44906	C	T	LOC108371250	FBgn0026314	Ugt35b
#NW_016157292.1	100999	C	T	LOC108358112	FBgn0027083	MetRS-m
#NW_016161350.1	1201	C	T	LOC108370487	FBgn0027585	CG8740
#NW_016159962.1	4644	A	G	LOC108367945	FBgn0027655	htt
#NW_016162104.1	20220	A	C	LOC108371490	FBgn0028476	Usp1
#NW_016165219.1	10275	G	A	LOC108374795	FBgn0029676	HIP-R
#NW_016161474.1	8123	A	G	LOC108370678	FBgn0029834	CG5937
#NW_016164066.1	27498	T	C	LOC108373722	FBgn0030114	CG17754
#NW_016160542.1	36830	C	G	LOC108369004	FBgn0030569	CG9411
#NW_016158069.1	86008	T	C	LOC108362310	FBgn0030674	CG8184
#NW_016157437.1	203823	A	G	LOC108359194	FBgn0031264	CG11835
#NW_016177062.1	2131	C	A	LOC108380668	FBgn0031874	CG13775
#NW_016157479.1	15679	T	A	LOC108359456	FBgn0032049	Bace
#NW_016158582.1	20462	T	A	LOC108364223	FBgn0032050	CG13096
#NW_016159772.1	47550	T	C	LOC108367443	FBgn0032538	CG16885
#NW_016163263.1	6860	A	G	LOC108372936	FBgn0032683	kon
#NW_016163483.1	30749	C	T	LOC108373131	FBgn0033072	CG17994
#NW_016157647.1	121260	A	C	LOC108360390	FBgn0033302	Cyp6a14
#NW_016157465.1	242261	T	C	LOC108359393	FBgn0033397	Cyp4p3
#NW_016157148.1	557808	T	A	LOC108371068	FBgn0033479	CG2292
#NW_016159670.1	45563	T	A	LOC108367222	FBgn0033624	CG12384
#NW_016159857.1	4499	T	C	LOC108367688	FBgn0033661	CG13185
#NW_016157519.1	211403	G	A	LOC108359689	FBgn0033668	exp
#NW_016158269.1	105420	G	C	LOC108363090	FBgn0033926	Arc1
#NW_016157647.1	93669	C	G	LOC108360382	FBgn0033978	Cyp6a23
#NW_016159265.1	40834	G	T	LOC108366182	FBgn0034284	CG14491
#NW_016158361.1	108230	C	G	LOC108363436	FBgn0034314	nopo
#NW_016158361.1	108253	T	G	LOC108363436	FBgn0034314	nopo
#NW_016157715.1	209587	C	G	LOC108360805	FBgn0034400	CG15099
#NW_016157730.1	95192	T	G	LOC108360863	FBgn0035077	CG9083
#NW_016169780.1	8263	C	T	LOC108377790	FBgn0035149	MED30
#NW_016158954.1	77122	C	A	LOC108365438	FBgn0035157	CG13894
#NW_016158954.1	77208	C	T	LOC108365438	FBgn0035157	CG13894
#NW_016159637.1	71700	G	A	LOC108367132	FBgn0035209	msd1
#NW_016159637.1	71823	T	C	LOC108367132	FBgn0035209	msd1
#NW_016164691.1	13951	G	A	LOC108374334	FBgn0035477	CG14982
#NW_016157598.1	126838	A	T	LOC108360066	FBgn0035807	CG7492
#NW_016158661.1	12334	G	C	LOC108364511	FBgn0036202	CG6024
#NW_016158326.1	86951	C	A	LOC108363315	FBgn0036337	AdenoK
#NW_016158326.1	86954	A	G	LOC108363315	FBgn0036337	AdenoK
#NW_016164204.1	15880	G	C	LOC108373892	FBgn0036398	upSET
#NW_016166040.1	14976	G	T	LOC108375452	FBgn0036494	Toll-6
#NW_016157309.1	357469	C	A	LOC108358248	FBgn0036849	CG14079
#NW_016163102.1	27304	C	T	LOC108372776	FBgn0037020	Pex14
#NW_016177188.1	6080	C	G	LOC108380706	FBgn0037720	CG8312
#NW_016157133.1	276227	G	A	LOC108368104	FBgn0037853	CG14696
#NW_016192047.1	1032	T	A	LOC108354317	FBgn0037877	CG6689
#NW_016157258.1	220837	T	A	LOC108356870	FBgn0037949	CG17360
#NW_016158379.1	101369	G	A	LOC108363522	FBgn0037993	dpr15
#NW_016164202.1	4453	T	G	LOC108373889	FBgn0038038	CG5167
#NW_016157964.1	105484	A	G	LOC108361910	FBgn0038542	TyrR
#NW_016159780.1	56950	C	T	LOC108367470	FBgn0038638	CG7702
#NW_016165373.1	18106	A	T	LOC108374939	FBgn0038676	CG6026
#NW_016165944.1	15257	A	G	LOC108375396	FBgn0038720	CG6231
#NW_016167943.1	5311	G	T	LOC108376722	FBgn0038846	CG5697
#NW_016181102.1	5913	C	G	LOC108381844	FBgn0038897	CG5849
#NW_016165596.1	5587	G	A	LOC108375106	FBgn0039038	CG6688
#NW_016159670.1	38132	T	C	LOC108367219	FBgn0039195	CG17782
#NW_016160714.1	29978	C	G	LOC108369365	FBgn0039223	CG5805
#NW_016170388.1	8970	C	T	LOC108378088	FBgn0039451	CG6420
#NW_016161342.1	41486	A	G	LOC108370469	FBgn0039461	CG5500
#NW_016163868.1	4825	A	G	LOC108373533	FBgn0039564	CG5527
#NW_016158836.1	20232	A	C	LOC108365032	FBgn0039602	CG1647
#NW_016193787.1	834	T	A	LOC108354520	FBgn0039714	Zip99C
#NW_016157357.1	113199	T	G	LOC108358597	FBgn0039883	RhoGAP100F
#NW_016157243.1	62683	C	T	LOC108354743	FBgn0040260	Ugt36Bc
#NW_016157120.1	395756	C	A	LOC108365214	FBgn0041245	Gr39b
#NW_016158232.1	61397	A	T	LOC108362933	FBgn0041622	Or69a
#NW_016158232.1	61641	C	A	LOC108362933	FBgn0041622	Or69a
#NW_016173076.1	3393	A	G	LOC108379271	FBgn0042111	CG18766
#NW_016159559.1	56923	A	G	LOC108366977	FBgn0045035	tefu
#NW_016158120.1	19771	T	G	LOC108362569	FBgn0050089	CG30089
#NW_016158183.1	19973	A	C	LOC108362783	FBgn0050361	mtt
#NW_016159612.1	64137	A	T	LOC108367075	FBgn0050375	CG30375
#NW_016159318.1	54843	T	C	LOC108366336	FBgn0050428	CG30428
#NW_016157784.1	24452	T	A	LOC108361127	FBgn0050428	CG30428
#NW_016157784.1	24474	T	A	LOC108361127	FBgn0050428	CG30428
#NW_016157784.1	24487	C	A	LOC108361127	FBgn0050428	CG30428
#NW_016157816.1	65846	T	G	LOC108361312	FBgn0051146	Nlg1
#NW_016157816.1	65855	T	C	LOC108361312	FBgn0051146	Nlg1
#NW_016163117.1	28799	C	A	LOC108372796	FBgn0052262	CG32262
#NW_016160072.1	6480	C	T	LOC108368171	FBgn0052549	CG32549
#NW_016159475.1	32888	T	C	LOC108366739	FBgn0086350	tef
#NW_016157641.1	192656	A	T	LOC108360334	FBgn0259174	Nedd4
#NW_016171372.1	9748	A	T	LOC108378525	FBgn0259193	Ir94d
#NW_016157217.1	36355	G	C	LOC108381194	FBgn0259193	Ir94d
#NW_016160620.1	15159	T	C	LOC108369195	FBgn0260399	gwl
#NW_016157437.1	222770	G	A	LOC108359198	FBgn0260933	rempA
#NW_016158818.1	17252	C	T	LOC108364980	FBgn0261286	Mat89Ba
#NW_016158818.1	17260	G	A	LOC108364980	FBgn0261286	Mat89Ba
#NW_016163298.1	9479	A	G	LOC108372967	FBgn0261799	dsx-c73A
#NW_016163298.1	9564	A	G	LOC108372967	FBgn0261799	dsx-c73A
#NW_016160148.1	12072	C	G	LOC108368308	FBgn0261823	Asx
#NW_016158645.1	24321	T	G	LOC108364450	FBgn0263241	Mocs1
#NW_016157096.1	391531	T	C	LOC108359014	FBgn0263490	mld
#NW_016159631.1	62549	C	T	LOC108367121	FBgn0264598	PsGEF
#NW_016159786.1	47405	T	G	LOC108367492	FBgn0266801	CG45263
#NW_016158155.1	66037	A	C	LOC108362690	FBgn0266916	Cenp-C
#NW_016188384.1	2448	G	A	LOC108383188	FBgn0267347	squ
#NW_016157501.1	31230	C	A	LOC108359597	FBgn0267428	CG45781
#NW_016160676.1	24336	A	G	LOC108369296	FBgn0267433	kl-5
#NW_016157417.1	94222	A	T	LOC108359061	FBgn0283521	lola
#NW_016166855.1	3307	C	A	LOC108376012	<NA>	<NA>
  #NW_016157134.1	226544	C	G	LOC108368207	<NA>	<NA>
  #NW_016161016.1	3148	C	G	LOC108369965	<NA>	<NA>
  #NW_016163577.1	20978	G	T	LOC108373235	<NA>	<NA>
  #NW_016159508.1	81384	A	T	LOC108366854	<NA>	<NA>
  #NW_016174058.1	9320	C	A	LOC108379672	<NA>	<NA>
  #NW_016190682.1	3320	G	C	LOC108354134	<NA>	<NA>
  #NW_016162394.1	20179	C	A	LOC108371930	<NA>	<NA>
  #NW_016164362.1	26128	T	G	LOC108374055	<NA>	<NA>
  #NW_016161867.1	25787	C	T	LOC108371210	<NA>	<NA>
  #NW_016157087.1	698985	G	A	LOC108368022	<NA>	<NA>
  #NW_016157087.1	700045	T	G	LOC108368022	<NA>	<NA>
  #NW_016198463.1	1856	A	C	LOC108355010	<NA>	<NA>
  #NW_016160333.1	53642	C	A	LOC108368676	<NA>	<NA>
  #NW_016165557.1	587	C	A	LOC108375071	<NA>	<NA>
  #NW_016165557.1	18129	A	C	LOC108375074	<NA>	<NA>
  #NW_016162420.1	38575	T	C	LOC108371963	<NA>	<NA>
  #NW_016157424.1	137056	T	G	LOC108359116	<NA>	<NA>
  #NW_016157215.1	143647	C	A	LOC108381053	<NA>	<NA>
  #NW_016239793.1	75	G	A	LOC108357704	<NA>	<NA>
  #NW_016158092.1	128334	C	A	LOC108362426	<NA>	<NA>
  #NW_016157096.1	240457	G	T	LOC108358971	<NA>	<NA>
  #NW_016157423.1	31616	T	C	LOC108359103	<NA>	<NA>
  #NW_016157632.1	173434	A	C	LOC108360267	<NA>	<NA>
  #NW_016162498.1	5136	G	C	LOC108372068	<NA>	<NA>
  #NW_016157658.1	26934	G	A	LOC108360452	<NA>	<NA>
  #NW_016157400.1	287332	G	T	LOC108358937	<NA>	<NA>
  #NW_016198888.1	410	T	C	LOC108355059	<NA>	<NA>
  #NW_016158399.1	2993	C	T	LOC108363609	<NA>	<NA>
  #NW_016174778.1	9466	G	A	LOC108379926	<NA>	<NA>
  #NW_016157358.1	293077	C	T	LOC108358604	<NA>	<NA>
  #NW_016157358.1	293131	A	G	LOC108358604	<NA>	<NA>
  #NW_016157358.1	293134	A	T	LOC108358604	<NA>	<NA>
  #NW_016159037.1	18538	C	T	LOC108365669	<NA>	<NA>
  #NW_016158757.1	78983	A	G	LOC108364793	<NA>	<NA>
  #NW_016158231.1	134817	G	A	LOC108362931	<NA>	<NA>
  #NW_016158231.1	134819	A	C	LOC108362931	<NA>	<NA>
  #NW_016158231.1	134197	G	A	LOC108362929	<NA>	<NA>
  #NW_016158231.1	134266	A	G	LOC108362929	<NA>	<NA>
  #NW_016158231.1	134315	C	A	LOC108362929	<NA>	<NA>


#ways to focus this down further

##################################################################################
#looking for things with a biologically different maf value between apple and haw#
##################################################################################
#say 0.5
nrow(query2) #7072
#dropped us in half
#splicing and exons & prime regions


query2 %>% select(.,Flybase_FBgn) %>% unique() %>% nrow() #2498
#david
#sig membrane/tranmembrane, immunoglobulin domain, DNA-binding, Homeobox, transcription, EGF, Ion transport, Flycoprotein, Metal-binding (Zinc finger), SH3 
#Lipoprotein, Kinase, Pleckstrin, Cell Junction/synapses, vision/sensory transduction, G-protein recptor/signalling

#homeobox embryonic development?

test <-query2 %>% filter(., effect!="start_lost" & effect!="intron_variant" & effect!="intragenic_variant" & effect!="upstream_gene_variant" & effect!="downstream_gene_variant" & effect!="intergenic_region")
nrow(test) #338

test%>% select(.,Flybase_FBgn) %>% unique()
#252 genes
#david nothing sig

query2 %>% filter(., effect=="missense_variant") 
#40 sites
query2 %>% filter(., effect=="missense_variant") %>%  select(.,Flybase_FBgn) %>% unique() %>% nrow()
#33
query2 %>% filter(., effect=="missense_variant") %>% group_by(loc) %>% dplyr::summarise( length_loc=length(loc)) %>% arrange(length_loc)  %>% print(n = Inf) 

query2 %>% filter(., loc=="LOC108382942") #2 Uro

query2 %>% filter(., loc=="LOC108363436") #2 in nopo #embryonic/cell death

query2 %>% filter(., effect=="missense_variant") %>% select(.,scaffold, position,ref,alt,loc,Flybase_FBgn,Flybase_gene_symbol,urbana_appleave_Maj,urbana_hawave_Maj)

#gene list
#scaffold	position	ref	alt	loc	Flybase_FBgn	Flybase_gene_symbol	urbana_appleave_Maj	urbana_hawave_Maj
#NW_016175384.1	9645	G	C	LOC108380127	FBgn0002931	net	0.849183	0.332999
#NW_016186674.1	3258	C	G	LOC108382942	FBgn0003961	Uro	0.434652	1
#NW_016186674.1	3259	A	C	LOC108382942	FBgn0003961	Uro	0.4165	1
#NW_016157390.1	90320	A	G	LOC108358862	FBgn0004379	Klp67A	0.295045	0.825651
#NW_016157132.1	361572	A	T	LOC108367859	FBgn0013813	Dhc98D	0.857859	0.34059
#NW_016185194.1	97	C	T	LOC108382676	FBgn0014396	tim	1	0.23023
#NW_016161892.1	19766	G	C	LOC108371232	FBgn0025697	santa-maria	0.368157	0.875752
#NW_016161910.1	44906	C	T	LOC108371250	FBgn0026314	Ugt35b	0.806168	0.249499
#NW_016162104.1	20220	A	C	LOC108371490	FBgn0028476	Usp1	0.76243	0.221666
#NW_016158069.1	86008	T	C	LOC108362310	FBgn0030674	CG8184	0.796511	0.289052
#NW_016157479.1	15679	T	A	LOC108359456	FBgn0032049	Bace	0.849183	0.090089
#NW_016158582.1	20462	T	A	LOC108364223	FBgn0032050	CG13096	0.422923	0.968679
#NW_016157465.1	242261	T	C	LOC108359393	FBgn0033397	Cyp4p3	0.818819	0.307307
#NW_016159670.1	45563	T	A	LOC108367222	FBgn0033624	CG12384	1	0.380714
#NW_016159265.1	40834	G	T	LOC108366182	FBgn0034284	CG14491	0.191691	0.767978
#NW_016158361.1	108230	C	G	LOC108363436	FBgn0034314	nopo	0.78105	0.215648
#NW_016158361.1	108253	T	G	LOC108363436	FBgn0034314	nopo	0.707733	0.157209
#NW_016157730.1	95192	T	G	LOC108360863	FBgn0035077	CG9083	0.934202	0.388666
#NW_016158326.1	86951	C	A	LOC108363315	FBgn0036337	AdenoK	0.323176	0.824178
#NW_016157133.1	276227	G	A	LOC108368104	FBgn0037853	CG14696	0.907064	0.382117
#NW_016158379.1	101369	G	A	LOC108363522	FBgn0037993	dpr15	0.853648	0.293705
#NW_016165944.1	15257	A	G	LOC108375396	FBgn0038720	CG6231	0.917502	0.291249
#NW_016181102.1	5913	C	G	LOC108381844	FBgn0038897	CG5849	0.181181	0.889668
#NW_016163868.1	4825	A	G	LOC108373533	FBgn0039564	CG5527	0.834001	0.191691
#NW_016193787.1	834	T	A	LOC108354520	FBgn0039714	Zip99C	0.41362	0.957437
#NW_016159559.1	56923	A	G	LOC108366977	FBgn0045035	tefu	0.879547	0.313914
#NW_016158183.1	19973	A	C	LOC108362783	FBgn0050361	mtt	0.215648	0.818819
#NW_016159612.1	64137	A	T	LOC108367075	FBgn0050375	CG30375	0.199399	1
#NW_016159318.1	54843	T	C	LOC108366336	FBgn0050428	CG30428	0.938377	0.347521
#NW_016163117.1	28799	C	A	LOC108372796	FBgn0052262	CG32262	0.448172	0.969689
#NW_016157641.1	192656	A	T	LOC108360334	FBgn0259174	Nedd4	1	0.454454
#NW_016171372.1	9748	A	T	LOC108378525	FBgn0259193	Ir94d	0.386871	0.971531
#NW_016158818.1	17252	C	T	LOC108364980	FBgn0261286	Mat89Ba	0.342542	0.864365
#NW_016188384.1	2448	G	A	LOC108383188	FBgn0267347	squ	0.434652	0.945335
#NW_016166855.1	3307	C	A	LOC108376012	<NA>	<NA>	0.3998	1
#NW_016157134.1	226544	C	G	LOC108368207	<NA>	<NA>	1	0.482724
#NW_016159508.1	81384	A	T	LOC108366854	<NA>	<NA>	0.959252	0.285285
#NW_016198463.1	1856	A	C	LOC108355010	<NA>	<NA>	0.234764	0.915116
#NW_016162498.1	5136	G	C	LOC108372068	<NA>	<NA>	0.850701	0.305166
#NW_016158757.1	78983	A	G	LOC108364793	<NA>	<NA>	0.454454	0.955456

################################################
#things that are sig and in the right direction#
################################################

nrow(query2) #
colnames(query2)

#out.snp<-final_table[!duplicated(final_table[,'loc']),]
test1<-query2 %>% dplyr::select(.,matches('appleave_Maj|hawave_Maj'),scaffold,position,ref,alt,loc,effect,gene_id,Flybase_tophit,Flybase_FBgn,Flybase_gene_symbol,urbana_appleave_hawave_fisher_pvalue_adjust,urbana_appleave_Maj,urbana_hawave_Maj,urbana_appleearly_Maj,urbana_applelate_Maj,urbana_hawearly_Maj,urbana_hawlate_Maj) %>% mutate(appleave_hawave = .[,1]-.[,2]) %>% dplyr::select(.,appleave_hawave,loc,scaffold,position,ref,alt,loc,effect,gene_id,Flybase_tophit,Flybase_FBgn,Flybase_gene_symbol,urbana_appleave_hawave_fisher_pvalue_adjust,urbana_appleave_Maj,urbana_hawave_Maj,urbana_appleearly_Maj,urbana_applelate_Maj,urbana_hawearly_Maj,urbana_hawlate_Maj)
                                 
test2<-query2 %>% dplyr::select(.,matches('hawearly_Maj|hawlate_Maj')) %>% mutate(hawearly_hawlate = .[,1]-.[,2]) %>% dplyr::select(.,hawearly_hawlate)
test3<-query2 %>% dplyr::select(.,matches('appleearly_Maj|applelate_Maj')) %>% mutate(appleearly_applelate = .[,1]-.[,2]) %>% dplyr::select(.,appleearly_applelate)
test<-cbind(test1,test2,test3) 
test2<-test %>% mutate(signappleave_hawave=ifelse(appleave_hawave > 0,'Positive','Negative')) %>% mutate(signappleearly_applelate=ifelse(appleearly_applelate > 0,'Positive','Negative'))  %>% mutate(signhawearly_hawlate =ifelse(hawearly_hawlate > 0,'Positive','Negative')) 
sumtab<-test2 %>% group_by(signappleave_hawave,signappleearly_applelate,signhawearly_hawlate)  %>% 

test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive') #1876
test2 %>% dplyr::filter(.,signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative') #2628

#david

test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' | signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative') %>%  select(.,Flybase_FBgn) %>% unique() %>% nrow() #1631
table1<-test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' | signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative') %>%  select(.,Flybase_FBgn) %>% unique() 
write.table(table1,'/home/eddydowle/Documents/pomonella/sigandrightdirectionhostraces',sep='\t',quote=F,row.names=F)
#Membrane, Immunoglobulin, Homeobox, EGF-like domain, transcription DNA-templated, Pleckstrin, Scr homology, G-protein coupled receptor, Glycoprotein, Dblhomology domain
#nucleotide binding, fibroblast growh factor receptor signaling pathway, metal binding/Zinc, Kinase, sensory transduction/Vision, golgi apparatus, 



#and abs
test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & abs(appleave_hawave) > 0.5) #1009
test2 %>% dplyr::filter(.,signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5) #1274

test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & abs(appleave_hawave) > 0.5 | signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5) %>%  select(.,Flybase_FBgn) %>% unique() %>% nrow() #968
table1<-test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' | signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative') %>%  select(.,Flybase_FBgn) %>% unique() 
write.table(table1,'/home/eddydowle/Documents/pomonella/sigandrightdirectionhostracesabsdiffover05',sep='\t',quote=F,row.names=F)

#membrane, immunoglobulin, homeobox, EGF-like domain, transcription, Pleckstrin homology like domain, Src homoloy-3 domain, G-protein coupled receptor, Glycoprotein, Dbl homology domain, Nucleotide-binding, fibroblast growth factor receptor signalling pathway, Metal-binding (Zinc), Kinase, sensory transduction/vision, fibroblast growth factor receptor signaling pathwaygolgi apparatus, 
#not quite sig R3/R4 cell fate commitment, 
#specific genes

test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & effect=="missense_variant"| signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' &effect=="missense_variant") %>% group_by(loc) %>% dplyr::summarise( length_loc=length(loc)) %>% arrange(length_loc)  %>% print(n = Inf) 
#49 total
#2 missense mutations:
query2 %>% filter(., loc=="LOC108365438") #CG13894
query2 %>% filter(., loc=="LOC108363315") #AdenoK
query2 %>% filter(., loc=="LOC108362931") #NA nothing in gff really
query2 %>% filter(., loc=="LOC108362929") #NA nothing in gff really
query2 %>% filter(., loc=="LOC108361312") #Nlg1

test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & effect=="missense_variant" & abs(appleave_hawave) > 0.5| signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' &effect=="missense_variant" & abs(appleave_hawave) > 0.5) %>% group_by(loc) %>% dplyr::summarise( length_loc=length(loc)) %>% arrange(length_loc)  %>% print(n = Inf) #13 snps no genes >1 snp

query2 %>% filter(., loc=="LOC108354520") #Zip99C zinc ion homeostatis iron ion transport transmembrane transporter
query2 %>% filter(., loc=="LOC108360863") #CG9083 unknown
query2 %>% filter(., loc=="LOC108362310") #CG8184 dsRNA transport  
query2 %>% filter(., loc=="LOC108362783") # mtt adult feeding behaviour
query2 %>% filter(., loc=="LOC108363315") #AdenoK adenosine salvage
query2 %>% filter(., loc=="LOC108366182") #CG14491 unkown
query2 %>% filter(., loc=="LOC108366336") #CG30428 apoptotic
query2 %>% filter(., loc=="LOC108366977") #tefu cellular processes
query2 %>% filter(., loc=="LOC108367859") #Dhc98D
query2 %>% filter(., loc=="LOC108368104") #CG14696 unknown
query2 %>% filter(., loc=="LOC108372068") #NA
query2 %>% filter(., loc=="LOC108373533") # CG5527/Nep7 proteolysis
query2 %>% filter(., loc=="LOC108378525") #Ir94d detection chemical stimulus



#filter(., effect=="missense_variant")
test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & effect=="missense_variant"| signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' &effect=="missense_variant") %>% select(scaffold, position,ref,alt,loc,Flybase_FBgn,Flybase_gene_symbol,urbana_appleave_Maj,urbana_hawave_Maj) %>%  unique()




#
#scaffold	position	ref	alt	loc	Flybase_FBgn	Flybase_gene_symbol	urbana_appleave_Maj	urbana_hawave_Maj
#1	NW_016158600.1	14291	G	A	LOC108364274	FBgn0002441	l(3)mbt	0.919758	0.424091
#2	NW_016157205.1	369200	A	C	LOC108378908	FBgn0004198	ct	0.965216	0.563764
#3	NW_016158536.1	101885	T	G	LOC108364058	FBgn0011576	Cyp4d2	0.78781	0.316707
#4	NW_016157132.1	361572	A	T	LOC108367859	FBgn0013813	Dhc98D	0.857859	0.34059
#5	NW_016157364.1	95438	C	A	LOC108358647	FBgn0014903	CG14630	0.536659	1
#6	NW_016161350.1	1201	C	T	LOC108370487	FBgn0027585	CG8740	0.465046	0.943745
#7	NW_016160542.1	36830	C	G	LOC108369004	FBgn0030569	CG9411	0.444333	0.906473
#8	NW_016158069.1	86008	T	C	LOC108362310	FBgn0030674	CG8184	0.796511	0.289052
#9	NW_016159772.1	47550	T	C	LOC108367443	FBgn0032538	CG16885	0.442192	0.873087
#10	NW_016159265.1	40834	G	T	LOC108366182	FBgn0034284	CG14491	0.191691	0.767978
#11	NW_016157730.1	95192	T	G	LOC108360863	FBgn0035077	CG9083	0.934202	0.388666
#12	NW_016158954.1	77122	C	A	LOC108365438	FBgn0035157	CG13894	0.726644	1
#13	NW_016158954.1	77208	C	T	LOC108365438	FBgn0035157	CG13894	0.716852	0.974282
#14	NW_016164691.1	13951	G	A	LOC108374334	FBgn0035477	CG14982	1	0.667001
#15	NW_016158326.1	86951	C	A	LOC108363315	FBgn0036337	AdenoK	0.323176	0.824178
#16	NW_016158326.1	86954	A	G	LOC108363315	FBgn0036337	AdenoK	0.299599	0.793039
#17	NW_016166040.1	14976	G	T	LOC108375452	FBgn0036494	Toll-6	0.986267	0.714715
#18	NW_016157309.1	357469	C	A	LOC108358248	FBgn0036849	CG14079	0.963191	0.555667
#19	NW_016157133.1	276227	G	A	LOC108368104	FBgn0037853	CG14696	0.907064	0.382117
#20	NW_016157258.1	220837	T	A	LOC108356870	FBgn0037949	CG17360	0.468021	0.871113
#21	NW_016159780.1	56950	C	T	LOC108367470	FBgn0038638	CG7702	0.561098	1
#22	NW_016165373.1	18106	A	T	LOC108374939	FBgn0038676	CG6026	0.619286	1
#23	NW_016159670.1	38132	T	C	LOC108367219	FBgn0039195	CG17782	0.54175	1
#24	NW_016160714.1	29978	C	G	LOC108369365	FBgn0039223	CG5805	0.718385	1
#25	NW_016170388.1	8970	C	T	LOC108378088	FBgn0039451	CG6420	0.625251	1
#26	NW_016163868.1	4825	A	G	LOC108373533	FBgn0039564	CG5527	0.834001	0.191691
#27	NW_016193787.1	834	T	A	LOC108354520	FBgn0039714	Zip99C	0.41362	0.957437
#28	NW_016157120.1	395756	C	A	LOC108365214	FBgn0041245	Gr39b	0.949617	0.483839
#29	NW_016159559.1	56923	A	G	LOC108366977	FBgn0045035	tefu	0.879547	0.313914
#30	NW_016158120.1	19771	T	G	LOC108362569	FBgn0050089	CG30089	0.516161	1
#31	NW_016158183.1	19973	A	C	LOC108362783	FBgn0050361	mtt	0.215648	0.818819
#32	NW_016159318.1	54843	T	C	LOC108366336	FBgn0050428	CG30428	0.938377	0.347521
#33	NW_016157784.1	24452	T	A	LOC108361127	FBgn0050428	CG30428	0.425777	0.864365
#34	NW_016157816.1	65846	T	G	LOC108361312	FBgn0051146	Nlg1	0.535786	1
#35	NW_016157816.1	65855	T	C	LOC108361312	FBgn0051146	Nlg1	0.592778	1
#36	NW_016171372.1	9748	A	T	LOC108378525	FBgn0259193	Ir94d	0.386871	0.971531
#37	NW_016157217.1	36355	G	C	LOC108381194	FBgn0259193	Ir94d	0.677278	0.258777
#38	NW_016157437.1	222770	G	A	LOC108359198	FBgn0260933	rempA	1	0.57515
#39	NW_016158155.1	66037	A	C	LOC108362690	FBgn0266916	Cenp-C	0.658211	1
#40	NW_016160676.1	24336	A	G	LOC108369296	FBgn0267433	kl-5	0.483839	0.952124
#41	NW_016161016.1	3148	C	G	LOC108369965	<NA>	<NA>	0.646126	1
#42	NW_016163577.1	20978	G	T	LOC108373235	<NA>	<NA>	0.791628	1
#43	NW_016157087.1	700045	T	G	LOC108368022	<NA>	<NA>	0.919391	1
#44	NW_016162420.1	38575	T	C	LOC108371963	<NA>	<NA>	0.710098	1
#45	NW_016157096.1	240457	G	T	LOC108358971	<NA>	<NA>	0.874765	1
#46	NW_016162498.1	5136	G	C	LOC108372068	<NA>	<NA>	0.850701	0.305166
#47	NW_016157658.1	26934	G	A	LOC108360452	<NA>	<NA>	0.743344	0.912283
#48	NW_016157400.1	287332	G	T	LOC108358937	<NA>	<NA>	0.924129	0.816053
#49	NW_016198888.1	410	T	C	LOC108355059	<NA>	<NA>	0.369304	0.861835
#50	NW_016174778.1	9466	G	A	LOC108379926	<NA>	<NA>	0.914214	1
#51	NW_016158231.1	134817	G	A	LOC108362931	<NA>	<NA>	0.781385	1
#52	NW_016158231.1	134819	A	C	LOC108362931	<NA>	<NA>	0.775198	1
#53	NW_016158231.1	134197	G	A	LOC108362929	<NA>	<NA>	0.779629	1
#54	NW_016158231.1	134266	A	G	LOC108362929	<NA>	<NA>	0.773881	1

#????????????????????????????????????????????????????????????????#
#what is the shit that is significant but in the wrong direction?#
#????????????????????????????????????????????????????????????????#

test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' | signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive') %>%  select(.,Flybase_FBgn) %>% unique() %>% nrow() #1617
table1<-test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' | signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive') %>%  select(.,Flybase_FBgn) %>% unique() 
write.table(table1,'/home/eddydowle/Documents/pomonella/sigandwrongdirectionhostraces',sep='\t',quote=F,row.names=F)
#david
#membrane, immunoglobulin like fold, transcription,ion transport, homeobox, glycoprotin, EGF like domain Zinc, cell junction, vision sensory transduction, lipid biosynthesis, 

#and abs
test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5 | signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & abs(appleave_hawave) > 0.5) %>%  select(.,Flybase_FBgn) %>% unique() %>% nrow() #983
table1<-test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5 | signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive'& abs(appleave_hawave) > 0.5) %>%  select(.,Flybase_FBgn) %>% unique() 
write.table(table1,'/home/eddydowle/Documents/pomonella/sigandwrongdirectionhostracesabsdiffover05',sep='\t',quote=F,row.names=F)
#transmembrane helix. immuoglobulin like fold, transcription, LDLa, immunoglobulin/Fibronectin, Metal-binding/zinc, glycoprotein, 
#almost in there: nuclear hormone receoptor~with a hint at ecdyosene


test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & effect=="missense_variant"| signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & effect=="missense_variant") %>% group_by(loc) %>% dplyr::summarise( length_loc=length(loc)) %>% arrange(length_loc)  %>% print(n = Inf) 
#49 genes

#3 missense mutations:
query2 %>% filter(., loc=="LOC108367089") #eyg Eyegone (Eyg) is a Pax family transcription factor that acts as a transcriptional repressor. In eye development, Eyg promotes cell proliferation in the larval eye disc through activation of Jak/STAT pathway. Eyg also plays both positive and negative roles in head vortex development. 

#2 missense mutations:
query2 %>% filter(., loc=="LOC108372967") #dsx-c73A cutitcle development
query2 %>% filter(., loc=="LOC108367132") #msd1 microtubule nucleartion
query2 %>% filter(., loc=="LOC108362933") #Or69a detection of chemical stimulus involved in sensory perception of smell !weird!

test2 %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & effect=="missense_variant" & abs(appleave_hawave) > 0.5| signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & effect=="missense_variant" & abs(appleave_hawave) > 0.5) %>% group_by(loc) %>% dplyr::summarise( length_loc=length(loc)) %>% arrange(length_loc)  %>% print(n = Inf) 
#13 genes


query2 %>% filter(., loc=="LOC108355010") #NA
query2 %>% filter(., loc=="LOC108358862") #Klp67A
query2 %>% filter(., loc=="LOC108360334") #Nedd4 Nedd4 is an E3 ubiquitin ligase that negatively regulates the Notch signaling pathway and the genes comm, Amph. Nedd4 contributes to neuromuscular synaptogenesis, transverse tubule formation in muscles, and in muscle function
query2 %>% filter(., loc=="LOC108363436") #nopo No poles is a RING domain-containing E3 ubiquitin ligase that is essential for early embryogenesis. It positively regulates caspase-dependent cell death. [Date last reviewed: 2016-12-01] 
query2 %>% filter(., loc=="LOC108363522") #dpr15 sensory perception of chemical stimulus
query2 %>% filter(., loc=="LOC108364223") #CG13096 translation
query2 %>% filter(., loc=="LOC108364793") #NA
query2 %>% filter(., loc=="LOC108367075") #CG30375 proteolysis
query2 %>% filter(., loc=="LOC108368207") #NA (has frame shift mutation too)
query2 %>% filter(., loc=="LOC108376012") #NA
query2 %>% filter(., loc=="LOC108380127") #net Net is a basic helix-loop-helix protein that probably acts as a transcriptional repressor. During wing vein formation it is expressed in all interveins territories and antagonises EGFR activity.
query2 %>% filter(., loc=="LOC108382676") # tim Timeless (Tim) is a key component of the Tim-per complex, required for the production of circadian rhythms.
query2 %>% filter(., loc=="LOC108383188") #squ Squash acts in the piRNA pathway that responds to transposase activity in the germline.


###############################################
#things that also have a DE between host races#
###############################################

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


#using together (which is within)

sigwithinFDR005<-together %>% filter(.,FDR < 0.05 ) %>% select(.,gene_id) %>% unique() #2902 gene ID's

#just sig

table1<-sigwithinFDR005 %>% left_join(.,query2,by='gene_id') %>% select(.,Flybase_FBgn) %>% unique() #807
write.table(table1,'/home/eddydowle/Documents/pomonella/sigandDEgene',sep='\t',quote=F,row.names=F)

#david
#metal-binding, nucleotide binding, immunoglobulin like domain, transcription, EGF-like domain, GTP-binding, zinc finger


sigwithinFDR005 %>% left_join(.,query2,by='gene_id') %>% filter(., effect=="missense_variant") %>% select(.,scaffold, position,ref,alt,loc,Flybase_FBgn,Flybase_gene_symbol,urbana_appleave_Maj,urbana_hawave_Maj) %>% group_by(loc) %>% dplyr::summarise( length_loc=length(loc)) %>% arrange(length_loc)  %>% print(n = Inf) 
#24 loci
#2 snps:
query2 %>% filter(., loc=="LOC108364980") #Mat89Ba rRNA processing
query2 %>% filter(., loc=="LOC108365438") #CG13894 DNA binding

#and biologically meaningful


table1<-sigwithinFDR005 %>% left_join(.,test2,by='gene_id') %>% filter(.,  abs(appleave_hawave) > 0.5)  %>% select(.,Flybase_FBgn) %>% unique() #511
write.table(table1,'/home/eddydowle/Documents/pomonella/sigandDEgeneBioMeaningful05',sep='\t',quote=F,row.names=F)

#david
#immunoglobulin domain, Metal-binding/Zinc,Nucleotide binding, Membrane

sigwithinFDR005 %>% left_join(.,query2,by='gene_id') %>% filter(., effect=="missense_variant" &  abs(urbana_appleave_Maj-urbana_hawave_Maj) > 0.5) %>% select(.,scaffold, position,ref,alt,loc,Flybase_FBgn,Flybase_gene_symbol,urbana_appleave_Maj,urbana_hawave_Maj) %>% group_by(loc) %>% dplyr::summarise( length_loc=length(loc)) %>% arrange(length_loc)  %>% print(n = Inf) 

#all ones
#  1 LOC108354520          1
#2 LOC108358862          1
#3 LOC108364223          1
#4 LOC108364980          1
#5 LOC108366977          1

query2 %>% filter(., loc=="LOC108354520") #Zip99C  Zinc/iron regulated transporter-related protein 99C (Zip99C) is an iron transporter located on the Golgi and endoplasmic reticulum that moves iron from the cytoplasm into the secretory compartments.
query2 %>% filter(., loc=="LOC108358862") #Klp67A meosis/mitosis
query2 %>% filter(., loc=="LOC108364223") #CG13096 translation
query2 %>% filter(., loc=="LOC108364980") #Mat89Ba rRNA processing
query2 %>% filter(., loc=="LOC108366977") #tefu mitotic/hisotone not well described

#sig and correct direction

sigwithinFDR005 %>% left_join(.,test2,by='gene_id')  %>% dplyr::filter(.,signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative'  | signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' ) %>%  select(.,Flybase_FBgn) %>% unique() %>% nrow() #354
table1<-sigwithinFDR005 %>% left_join(.,test2,by='gene_id') %>% dplyr::filter(.,signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' | signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' ) %>%  select(.,Flybase_FBgn) %>% unique() 
write.table(table1,'/home/eddydowle/Documents/pomonella/sigandDEgeneAndCorrectDir',sep='\t',quote=F,row.names=F)

#david
#nothing sig

#specfic genes with missense mutation

sigwithinFDR005 %>% left_join(.,test2,by='gene_id')  %>% dplyr::filter(.,effect=="missense_variant" & signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative'  | effect=="missense_variant" & signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' ) %>% select(.,scaffold, position,ref,alt,loc,Flybase_FBgn,Flybase_gene_symbol,urbana_appleave_Maj,urbana_hawave_Maj) %>% group_by(loc) %>% dplyr::summarise( length_loc=length(loc)) %>% arrange(length_loc)  %>% print(n = Inf) 

#10 loci
#loc          length_loc
#  1 LOC108354520          1
#2 LOC108358248          1
#3 LOC108362569          1
#4 LOC108362690          1
#5 LOC108366977          1
#6 LOC108369365          1
#7 LOC108375452          1
#8 LOC108378088          1
#9 LOC108378908          1
#10 LOC108365438          2

query2 %>% filter(., loc=="LOC108354520") # Zip99C Zinc/iron regulated transporter-related protein 99C (Zip99C) is an iron transporter located on the Golgi and endoplasmic reticulum that moves iron from the cytoplasm into the secretory compartments.
query2 %>% filter(., loc=="LOC108358248") # CG14079 no data
query2 %>% filter(., loc=="LOC108362569") # CG30089 no data
query2 %>% filter(., loc=="LOC108362690") # Cenp-C Centromeric protein-C (Cenp-C) is an essential centromere protein. It binds to the centromere-specific cid nucleosomes and provides a binding site for the Mis12 kinetochore protein complex, which is recruited to the centromere at the start of mitotic and meiotic M phases. Cenp-C also binds cal1, a cid loading factor crucial for propagation of the epigenetic mark that specifies centromere identity during progression through the cell division cycle.
query2 %>% filter(., loc=="LOC108366977") # tefu mitotic/hisotone not well described
query2 %>% filter(., loc=="LOC108369365") # CG5805 mitochondrial transport
query2 %>% filter(., loc=="LOC108375452") # Toll-6 Toll-6 (Toll-6) encodes a member of the Toll-like receptor family. It has neurotrophin receptor activity and contributes to dendrite guidance. Genetic interaction of Toll-6 with 18w and Tollo suggests a role in convergent extension during early embryogenesis
query2 %>% filter(., loc=="LOC108378088") # CG6420 protein deubiquitination
query2 %>% filter(., loc=="LOC108378908") # ct The homeoprotein Cut functions as a transcriptional factor in many different cells such as wing disc, muscle, oocyte and sense organ cells. It is a regulator of type-specific neuronal identity in the peripheral nervous system. Cut is expressed at variable levels in the dendritic arborization (DA) neurons and these levels control the different dendritic morphologies specific for each class of DA neurons. 
query2 %>% filter(., loc=="LOC108365438") # CG13894 DNA binding


#and abs value >0.5

sigwithinFDR005 %>% left_join(.,test2,by='gene_id')  %>% dplyr::filter(.,signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5 | signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & abs(appleave_hawave) > 0.5) %>%  select(.,Flybase_FBgn) %>% unique() %>% nrow() #205
table1<-sigwithinFDR005 %>% left_join(.,test2,by='gene_id') %>% dplyr::filter(.,signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5 | signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive'& abs(appleave_hawave) > 0.5) %>%  select(.,Flybase_FBgn) %>% unique() 
write.table(table1,'/home/eddydowle/Documents/pomonella/sigandDEgeneAndCorrectDirabsdiffover05',sep='\t',quote=F,row.names=F)

#david
#nothing sig

#genes with missense mutation

sigwithinFDR005 %>% left_join(.,test2,by='gene_id')  %>% dplyr::filter(.,effect=="missense_variant" & signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5 | effect=="missense_variant" & signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & abs(appleave_hawave) > 0.5) %>% select(.,scaffold, position,ref,alt,loc,Flybase_FBgn,Flybase_gene_symbol,urbana_appleave_Maj,urbana_hawave_Maj) %>% group_by(loc) %>% dplyr::summarise( length_loc=length(loc)) %>% arrange(length_loc)  %>% print(n = Inf) 

#two loci
#loc          length_loc
#1 LOC108354520          1
#2 LOC108366977          1

query2 %>% filter(., loc=="LOC108354520") # Zip99C Zinc/iron regulated transporter-related protein 99C (Zip99C) is an iron transporter located on the Golgi and endoplasmic reticulum that moves iron from the cytoplasm into the secretory compartments.
query2 %>% filter(., loc=="LOC108366977") # tefu mitotic/hisotone not well described

#what about the other 205 which have stuff in 5prime etc
table1<-sigwithinFDR005 %>% left_join(.,test2,by='gene_id')  %>% dplyr::filter(.,signappleave_hawave=='Negative' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5 | signappleave_hawave=='Positive' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & abs(appleave_hawave) > 0.5) %>% unique() 

write.table(table1,'/home/eddydowle/Documents/pomonella/sigandDEgeneAndCorrectDirabsdiffover05fullinfo',sep='\t',quote=F,row.names=F)
#have table 



#what is going in opposite direction
#sig and wrong direction

sigwithinFDR005 %>% left_join(.,test2,by='gene_id')  %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative'  | signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' ) %>%  select(.,Flybase_FBgn) %>% unique() %>% nrow() #322
table1<-sigwithinFDR005 %>% left_join(.,test2,by='gene_id') %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' | signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' ) %>%  select(.,Flybase_FBgn) %>% unique() 
write.table(table1,'/home/eddydowle/Documents/pomonella/sigandDEgeneAndWrongDir',sep='\t',quote=F,row.names=F)

#sig and wrong direction and biological meaningful


sigwithinFDR005 %>% left_join(.,test2,by='gene_id')  %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5 | signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & abs(appleave_hawave) > 0.5) %>%  select(.,Flybase_FBgn) %>% unique() %>% nrow() #195
table1<-sigwithinFDR005 %>% left_join(.,test2,by='gene_id') %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5 | signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive'& abs(appleave_hawave) > 0.5) %>%  select(.,Flybase_FBgn) %>% unique() 
write.table(table1,'/home/eddydowle/Documents/pomonella/sigandDEgeneAndWrongDirabsdiffover05',sep='\t',quote=F,row.names=F)

#output whole table for 195
table1<-sigwithinFDR005 %>% left_join(.,test2,by='gene_id')  %>% dplyr::filter(.,signappleave_hawave=='Positive' & signappleearly_applelate=='Negative' & signhawearly_hawlate=='Negative' & abs(appleave_hawave) > 0.5 | signappleave_hawave=='Negative' & signappleearly_applelate=='Positive' & signhawearly_hawlate=='Positive' & abs(appleave_hawave) > 0.5) %>% unique() 

write.table(table1,'/home/eddydowle/Documents/pomonella/sigandDEgeneAndWrongDirabsdiffover05fullinfo',sep='\t',quote=F,row.names=F)


