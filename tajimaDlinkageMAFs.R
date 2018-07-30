#want to build in tajima D and ldx and radlinkage and then do heatmaps and look at average
#statistics per block

#ejd 2018

library(pheatmap)
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
library(RDAVIDWebService)



con <- dbConnect(MySQL(),user="raglandlab", password="pomonella",dbname="PomUrbanaGrant", host="localhost")

apple<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/prelimResults/Table.apple.JuneEJD",header=T,row.names=1,sep="\t",stringsAsFactors = F)
haw<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/prelimResults/Table.haw.JuneEJD",header=T,row.names=1,sep="\t",stringsAsFactors = F)

#to start with lets just work on the RNAtable
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


#pathways
#wnt<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/wntSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
#tor<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/torSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
#insulin<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/insulinSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)


wnt<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/wntSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
tor<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/torSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
insulin<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/insulinSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)

#https://stackoverflow.com/questions/48565661/dynamically-adding-and-removing-objects-in-response-to-selectinput-in-shiny


#flybase_symbol_ID<-read.csv("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/fbgn_annotation_ID_fb_2017_06.tsv",header=T,row.names=NULL,sep="\t",stringsAsFactors = F,strip.white=T)
flybase_symbol_ID<-read.csv("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/fbgn_annotation_ID_fb_2017_06.tsv",header=T,row.names=NULL,sep="\t",stringsAsFactors = F,strip.white=T)

my_palette <- colorRampPalette(c("darkturquoise", "darkgoldenrod2"))(n = 1000)


snpeff_effects<-c('3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant', '5_prime_UTR_variant', 'downstream_gene_variant', 'initiator_codon_variant', 'initiator_codon_variant&non_canonical_start_codon', 'initiator_codon_variant&splice_region_variant', 'intergenic_region', 'intragenic_variant', 'intron_variant', 'missense_variant', 'missense_variant&splice_region_variant', 'non_coding_transcript_exon_variant', 'non_coding_transcript_variant', 'splice_acceptor_variant&intron_variant', 'splice_acceptor_variant&splice_donor_variant&intron_variant', 'splice_donor_variant&intron_variant', 'splice_region_variant', 'splice_region_variant&initiator_codon_variant&non_canonical_start_codon', 'splice_region_variant&intron_variant', 'splice_region_variant&non_coding_transcript_exon_variant', 'splice_region_variant&stop_retained_variant', 'splice_region_variant&synonymous_variant', 'start_lost', 'start_lost&splice_region_variant', 'stop_gained', 'stop_gained&splice_region_variant', 'stop_lost', 'stop_lost&splice_region_variant', 'stop_retained_variant', 'synonymous_variant', 'upstream_gene_variant')
snpeff_effects
table_options<-c('scaffold','position','ref','alt','appleave_MAF','appleearly_MAF','applelate_MAF','hawave_MAF','hawearly_MAF','hawlate_MAF','protein_changes','impact'
                 ,'appleave_hawave_fisher_score','appleave_hawave_fisher_pvalue','appleave_hawave_fisher_pvalue_adjust','appleearly_applelate_fisher_score','appleearly_applelate_fisher_pvalue','appleearly_applelate_fisher_pvalue_adjust','hawearly_hawlate_fisher_score','hawearly_hawlate_fisher_pvalue','hawearly_hawlate_fisher_pvalue_adjust',"LDxApple","LDxHaw")
#link to mysql
#slider for ldx, fisher values type of snp "protein change" "upstream" etc

#thinking about:
#do network plots and then change highlighting?? with radio buttons??
#pathway package in R kegg pathways uses flybase IDs I think

#link to mysql
#slider for ldx, fisher values type of snp "protein change" "upstream" etc
#could do graph of gene projectory and then along side a table of SNPS?
#?column
input <- as.data.frame(t(as.data.frame(c('urbana',22,9),ncol=1,byrow=F)))
colnames(input) <- c("population","Low","Middle")
input$population
cleanup<-function(strsql){
  strsql<-strsql[strsql != "LDxApple"]
  strsql<-strsql[strsql != "LDxHaw"]
  strsql<-gsub("scaffold", paste("snpFisher",input$population,'.',"scaffold",sep=''), strsql)
  strsql<-gsub("position", paste("snpFisher",input$population,'.',"position",sep=''), strsql)
  strsql<-gsub("ref", paste("snpFisher",input$population,'.',"ref",sep=''), strsql)
  strsql<-gsub("alt", paste("snpFisher",input$population,'.',"alt",sep=''), strsql)
  strsql<-gsub('appleave_MAF',paste(input$population,'poolmaf.',input$population,'_appleave_Min',sep=''), strsql)
  strsql<-gsub('appleearly_MAF',paste(input$population,'poolmaf.',input$population,'_appleearly_Min',sep=''), strsql)
  strsql<-gsub('applelate_MAF',paste(input$population,'poolmaf.',input$population,'_applelate_Min',sep=''), strsql)
  strsql<-gsub('hawave_MAF',paste(input$population,'poolmaf.',input$population,'_hawave_Min',sep=''), strsql)
  strsql<-gsub('hawearly_MAF',paste(input$population,'poolmaf.',input$population,'_hawearly_Min',sep=''), strsql)
  strsql<-gsub('hawlate_MAF',paste(input$population,'poolmaf.',input$population,'_hawlate_Min',sep=''), strsql)
  strsql<-gsub("protein_changes", "annotation.protein_changes", strsql)
  strsql<-gsub("impact", "annotation.impact",strsql)
  strsql<-gsub('appleave_hawave_fisher_score', paste("snpFisher",input$population,'.',input$population,'_appleave_hawave_fisher_score',sep=''),strsql)
  strsql<-gsub('appleave_hawave_fisher_pvalue', paste("snpFisher",input$population,'.',input$population,'_appleave_hawave_fisher_pvalue',sep=''),strsql)
  strsql<-gsub('appleearly_applelate_fisher_score', paste("snpFisher",input$population,'.',input$population,'_appleearly_applelate_fisher_score',sep=''),strsql)
  strsql<-gsub('appleearly_applelate_fisher_pvalue', paste("snpFisher",input$population,'.',input$population,'_appleearly_applelate_fisher_pvalue',sep=''),strsql)
  strsql<-gsub('appleearly_applelate_fisher_pvalue_adjust', paste("snpFisher",input$population,'.',input$population,'_appleearly_applelate_fisher_pvalue_adjust',sep=''),strsql)
  strsql<-gsub('hawearly_hawlate_fisher_score', paste("snpFisher",input$population,'.',input$population,'_hawearly_hawlate_fisher_score',sep=''),strsql)
  strsql<-gsub('hawearly_hawlate_fisher_pvalue', paste("snpFisher",input$population,'.',input$population,'_hawearly_hawlate_fisher_pvalue',sep=''),strsql)
  strsql<-gsub('hawearly_hawlate_fisher_pvalue_adjust', paste("snpFisher",input$population,'.',input$population,'_hawearly_hawlate_fisher_pvalue_adjust',sep=''),strsql)
  return(strsql)
}


cleanupldx<-function(strsqlb){
  strsqlb<-gsub("LDxApple", paste(input$population,"poolLDx.MLE_AppleEarlyAppleLate_",input$population,sep=''), strsqlb)
  strsqlb<-gsub("LDxHaw", paste(input$population,"poolLDx.MLE_HawEarlyHawLate_",input$population,sep=''), strsqlb)
  strsqlb<-strsqlb[strsqlb != "LDxHaw"]
  strsqlb<-strsqlb[strsqlb != "scaffold"]
  strsqlb<-strsqlb[strsqlb != "position"]
  strsqlb<-strsqlb[strsqlb != "ref"]
  strsqlb<-strsqlb[strsqlb != "alt"]
  strsqlb<-strsqlb[strsqlb != 'appleave_MAF']
  strsqlb<-strsqlb[strsqlb != 'appleearly_MAF']
  strsqlb<-strsqlb[strsqlb != 'applelate_MAF']
  strsqlb<-strsqlb[strsqlb != 'hawave_MAF']
  strsqlb<-strsqlb[strsqlb != 'hawearly_MAF']
  strsqlb<-strsqlb[strsqlb != 'hawlate_MAF']
  strsqlb<-strsqlb[strsqlb != "protein_changes"]
  strsqlb<-strsqlb[strsqlb != "impact"]
  strsqlb<-strsqlb[strsqlb != 'appleave_hawave_fisher_score']
  strsqlb<-strsqlb[strsqlb != 'appleave_hawave_fisher_pvalue']
  strsqlb<-strsqlb[strsqlb != 'appleearly_applelate_fisher_score']
  strsqlb<-strsqlb[strsqlb != 'appleearly_applelate_fisher_pvalue']
  strsqlb<-strsqlb[strsqlb != 'appleearly_applelate_fisher_pvalue_adjust']
  strsqlb<-strsqlb[strsqlb != 'hawearly_hawlate_fisher_score']
  strsqlb<-strsqlb[strsqlb != 'hawearly_hawlate_fisher_pvalue']
  strsqlb<-strsqlb[strsqlb != 'hawearly_hawlate_fisher_pvalue_adjust']
  return(strsqlb)
}

table_options<-c('scaffold','position','ref','alt','appleave_MAF','appleearly_MAF','applelate_MAF','hawave_MAF','hawearly_MAF','hawlate_MAF','appleave_hawave_fisher_pvalue_adjust','appleearly_applelate_fisher_pvalue_adjust','hawearly_hawlate_fisher_pvalue_adjust',"MLE_AppleEarlyAppleLate_urbana","MLE_HawEarlyHawLate_urbana","LDgr","gene_id","loc","Flybase_gene_symbol","Flybase_FBgn","effect","Urbana_apple_num_gene_TajD","Urbana_haw_num_gene_TajD")

ui<-fluidPage(
  conditionalPanel(condition="input.conditionedPanels==1",
                   fluidRow(column(3,
                                   radioButtons('population', label = 'Population choice for SNP',
                                                choices = c("urbana","grant"),selected='urbana'),
                                   sliderInput("integer", "Absolute MAF difference:",
                                               min = 0, max = 1,
                                               value = 0.3)),
                        column(4,offset=1,
                               checkboxGroupButtons('comparisons', label='MAF comparisons',
                                            choices=c('AppleAverage & HawAverage','AppleEarly & AppleLate','HawEarly & HawLate'),
                                            selected=c('AppleAverage & HawAverage','HawEarly & HawLate')), #AppleAverage & HawAverage','HawEarly & HawLate'AppleEarly & AppleLate', 'HawEarly & HawLate'
                               sliderInput("pvalue", "Adjusted pvalue (is true for populations in comparisons):",
                                           min = 0, max = 1,
                                           value = 0.05)),
                        column(4,
                                   pickerInput(inputId = "tableoptions", 
                                               label = "Select/deselect raw table output", 
                                               choices = table_options, options = list(`max-options` = 10),multiple = TRUE),
                radioButtons('tabletype', label = 'Table Type',
                             choices = c("Gene Summary","Raw SNPs"),selected='Raw SNPs')))),
  conditionalPanel(condition="input.conditionedPanels==2",
                   fluidRow(column(3,
                                   radioButtons('population2', label = 'Population choice for SNP (only urbana works)',
                                                choices = c("urbana","grant"),selected='urbana'),
                                   sliderInput("tajimaDApple", "Absolute TajimaD values Apple:",
                                               min = 0, max = 10,
                                               value = 2)),
                            sliderInput("tajimaDHaw", "Absolute TajimaD values Haw:",
                                        min = 0, max = 10,
                                        value = 2)),
                   column(4, offset=1,
                                   radioButtons('sign', label = 'Sign of TajimaD values',
                                                choices = c("Both","positive","negative"),selected='Both'),
                                   radioButtons('tabletype2', label = "Summary, Raw or David Table",
                                                choices = c ("David","Raw","Summary Gene"),selected='David')),
                             #      uiOutput(outputId = "gene2")),
                            column(4,
                                   sliderInput("pvalue2", "Adjusted pvalue AppleAve and HawAve:",
                                               min = 0, max = 1,
                                               value = 0.5),
          #                         uiOutput(outputId = "module_selection"),
           #                        uiOutput(outputId = "gene_module"),
                                   uiOutput(outputId = "gene_module_enrichment"))),

  mainPanel(
    h4("SNP MAF TajimaD pvalue"),
    tabsetPanel(
      tabPanel(
        h4("Tables MAF pvalues"),
        plotOutput(outputId = "main_plot") ,
        tableOutput("SNP_table"),
        value=1 
      ),
      tabPanel(
        h4("TajimaD"),
        tableOutput("David_table"),
        value=2
      )
       , id='conditionedPanels'
    )))



server<-function(input,output) {
  
  mysqlcall <- reactive({
    strsql<-input$comparisons
    pop=input$population
    if ('AppleAverage & HawAverage' %in% strsql & 'HawEarly & HawLate' %in% strsql & !('AppleEarly & AppleLate' %in% strsql)){
      sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.scaffold, testpoolmaf.position, testpoolmaf.ref, testpoolmaf.alt, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, feature_alias.Flybase_gene_symbol, feature_alias.Flybase_FBgn, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleave_hawave_fisher_pvalue_adjust < ", input$pvalue, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", input$pvalue, " AND ABS(testpoolmaf.test_appleave_Maj - testpoolmaf.test_hawave_Maj) > ", input$integer, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", input$integer, ";")
      sql<-gsub('test',input$population,sql)
      sql<-paste(sql,collapse='')
      query <- dbGetQuery(con, sql)
      final_table<-query
      return(final_table)
    }
    else if ('AppleAverage & HawAverage' %in% strsql & 'AppleEarly & AppleLate' %in% strsql & !('HawEarly & HawLate' %in% strsql)){
      #    strsql<-input$comparisons
      #    pop=input$population
      
      sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.scaffold, testpoolmaf.position, testpoolmaf.ref, testpoolmaf.alt, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, feature_alias.Flybase_gene_symbol, feature_alias.Flybase_FBgn, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE  snpFishertest.test_appleave_hawave_fisher_pvalue_adjust < ", input$pvalue, " AND snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust < ", input$pvalue, " AND ABS (testpoolmaf.test_appleearly_Maj - testpoolmaf.test_applelate_Maj) > ", input$integer, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", input$integer, ";")
      sql<-gsub('test',input$population,sql)
      sql<-paste(sql,collapse='')
      query <- dbGetQuery(con, sql)
      final_table<-query
      return(final_table)
    }
    else if ('AppleEarly & AppleLate' %in% strsql & 'HawEarly & HawLate' %in% strsql & !('AppleAverage & HawAverage' %in% strsql)){
      #   strsql<-input$comparisons
      #   pop=input$population
       sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.scaffold, testpoolmaf.position, testpoolmaf.ref, testpoolmaf.alt, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, feature_alias.Flybase_gene_symbol, feature_alias.Flybase_FBgn, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust < ", input$pvalue, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", input$pvalue, " AND ABS(testpoolmaf.test_appleearly_Maj - testpoolmaf.test_applelate_Maj) > ", input$integer, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", input$integer, ";")
      sql<-gsub('test',input$population,sql)
      sql<-paste(sql,collapse='')
      query <- dbGetQuery(con, sql)
      final_table<-query
      return(final_table)
    }
    else {
      #   strsql<-input$comparisons
      #   pop=input$population
      sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.scaffold, testpoolmaf.position, testpoolmaf.ref, testpoolmaf.alt, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, feature_alias.Flybase_gene_symbol, feature_alias.Flybase_FBgn, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleave_hawave_fisher_pvalue_adjust < ", input$pvalue, " AND snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust < ", input$pvalue, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", input$pvalue, " AND ABS(testpoolmaf.test_appleave_Maj - testpoolmaf.test_hawave_Maj) > ", input$integer, " AND ABS (testpoolmaf.test_appleearly_Maj - testpoolmaf.test_applelate_Maj) > ", input$integer, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", input$integer, ";")
      sql<-gsub('test',input$population,sql)
      sql<-paste(sql,collapse='')
      query <- dbGetQuery(con, sql)
      final_table<-query
      return(final_table)
    }
  })
 
  #build an escape for when there is no rows in the table
  
   output$main_plot<-renderPlot({
    strsql<-input$comparisons
    pop=input$population
    if ('AppleAverage & HawAverage' %in% strsql & 'HawEarly & HawLate' %in% strsql & !('AppleEarly & AppleLate' %in% strsql)){
     final_table<-mysqlcall()
      out.snp<-final_table[!duplicated(final_table[,'snpId']),]
      test1<-out.snp %>% dplyr::select(.,matches('appleave_Maj|hawave_Maj')) %>% mutate(appleave_hawave = .[,1]-.[,2]) %>% dplyr::select(.,appleave_hawave)
      test2<-out.snp %>% dplyr::select(.,matches('hawearly_Maj|hawlate_Maj')) %>% mutate(hawearly_hawlate = .[,1]-.[,2]) %>% dplyr::select(.,hawearly_hawlate)
      test<-cbind(test1,test2) %>% as.matrix
    }

     else if ('AppleAverage & HawAverage' %in% strsql & 'AppleEarly & AppleLate' %in% strsql & !('HawEarly & HawLate' %in% strsql)){
      final_table<-mysqlcall()
      out.snp<-final_table[!duplicated(final_table[,'snpId']),]
      test1<-out.snp %>% dplyr::select(.,matches('appleave_Maj|hawave_Maj')) %>% mutate(appleave_hawave = .[,1]-.[,2]) %>% dplyr::select(.,appleave_hawave)
      test2<-out.snp %>% dplyr::select(.,matches('appleearly_Maj|applelate_Maj')) %>% mutate(appleearly_applelate = .[,1]-.[,2]) %>% dplyr::select(.,appleearly_applelate)
      test<-cbind(test1,test2) %>% as.matrix
       }
    
    else if ('AppleEarly & AppleLate' %in% strsql & 'HawEarly & HawLate' %in% strsql & !('AppleAverage & HawAverage' %in% strsql)){
      final_table<-mysqlcall()
      out.snp<-final_table[!duplicated(final_table[,'snpId']),]
      test1<-out.snp %>% dplyr::select(.,matches('appleearly_Maj|applelate_Maj')) %>% mutate(appleearly_applelate = .[,1]-.[,2]) %>% dplyr::select(.,appleearly_applelate)
      test2<-out.snp %>% dplyr::select(.,matches('hawearly_Maj|hawlate_Maj')) %>% mutate(hawearly_hawlate = .[,1]-.[,2]) %>% dplyr::select(.,hawearly_hawlate)
      test<-cbind(test1,test2) %>% as.matrix
       }

        else {
      final_table<-mysqlcall()
      out.snp<-final_table[!duplicated(final_table[,'snpId']),]
      test1<-out.snp %>% dplyr::select(.,matches('appleave_Maj|hawave_Maj')) %>% mutate(appleave_hawave = .[,1]-.[,2]) %>% dplyr::select(.,appleave_hawave)
      test2<-out.snp %>% dplyr::select(.,matches('hawearly_Maj|hawlate_Maj')) %>% mutate(hawearly_hawlate = .[,1]-.[,2]) %>% dplyr::select(.,hawearly_hawlate)
      test3<-out.snp %>% dplyr::select(.,matches('appleearly_Maj|applelate_Maj')) %>% mutate(appleearly_applelate = .[,1]-.[,2]) %>% dplyr::select(.,appleearly_applelate)
      test<-cbind(test1,test2,test3) %>% as.matrix
      
       }
 #   library(pheatmap)
    #?pheatmap
    pheatmap(test,col = my_palette,treeheight_row=0,treeheight_col=0)
    
  })
  
   #output summary table
   
   #output raw table 
  output$SNP_table<-renderTable({
    final_table<-mysqlcall()
    tabtype<-input$tabletype
    if (tabtype=='Raw SNPs') {
    out.snp<-final_table[!duplicated(final_table),]
    tablepicks<-input$tableoptions
    tablepicks<-c("snpId",tablepicks)
    tablepicks<-gsub("appleave_MAF", paste(input$population,"_appleave_MAF",input$population,sep=''),tablepicks)
    tablepicks<-gsub("appleearly_MAF", paste(input$population,"_appleearly_MAF",input$population,sep=''),tablepicks)
    tablepicks<-gsub("applelate_MAF", paste(input$population,"_applelate_MAF",input$population,sep=''),tablepicks)
    tablepicks<-gsub("hawave_MAF", paste(input$population,"_hawave_MAF",input$population,sep=''),tablepicks)
    tablepicks<-gsub("hawearly_MAF", paste(input$population,"_hawearly_MAF",input$population,sep=''),tablepicks)
    tablepicks<-gsub("hawlate_MAF", paste(input$population,"_hawlate_MAF",input$population,sep=''),tablepicks)
   tableout<- out.snp[, tablepicks]
    }
    
    #table_options<-c('scaffold','position','ref','alt','appleave_MAF','appleearly_MAF','applelate_MAF','hawave_MAF','hawearly_MAF','hawlate_MAF','appleave_hawave_fisher_pvalue_adjust','appleearly_applelate_fisher_pvalue_adjust','hawearly_hawlate_fisher_pvalue_adjust',"MLE_AppleEarlyAppleLate_urbana","MLE_HawEarlyHawLate_urbana","LDgr","gene_id","loc","Flybase_gene_symbol","Flybase_FBgn","effect","Urbana_apple_num_gene_TajD","Urbana_haw_num_gene_TajD")                             choices = c("Gene Summary","Raw SNPs"),selected='Raw SNPs')))),


    else if (tabtype=='Gene Summary') {
      
      strsql<-input$comparisons
      pop=input$population
      if ('AppleAverage & HawAverage' %in% strsql & 'HawEarly & HawLate' %in% strsql & !('AppleEarly & AppleLate' %in% strsql)){
        out.snp<-final_table[!duplicated(final_table[,'loc']),]
        test1<-out.snp %>% dplyr::select(.,matches('appleave_Maj|hawave_Maj'),loc,Urbana_apple_num_gene_TajD,Urbana_haw_num_gene_TajD,LDgr,MLE_AppleEarlyAppleLate_urbana,MLE_HawEarlyHawLate_urbana) %>% mutate(appleave_hawave = .[,1]-.[,2]) %>% dplyr::select(.,appleave_hawave,loc,Urbana_apple_num_gene_TajD,Urbana_haw_num_gene_TajD,LDgr,MLE_AppleEarlyAppleLate_urbana,MLE_HawEarlyHawLate_urbana)
        test2<-out.snp %>% dplyr::select(.,matches('hawearly_Maj|hawlate_Maj')) %>% mutate(hawearly_hawlate = .[,1]-.[,2]) %>% dplyr::select(.,hawearly_hawlate)
        test<-cbind(test1,test2) 
        test2<-test %>% mutate(signappleave_hawave=ifelse(appleave_hawave > 0,'Positive','Negative'))  %>% mutate(signhawearly_hawlate =ifelse(hawearly_hawlate > 0,'Positive','Negative'))
        sumtab<-test2 %>% group_by(signappleave_hawave,signhawearly_hawlate)  %>% summarise(loc_count=n_distinct(loc,na.rm = TRUE), Urbana_apple_num_gene_TajD_mean=mean(Urbana_apple_num_gene_TajD,na.rm = TRUE),Urbana_haw_num_gene_TajD_mean=mean(Urbana_haw_num_gene_TajD,na.rm = TRUE), MLE_AppleEarlyAppleLate_urbana_mean= mean(MLE_AppleEarlyAppleLate_urbana,na.rm = TRUE), MLE_HawEarlyHawLate_urbana_mean=mean(MLE_HawEarlyHawLate_urbana,na.rm = TRUE),LDgr_countHigh=sum(LDgr %in% c("H")),LDgr_countMed=sum(LDgr %in% c("M")),LDgr_countLow=sum(LDgr %in% c("L")))
      }
      else if ('AppleAverage & HawAverage' %in% strsql & 'AppleEarly & AppleLate' %in% strsql & !('HawEarly & HawLate' %in% strsql)){
        final_table<-mysqlcall()
        out.snp<-final_table[!duplicated(final_table[,'loc']),]
        test1<-out.snp %>% dplyr::select(.,matches('appleave_Maj|hawave_Maj'),loc,Urbana_apple_num_gene_TajD,Urbana_haw_num_gene_TajD,LDgr,MLE_AppleEarlyAppleLate_urbana,MLE_HawEarlyHawLate_urbana) %>% mutate(appleave_hawave = .[,1]-.[,2]) %>% dplyr::select(.,appleave_hawave,loc,Urbana_apple_num_gene_TajD,Urbana_haw_num_gene_TajD,LDgr,MLE_AppleEarlyAppleLate_urbana,MLE_HawEarlyHawLate_urbana)
        test2<-out.snp %>% dplyr::select(.,matches('appleearly_Maj|applelate_Maj')) %>% mutate(appleearly_applelate = .[,1]-.[,2]) %>% dplyr::select(.,appleearly_applelate)
        test<-cbind(test1,test2) 
        test2<-test %>% mutate(signappleave_hawave=ifelse(appleave_hawave > 0,'Positive','Negative'))  %>% mutate(signappleearly_applelate =ifelse(appleearly_applelate > 0,'Positive','Negative'))
        sumtab<-test2 %>% group_by(signappleave_hawave,signappleearly_applelate)  %>% summarise(loc_count=n_distinct(loc,na.rm = TRUE), Urbana_apple_num_gene_TajD_mean=mean(Urbana_apple_num_gene_TajD,na.rm = TRUE),Urbana_haw_num_gene_TajD_mean=mean(Urbana_haw_num_gene_TajD,na.rm = TRUE),MLE_AppleEarlyAppleLate_urbana_mean= mean(MLE_AppleEarlyAppleLate_urbana,na.rm = TRUE),MLE_HawEarlyHawLate_urbana_mean=mean(MLE_HawEarlyHawLate_urbana,na.rm = TRUE),LDgr_countHigh=sum(LDgr %in% c("H")),LDgr_countMed=sum(LDgr %in% c("M")),LDgr_countLow=sum(LDgr %in% c("L")))
        
      }
      
      else if ('AppleEarly & AppleLate' %in% strsql & 'HawEarly & HawLate' %in% strsql & !('AppleAverage & HawAverage' %in% strsql)){
        final_table<-mysqlcall()
        out.snp<-final_table[!duplicated(final_table[,'loc']),]
        test1<-out.snp %>% dplyr::select(.,matches('appleearly_Maj|applelate_Maj'),loc,Urbana_apple_num_gene_TajD,Urbana_haw_num_gene_TajD,LDgr,MLE_AppleEarlyAppleLate_urbana,MLE_HawEarlyHawLate_urbana) %>% mutate(appleearly_applelate = .[,1]-.[,2]) %>% dplyr::select(.,appleearly_applelate,loc,Urbana_apple_num_gene_TajD,Urbana_haw_num_gene_TajD,LDgr,MLE_AppleEarlyAppleLate_urbana,MLE_HawEarlyHawLate_urbana)
        test2<-out.snp %>% dplyr::select(.,matches('hawearly_Maj|hawlate_Maj')) %>% mutate(hawearly_hawlate = .[,1]-.[,2]) %>% dplyr::select(.,hawearly_hawlate)
        test<-cbind(test1,test2) 
        test2<-test %>% mutate(signappleearly_applelate=ifelse(appleearly_applelate > 0,'Positive','Negative'))  %>% mutate(signhawearly_hawlate =ifelse(hawearly_hawlate > 0,'Positive','Negative'))
        sumtab<-test2 %>% group_by(signappleearly_applelate,signhawearly_hawlate)  %>% summarise(loc_count=n_distinct(loc,na.rm = TRUE), Urbana_apple_num_gene_TajD_mean=mean(Urbana_apple_num_gene_TajD,na.rm = TRUE),Urbana_haw_num_gene_TajD_mean=mean(Urbana_haw_num_gene_TajD,na.rm = TRUE),MLE_AppleEarlyAppleLate_urbana_mean= mean(MLE_AppleEarlyAppleLate_urbana,na.rm = TRUE),MLE_HawEarlyHawLate_urbana_mean=mean(MLE_HawEarlyHawLate_urbana,na.rm = TRUE),LDgr_countHigh=sum(LDgr %in% c("H")),LDgr_countMed=sum(LDgr %in% c("M")),LDgr_countLow=sum(LDgr %in% c("L")))
        
      }
      
      else {
        final_table<-mysqlcall()
        out.snp<-final_table[!duplicated(final_table[,'loc']),]
        test1<-out.snp %>% dplyr::select(.,matches('appleave_Maj|hawave_Maj'),loc,Urbana_apple_num_gene_TajD,Urbana_haw_num_gene_TajD,LDgr,MLE_AppleEarlyAppleLate_urbana,MLE_HawEarlyHawLate_urbana) %>% mutate(appleave_hawave = .[,1]-.[,2]) %>% dplyr::select(.,appleave_hawave,loc,Urbana_apple_num_gene_TajD,Urbana_haw_num_gene_TajD,LDgr,MLE_AppleEarlyAppleLate_urbana,MLE_HawEarlyHawLate_urbana)
        test2<-out.snp %>% dplyr::select(.,matches('hawearly_Maj|hawlate_Maj')) %>% mutate(hawearly_hawlate = .[,1]-.[,2]) %>% dplyr::select(.,hawearly_hawlate)
        test3<-out.snp %>% dplyr::select(.,matches('appleearly_Maj|applelate_Maj')) %>% mutate(appleearly_applelate = .[,1]-.[,2]) %>% dplyr::select(.,appleearly_applelate)
        test<-cbind(test1,test2,test3) 
        test2<-test %>% mutate(signappleave_hawave=ifelse(appleave_hawave > 0,'Positive','Negative')) %>% mutate(signappleearly_applelate=ifelse(appleearly_applelate > 0,'Positive','Negative'))  %>% mutate(signhawearly_hawlate =ifelse(hawearly_hawlate > 0,'Positive','Negative')) 
        sumtab<-test2 %>% group_by(signappleave_hawave,signappleearly_applelate,signhawearly_hawlate)  %>% summarise(loc_count=n_distinct(loc,na.rm = TRUE), Urbana_apple_num_gene_TajD_mean=mean(Urbana_apple_num_gene_TajD,na.rm = TRUE),Urbana_haw_num_gene_TajD_mean=mean(Urbana_haw_num_gene_TajD,na.rm = TRUE),MLE_AppleEarlyAppleLate_urbana_mean= mean(MLE_AppleEarlyAppleLate_urbana,na.rm = TRUE),MLE_HawEarlyHawLate_urbana_mean=mean(MLE_HawEarlyHawLate_urbana,na.rm = TRUE),LDgr_countHigh=sum(LDgr %in% c("H")),LDgr_countMed=sum(LDgr %in% c("M")),LDgr_countLow=sum(LDgr %in% c("L")))
        
      }
      
      tableout<-sumtab
    }
    tableout
     })
  
  FuncAnnotClust <- reactive({
      if(input$sign=='Both'){
        sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.scaffold, testpoolmaf.position, testpoolmaf.ref, testpoolmaf.alt, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, feature_alias.Flybase_gene_symbol, feature_alias.Flybase_FBgn, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleave_hawave_fisher_pvalue_adjust < ", input$pvalue2, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", input$pvalue2, " AND ABS (urbanapoolTajDgene.Urbana_apple_num_gene_TajD) > ", input$tajimaDApple, " AND ABS (urbanapoolTajDgene.Urbana_haw_num_gene_TajD) >", input$tajimaDHaw, ";")
        sql<-gsub('test',"urbana",sql)
        sql<-paste(sql,collapse='')
        query2 <- dbGetQuery(con, sql)
        final_table2<-query2 
        #  final_table2<-mysqlcalltab2()
 #     df6<-final_table2 %>% filter (., abs(Urbana_apple_num_gene_TajD) > input$tajimaD, abs(Urbana_haw_num_gene_TajD) > input$tajimaD) %>% distinct(.,loc,.keep_all = TRUE)
 #     test1<-DAVIDWebService$new(email='eddy.dowle@otago.ac.nz',url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
#      FG <- addList(test1, final_table2$Flybase_FBgn, idType="FLYBASE_GENE_ID", listName="isClass", listType="Gene")
#      FuncAnnotClust1 <- getClusterReport(test1)
#      FuncAnnotClust1
#      return(FuncAnnotClust1)
    }
    else if(input$sign=='positive'){
    #  final_table2<-mysqlcalltab2()
      sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.scaffold, testpoolmaf.position, testpoolmaf.ref, testpoolmaf.alt, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, feature_alias.Flybase_gene_symbol, feature_alias.Flybase_FBgn, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleave_hawave_fisher_pvalue_adjust < ", input$pvalue2, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", input$pvalue2, " AND urbanapoolTajDgene.Urbana_apple_num_gene_TajD > ", input$tajimaDApple, " AND urbanapoolTajDgene.Urbana_haw_num_gene_TajD >", input$tajimaDHaw, ";")
      sql<-gsub('test',"urbana",sql)
      sql<-paste(sql,collapse='')
      query2 <- dbGetQuery(con, sql)
      final_table2<-query2 #%>% distinct(.,loc,.keep_all = TRUE)
  #   df6<-final_table2 %>% filter (., Urbana_apple_num_gene_TajD > input$tajimaD, Urbana_haw_num_gene_TajD > input$tajimaD)  %>% distinct(.,loc,.keep_all = TRUE)
   #   test2<-DAVIDWebService$new(email='eddy.dowle@otago.ac.nz',url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  #    FG <- addList(test2, final_table2$Flybase_FBgn, idType="FLYBASE_GENE_ID", listName="isClass", listType="Gene")
  #    FuncAnnotClust1 <- getClusterReport(test2)
      #     FuncAnnotClust1
   #   return(FuncAnnotClust1)
    }
    else if(input$sign=='negative'){
      sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.scaffold, testpoolmaf.position, testpoolmaf.ref, testpoolmaf.alt, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, feature_alias.Flybase_gene_symbol, feature_alias.Flybase_FBgn, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleave_hawave_fisher_pvalue_adjust < ", input$pvalue2, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", input$pvalue2, " AND urbanapoolTajDgene.Urbana_apple_num_gene_TajD < -", input$tajimaDApple, " AND urbanapoolTajDgene.Urbana_haw_num_gene_TajD < -", input$tajimaDHaw, ";")
      sql<-gsub('test',"urbana",sql)
      sql<-paste(sql,collapse='')
      query2 <- dbGetQuery(con, sql)
      final_table2<-query2 #%>% distinct(.,loc,.keep_all = TRUE)
    }
    
    if (input$tabletype2=="David"){
      final_table2<-final_table2 %>% distinct(.,loc,.keep_all = TRUE)
#      df6<-final_table2 %>% filter (., Urbana_apple_num_gene_TajD < -input$tajimaD, Urbana_haw_num_gene_TajD < -input$tajimaD) %>% distinct(.,loc,.keep_all = TRUE)
      test2<-DAVIDWebService$new(email='eddy.dowle@otago.ac.nz',url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
      FG <- addList(test2, final_table2$Flybase_FBgn, idType="FLYBASE_GENE_ID", listName="isClass", listType="Gene")
      FuncAnnotClust1 <- getClusterReport(test2)
      #     FuncAnnotClust1
      return(FuncAnnotClust1)
  }
 else (
   return(final_table2))
  })
  
  output$gene_module_enrichment <- renderUI({
    if (input$tabletype2=="David"){
    FuncAnnotClust<-FuncAnnotClust()
    choicesenrich=1:nrow(summary(FuncAnnotClust))
    selectizeInput(inputId = "table_var2",
                   label= "David Cluster (switch to raw/summary for pathways)",
                   #                  choices=together.across %>% distinct(.,gene_id), 
                   choices=choicesenrich,
                   options=list(maxOptions=3000))
    }
    else {
      selectizeInput(inputId = "table_var3",
                     label= "Pathway",
                     #                  choices=together.across %>% distinct(.,gene_id), 
                     choices=c("All","insulin","wnt","tor"),
                     selected='All')
      }
  })
#  radioButtons('tabletype', label = "Raw or David Table",
 #              choices = c ("David","Raw"),selected='David')),

  output$David_table<-renderTable({
    if(input$tabletype2=="David"){
       clus<-as.integer(input$table_var2)
      FuncAnnotClust<-FuncAnnotClust()
      final_table<-members(FuncAnnotClust)[[clus]] %>% dplyr::select(-one_of("Genes"))
    }
    else if  (input$tabletype2=="Raw"){
      final_table<-FuncAnnotClust()
      if (input$table_var3=='All')
      {
        final_table<-final_table
      }
      else{
        test<-input$table_var3
        final_table<-final_table %>% filter(., Flybase_FBgn %in% get(test)$flyid)
      }
    }
    else if  (input$tabletype2=="Summary Gene"){
        final_table<-FuncAnnotClust()
        if (input$table_var3=='All')
        {
        out.snp<-final_table[!duplicated(final_table[,'loc']),]
        test1<-out.snp %>% dplyr::select(.,loc,Flybase_gene_symbol,snpId,Flybase_FBgn,Urbana_apple_num_gene_TajD,Urbana_haw_num_gene_TajD,LDgr,MLE_AppleEarlyAppleLate_urbana,MLE_HawEarlyHawLate_urbana) %>% group_by(loc,Flybase_gene_symbol,Flybase_FBgn)  %>% summarise(snp_count=n_distinct(snpId,na.rm = TRUE), Urbana_apple_num_gene_TajD_mean=mean(Urbana_apple_num_gene_TajD,na.rm = TRUE),Urbana_haw_num_gene_TajD_mean=mean(Urbana_haw_num_gene_TajD,na.rm = TRUE),MLE_AppleEarlyAppleLate_urbana_mean= mean(MLE_AppleEarlyAppleLate_urbana,na.rm = TRUE),MLE_HawEarlyHawLate_urbana_mean=mean(MLE_HawEarlyHawLate_urbana,na.rm = TRUE),LDgr_countHigh=sum(LDgr %in% c("H")),LDgr_countMed=sum(LDgr %in% c("M")),LDgr_countLow=sum(LDgr %in% c("L")))
        final_table<-test1
        }
        else{
          test<-input$table_var3
          out.snp<-final_table %>% filter(., Flybase_FBgn %in% get(test)$flyid) %>% distinct(.,loc,Flybase_gene_symbol,Flybase_FBgn,.keep_all = TRUE)
          test1<-out.snp %>% dplyr::select(.,loc,Flybase_gene_symbol,snpId,Flybase_FBgn,Urbana_apple_num_gene_TajD,Urbana_haw_num_gene_TajD,LDgr,MLE_AppleEarlyAppleLate_urbana,MLE_HawEarlyHawLate_urbana) %>% group_by(loc,Flybase_gene_symbol,Flybase_FBgn)  %>% summarise(snp_count=n_distinct(snpId,na.rm = TRUE), Urbana_apple_num_gene_TajD_mean=mean(Urbana_apple_num_gene_TajD,na.rm = TRUE),Urbana_haw_num_gene_TajD_mean=mean(Urbana_haw_num_gene_TajD,na.rm = TRUE),MLE_AppleEarlyAppleLate_urbana_mean= mean(MLE_AppleEarlyAppleLate_urbana,na.rm = TRUE),MLE_HawEarlyHawLate_urbana_mean=mean(MLE_HawEarlyHawLate_urbana,na.rm = TRUE),LDgr_countHigh=sum(LDgr %in% c("H")),LDgr_countMed=sum(LDgr %in% c("M")),LDgr_countLow=sum(LDgr %in% c("L")))
          final_table<-test1
          
          }
         }
        # final_table
    #  final_table<-head(together)
    final_table
  })
  
 
  
}

shinyApp(ui=ui,server=server)












#ok so lets look at joining the tajimaD estiamte in to the mysql return

#ok so for tajimas d I have built the gtf file into the mysql
#this means we can get the gene start and stop position and then average the tajima's D estiamte 
#across the 'gene region' ??????

test<-0.05
test2<-0.5

sql <- c( "SELECT annotation.snpId, annotation.loc, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana,RAD_linkage_matchPool.LDgr, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj FROM annotation, snpFishertest, testpoolmaf WHERE snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust < ", 0.05, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", 0.05, " AND ABS(testpoolmaf.test_appleearly_Maj - testpoolmaf.test_applelate_Maj) > ", 0.5, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", 0.5, " AND testpoolmaf.snpId=annotation.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN RAD_linkage_matchPool ON annotation.snpId=RAD_linkage_matchPool.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id;")

library(shiny)
library(pool)
library(dplyr)
library(RMySQL)

con <- dbConnect(MySQL(),user="raglandlab", password="pomonella",dbname="PomUrbanaGrant", host="localhost")

sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD, feature_alias.loc, urbanapoolTajDgene.gene_id, snpFishertest.snpId, urbanapoolLDx.snpId, annotation.snpId, RAD_linkage_matchPool.snpId FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust < ", 0.05, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", 0.05, " AND ABS(testpoolmaf.test_appleearly_Maj - testpoolmaf.test_applelate_Maj) > ", 0.5, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", 0.5, ";")

sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, feature_alias.Flybase_gene_symbol, feature_alias.Flybase_FBgn, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust < ", 0.05, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", 0.05, " AND ABS(testpoolmaf.test_appleearly_Maj - testpoolmaf.test_applelate_Maj) > ", 0.5, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", 0.5, ";")


sql<-gsub('test','urbana',sql)

strsql<-paste(sql,collapse='')
#gsub('test','urbana',sql)
strsql
query <- dbGetQuery(con, strsql)
final_table<-query
head(final_table)
nrow(final_table)


sql <- c( "SELECT annotation.snpId, annotation.loc, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana,RAD_linkage_matchPool.LDgr, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj FROM annotation, snpFishertest, testpoolmaf LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id AND snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust < ", 0.05, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", 0.05, " AND ABS(testpoolmaf.test_appleearly_Maj - testpoolmaf.test_applelate_Maj) > ", 0.5, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", 0.5, " AND testpoolmaf.snpId=annotation.snpId;")

sql<-c("SELECT * FROM annotation.snpId LEFT JOIN  feature_alias ON annotation.loc=feature_alias.loc AND annotation.loc=LOC108374498")


sql <- c( "SELECT testpoolmaf.snpId, annotation.snpId, annotation.loc, urbanapoolLDx.snpId, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, feature_alias.loc,feature_alias.gene_id, RAD_linkage_matchPool.snpId,RAD_linkage_matchPool.LDgr,urbanapoolTajDgene.gene_id, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD, snpFishertest.snpId, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust < ", 0.05, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", 0.05, " AND ABS(testpoolmaf.test_appleearly_Maj - testpoolmaf.test_applelate_Maj) > ", 0.5, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", 0.5, " AND testpoolmaf.snpId=annotation.snpId;")

LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id 

 feature_alias.loc,feature_alias.gene_id, urbanapoolTajDgene.gene_id, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD,

 sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust < ", 0.05, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", 0.05, " AND ABS(testpoolmaf.test_appleearly_Maj - testpoolmaf.test_applelate_Maj) > ", 0.5, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", 0.5, ";")
 
 sql <- c( "SELECT testpoolmaf.snpId, testpoolmaf.scaffold, testpoolmaf.position, testpoolmaf.ref, testpoolmaf.alt, testpoolmaf.test_appleave_Maj, testpoolmaf.test_hawave_Maj, testpoolmaf.test_appleearly_Maj,testpoolmaf.test_applelate_Maj, testpoolmaf.test_hawearly_Maj, testpoolmaf.test_hawlate_Maj, snpFishertest.test_appleave_hawave_fisher_pvalue_adjust, snpFishertest.test_appleearly_applelate_fisher_pvalue_adjust, snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana, RAD_linkage_matchPool.LDgr, feature_alias.gene_id, annotation.loc, feature_alias.Flybase_gene_symbol, feature_alias.Flybase_FBgn, annotation.effect, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD FROM testpoolmaf LEFT JOIN annotation ON testpoolmaf.snpId=annotation.snpId LEFT JOIN snpFishertest ON testpoolmaf.snpId=snpFishertest.snpId LEFT JOIN urbanapoolLDx ON testpoolmaf.snpId=urbanapoolLDx.snpId LEFT JOIN feature_alias ON annotation.loc=feature_alias.loc LEFT JOIN urbanapoolTajDgene ON feature_alias.gene_id=urbanapoolTajDgene.gene_id LEFT JOIN RAD_linkage_matchPool ON testpoolmaf.snpId=RAD_linkage_matchPool.snpId WHERE snpFishertest.test_appleave_hawave_fisher_pvalue_adjust < ", 0.05, " AND snpFishertest.test_hawearly_hawlate_fisher_pvalue_adjust < ", 0.05, " AND ABS(testpoolmaf.test_appleave_Maj - testpoolmaf.test_hawave_Maj) > ", 0.5, " AND ABS (testpoolmaf.test_hawearly_Maj - testpoolmaf.test_hawlate_Maj) > ", 0.5, ";")
 sql<-gsub('test','urbana',sql)
 sql<-paste(sql,collapse='')
 query <- dbGetQuery(con, sql)
 final_table<-query
 
sql<-gsub('test','urbana',sql)

strsql<-paste(sql,collapse='')
#gsub('test','urbana',sql)
strsql
query <- dbGetQuery(con, strsql)
final_table2<-query

head(final_table2)
nrow(final_table2)


SELECT * FROM annotation LEFT JOIN  feature_alias ON annotation.loc=feature_alias.loc AND annotation.loc=LOC108374498
#so we want a heatmap of mafs

out=NULL
for(i in 1:length(unique(together$gene_id))){
  row<-(unique(together$gene_id[i]))
  query <- dbGetQuery(con, sql)
  out=rbind(out,query)
}
out

out.snp<-out[!duplicated(out[,'snpId']),]

test<-as.matrix(cbind((out.snp$urbana_hawearly_Maj- out.snp$urbana_hawlate_Maj),(out.snp$urbana_appleave_Maj- out.snp$urbana_hawave_Maj)))
my_palette <- colorRampPalette(c("darkturquoise", "darkgoldenrod2"))(n = 1000)
library(pheatmap)
#?pheatmap
pheatmap(test,col = my_palette,treeheight_row=0,treeheight_col=0)





#pheatmap(test,col = my_palette,treeheight_row=0,treeheight_col=0, cluster_cols=FALSE, cluster_rows=FALSE)

#so some of them are in the wrong direction than others which genes are in which direction
out.snp$hawearlylate<- out.snp$hawearly_MAF- out.snp$hawlate_MAF
out.snp$applehawave<-out.snp$appleave_MAF- out.snp$hawave_MAF

out.snp.samesign<-out.snp[sign(out.snp$hawearlylate)==sign(out.snp$applehawave),]
test<-as.matrix(cbind(out.snp.samesign$hawearlylate,out.snp.samesign$applehawave))
pheatmap(test,col = my_palette,treeheight_row=0,treeheight_col=0)

out$hawearlylate<- out$hawearly_MAF- out$hawlate_MAF
out$applehawave<-out$appleave_MAF- out$hawave_MAF

out$sign<-sign(out$hawearlylate)==sign(out$applehawave)



"SELECT annotation.snpId, annotation.loc, urbanapoolLDx.MLE_AppleEarlyAppleLate_urbana, urbanapoolLDx.MLE_HawEarlyHawLate_urbana,RAD_linkage_matchPool.LDgr, urbanapoolTajDgene.Urbana_apple_num_gene_TajD, urbanapoolTajDgene.Urbana_haw_num_gene_TajD, snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust, snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust, snpFisherurbana.urbana_appleearly_applelate_fisher_pvalue_adjust, urbanapoolmaf.urbana_appleave_Maj, urbanapoolmaf.urbana_hawave_Maj, urbanapoolmaf.urbana_hawearly_Maj, urbanapoolmaf.urbana_hawlate_Maj, urbanapoolmaf.urbana_appleearly_Maj,urbanapoolmaf.urbana_applelate_Maj WHERE snpFisherurbana.urbana_appleave_hawave_fisher_pvalue_adjust < 0.05 AND snpFisherurbana.urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05 AND ABS(urbanapoolmaf.urbana_appleave_Maj - urbanapoolmaf.urbana_hawave_Maj) > 0.5 AND ABS (urbanapoolmaf.urbana_hawearly_Maj - urbanapoolmaf.urbana_hawlate_Maj) > 0.5;"



library(RMySQL)
con <- dbConnect(MySQL(),user="raglandlab", password="pomonella",dbname="PomUrbanaPoolseqSNP", host="localhost")

sql <- sprintf("SELECT feature_alias.gene_id, annotation.loc, annotation.effect,snpFisher.appleave_hawave_fisher_pvalue FROM feature_alias, annotation, snpFisher WHERE gene_id='%s' AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisher.snpId AND snpFisher.appleave_hawave_fisher_pvalue < 0.05;",gene)
query <- dbGetQuery(con, sql)

#to get those with a MAF difference of >0.55
sql<-("SELECT poolmaf2.snpId FROM poolmaf2 WHERE ABS(poolmaf2.appleearly_MAF - poolmaf2.applelate_MAF) > 0.55 AND ABS(poolmaf2.hawearly_MAF - poolmaf2.hawlate_MAF) > 0.55")
query <- dbGetQuery(con, sql)

dbGetQuery(con, "select * from poolmaf2 where snpId=251970;")

#create a table and figure
sql<-("SELECT * FROM poolmaf2 WHERE ABS(poolmaf2.appleearly_MAF - poolmaf2.applelate_MAF) > 0.55 AND ABS(poolmaf2.hawearly_MAF - poolmaf2.hawlate_MAF) > 0.55")
query <- dbGetQuery(con, sql)
