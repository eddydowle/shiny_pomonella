#trying to build poolseq mysql stuff into a shiny app for figures to make it more interactive and such
#firstly I've just started with just the RNAseq to try and make it simplier so I can get the syntax for shiny down first

#This is the App that Im going to use for the mysql poolseq data
#currently Im not there yet. Next on the list
#
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
con <- dbConnect(MySQL(),user="raglandlab", password="pomonella",dbname="PomUrbanaPoolseqSNP", host="localhost")


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

snpeff_effects<-c('3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant', '5_prime_UTR_variant', 'downstream_gene_variant', 'initiator_codon_variant', 'initiator_codon_variant&non_canonical_start_codon', 'initiator_codon_variant&splice_region_variant', 'intergenic_region', 'intragenic_variant', 'intron_variant', 'missense_variant', 'missense_variant&splice_region_variant', 'non_coding_transcript_exon_variant', 'non_coding_transcript_variant', 'splice_acceptor_variant&intron_variant', 'splice_acceptor_variant&splice_donor_variant&intron_variant', 'splice_donor_variant&intron_variant', 'splice_region_variant', 'splice_region_variant&initiator_codon_variant&non_canonical_start_codon', 'splice_region_variant&intron_variant', 'splice_region_variant&non_coding_transcript_exon_variant', 'splice_region_variant&stop_retained_variant', 'splice_region_variant&synonymous_variant', 'start_lost', 'start_lost&splice_region_variant', 'stop_gained', 'stop_gained&splice_region_variant', 'stop_lost', 'stop_lost&splice_region_variant', 'stop_retained_variant', 'synonymous_variant', 'upstream_gene_variant')
snpeff_effects
#link to mysql
#slider for ldx, fisher values type of snp "protein change" "upstream" etc

#thinking about:
#do network plots and then change highlighting?? with radio buttons??
#pathway package in R kegg pathways uses flybase IDs I think

#link to mysql
#slider for ldx, fisher values type of snp "protein change" "upstream" etc
#could do graph of gene projectory and then along side a table of SNPS?
#?column

ui<-fluidPage(
  conditionalPanel(condition="input.conditionedPanels==1",
                   fluidRow(column(3,
                   sliderInput("integer", "FDR:",
                               min = 0, max = 1,
                               value = 0.05)),
#                   fluidRow(
                    column(4, offset=1,
                   radioButtons('within', label = 'pathway_within',
                                choices = c("all","wnt","insulin","tor"),selected='all'),
 #                  fluidRow(
                   uiOutput(outputId = "gene")),
  #                 fuildRow(
                  column(4,
                   sliderInput("integerFisher", "Max Fisher value SNP Apple Ave Haw Ave:",
                               min = 0, max = 1,
                               value = 0.05),
 #                  fuildRow(
                   sliderInput("integerLDX", "Min LDx value SNP Apple Early Late:",
                               min = 0, max = 1,
                               value = 0.05),
#                   fuildRow(
                   sliderInput("integerLDXhaw", "Min LDx value SNP Haw Early Late:",
                               min = 0, max = 1,
                               value = 0.05),
                  pickerInput(inputId = "annotation", 
                                     label = "Select/deselect annotations", 
                                    choices = snpeff_effects, options = list(`max-options` = 10), #, options = list(`max-options` = 4)
                                   multiple = TRUE)))),
            #         )),),
                                      #,
# I think this will be more important for just looking at poolseq for now its both RNAseq and pool and that seems to limit the SNPS enough
                  
    #  column(4,
     #  pickerInput(inputId = "annotation", 
      #                         label = "Select/deselect annotations", 
       #                        choices = snpeff_effects, options = list(`actions-box` = TRUE), 
        #                       multiple = TRUE)
         #         )),
  conditionalPanel(condition="input.conditionedPanels==2",
                   fluidRow(column(3,
                   sliderInput("integer2", "FDR:",
                               min = 0, max = 1,
                               value = 0.05)),
                   column(4, offset=1,
                   radioButtons('across', label = 'pathway_across',
                                choices = c("all","wnt","insulin","tor"),selected='all'),
                   uiOutput(outputId = "gene2")),
                   column(4,
                   sliderInput("integerFisher2", "Max Fisher value SNP Apple Ave Haw Ave:",
                               min = 0, max = 1,
                               value = 0.05),
                   sliderInput("integerLDX2", "Min LDx value SNP Apple Early Late:",
                               min = 0, max = 1,
                               value = 0.05),
                   sliderInput("integerLDX2haw", "Min LDx value SNP Haw Early Late:",
                               min = 0, max = 1,
                               value = 0.05),
                   pickerInput(inputId = "annotation2", 
                               label = "Select/deselect annotations", 
                               choices = snpeff_effects, options = list(`max-options` = 10), 
                               multiple = TRUE)))),
  
  
  
  
  
  mainPanel(
    h4("Significant Across Months"),
    tabsetPanel(
      tabPanel(
        h4("Pathways DE within months FDR"),
        plotOutput(outputId = "main_plot") ,#,width = '800px', height = '300px'
        tableOutput("SNP_table"),
        #define the value of this tab
        value=1 #value relates to 'id=' in tabsetpanel. Except to get a value and relates that to'id' in tabsetPanel
      ),
      #add a tab for the boxplot
      tabPanel(
        h4("Pathways DE across months FDR"),
        plotOutput(outputId = "second_plot") ,
        tableOutput("SNP_table2"),
        #define the value of this tab (not necessary here-just example)
        value=2
      )
      #name the panel to correspond ot the condition defined in the sidebarpanel
      , id='conditionedPanels'
    )
  )
)

#create the function to be excuted by shiny
server<-function(input,output) {
  # output$menuitem<-renderMenu({
  
  #})
  output$gene <- renderUI({
    df<-together
    #test<-'wnt'
    if (input$within=='wnt'){
      choicesare<-df %>% filter(., flybase %in% wnt$flyid)%>% filter(.,FDR < input$integer) %>% distinct(.,gene_id,.keep_all = TRUE) %>% select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase)) }
    if(input$within=='tor'){
      choicesare<-df %>% filter(., flybase %in% tor$flyid)%>% filter(.,FDR < input$integer) %>% distinct(.,gene_id,.keep_all = TRUE) %>% select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase)) }
    if(input$within=='insulin'){
      choicesare<-df %>% filter(., flybase %in% insulin$flyid)%>% filter(.,FDR < input$integer) %>% distinct(.,gene_id,.keep_all = TRUE) %>% select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase)) }
    if(input$within=='all'){
      choicesare<-df %>% filter(.,FDR < input$integer) %>% distinct(.,gene_id,.keep_all = TRUE) %>% select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase)) }
    
    choicesare<-rbind(paste("All_Genes","Overlayed"),choicesare)
    
    selectizeInput(inputId = "plot_var",
                   label= "variable to plot",
                   #                  choices=together.across %>% distinct(.,gene_id), 
                   choices=choicesare,
                   selected='gene10053',
                   options=list(maxOptions=3000))
  })
  
  output$gene2 <- renderUI({
    df2<-together.across
    #test<-'wnt'
    #I think this should bring up things that are significant in one and both??
    if (input$across=='wnt'){
      choicesare2<-df2 %>% filter(., flybase %in% wnt$flyid)%>% filter(.,FDR < input$integer2) %>% distinct(.,gene_id,.keep_all = TRUE) %>% select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase))
    }
    if(input$across=='tor'){
      choicesare2<-df2 %>% filter(., flybase %in% tor$flyid)%>% filter(.,FDR < input$integer2) %>% distinct(.,gene_id,.keep_all = TRUE) %>% select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase)) }
    if(input$across=='insulin'){
      choicesare2<-df2 %>% filter(., flybase %in% insulin$flyid)%>% filter(.,FDR < input$integer2) %>% distinct(.,gene_id,.keep_all = TRUE) %>% select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase)) }
    if(input$across=='all'){
      choicesare2<-df2 %>% filter(.,FDR < input$integer2) %>% distinct(.,gene_id,.keep_all = TRUE) %>% select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase)) }
    
    choicesare2<-rbind(paste("All_Genes","Overlayed"),choicesare2)
    
    selectizeInput(inputId = "plot_var2",
                   label= "variable to plot",
                   #                  choices=together.across %>% distinct(.,gene_id), 
                   choices=choicesare2,
                   selected='gene10053',
                   options=list(maxOptions=3000))
  })
  
  
  output$main_plot <-renderPlot({
    if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$within=='all') {
      together.test<-together %>% filter(.,FDR < input$integer)
      title<-"All_Genes"
      linesize<-1}
    else if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$within=="insulin") {
      together.test<-together%>% filter(., flybase %in% insulin$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes" 
      linesize<-1}
    else if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$within=="wnt") {
      together.test<-together%>% filter(., flybase %in% wnt$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"
      linesize<-1}
    else if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$within=="tor") {
      together.test<-together%>% filter(., flybase %in% tor$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"
      linesize<-1}
    else {
      together.test<-together %>% filter(.,gene_id==strsplit(input$plot_var," +")[[1]][1])
      title<-paste(together.test[1,2],(flybase_symbol_ID %>% filter(.,flybase==together.test[1,2]))[1,1],": FDR",format(round(unique(together.test$FDR),3),nsmall=3),collapse=" ")
      linesize<-3}
    #specify the plot that will be generated
    ggplot(together.test, aes(variable, value,group = interaction(gene_id,pop) ,colour=pop)) +
      geom_line(size=linesize) +
      ggtitle("Within Months",title) +
      scale_x_discrete("Time series from 2M",labels=c('2M','3M','4M','5M', '6M')) + #discrete
      scale_y_continuous("Log Fold Change",limits = c(-3.5, 3.5)) + #continuous
      theme_bw(base_size=20) +
      scale_colour_manual(values=c('red','blue'))
  })
  
  output$second_plot <-renderPlot({
    if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes" && input$across=='all') {
      together.across.test<-together.across %>% filter(.,FDR < input$integer)
      title<-"All_Genes"
      linesize<-1}
    else if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes" && input$across=="insulin") {
      together.across.test<-together.across%>% filter(., flybase %in% insulin$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"
      linesize<-1}
    else if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes" && input$across=="wnt") {
      together.across.test<-together.across%>% filter(., flybase %in% wnt$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"
      linesize<-1}
    else if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes" && input$across=="tor") {
      together.across.test<-together.across%>% filter(., flybase %in% tor$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"
      linesize<-1}
    else {
      together.across.test<-together.across %>% filter(.,gene_id==strsplit(input$plot_var2," +")[[1]][1])
      title<-paste(c(together.across.test[1,2],(flybase_symbol_ID %>% filter(.,flybase==together.across.test[1,2]))[1,1],": FDR apple and haw", format(round(unique(together.across.test$FDR),3),nsmall=3)),collapse=" ")
      linesize<-3}
    
    #specify the plot that will be generated
    #    together.across.test<-together.across %>% filter(.,gene_id==strsplit(input$plot_var2," +")[[1]][1])
    #    title<-paste(c(together.across.test[1,2],(flybase_symbol_ID %>% filter(.,flybase==together.across.test[1,2]))[1,1],": FDR apple and haw", format(round(unique(together.across.test$FDR),3),nsmall=3)),collapse=" ")
    ggplot(together.across.test, aes(variable, value,group = interaction(gene_id,pop) ,colour=pop)) +
      geom_line(size=linesize) +
      ggtitle("Across Months",title) +
      scale_x_discrete("Time series from 2M",labels=c('2M','3M','4M','5M', '6M')) + #discrete
      scale_y_continuous("Log Fold Change",limits = c(-3.5, 3.5)) + #continuous
      theme_bw(base_size=20) +
      scale_colour_manual(values=c('red','blue'))
  })
#ok so we want a table of the snps  
  
#not sure what to do with ldx not every snp has LDx value might have to do a seperate table for ldx
  output$SNP_table<-renderTable({
    if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes") {
      final_table<-'Select single gene'
    }
    else {
      final_table<-'This is single gene'
      transcript<-strsplit(input$plot_var," +")[[1]][1]
#      final_table<-transcript
      sql <- sprintf("SELECT annotation2.snpId, snpFisher2.scaffold, snpFisher2.position,annotation2.effect,annotation2.protein_changes,snpFisher2.appleave_hawave_fisher_pvalue,snpFisher2.appleearly_applelate_fisher_pvalue,snpFisher2.hawearly_hawlate_fisher_pvalue FROM  annotation2, snpFisher2, feature_alias WHERE feature_alias.gene_id='%s' AND feature_alias.loc = annotation2.loc AND annotation2.snpId=snpFisher2.snpId AND  snpFisher2.appleave_hawave_fisher_pvalue < '%f';",transcript,input$integerFisher)
      query <- dbGetQuery(con, sql)
      sql <- sprintf("SELECT annotation2.snpId, poolLDx.MLE_AppleEarlyAppleLate, poolLDx.MLE_HawEarlyHawLate FROM  annotation2, poolLDx, feature_alias WHERE feature_alias.gene_id='%s' AND feature_alias.loc = annotation2.loc AND annotation2.snpId=poolLDx.snpId AND  (poolLDx.MLE_AppleEarlyAppleLate > '%f' OR poolLDx.MLE_HawEarlyHawLate > '%f');",transcript,input$integerLDX,input$integerLDXhaw)
      query2 <- dbGetQuery(con, sql)
      final_table<-left_join(query,query2,by='snpId') %>% distinct(.)  %>% filter(effect %in% input$annotation) %>% distinct(., snpId,effect, .keep_all = TRUE) #must filter this way to make this work distinct for multiple missense variations. where there is isoforms the missesnse variation will be slighly different number and now come up under the first distinct
      }
    # final_table
    #  final_table<-head(together)
    final_table
  })
  
  output$SNP_table2<-renderTable({
    if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes") {
      final_table<-'Select single gene'
      }
    else {
#    final_table<-'This is single gene'
    transcript<-strsplit(input$plot_var2," +")[[1]][1]
    #      final_table<-transcript
    sqlb <- sprintf("SELECT annotation2.snpId, snpFisher2.scaffold, snpFisher2.position, annotation2.effect,annotation2.protein_changes,snpFisher2.appleave_hawave_fisher_pvalue,snpFisher2.appleearly_applelate_fisher_pvalue,snpFisher2.hawearly_hawlate_fisher_pvalue FROM  annotation2, snpFisher2, feature_alias WHERE feature_alias.gene_id='%s' AND feature_alias.loc = annotation2.loc AND annotation2.snpId=snpFisher2.snpId AND snpFisher2.appleave_hawave_fisher_pvalue <  '%f';",transcript,input$integerFisher2)
    queryb <- dbGetQuery(con, sqlb)
    sqlb <- sprintf("SELECT annotation2.snpId, poolLDx.MLE_AppleEarlyAppleLate, poolLDx.MLE_HawEarlyHawLate FROM  annotation2, poolLDx, feature_alias WHERE feature_alias.gene_id='%s' AND feature_alias.loc = annotation2.loc AND annotation2.snpId=poolLDx.snpId AND  (poolLDx.MLE_AppleEarlyAppleLate > '%f' OR poolLDx.MLE_HawEarlyHawLate > '%f');",transcript,input$integerLDX2,input$integerLDX2haw)
    queryb2 <- dbGetQuery(con, sqlb)
    final_table<-left_join(queryb,queryb2,by='snpId') %>% distinct(.)  %>% filter(effect %in% input$annotation2)  %>% distinct(., snpId,effect, .keep_all = TRUE )#must filter this way to make this work
    }
        # final_table
    #  final_table<-head(together)
    final_table
  })
  

}

#run the application
shinyApp(ui=ui,server=server)

#click run app botton on rstudio



#going to try and switch to pool and dbplyr for connecting to the mysql server as apparently more efficient for shiny?
#Im going to switch to dbplyr but not really because I can already do this in mysql language and I cant be arsed reworking in dbplyr
#so this:
#con <- dbConnect(MySQL(),user="raglandlab", password="pomonella",dbname="PomUrbanaPoolseqSNP", host="localhost")
#query <- dbGetQuery(con,"SELECT feature_alias.gene_id, annotation.loc, annotation.effect,snpFisher.appleave_hawave_fisher_pvalue FROM feature_alias, annotation, snpFisher WHERE gene_id='gene10083' AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisher.snpId AND snpFisher.appleave_hawave_fisher_pvalue < 0.05;")

#becomes:
#my_db <- dbPool(
#  RMySQL::MySQL(), 
#  dbname = "PomUrbanaPoolseqSNP",
#  host = "localhost",
#  username = "raglandlab",
#  password = "pomonella"
#)

#x<-tbl(my_db, sql("SELECT feature_alias.gene_id, annotation.loc, annotation.effect,snpFisher.appleave_hawave_fisher_pvalue FROM feature_alias, annotation, snpFisher WHERE gene_id='gene10083' AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisher.snpId AND snpFisher.appleave_hawave_fisher_pvalue < 0.05")) %>% explain()
#i think to do it properly  in dbplyr you would have to do joins
#https://stackoverflow.com/questions/39864427/dplyr-sql-joins
#Im too lazy to relearn how to call mysql from R so Im sticking with the mysql syntax but I think it would be better for new people to learn the dbplyr syntax so that its less learning and fits with tidyverse
