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
table_options<-c('scaffold','position','ref','alt','appleave_MAF','appleearly_MAF','applelate_MAF','hawave_MAF','hawearly_MAF','hawlate_MAF','protein_changes','impact'
                 ,'appleave_hawave_fisher_score','appleave_hawave_fisher_pvalue','appleearly_applelate_fisher_score','appleearly_applelate_fisher_pvalue','hawearly_hawlate_fisher_score',
                 'hawearly_hawlate_fisher_pvalue',"LDxApple","LDxHaw")
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
                   radioButtons('population', label = 'population choice for snps',
                              choices = c("urbana","grant"),selected='urbana'),
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
 #                  sliderInput("integerLDX", "Min LDx value SNP Apple Early Late:",
#                               min = 0, max = 1,
#                               value = 0.05),
#                   fuildRow(
#                   sliderInput("integerLDXhaw", "Min LDx value SNP Haw Early Late:",
#                               min = 0, max = 1,
#                               value = 0.05),
              pickerInput(inputId = "tableoptions", 
                      label = "Select/deselect table output", 
                      choices = table_options, options = list(`max-options` = 10), 
                      multiple = TRUE),
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
                   radioButtons('population2', label = 'population choice for snps',
                               choices = c("urbana","grant"),selected='urbana'),
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
               #    sliderInput("integerLDX2", "Min LDx value SNP Apple Early Late:",
                #               min = 0, max = 1,
                #               value = 0.05),
                #   sliderInput("integerLDX2haw", "Min LDx value SNP Haw Early Late:",
                #               min = 0, max = 1,
                #               value = 0.05),
               pickerInput(inputId = "tableoptions2", 
                           label = "Select/deselect table output", 
                           choices = table_options, options = list(`max-options` = 10), 
                           multiple = TRUE),
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
#Im making snpId, effect and ave pvalue compulsary collumns
#snpId because then I can join to LDx
#pvalue because of the slider
#effect because of annotation drop down memu
#not sure what to do with ldx not every snp has LDx value might have to do a seperate table for ldx
  output$SNP_table<-renderTable({
    if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes") {
      final_table1<-'Select single gene'
    }
    else {
      strsql<-input$tableoptions
      final_table<-'This is single gene'
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
      strsql<-gsub('appleearly_applelate_fisher_score', paste("snpFisher",input$population,'.',input$population,'_appleearly_applelate_fisher_score',sep=''),strsql)
      strsql<-gsub('appleearly_applelate_fisher_pvalue', paste("snpFisher",input$population,'.',input$population,'_appleearly_applelate_fisher_pvalue',sep=''),strsql)
      strsql<-gsub('hawearly_hawlate_fisher_score', paste("snpFisher",input$population,'.',input$population,'_hawearly_hawlate_fisher_score',sep=''),strsql)
      strsql<-gsub('hawearly_hawlate_fisher_pvalue', paste("snpFisher",input$population,'.',input$population,'_hawearly_hawlate_fisher_pvalue',sep=''),strsql)
      strsql<-c('SELECT annotation.snpId',strsql,paste("snpFisher",input$population,'.',input$population,"_appleave_hawave_fisher_pvalue",sep=''), "annotation.effect FROM annotation", paste('snpFisher',input$population,sep=''), "feature_alias", paste(input$population,"poolmaf WHERE feature_alias.gene_id='%s' AND feature_alias.loc = annotation.loc AND annotation.snpId=snpFisher", input$population,'.snpId AND ','annotation.snpId=',input$population,'poolmaf.snpId AND snpFisher',input$population,'.',input$population,"_appleave_hawave_fisher_pvalue <  '%f';",sep=''))
      strsql<-paste(strsql,collapse=', ')
      transcript<-strsplit(input$plot_var," +")[[1]][1]
      sql <- sprintf(strsql,transcript,input$integerFisher)
      query <- dbGetQuery(con, sql)
      final_table<-query
      
      #because Im having trouble doing the join in mysql running the LDx seperatly
      if ('LDxApple' %in% input$tableoptions | 'LDxHaw' %in% input$tableoptions){
        strsqlb<-input$tableoptions
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
        strsqlb<-strsqlb[strsqlb != 'appleearly_applelate_fisher_score']
        strsqlb<-strsqlb[strsqlb != 'appleearly_applelate_fisher_pvalue']
        strsqlb<-strsqlb[strsqlb != 'hawearly_hawlate_fisher_score']
        strsqlb<-strsqlb[strsqlb != 'hawearly_hawlate_fisher_pvalue']
        strsqlb<-c('SELECT annotation.snpId',strsqlb,paste("FROM annotation, ",input$population,"poolLDx, feature_alias WHERE feature_alias.gene_id='%s' AND feature_alias.loc = annotation.loc AND annotation.snpId=",input$population,"poolLDx.snpId;",sep=''))
         strsqlb<-paste(strsqlb,collapse=', ')
        strsqlb<-gsub(', FROM ',' FROM ',strsqlb)
        transcript<-strsplit(input$plot_var," +")[[1]][1]
        sqlb <- sprintf(strsqlb,transcript)
        query2 <- dbGetQuery(con, sqlb)
        final_table1<-left_join(final_table,query2,by='snpId') %>% distinct(.)  %>% filter(effect %in% input$annotation) %>% distinct(., snpId,effect, .keep_all = TRUE) #must filter this way to make this work distinct for multiple missense variations. where there is isoforms the missesnse variation will be slighly different number and now come up under the first distinct
        #remove the last distinct if you want this information.
      }
      else{
        final_table1<-final_table %>% distinct(.)  %>% filter(effect %in% input$annotation) %>% distinct(., snpId,effect, .keep_all = TRUE)
      }  
      
    }

    final_table1
  })
  
  output$SNP_table2<-renderTable({
#    if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes") {
#      final_table<-'Select single gene'
#      }
#    else {
#    transcript<-strsplit(input$plot_var2," +")[[1]][1]
#    sqlb <- sprintf("SELECT annotation.snpId, snpFisher.scaffold, snpFisher.position, annotation.effect,annotation.protein_changes,snpFisher.appleave_hawave_fisher_pvalue,snpFisher.appleearly_applelate_fisher_pvalue,snpFisher.hawearly_hawlate_fisher_pvalue FROM  annotation, snpFisher, feature_alias WHERE feature_alias.gene_id='%s' AND feature_alias.loc = annotation.loc AND annotation.snpId=snpFisher.snpId AND snpFisher.appleave_hawave_fisher_pvalue <  '%f';",transcript,input$integerFisher2)
 #   queryb <- dbGetQuery(con, sqlb)
#    sqlb <- sprintf("SELECT annotation.snpId, poolLDx.MLE_AppleEarlyAppleLate, poolLDx.MLE_HawEarlyHawLate FROM  annotation, poolLDx, feature_alias WHERE feature_alias.gene_id='%s' AND feature_alias.loc = annotation.loc AND annotation.snpId=poolLDx.snpId AND  (poolLDx.MLE_AppleEarlyAppleLate > '%f' OR poolLDx.MLE_HawEarlyHawLate > '%f');",transcript,input$integerLDX2,input$integerLDX2haw)
 #   queryb2 <- dbGetQuery(con, sqlb)
  #  final_table<-left_join(queryb,queryb2,by='snpId') %>% distinct(.)  %>% filter(effect %in% input$annotation)  %>% distinct(., snpId,effect, .keep_all = TRUE )#must filter this way to make this work
#    }
#    final_table
#  })
    if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes") {
      final_table2<-'Select single gene'
    }
    else {
      strsql2<-input$tableoptions2
      final_tableb<-'This is single gene'
      strsql2<-strsql2[strsql2 != "LDxApple"]
      strsql2<-strsql2[strsql2 != "LDxHaw"]
      strsql2<-gsub("scaffold", paste("snpFisher",input$population2,'.',"scaffold",sep=''), strsql2)
      strsql2<-gsub("position", paste("snpFisher",input$population2,'.',"position",sep=''), strsql2)
      strsql2<-gsub("ref", paste("snpFisher",input$population2,'.',"ref",sep=''), strsql2)
      strsql2<-gsub("alt", paste("snpFisher",input$population2,'.',"alt",sep=''), strsql2)
      strsql2<-gsub('appleave_MAF',paste(input$population2,'poolmaf.',input$population2,'_appleave_Min',sep=''), strsql2)
      strsql2<-gsub('appleearly_MAF',paste(input$population2,'poolmaf.',input$population2,'_appleearly_Min',sep=''), strsql2)
      strsql2<-gsub('applelate_MAF',paste(input$population2,'poolmaf.',input$population2,'_applelate_Min',sep=''), strsql2)
      strsql2<-gsub('hawave_MAF',paste(input$population2,'poolmaf.',input$population2,'_hawave_Min',sep=''), strsql2)
      strsql2<-gsub('hawearly_MAF',paste(input$population2,'poolmaf.',input$population2,'_hawearly_Min',sep=''), strsql2)
      strsql2<-gsub('hawlate_MAF',paste(input$population2,'poolmaf.',input$population2,'_hawlate_Min',sep=''), strsql2)
      strsql2<-gsub("protein_changes", "annotation.protein_changes", strsql2)
      strsql2<-gsub("impact", "annotation.impact",strsql2)
      strsql2<-gsub('appleave_hawave_fisher_score', paste("snpFisher",input$population2,'.',input$population2,'_appleave_hawave_fisher_score',sep=''),strsql2)
      strsql2<-gsub('appleearly_applelate_fisher_score', paste("snpFisher",input$population2,'.',input$population2,'_appleearly_applelate_fisher_score',sep=''),strsql2)
      strsql2<-gsub('appleearly_applelate_fisher_pvalue', paste("snpFisher",input$population2,'.',input$population2,'_appleearly_applelate_fisher_pvalue',sep=''),strsql2)
      strsql2<-gsub('hawearly_hawlate_fisher_score', paste("snpFisher",input$population2,'.',input$population2,'_hawearly_hawlate_fisher_score',sep=''),strsql2)
      strsql2<-gsub('hawearly_hawlate_fisher_pvalue', paste("snpFisher",input$population2,'.',input$population2,'_hawearly_hawlate_fisher_pvalue',sep=''),strsql2)
      strsql2<-c('SELECT annotation.snpId',strsql2,paste("snpFisher",input$population2,'.',input$population2,"_appleave_hawave_fisher_pvalue",sep=''), "annotation.effect FROM annotation", paste('snpFisher',input$population2,sep=''), "feature_alias", paste(input$population2,"poolmaf WHERE feature_alias.gene_id='%s' AND feature_alias.loc = annotation.loc AND annotation.snpId=snpFisher", input$population2,'.snpId AND ','annotation.snpId=',input$population2,'poolmaf.snpId AND snpFisher',input$population2,'.',input$population2,"_appleave_hawave_fisher_pvalue <  '%f';",sep=''))
      strsql2<-paste(strsql2,collapse=', ')
      transcript<-strsplit(input$plot_var2," +")[[1]][1]
      sql2 <- sprintf(strsql2,transcript,input$integerFisher2)
      queryb <- dbGetQuery(con, sql2)
      final_table2<-queryb
      
      #because Im having trouble doing the join in mysql running the LDx seperatly
      if ('LDxApple' %in% input$tableoptions2 | 'LDxHaw' %in% input$tableoptions2){
        strsql2b<-input$tableoptions2
        strsql2b<-gsub("LDxApple", paste(input$population2,"poolLDx.MLE_AppleEarlyAppleLate_",input$population2,sep=''), strsql2b)
        strsql2b<-gsub("LDxHaw", paste(input$population2,"poolLDx.MLE_HawEarlyHawLate_",input$population2,sep=''), strsql2b)
        strsql2b<-strsql2b[strsql2b != "LDxHaw"]
        strsql2b<-strsql2b[strsql2b != "scaffold"]
        strsql2b<-strsql2b[strsql2b != "position"]
        strsql2b<-strsql2b[strsql2b != "ref"]
        strsql2b<-strsql2b[strsql2b != "alt"]
        strsql2b<-strsql2b[strsql2b != 'appleave_MAF']
        strsql2b<-strsql2b[strsql2b != 'appleearly_MAF']
        strsql2b<-strsql2b[strsql2b != 'applelate_MAF']
        strsql2b<-strsql2b[strsql2b != 'hawave_MAF']
        strsql2b<-strsql2b[strsql2b != 'hawearly_MAF']
        strsql2b<-strsql2b[strsql2b != 'hawlate_MAF']
        strsql2b<-strsql2b[strsql2b != "protein_changes"]
        strsql2b<-strsql2b[strsql2b != "impact"]
        strsql2b<-strsql2b[strsql2b != 'appleave_hawave_fisher_score']
        strsql2b<-strsql2b[strsql2b != 'appleearly_applelate_fisher_score']
        strsql2b<-strsql2b[strsql2b != 'appleearly_applelate_fisher_pvalue']
        strsql2b<-strsql2b[strsql2b != 'hawearly_hawlate_fisher_score']
        strsql2b<-strsql2b[strsql2b != 'hawearly_hawlate_fisher_pvalue']
        strsql2b<-c('SELECT annotation.snpId',strsql2b,paste("FROM annotation, ",input$population2,"poolLDx, feature_alias WHERE feature_alias.gene_id='%s' AND feature_alias.loc = annotation.loc AND annotation.snpId=",input$population2,"poolLDx.snpId;",sep=''))
        strsql2b<-paste(strsql2b,collapse=', ')
        strsql2b<-gsub(', FROM ',' FROM ',strsql2b)
        transcript<-strsplit(input$plot_var2," +")[[1]][1]
        sql2b <- sprintf(strsql2b,transcript)
        query2b <- dbGetQuery(con, sql2b)
        final_table2<-left_join(final_table2,query2b,by='snpId') %>% distinct(.)  %>% filter(effect %in% input$annotation2) %>% distinct(., snpId,effect, .keep_all = TRUE) #must filter this way to make this work distinct for multiple missense variations. where there is isoforms the missesnse variation will be slighly different number and now come up under the first distinct
        #remove the last distinct if you want this information.
      }
      else{
        final_table2<-final_table2 %>% distinct(.)  %>% filter(effect %in% input$annotation2) %>% distinct(., snpId,effect, .keep_all = TRUE)
      }  
      
    }
    
    final_table2
  })
  

}

#run the application
shinyApp(ui=ui,server=server)

#click run app botton on rstudio



#going to try and switch to pool and dbplyr for connecting to the mysql server as apparently more efficient for shiny?
#Im going to switch to dbplyr but not really because I can already do this in mysql language and I cant be arsed reworking in dbplyr
#so this:
#con <- dbConnect(MySQL(),user="raglandlab", password="pomonella",dbname="PomUrbanaGrant", host="localhost")
#query <- dbGetQuery(con,"SELECT feature_alias.gene_id, snpFisherurbana.ref, snpFisherurbana.alt, annotation.loc, annotation.effect,snpFisherurbana.urbana_appleave_hawave_fisher_pvalue FROM feature_alias, annotation, snpFisherurbana WHERE (char_length(snpFisherurbana.ref) > 1 OR char_length(snpFisherurbana.alt) > 1) AND gene_id='gene10083' AND feature_alias.loc = annotation.loc AND annotation.snpId = snpFisherurbana.snpId AND snpFisherurbana.urbana_appleave_hawave_fisher_pvalue < 0.05;")

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





