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

#to start with lets just work on the RNAtable
#apple<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/prelimResults/Table.apple.JuneEJD",header=T,row.names=1,sep="\t",stringsAsFactors = F)
#haw<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/prelimResults/Table.haw.JuneEJD",header=T,row.names=1,sep="\t",stringsAsFactors = F)
#my connection through uni and the vpn is super slow so Im going to start with running on my machine so I can work out shiny 
apple<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/Table.apple.JuneEJD",header=T,row.names=1,sep="\t",stringsAsFactors = F)
colnames(apple)
WithinMonthsDE.apple<-data.frame(apple[c(1:5,44)])
haw<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/Table.haw.JuneEJD",header=T,row.names=1,sep="\t",stringsAsFactors = F)
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
wnt<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/wntSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
tor<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/torSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
insulin<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/insulinSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)


#wnt<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/wntSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
#tor<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/torSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
#insulin<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/insulinSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)

#https://stackoverflow.com/questions/48565661/dynamically-adding-and-removing-objects-in-response-to-selectinput-in-shiny


flybase_symbol_ID<-read.csv("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/fbgn_annotation_ID_fb_2017_06.tsv",header=T,row.names=NULL,sep="\t",stringsAsFactors = F,strip.white=T)
#flybase_symbol_ID<-read.csv("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/fbgn_annotation_ID_fb_2017_06.tsv",header=T,row.names=NULL,sep="\t",stringsAsFactors = F,strip.white=T)


#link to mysql
#slider for ldx, fisher values type of snp "protein change" "upstream" etc

#thinking about:
#do network plots and then change highlighting?? with radio buttons??
#pathway package in R kegg pathways uses flybase IDs I think

#link to mysql
#slider for ldx, fisher values type of snp "protein change" "upstream" etc
#could do graph of gene projectory and then along side a table of SNPS?


ui<-fluidPage(
  conditionalPanel(condition="input.conditionedPanels==1",
                   sliderInput("integer", "FDR:",
                               min = 0, max = 1,
                               value = 0.05),
                   radioButtons('within', label = 'pathway_within',
                                choices = c("all","wnt","insulin","tor"),selected='all'),
                   uiOutput(outputId = "gene")),
  conditionalPanel(condition="input.conditionedPanels==2",
                   sliderInput("integer2", "FDR:",
                               min = 0, max = 1,
                               value = 0.05),
                   radioButtons('across', label = 'pathway_across',
                                choices = c("all","wnt","insulin","tor"),selected='all'),
                   uiOutput(outputId = "gene2")),
  mainPanel(
    h4("Significant Across Months"),
    tabsetPanel(
      tabPanel(
        h4("Pathways DE within months FDR"),
        plotOutput(outputId = "main_plot") ,
        #define the value of this tab
        value=1 #value relates to 'id=' in tabsetpanel. Except to get a value and relates that to'id' in tabsetPanel
      ),
      #add a tab for the boxplot
      tabPanel(
        h4("Pathways DE across months FDR"),
        plotOutput(outputId = "second_plot") ,
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
      title<-"All_Genes"}
    else if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$within=="insulin") {
      together.test<-together%>% filter(., flybase %in% insulin$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"}
    else if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$within=="wnt") {
      together.test<-together%>% filter(., flybase %in% wnt$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"}
    else if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$within=="tor") {
      together.test<-together%>% filter(., flybase %in% tor$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"}
    else {
      together.test<-together %>% filter(.,gene_id==strsplit(input$plot_var," +")[[1]][1])
    title<-paste(together.test[1,2],(flybase_symbol_ID %>% filter(.,flybase==together.test[1,2]))[1,1],": FDR",format(round(unique(together.test$FDR),3),nsmall=3),collapse=" ")}
    #specify the plot that will be generated
    ggplot(together.test, aes(variable, value,group = interaction(gene_id,pop) ,colour=pop)) +
      geom_line() +
      ggtitle("Within Months",title) +
      scale_x_discrete("Time series from 2M",labels=c('2M','3M','4M','5M', '6M')) + #discrete
      scale_y_continuous("Log Fold Change",limits = c(-3.5, 3.5)) + #continuous
      theme_bw()
  })
  
  output$second_plot <-renderPlot({
    if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes" && input$across=='all') {
      together.across.test<-together.across %>% filter(.,FDR < input$integer)
      title<-"All_Genes"}
    else if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes" && input$across=="insulin") {
      together.across.test<-together.across%>% filter(., flybase %in% insulin$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"}
    else if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes" && input$across=="wnt") {
      together.across.test<-together.across%>% filter(., flybase %in% wnt$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"}
    else if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes" && input$across=="tor") {
      together.across.test<-together.across%>% filter(., flybase %in% tor$flyid)%>% filter(.,FDR < input$integer2)
      title<-"All_Genes"}
    else {
      together.across.test<-together.across %>% filter(.,gene_id==strsplit(input$plot_var2," +")[[1]][1])
      title<-paste(c(together.across.test[1,2],(flybase_symbol_ID %>% filter(.,flybase==together.across.test[1,2]))[1,1],": FDR apple and haw", format(round(unique(together.across.test$FDR),3),nsmall=3)),collapse=" ")}
    
    #specify the plot that will be generated
#    together.across.test<-together.across %>% filter(.,gene_id==strsplit(input$plot_var2," +")[[1]][1])
#    title<-paste(c(together.across.test[1,2],(flybase_symbol_ID %>% filter(.,flybase==together.across.test[1,2]))[1,1],": FDR apple and haw", format(round(unique(together.across.test$FDR),3),nsmall=3)),collapse=" ")
    ggplot(together.across.test, aes(variable, value,group = interaction(gene_id,pop) ,colour=pop)) +
      geom_line() +
      ggtitle("Across Months",title) +
      scale_x_discrete("Time series from 2M",labels=c('2M','3M','4M','5M', '6M')) + #discrete
      scale_y_continuous("Log Fold Change",limits = c(-3.5, 3.5)) + #continuous
      theme_bw()
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
