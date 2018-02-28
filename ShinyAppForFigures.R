#trying to build poolseq mysql stuff into a shiny app for figures to make it more interactive and such
#firstly I've just started with just the RNAseq to try and make it simplier so I can get the syntax for shiny down first


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

#to start with lets just work on the RNAtable
#WithinMonthsDE<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/PomModsDEwithinMonth.csv",header=T,row.names=1,sep=",",stringsAsFactors = F)
#my connection through uni and the vpn is super slow so Im going to start with running on my machine so I can work out shiny 
WithinMonthsDE<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/PomModsDEwithinMonth.csv",header=T,row.names=1,sep=",",stringsAsFactors = F)
WithinMonthsDE.haw<-data.frame(WithinMonthsDE[1:6])
WithinMonthsDE.apple<-data.frame(WithinMonthsDE[c(1,2,7:11)])
WithinMonthsDE.haw$month2vs2_haw_logFC.x<-0
WithinMonthsDE.haw<-WithinMonthsDE.haw[,c(1,2,7,3,4,5,6)]
WithinMonthsDE.apple<-WithinMonthsDE.apple[,c(1,2,3,4,5,6,7)]
colnames(WithinMonthsDE.haw) <- c("gene_id","flybase", "2M","3M","4M","5M","6M")
colnames(WithinMonthsDE.apple) <- c("gene_id", "flybase", "2M","3M","4M","5M","6M")
WithinMonthsDE.haw.melt <- melt(WithinMonthsDE.haw,  c("gene_id","flybase"))
WithinMonthsDE.apple.melt <- melt(WithinMonthsDE.apple, c("gene_id","flybase"))
WithinMonthsDE.haw.melt$pop<-"haw"
WithinMonthsDE.apple.melt$pop<-"apple"
#together<-rbind(WithinMonthsDE.haw.melt,WithinMonthsDE.apple.melt)
together<-rbind(WithinMonthsDE.apple.melt,WithinMonthsDE.haw.melt)

#Across months:
#WithinMonthsDE<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/PomModsDEacrossMonth.csv",header=T,row.names=1,sep=",",stringsAsFactors = F)
#my connection through uni and the vpn is super slow so Im going to start with running on my machine so I can work out shiny 
AcrossMonthsDE<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/PomModsDEacrossMonth.csv",header=T,row.names=1,sep=",",stringsAsFactors = F)
colnames(AcrossMonthsDE)
AcrossMonthsDE.haw<-data.frame(AcrossMonthsDE[c(1,12,2:6)])
AcrossMonthsDE.apple<-data.frame(AcrossMonthsDE[c(1,12,2,7:11)])
AcrossMonthsDE.haw$month2vs2_haw_logFC.x<-0
AcrossMonthsDE.haw<-AcrossMonthsDE.haw[,c(1,2,3,8,4,5,6,7)]
AcrossMonthsDE.apple<-AcrossMonthsDE.apple[,c(1,2,3,4,5,6,7,8)]
colnames(AcrossMonthsDE.haw) <- c("gene_id","module","flybase", "2M","3M","4M","5M","6M")
colnames(AcrossMonthsDE.apple) <- c("gene_id", "module","flybase", "2M","3M","4M","5M","6M")
AcrossMonthsDE.haw.melt <- melt(AcrossMonthsDE.haw,  c("gene_id","flybase","module"))
AcrossMonthsDE.apple.melt <- melt(AcrossMonthsDE.apple, c("gene_id","flybase","module"))
AcrossMonthsDE.haw.melt$pop<-"haw"
AcrossMonthsDE.apple.melt$pop<-"apple"
#together<-rbind(WithinMonthsDE.haw.melt,WithinMonthsDE.apple.melt)
together.across<-rbind(AcrossMonthsDE.apple.melt,AcrossMonthsDE.haw.melt)


#pathways
wnt<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/wntSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
tor<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/torSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
insulin<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/insulinSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)


together.across %>% distinct(., module)
#https://stackoverflow.com/questions/48565661/dynamically-adding-and-removing-objects-in-response-to-selectinput-in-shiny

#thinking about:
#put a slider for pvalue
#put a table in the modules side that has enriched pathways
#do network plots and then change highlighting?? with radio buttons??
#pathway package in R kegg pathways uses flybase IDs I think


#link to mysql
#slider for ldx, fisher values type of snp "protein change" "upstream" etc

ui<-fluidPage(
  conditionalPanel(condition="input.conditionedPanels==1",
                   radioButtons('across', label = 'pathway_across',
                                choices = c("all","wnt","insulin","tor"),selected='all'),
                   uiOutput(outputId = "gene")),
  conditionalPanel(condition="input.conditionedPanels==2",
                   radioButtons('across_modules', label = 'Modules',
                                choices = c("1","2","3","4","5","6","7"),selected="1"),
                   uiOutput(outputId = "gene_module"))
  
,
  mainPanel(
        h4("Significant Across Months"),
        tabsetPanel(
          tabPanel(
            h4("Pathways"),
            plotOutput(outputId = "main_plot") ,
            #define the value of this tab
            value=1 #value relates to 'id=' in tabsetpanel. Except to get a value and relates that to'id' in tabsetPanel
          ),
          #add a tab for the boxplot
          tabPanel(
            h4("Modules"),
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
    df<-together.across
    #test<-'wnt'
    if (input$across=='wnt'){
      choicesare<-df %>% filter(., flybase %in% wnt$flyid) %>% distinct(.,gene_id)}
    if(input$across=='tor'){
    choicesare<-df %>% filter(., flybase %in% tor$flyid) %>% distinct(.,gene_id)}
    if(input$across=='insulin'){
    choicesare<-df %>% filter(., flybase %in% insulin$flyid) %>% distinct(.,gene_id)}
    if(input$across=='all'){
    choicesare<-df %>% distinct(.,gene_id)}
 
       selectizeInput(inputId = "plot_var",
                   label= "variable to plot",
 #                  choices=together.across %>% distinct(.,gene_id), 
                    choices=choicesare,
                   selected='gene10053',
                   options=list(maxOptions=3000))
    
  })
  output$gene_module <- renderUI({
    df<-together.across
    
    #test<-'wnt'
    selectizeInput(inputId = "plot_var2",
                   label= "variable to plot",
                   #                  choices=together.across %>% distinct(.,gene_id), 
                   choices=together.across %>% filter(., module %in% input$across_modules) %>% distinct(.,gene_id),
                   options=list(maxOptions=3000))
  })

  

  output$main_plot <-renderPlot({
    #specify the plot that will be generated
    together.across.test<-together.across %>% filter(.,gene_id==input$plot_var)
    title<-together.across.test[1,2]
    ggplot(together.across.test, aes(variable, value,group = interaction(gene_id,pop) ,colour=pop)) +
      geom_line() +
      ggtitle("Across Months",title) +
      scale_x_discrete("Time series from 2M",labels=c('2M','3M','4M','5M', '6M')) + #discrete
      scale_y_continuous("Log Fold Change",limits = c(-3.5, 3.5)) + #continuous
      theme_bw()
  })
  
  output$second_plot <-renderPlot({
    #specify the plot that will be generated
    together.across.test<-together.across %>% filter(.,gene_id==input$plot_var2)
    title<-together.across.test[1,2]
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




