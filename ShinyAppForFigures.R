#trying to build poolseq mysql stuff into a shiny app for figures to make it more interactive and such
#firstly I've just started with just the RNAseq to try and make it simplier so I can get the syntax for shiny down first

#source("https://bioconductor.org/biocLite.R")
#biocLite("RDAVIDWebService")
#

#You will need to change the directories for the files that I am using
#I have them set to my computer as for some reason at the uni I am at the VPN is horible just flick the comments over
#other than that it should run

#this one loads very slowly especially for the module side so just give it a minute to sort its stuff out.

####THESE MUST BE LOADED SLOWLY IN ORDER!!!####
####Otherwise shit goes wrong for some reason####
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
library(RDAVIDWebService)
#this uses my profile for DAVID feel free to build your own
#install.packages('BACA')
#library(BACA)

#to start with lets just work on the RNAtable
#WithinMonthsDE<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/PomModsDEwithinMonth.csv",header=T,row.names=1,sep=",",stringsAsFactors = F)
#my connection through uni and the vpn is super slow so Im going to start with running on my machine so I can work out shiny 
WithinMonthsDE<-read.table("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/PomModsDEwithinMonth.csv",header=T,row.names=1,sep=",",stringsAsFactors = F)
WithinMonthsDE.haw<-data.frame(WithinMonthsDE[c(1:6,12)])
WithinMonthsDE.apple<-data.frame(WithinMonthsDE[c(1,2,7:12)])
WithinMonthsDE.haw$month2vs2_haw_logFC.x<-0
WithinMonthsDE.haw<-WithinMonthsDE.haw[,c(1,7,2,8,3,4,5,6)]
WithinMonthsDE.apple<-WithinMonthsDE.apple[,c(1,8,2,3,4,5,6,7)]
colnames(WithinMonthsDE.haw) <- c("gene_id","module","flybase", "2M","3M","4M","5M","6M")
colnames(WithinMonthsDE.apple) <- c("gene_id","module","flybase", "2M","3M","4M","5M","6M")
WithinMonthsDE.haw.melt <- melt(WithinMonthsDE.haw,  c("gene_id","flybase","module"))
WithinMonthsDE.apple.melt <- melt(WithinMonthsDE.apple, c("gene_id","flybase","module"))
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

#wnt<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/wntSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
#tor<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/torSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)
#insulin<-read.table("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/insulinSignalingFlybase.txt",header=T,row.names=NULL,sep="\t",stringsAsFactors = F)

#https://stackoverflow.com/questions/48565661/dynamically-adding-and-removing-objects-in-response-to-selectinput-in-shiny


flybase_symbol_ID<-read.csv("/Users/edwinadowle/Documents/Cerasi/Pomonella/RSEMresults/AppleBackToHaw2MremoveBadHaw4M/shiny/fbgn_annotation_ID_fb_2017_06.tsv",header=T,row.names=NULL,sep="\t",stringsAsFactors = F,strip.white=T)
#flybase_symbol_ID<-read.csv("/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/RNAseqToMysql/fbgn_annotation_ID_fb_2017_06.tsv",header=T,row.names=NULL,sep="\t",stringsAsFactors = F,strip.white=T)

flybase_symbol_ID<-flybase_symbol_ID[c(1,3)]



ui<-fluidPage(
  conditionalPanel(condition="input.conditionedPanels==1",
                   radioButtons('choose_across_between', label="Dataset_choice",
                                choices=c('Across','Between'),selected='Across'),
                   radioButtons('across', label = 'pathway',
                                choices = c("all","wnt","insulin","tor"),selected='all'),
                   uiOutput(outputId = "gene")),
  conditionalPanel(condition="input.conditionedPanels==2",
                   radioButtons('choose_across_between2', label="Dataset_choice",
                                choices=c('Across','Between'),selected='Across'),
                   uiOutput(outputId = "module_selection"),
                   uiOutput(outputId = "gene_module"),
                   uiOutput(outputId = "gene_module_enrichment"))
  
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
            tableOutput("David_table"),
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
  FuncAnnotClust <- reactive({
    if(input$choose_across_between2=='Across'){
      df6<-together.across
      df6<-df6 %>% filter(., module %in% input$across_modules) %>% na.omit(object, cols=flybase)
      test1<-DAVIDWebService$new(email='eddy.dowle@otago.ac.nz',url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
      FG <- addList(test1, df6$flybase, idType="FLYBASE_GENE_ID", listName="isClass", listType="Gene")
      FuncAnnotClust1 <- getClusterReport(test1)
      FuncAnnotClust1
      return(FuncAnnotClust1)
    }
    if(input$choose_across_between2=='Between'){
      df6<-together
      df6<-together %>% filter(., module %in% input$across_modules) %>% na.omit(object, cols=flybase)
      test2<-DAVIDWebService$new(email='eddy.dowle@otago.ac.nz',url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
      FG <- addList(test2, df6$flybase, idType="FLYBASE_GENE_ID", listName="isClass", listType="Gene")
      FuncAnnotClust1 <- getClusterReport(test2)
#     FuncAnnotClust1
      return(FuncAnnotClust1)
   }
  #  FuncAnnotClust1
 # })
  })
  output$gene <- renderUI({
    if (input$choose_across_between=='Across'){
    df<-together.across}
    if(input$choose_across_between=='Between'){
      df<-together}
    #test<-'wnt'
    if (input$across=='wnt'){
      choicesare<-df %>% filter(., flybase %in% wnt$flyid) %>% distinct(.,gene_id,.keep_all = TRUE) %>% dplyr::select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% dplyr::select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase))}
    if(input$across=='tor'){
    choicesare<-df %>% filter(., flybase %in% tor$flyid) %>% distinct(.,gene_id,.keep_all = TRUE) %>% dplyr::select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% dplyr::select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase))}
    if(input$across=='insulin'){
    choicesare<-df %>% filter(., flybase %in% insulin$flyid) %>% distinct(.,gene_id,.keep_all = TRUE) %>% dplyr::select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% dplyr::select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase))}
    if(input$across=='all'){
    choicesare<-df %>% distinct(.,gene_id,.keep_all = TRUE) %>% dplyr::select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% dplyr::select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase))}
 
    choicesare<-rbind(paste("All_Genes","Overlayed"),choicesare)
       selectizeInput(inputId = "plot_var",
                   label= "Gene to plot",
 #                  choices=together.across %>% distinct(.,gene_id), 
                    choices=choicesare,
                   selected='gene10053',
                   options=list(maxOptions=3000))
    
  })
  
output$module_selection <-renderUI({
  if (input$choose_across_between2=='Between'){
    optionsare = c("1","2","3","4","5","6","7","8")
    }
  if(input$choose_across_between2=='Across'){
    optionsare = c("1","2","3","4","5","6","7")
    }
  radioButtons('across_modules', label = 'Modules',
               choices = optionsare,selected="1")
  
  
})

  output$gene_module <- renderUI({
    if(input$choose_across_between2=='Across'){
      df2<-together.across
      choicesgenes=df2 %>% filter(., module %in% input$across_modules) %>% distinct(.,gene_id,.keep_all = TRUE) %>% dplyr::select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% dplyr::select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase))
    }
    if(input$choose_across_between2=='Between'){
      df2<-together
      choicesgenes=df2 %>% filter(., module %in% input$across_modules) %>% distinct(.,gene_id,.keep_all = TRUE) %>% dplyr::select(.,gene_id,flybase) %>% left_join(.,flybase_symbol_ID,by="flybase") %>% dplyr::select(.,gene_id,gene_symbol,flybase) %>% transmute(.,choice=paste(gene_id,gene_symbol,flybase))
    }
    choicesgenes<-rbind(paste("All_Genes","Overlayed"),choicesgenes)
    
    selectizeInput(inputId = "plot_var2",
                   label= "variable to plot",
                   #                  choices=together.across %>% distinct(.,gene_id), 
                   choices=choicesgenes,
                   options=list(maxOptions=3000))
  })

  output$gene_module_enrichment <- renderUI({
 #   if(input$choose_across_between2=='Across'){
#      df6<-together.across
#      df6<-df6 %>% filter(., module %in% input$across_modules) %>% na.omit(object, cols=flybase)
#      test1<-DAVIDWebService$new(email='eddy.dowle@otago.ac.nz',url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
#      FG <- addList(test1, df6$flybase, idType="FLYBASE_GENE_ID", listName="isClass", listType="Gene")
 #     FuncAnnotClust <- getClusterReport(test1)
#      choicesenrich=1:nrow(summary(FuncAnnotClust))
 #   }
  #  if(input$choose_across_between2=='Between'){
   #   df6<-together
    #  df6<-together %>% filter(., module %in% input$across_modules) %>% na.omit(object, cols=flybase)
     # test2<-DAVIDWebService$new(email='eddy.dowle@otago.ac.nz',url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
    #  FG <- addList(test2, df6$flybase, idType="FLYBASE_GENE_ID", listName="isClass", listType="Gene")
    #  FuncAnnotClust <- getClusterReport(test2)
    #  choicesenrich=1:nrow(summary(FuncAnnotClust))
  #  }
    FuncAnnotClust<-FuncAnnotClust()
    choicesenrich=1:nrow(summary(FuncAnnotClust))
    selectizeInput(inputId = "table_var2",
                   label= "David Cluster",
                   #                  choices=together.across %>% distinct(.,gene_id), 
                   choices=choicesenrich,
                   options=list(maxOptions=3000))
  })
  
  output$David_table<-renderTable({
  #  if(input$choose_across_between2=='Across'){
  #   df7<-together.across %>% filter(., module %in% input$across_modules) %>% na.omit(object, cols=flybase)
  #    test3<-DAVIDWebService$new(email='eddy.dowle@otago.ac.nz',url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  #    FG2 <- addList(test3, df7$flybase, idType="FLYBASE_GENE_ID", listName="isClass", listType="Gene")
  #    FuncAnnotClust <- getClusterReport(test3)
   # }
  #  if(input$choose_across_between2=='Between'){
  #    df7<-together %>% filter(., module %in% input$across_modules) %>% na.omit(object, cols=flybase)
  #    test4<-DAVIDWebService$new(email='eddy.dowle@otago.ac.nz',url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  #   FG2 <- addList(test4, df7$flybase, idType="FLYBASE_GENE_ID", listName="isClass", listType="Gene")
  #    FuncAnnotClust <- getClusterReport(test4)
  #  }
clus<-as.integer(input$table_var2)
FuncAnnotClust<-FuncAnnotClust()
final_table<-members(FuncAnnotClust)[[clus]] %>% dplyr::select(-one_of("Genes"))
# final_table
  #  final_table<-head(together)
final_table
     })

  output$main_plot <-renderPlot({
    #specify the plot that will be generated
    if (input$choose_across_between=='Across'){
      df3<-together.across}
    if (input$choose_across_between=='Between')
      {df3<-together}
    if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$across=='all') {
      together.across.test<-df3 
      title<-"All_Genes"
      linesize<-1}
    else if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$across=="insulin") {
      together.across.test<-df3  %>% filter(.,flybase %in% insulin$flyid) 
      title<-"All_Genes"
      linesize<-1}
    else if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$across=='wnt' ) {
      together.across.test<-df3 %>% filter(.,flybase %in% wnt$flyid) 
      title<-"All_Genes"
      linesize<-1}
    else if(strsplit(input$plot_var," +")[[1]][1]=="All_Genes" && input$across=="tor" ) {
      together.across.test<-df3 %>% filter(.,flybase %in% tor$flyid) 
      title<-"All_Genes" 
      linesize<-1}
     else
     {together.across.test<-df3 %>% filter(.,gene_id==strsplit(input$plot_var," +")[[1]][1])
     title<- together.across.test[1,2]
     linesize<-3
     }
#    title<-together.across.test[1,2]
    ggplot(together.across.test, aes(variable, value,group = interaction(gene_id,pop) ,colour=pop)) +
      geom_line(size=linesize) +
      ggtitle("Across Months",title) +
      scale_x_discrete("Time series from 2M",labels=c('2M','3M','4M','5M', '6M')) + #discrete
      scale_y_continuous("Log Fold Change",limits = c(-4, 4)) + #continuous
      theme_bw(base_size=20) +
      scale_colour_manual(values=c('red','blue'))
  })
  
  output$second_plot <-renderPlot({
    #specify the plot that will be generated
    if (input$choose_across_between2=='Across'){
      df4<-together.across}
    if(input$choose_across_between2=='Between'){
      df4<-together}
    if(strsplit(input$plot_var2," +")[[1]][1]=="All_Genes")
    {
      together.across.test<-df4 %>% filter(.,module %in% input$across_modules)
      title<-"All_Genes"
      linesize<-1}
     else
     {together.across.test<-df4 %>% filter(.,gene_id==strsplit(input$plot_var2," +")[[1]][1])
     title<-together.across.test[1,2]
     linesize<-3}
    ggplot(together.across.test, aes(variable, value,group = interaction(gene_id,pop) ,colour=pop)) +
      geom_line(size=linesize) +
      ggtitle("Across Months",title) +
      scale_x_discrete("Time series from 2M",labels=c('2M','3M','4M','5M', '6M')) + #discrete
      scale_y_continuous("Log Fold Change") + #continuous
      coord_cartesian(ylim = c(-4, 4)) +
      theme_bw(base_size=20) +
      scale_colour_manual(values=c('red','blue'))
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
test<-'Across'
FuncAnnotClust3 <- 
    if(test=='Across'){
      return(FuncAnnotClust1)
    }
    if(input$choose_across_between2=='Between'){
      return(FuncAnnotClust1)
    }

