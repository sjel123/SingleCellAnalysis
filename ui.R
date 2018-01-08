
if(Sys.info()["nodename"] == "amre1al709")
  .libPaths("/vol00/shiny-server/ADME/RLib")

library(shiny)
library(DT) 
library(sqldf)
library(plotly)
library(Biobase)
library(knitr)
library(shinycssloaders)
# library(rmarkdown)
#source("Factors.R")

#Load Data
################################################
#Load Data.  This is run only once when the server.R is called
#df3 <<- read.csv("Data/data.csv", nrows=1)
Metadata <- read.table("/home/sjelinsk/projects/AMP_Lupus/Data/metadata.csv", 
                       header=TRUE, sep=",")
##End Load Data     
options(spinner.color.background="#F5F5F5")
##################################################
ui <- fluidPage(
  mainPanel(  
    #img(src="532306923_640.jpg", align="right",height = 72, width = 128),
    fluidRow(
      img(src = "CSI_3.jpg", align="right", height = 72, width = 336))
    
      #img(src = "Picture1.png", align="left", height = 72, width=336))
    ),#'End Main Panle'
    tabsetPanel(  
      tabPanel("Graphs",
        selectizeInput("e1", "Gene1", choices=NULL),
        selectizeInput("e2", "Gene2", choices=NULL),
        
        titlePanel("Expression of AMP Phase I Data Lupus Kidney/Skin"),
        #fluidRow(checkboxGroupInput('show_vars1', 'Columns in Table to show:', choices = NULL, selected="genename", inline=TRUE)),
         fluidRow(radioButtons("action_selectiontype", "Disease",
                       choices = c("All", "Kidney", "Skin"),
                       selected = "All", inline = TRUE)),
        
         fluidRow(radioButtons("action_selectiontype1", "Cell",
                              choices = c("All", unique(as.character(Metadata$SUBJECT_ID))),
                                  selected = "All", inline = TRUE)),

    fluidRow(
      column(6, #h3(verbatimTextOutput("Gene1")),
             withSpinner(plotOutput("main_plot",  height = "300px"))),
      column(6, #h3(verbatimTextOutput("Gene2")),
             withSpinner(plotOutput("main_plot2",  height = "300px"),type=4))
               ),#End Fluid Row
    fluidRow(
      column(6, h3("Expression AMP Phase I Data"),
             withSpinner(plotOutput("main_plot3",  height = "300px"), type=5, size=2)),   
      column(6, h3("Expression AMP Phase I Data"),
             withSpinner(plotOutput("main_plot4",  height = "300px"),type=6))
      )#End Fluid Row  
    
    ) # end tabPanel
    

    )#End Tabset Panel
  
  #uiOutput("DataTable"),
	

)#End UI






