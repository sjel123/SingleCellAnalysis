
if(Sys.info()["nodename"] == "amre1al709")
  .libPaths("/vol00/shiny-server/ADME/RLib")

library(shiny)
library(DT)   
library(lubridate)
##library(ggplot2)
library(plotly)
library(Biobase)
#setwd("X:\\R Scripts\\RA\\RAData\\ShinyApp\\")
#source("Factors.R")

################################################
#Load Data.  This is run only once when the server.R is called

#res <- read.table("Data/res.txt", header=TRUE, sep="\t")
mydata1 <- read.table("/home/sjelinsk/projects/AMP_Lupus/Data/ExprsData.txt", 
                      header=TRUE, sep="\t")
#   numZero <- data.frame(apply(mydata1,1,function(x) sum(x==0)))
#     numZero$Gene <- rownames(numZero)
#     numZero[order(numZero[,1]),][1:10,]

# z <- t(apply(mydata1[,],1,function(x) t(aggregate(as.numeric(x)~meta$TISSUE,FUN=mean))[2,]))
# z <- data.frame(z)
#   z$X1 <- as.numeric(as.character(z$X1))
#   z$X2 <- as.numeric(as.character(z$X2))
# z$diff <- z$X1/z$X2
# head(z[z$X2>.2 & z$diff<.5,])  


pca_rotation <- read.table("/home/sjelinsk/projects/AMP_Lupus/Data/TNSE.txt", 
                      header=TRUE, sep="\t")
  names(pca_rotation)[1:2]<- c("PC1","PC2")
Metadata <- read.table("/home/sjelinsk/projects/AMP_Lupus/Data/metadata.csv", 
                      header=TRUE, sep=",")

meta <- Metadata[Metadata$SAMPLE_ID %in% names(mydata1),]
  meta <- meta[match(row.names(pca_rotation),meta$SAMPLE_ID),]
table(meta$SAMPLE_ID==row.names(pca_rotation))
table(meta$SAMPLE_ID==colnames(mydata1))
table(colnames(mydata1)==row.names(pca_rotation))

res <- as.data.frame(rownames(mydata1))
names(res) <- "geneName"
res$Gene <- res$geneName
#res[,3:13] <- signif(res[,3:13],2)


##End Load Data  
##################################################

shinyServer <- function(input, output, session) {  
    
  Title <- reactive({
    s = as.numeric(rownames(res[RowValue(),Names()])[input$x1_rows_selected])
    if (is.null(s)) s=2  
    if (length(s)==0) s=2 
    res[s,2]
  })
  
  updateSelectizeInput(session, "e1", choices=res[,1], server=TRUE, selected = "GPX3")
  updateSelectizeInput(session, "e2", choices=res[,1], server=TRUE, selected = "KRT10")
  
   output$Gene1 <- renderText(input$e1)
   output$Gene2 <- renderText(input$e2)
  
  myGene1Inputs <- reactive({
    myData1 <- input$'e1'
  })
  myGene2Inputs <- reactive({
    myData2 <- input$'e2'
  })
    
  observe({
    choices2 <-  names(res)
    updateCheckboxGroupInput(session, "show_vars1",
                             label = "Select Columns to Include",
                             choices = choices2,
                             selected = choices2, 
                             inline = TRUE)
  })
  

  #Used to remove columns rows from data table() based on categories selected from show_vars1 checkbox on client
  Names <- reactive({
    NameLabels <- (input$'show_vars1')
    #print(sprintf("NameLables %s", NameLabels))
    NameLabels
  }) 
  
   ## used to select rows based on slider inputs
   RowValue <- reactive({
     Value <- (row.names(unique(res[,Names()])))
    Value
   })

  SW <- reactive({ 
    s = NULL
    s = which(res[,1] == myGene1Inputs() )
    print(sprintf("Names() %s", Names()))
    print(sprintf("input$x1_rows_selected %s", input$x1_rows_selected))
      if (is.null(s)) s=2  
      if (length(s)==0) s=2 
    print(sprintf("row selected %s", s))
    s
  })
  
  SW1 <- reactive({ 
    s = NULL
    s = which(res[,1] == myGene2Inputs() )
    #s = as.numeric(rownames(res[RowValue(),]))[input$x2_rows_selected]
    #print(sprintf("RowValue() %s", RowValue()))
      if (is.null(s)) s=2  
      if (length(s)==0) s=2 
    s
  })
  
    myDatapca <- reactive({
        myData1 <- input$'action_selectiontype'
        myData2 <- input$'action_selectiontype1'
          pca <- pca_rotation
            if(myData1 != "All"){
              pca <- pca_rotation[which(meta$TISSUE == myData1),]
             }
            if(myData2 != "All"){
              pca <- pca_rotation[which(meta$SUBJECT_ID == myData2),]
            }
            if(myData1 != "All"&myData2 != "All"){
              pca <- pca_rotation#[which(meta$TISSUE == myData1),]# & meta$SUBJECT_ID == myData2),]
            }
        if(myData1 == "All"){
          pca <- pca_rotation#[which(meta$TISSUE == myData1),]# & meta$SUBJECT_ID == myData2),]
        }
          print(dim(pca))
          print(myData1)
          pca
    })

    myDatacpm <- reactive({
      myData1 <- input$'action_selectiontype'
      myData2 <- input$'action_selectiontype1'
        log2cpm <- mydata1
          if(myData1 != "All"){
            cpm <- log2cpm[,which(meta$TISSUE == myData1)]
          }
          if(myData2 != "All"){
            cpm <- log2cpm[,which(meta$SUBJECT_ID == myData2)]
          }
          if(myData1 != "All"&myData2 != "All"){
            cpm <- log2cpm[,which(meta$TISSUE == myData1)]# & meta$SUBJECT_ID == myData2)]
          }
      if(myData1 == "All"){
        cpm <- log2cpm# & meta$SUBJECT_ID == myData2)]
      }
      print(sprintf("dim(CPM)%s", dim(cpm)))
      print(myData1)
      cpm
    })
      
  
  myDataMetaData <- reactive({
    myData1 <- input$'action_selectiontype'
    myData2 <- input$'action_selectiontype1'
    cpm <- meta
    if(myData1 != "All"){
      cpm <- meta[which(meta$disease == myData1),]
    }
    if(myData2 != "All"){
      cpm <- meta[which(meta$type == myData2),]
    }
    if(myData1 != "All"&myData2 != "All"){
      cpm <- meta[which(meta$disease == myData1 & meta$type == myData2),]
    }
    print(sprintf("dim(CPM)%s", dim(cpm)))
    print(myData1)
    cpm
  })
                        
  output$main_plot <-  renderPlot({
    generow <- as.numeric(SW())
    generow1 <- generow#which(generow==row.names(exprs(gset)))
    print(sprintf("main plot SW() %s", generow))
    Data1 <-myDatapca()
    #print(str(Data1))
    Data2 <- myDatacpm()
      x <- as.numeric(Data2[generow,])
    print(length(Data2[generow,]))
    scale01 <- (x-min(x))/(max(x)-min(x))
      scale01 <- scale01+.1
      scale01[scale01>1] <- 1
   print(ggplot(Data1, aes(PC1,PC2, color = as.numeric(Data2[generow,])))+
           geom_point(alpha=scale01)+
           ggtitle(myGene1Inputs())+
           scale_colour_gradient2(high="red", mid="blue", low="grey50", midpoint=4)+
           labs(color=NULL)
         )
   #print(ggplot(Data1, aes(PC1,PC2))+geom_point())
   print("Complete1")
  })
  
  output$main_plot2 <-  renderPlot({
    generow <- as.numeric(SW1())
    generow1 <- generow#which(generow==row.names(exprs(gset)))
    print(sprintf("main plot SW() %s", generow))
    Data1 <-myDatapca()
    #print(str(Data1))
    Data2 <- myDatacpm()
    print(length(Data2[generow,]))
    x <- as.numeric(Data2[generow,])
    scale01 <- (x-min(x))/(max(x)-min(x))
      scale01 <- scale01+.1
      scale01[scale01>1] <- 1
    print(ggplot(Data1, aes(PC1,PC2, color = as.numeric(Data2[generow,])))+geom_point(alpha=scale01)+
            ggtitle(myGene2Inputs())+ 
            scale_colour_gradient2(high="red", mid="blue", low="grey50", midpoint=4)+
            labs(color=NULL)
    )
    print("Complete2")
  })
  
  output$main_plot3 <-  renderPlot({
    generow1 <-  as.numeric(SW())
    generow2 <-  as.numeric(SW1())
 
    Data1 <- myDatacpm()[generow1,]
    Data2 <- myDatacpm()[generow2,]
    Data3 <- data.frame(t(Data1), t(Data2))
      names(Data3) <- c("X", "Y")
        DoubleZero <- sum(Data3$X<1 & Data3$Y <1)
        DoublePos <- sum(Data3$X>3 & Data3$Y >3) 
        XPos <- sum(Data3$X>3 & Data3$Y <3) 
        YPos <- sum(Data3$X<3 & Data3$Y >3) 
    colorM = as.character(unlist(myDataMetaData()[,4]))
    Data3$colorM = myDataMetaData()[,4]
    
    print(ggplot(Data3, aes(X,Y, color=colorM))+geom_point(size=3) +
          labs(x=myGene1Inputs(),y= myGene2Inputs())+ 
            theme(legend.text=element_text(size=14))+
            guides(colorM=guide_legend(override.aes=list(size=30)))+
            annotate("text", x=3, y=3, label=DoublePos)+
            annotate("text", x=1, y=1, label=DoubleZero)+
            annotate("text", x=1, y=3, label=YPos)+
            annotate("text", x=3, y=1, label=XPos)+
            annotate("text", x=2, y=2, label=XPos+DoublePos+DoubleZero+YPos)
    )
    print("Complete3")
  })
  
  output$main_plot4 <-  renderPlot({
    generow <- as.numeric(SW1())
    generow1 <- generow#which(generow==row.names(exprs(gset)))
    print(sprintf("main plot SW() %s", generow))
    Data1 <-pca_rotation
    print(ggplot(Data1, aes(PC1,PC2, color = as.character(unlist(meta[,4]))))+
            geom_point(size=3)+
            theme(legend.text=element_text(size=14))+
            labs(color=NULL)
    )
    print("Complete2")
  })
  
# #############Variable Panel
# # Combine the selected variables into a new data frame
# selectedData <- reactive({
#   pData(gset)[, c(input$xcol, input$ycol, "Groups", "visit", "pat_id")]
# })
# 
# 
# 
# 
# output$text1 <- renderTable(head(selectedData()))
# #data1 <- pData(gset)[c("Groups", "CRP", "visit", "pat_id")]
# #data1 <- pData(gset)[c("cohort", "CRP", "Groups", "pat_id")]
# output$text2 <- renderPrint({
#  data1 <- selectedData()
#  #if (class(data1[,2])=='character'){data1[,2]<- (factor(data1[,2]))}
#  data1a <- data1[grep("Aa|Ab", data1$Groups),]
#  data1b <- data1[grep("Ba|Bb", data1$Groups),]
#  data1c <- data1[grep("Ca|Cb", data1$Groups),]
#  data1hc <- data1[grep("HCa|HCb", data1$Groups),]
#  data1con <- data1[grep("HCa|Aa", data1$Groups),]
#  data1con2 <- data1[grep("HCa|Ba", data1$Groups),]
#  
#  GroupList <- list(data1a, data1b, data1c, data1hc, data1con, data1con2)
#     dataTabRes11 <- data.frame() 
#     dataTabRes1a <- data.frame()
#       for (i in GroupList){
#         #Names1 <- names(summary(lm(as.numeric(as.character(i[,2]))~i$Groups+(i$pat_id)))$coefficients[2,])
#         x <- (i[,2])
#         x[is.na(x)]<- 0
#         if(length(unique(i$pat_id))!=nrow(i)){  
#             dataTabRes11 <- rbind(dataTabRes11, 
#                                 summary(lm(x~i$Groups+(i$pat_id)))$coefficients[2,])
#             dataTabRes1a <- rbind(dataTabRes1a,
#                                 summary(lm(x~i[,1]+(i$pat_id)))$coefficients[2,])
#             }
#           else { dataTabRes11 <- rbind(dataTabRes11, 
#                                 summary(lm(x~i$Groups))$coefficients[2,])
#                  dataTabRes1a <- rbind(dataTabRes1a, 
#                                 summary(lm(x~i[,1]))$coefficients[2,])
#                 }
#         }
#       names(dataTabRes11) <- Names1
#       names(dataTabRes1a) <- Names1 
#       row.names(dataTabRes11) <- c("MTX", "TNF", "CohortC", "HealthyCon", "EarlyRA/Con","EstRA/Con" )
#       row.names(dataTabRes1a) <- row.names(dataTabRes11)
# 
# #    dataTabRes1a <- summary(lm(as.numeric(as.character(data1a[,2]))~data1a[,1]+(data1a$pat_id)))$coefficients[2,]
# #    dataTabRes2a <- summary(lm(as.numeric(as.character(data1b[,2]))~data1b[,1]+(data1b$pat_id)))$coefficients[2,]
# #    dataTabRes3a <- summary(lm(as.numeric(as.character(data1c[,2]))~data1c[,1]+(data1c$pat_id)))$coefficients[2,]
# #    dataTabRes4a <- summary(lm(as.numeric(as.character(data1hc[,2]))~data1hc[,1]+(data1hc$pat_id)))$coefficients[2,]
# #  
# #  dataTabRes <- rbind(dataTabRes1, dataTabRes2, dataTabRes3, dataTabRes4)
# #  dataTabResa <- rbind(dataTabRes1a, dataTabRes2a, dataTabRes3a, dataTabRes4a)
# #  row.names(dataTabRes) <- c("MTX", "TNF", "CohortC", "HealthyCon")
# #  row.names(dataTabResa) <- c("MTX", "TNF", "CohortC", "HealthyCon")
#  list(dataTabRes11, dataTabRes1a)
#  })
# 
# output$plot1 <- renderPlot({
#   Data <- selectedData()
#   ggplot(Data, aes(x=selectedData()[,1], y=Data[,2], group=Data[,1], 
#          color=Data[,1]))+ geom_boxplot(outlier.shape=NA)+ geom_jitter(width = .2, height=0)+
#         labs(list(x=colnames(Data[1]), y= colnames(Data[2]),  title=colnames(Data[2])))
#    
# })
# output$plot2 <- renderPlot({
#   Data <- selectedData()
#   ggplot(Data, aes(y=selectedData()[,2], x=Groups, group=Groups, 
#                    color=Data[,1]))+ geom_boxplot(outlier.shape=NA)+ 
#                   geom_jitter(width = .2, height=0)+
#                   labs(list(x="Cohorts", y= colnames(Data[2]), title=colnames(Data[2])))
#   
# })

}


