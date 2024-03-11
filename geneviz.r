library(shiny)
library(pheatmap)
library(plotly)
library(tidyverse)
library(ggpubr)
library(bslib)


#UPLOAD AND PROCESS DATA
cleaned_data = read.csv("cleaned_data.csv", row.names = 1)
full_exp = read.csv("datafull_w_demograph.csv", row.names = 1)
full_exp$gender <- ifelse(full_exp$gender == 1, "male", "female")
matrix = t(cleaned_data)
corr = cor(matrix)

cleaned_dataT = data.frame(t(cleaned_data))
cleaned_dataT$mean = rowMeans(cleaned_dataT)

pca  =  prcomp(matrix)

compute_counts <- function(dataset, grouping_vars) {
  counts <- dataset %>%
    count(across(all_of(grouping_vars)), name = "Count")
  
  return(counts)
}

# Specify the variables for grouping
grouping_variables <- c("age", "gender", "hardyscale")

# Compute counts
result_counts <- compute_counts(full_exp, grouping_variables)[-c(1),]


#Define UI
ui =  navbarPage("Visualizations",
      #theme = bs_theme(version = 5, bootswatch = "cerulean" ),
      tabPanel("Overview",
               h3("These visualizations are the summary of the gene expressions in human brain frontal cortex tissues."),
               h3(" The goal of these graph is to make scientists easily notice the trends and outliers in this dataset, 
                  and, therefore, gain a direction in further investigation."),
              h3("The questions can be answered through these visualtions:"),
              h4("- What are the demographics of the samples"),
              h4("- Solely based on expression levels, what are the correlation between the genes"),
              h4("- Is there any difference between male and female samples?"),
              h4("- Is there any seggragation between samples of different age group?")),
      
      tabPanel("Demographic Summary",
             h2("Summary of demographics information of the samples"),
             plotOutput("Balloon_Plot"),
             h5("Note: Hardy Scale is Death classification based on the 4-point scale:
1) Violent and fast death Deaths due to accident, blunt force trauma or suicide, terminal phase estimated at < 10 min. 
2) Fast death of natural causes Sudden unexpected deaths of people who had been reasonably healthy, after a terminal phase estimated at < 1 hr (with sudden death from a myocardial infarction as a model cause of death for this category) 
3) Intermediate death Death after a terminal phase of 1 to 24 hrs (not classifiable as 2 or 4); patients who were ill but death was unexpected 
4) Slow death Death after a long illness, with a terminal phase longer than 1 day (commonly cancer or chronic pulmonary disease); deaths that are not unexpected 
0) Ventilator Case All cases on a ventilator immediately before death.")
    ),
    
    tabPanel("Expression Correlation",
             h2("Gene expression correlation"),
             plotlyOutput(outputId = "corr"),
             h4("Axes: GeneID. Select region to zoom in. Hover over to view datapoint")
    ),
    tabPanel("Gender groups",
             h4("Below is a scatterplot showing the expressions of 3 genes among male and female samples with confidence ellipses"),
             canvasXpressOutput("threeDPlot"),
             h4("To hide/unhide a group, click on its icon on the right")
    ),
    tabPanel("Age Groups",
             h4("Principal component analysis (PCA) reduces the dimensionality of multivariate data, 
                to two dimensions that can be visualized graphically with minimal loss of information. 
                Overlapping confidence ellipses showed that there is little seggration between samples of different ages"),
             plotlyOutput("PCA"),
             h4("To hide/unhide a group, click on its icon on the right")
    )
  )

  

 

#Define server logic

server  = function (input, output) {
  output$corr = renderPlotly({
    heatmaply_cor(
      corr,
      colors = GnBu,
      k_col = 1, 
      k_row = 1,
      fontsize_row = 8,
      fontsize_col = 8,
      show_dendrogram = c(FALSE, FALSE),
      showticklabels = c(FALSE, FALSE),
    )
  })
  
  output$threeDPlot = renderCanvasXpress ({
    canvasXpress(
      data=matrix,
      varAnnot=as.data.frame(full_exp$gender, row.names=rownames(full_exp)),
      axisTickScaleFontFactor=0.6,
      axisTitleScaleFontFactor=0.6,
      grid = FALSE,
      ellipseBy="full_exp$gender",
      colorBy="full_exp$gender",
      graphType="Scatter3D",
      title="Gender groups",
      xAxisTitle = "log(Gene1)",
      yAxisTitle = "log(Gene2)",
      zAxisTitle = "log(Gene3)",
      xAxisMinorTicks         = FALSE,
      yAxisMinorTicks         = FALSE,
      showLoessFit = FALSE)
  })
  
  output$PCA = renderPlotly({ fviz_pca_ind(pca,
                                           
                                           geom = "point",
                                           
                                           addEllipses = TRUE,
                                           
                                           label = "none",
                                           xlab = "First Main Dimension",
                                           ylab = "Second Main Dimension",
                                           habillage=full_exp$age,
                                           
                                           
                                           title=NULL)
  })
  output$Balloon_Plot = renderPlot({
    ggplot(result_counts, aes(x =hardyscale, y = age)) +
      geom_point( aes(size=Count, color = Count))+facet_wrap(~gender)+
      theme_bw()+scale_size(range = c(10, 20))+
      xlab("Hardy Scale") +
      ylab("Age") +scale_color_distiller(palette = "PuBuGn", direction = 1)
  })
}


#Run the app
shinyApp(ui = ui, server = server)