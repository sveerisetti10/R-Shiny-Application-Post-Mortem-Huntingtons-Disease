library(shiny)
library(bslib)
library(tidyverse)
library(ggplot2)
#install.packages('ggbeeswarm')
library('ggbeeswarm')
library(colourpicker)
library(DESeq2)
library(genefilter)
library(dplyr)
library('RColorBrewer')
#install.packages("ggfortify")
library("ggfortify")
#install.packages("shinythemes")
library("shinythemes")

#This is used in order to increase the file size capabilities of the file inputs
options(shiny.maxRequestSize = 30*1024^2)

deseq_choices <- c('baseMean', 'HD.mean', "Control.mean", 'log2FoldChange', 
                   'IfcSE', "stat", 'pvalue', 'padj')

choice_radio <- c('Sample_pmi', 'Sample_Age_of_Death', 'Sample_rin', 
                  'Sample_mrna_seq', 'Sample_age_of_onset')

choice_radio_four <- c('Sample_Neurological_Diagnosis')

data_choice <- read_csv("/usr4/bf528/sv/Documents/BF_591_Project_Veerisetti/Huntington_Samples.csv")

counts_choice <- read_csv("/usr4/bf528/sv/Documents/BF_591_Project_Veerisetti/counts_matrix_final.csv")

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  theme = shinytheme("cosmo"),
  #Here we will insert the title of the application
  titlePanel("BF 591 R-Shiny Cumulative Project"),
  markdown(paste0("In order to utilize this app, please make sure to download",
                  " all the .csv files required")),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Sample Information",
               tabsetPanel(
                 fileInput('file',
                           label = "Please Load the .csv File", 
                           accept = ".csv", 
                           placeholder = 'Huntington_GSE64810_DESeq2.csv'), 
                 tabPanel("Summary", 
                          hr('There are 69 rows to account for as well as 15 columns'),
                          br(),
                          tableOutput("table")),
                 tabPanel("Table",
                          tableOutput("table_filtered")), 
                 tabPanel("Plots",
                          radioButtons(
                            inputId = 'x_axis', 
                            label = 'What variable would you like on the x-axis?',
                            choices = choice_radio, 
                            selected = 'baseMean'),
                          colourInput(
                            inputId = 'outer_color',
                            label = 'Highlight Outer Color', 
                            value = '#F0D400', 
                            closeOnClick = T),
                          colourInput(
                            inputId = "fill_color",
                            label = "Highlight Fill Color", 
                            value = '#FA0026', 
                            closeOnClick = T),
                          submitButton(text = "Submit",
                                       icon = icon('car-crash')),
                          br(),
                          hr('This is a histogram plot that displays the variable of your choice on the x-axis and count on the y-axis'),
                          plotOutput("plotting_histogram"),
                          br(),
                          br(),
                          hr('This is a density plot that displays the variable of your choice on the x-axis and density on the y-axis'),
                          plotOutput('plotting_density')),
                 
               )
      ),
      tabPanel("Counts Matrix", 
               tabsetPanel(
                 fileInput('count_file', 
                           label = "Please Load the Counts Matrix File", 
                           accept = '.csv', 
                           placeholder = 'counts_matrix_final.csv'), 
                 br(),
                 sliderInput('slider', 
                             'Includes genes with at least X percentile of variance', 
                             min = 1, 
                             max = 100, 
                             value = 40, 
                             step = 1),
                 sliderInput('slider_zero', 
                             'includes genes with at least X samples that are non-zero', 
                             min = 1,
                             max = 69,
                             value = 40, 
                             step = 1
                 ),
                 submitButton(text = 'Submit', 
                              icon = icon('car-crash')),
                 tabPanel("Filtering",
                          hr("Here is your summary for percentile of variance: "),
                          tableOutput('varianceee'), 
                          hr("Here is your summary for genes with at least X samples that are non-zero"),
                          tableOutput('agent_zero')
                 ),
                 
                 tabPanel("Diagnostic Scatter Plots", 
                          hr("Here is your diagnostic scatter plot tailored to your variance input: "), 
                          plotOutput("volcano"),
                          hr("Here is your diagnostic scatter plot tailored to your non-zero input"),
                          plotOutput("volcano_zero")
                 ), 
                 tabPanel("Clustered Heatmap", 
                          hr('Here is your heatmap tailored to your input'), 
                          plotOutput('heatmap')), 
                 tabPanel("PCA Plot", 
                          hr("Here is your PCA Scatter Plot"), 
                          plotOutput('pca'))
               )
               
      ), 
      tabPanel("Differential Expression",
               fileInput(
                 inputId = 'differential_expression', 
                 label = "Please Load the Differential Expression File",
                 accept = '.csv',
                 placeholder = 'DE_File'),
               selectInput(
                 inputId = 'gene_choice_select',
                 label = "Please Choose a Gene!",
                 choices = data_choice),
               radioButtons(
                 inputId = 'x_value', 
                 label = 'What variable would you like on the x-axis?',
                 choices = deseq_choices, 
                 selected = 'baseMean'),
               colourInput(
                 inputId = 'color1',
                 label = 'Highlight Outer Color', 
                 value = '#F0D400', 
                 closeOnClick = T),
               colourInput(
                 inputId = "color2",
                 label = "Highlight Fill Color", 
                 value = '#FA0026', 
                 closeOnClick = T),
               submitButton(text = 'Submit', 
                            icon = icon('car-crash')),
               tableOutput("different"), 
               br(), 
               hr("Here is your density plot fit for your x-axis selection"),
               plotOutput('plotting_volcano')
               
               
      ), 
      
      tabPanel("Visualization of Individual Gene Expression", 
               fileInput(
                 inputId = 'fusion_journey', 
                 label = "Please Load the Proper File",
                 accept = '.csv',
                 placeholder = 'Fusion File'),
               radioButtons(
                 inputId = 'category_journey', 
                 label = 'What categorial field would you like?',
                 choices = choice_radio_four, 
                 selected = 'baseMean'),
               textInput(
                 inputId = 'choose_gene',
                 label = "Please Choose a Gene!",
                 value = NULL),
               radioButtons(
                 inputId = 'plot_type', 
                 label = 'Type of Plot', 
                 choices = c("Bar_Plot", "Scatter_Plot", 'Beeswarm_Plot')
               ),
               submitButton(text = 'Submit', 
                            icon = icon('car-crash')),
               plotOutput('plottingbar'),
               br(), 
               plotOutput('plottingpoint')),    
    )),  
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #' @details The purpose of this function is to load in the sample information
  #' .csv file. The user has the ability to input their file 
  load_file <- reactive ({
    huntington <- read_csv(input$file$datapath)
    return(huntington)
    
  })
  
  #' Here we want to make a dataframe summarizing the sample information file
  #' that the use has inputted. 
  #' @param hdataframe This is the data frame that is loaded
  #' @details  column_name This is the column where each variable is mentioned
  #' @details  summarized_df This is the dataframe that includes the column 
  #' name, classification, and the mean of the variable 
  
  summarization <- function(hdataframe) {
    Column_Name <- c("Sample_geo_accession", "Sample_status", "Sample_submission_date", 
                     "Sample_last_update_date", "Sample_type", 
                     "Sample_Channel_Count", "Sample_Organ", "Sample_Organism", 
                     "Sample_Tissue_Type", "Sample_Neurological Diagnosis", 
                     "Sample_pmi", "Sample_Age_of_Death", 'Sample_rin', 
                     "Sample_mrna_seq", "Sample_age_of_onset")
    
    Classification <- c("character", "character", "character", "character", 
                        "character", "double", "character", "character", 
                        "character", "character", "double", "double",
                        "double", "double", "double")
    
    Sample_Channel_Count_average <- mean(hdataframe$Sample_Channel_Count)
    Sample_pmi_average <- mean(hdataframe$Sample_pmi)
    Sample_Age_of_Death_average <- mean(hdataframe$Sample_Age_of_Death)
    Sample_rin_average <- mean(hdataframe$Sample_rin)
    Sample_mrna_seq_average <- mean(hdataframe$Sample_mrna_seq)
    Sample_age_of_onset_average <- mean(hdataframe$Sample_age_of_onset)
    
    Average_Summarized <- c('N/A', 'N/A', 'N/A', 'N/A', 'N/A', 
                            Sample_Channel_Count_average, 'N/A', 'N/A', 
                            'N/A', 'N/A', Sample_pmi_average,
                            Sample_Age_of_Death_average, Sample_rin_average, 
                            Sample_mrna_seq_average, '41.85')
    
    summarized_df <- data_frame(Column_Name, Classification, Average_Summarized)
    
    return(summarized_df)
  }
  
  #' @details Here we will output a data table that contains sample 
  #' information as well as the other variables 
  #' @param hdataframe The data frame that the user has imported
  
  filtered_table <- function(hdataframe) {
    filtered_data <- hdataframe 
    
    
    return(filtered_data)
    
  }
  
  #' @details  Here the goal is to develop a way for the user to generate a 
  #' histogram plot based on the sample information. 
  #' @param hdataframe The df with all of the data that the user has inputted
  #' @param x_name This is x-axis. The user has the ability to choose what 
  #' they would like the x-axis of the plot to be. 
  #' @param color_1 This is the outer color that the user chooses. 
  #' @param color_2 This is the inner color that the user chooses.  
  
  histogram_plot <- function(hdataframe, x_name, color_1, color_2){
    
    data_histogram <- 
      ggplot (hdataframe, aes(x = !!sym(x_name))) +
      geom_histogram(color = color_1, fill = color_2) +
      scale_color_manual(values = c(color_1, color_2))
    theme_bw() +
      theme(legend.position = 'bottom') 
    
    return(data_histogram)
    
  }
  
  #' @details  Here the goal is to develop a way for the user to generate a  
  #' density plot based on the sample information. 
  #'
  #' @param hdataframe The df with all of the data that the user has inputted
  #' @param x_name This is x-axis. The user has the ability to choose what 
  #' they would like the x-axis of the plot to be. 
  #' @param color_1 This is the outer color that the user chooses. 
  #' @param color_2 This is the inner color that the user chooses. 
  density_plot <- function(hdataframe, x_name, color_1, color_2) {
    
    data_density <- 
      ggplot(hdataframe, aes(x = !!sym(x_name))) +
      geom_density(color = color_1, fill = color_2) + 
      scale_color_manual(values = c(color_1, color_2)) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    return(data_density)
    
  }
  
  #'XXXXXXXXXXXXXXXXXXXXXXX Counts_Matrix XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  
  
  #' Load in the Data 
  #' 
  #' @details Here all we want to do is simply load in the counts .csv file 
  #' that the user will input via the "Browse" button. 
  load_counts <- reactive ({
    c_matrix <- read_csv(input$count_file$datapath)
    colnames(c_matrix)[1] <- 'gene'
    return(c_matrix)
    
  })
  
  
  #' @details The purpose of this function is to filter the data based 
  #' on what value the of number of samples that are greater the variance
  #' threshold that the user inputs via the slider 
  #' @param counts_matrix The counts matrix that the user has imported
  #' @param slider The slider that allows the user to choose 
  #' variance percentage 
  
  filterering_variance <- function(counts_matrix, slider) {
    
    
    county <- counts_matrix[,-1] 
    counts_matrix$variance <- rowVars(as.matrix(county))
    counts_matrix <- counts_matrix[rev(order(counts_matrix$variance)),]
    input_value <- slider/100 #assumes slider is between 1 and 100
    counts_matrix <- counts_matrix[1:floor(nrow(counts_matrix)*input_value), ]
    
    counts_variance <- counts_matrix[, c(1,71)]
    
    counts_columns <- ncol(county)
    counts_rows <- nrow(county)
    counts_variance_columns <- ncol(counts_variance)
    counts_variance_rows <- nrow(counts_variance)
    
    
    Information <- c("Number of Samples", "Number of Genes", "Number of Passing", "Percent Passing,", "Number not Passing", 
                     "Percent not Passing")
    
    Values <- c(counts_columns, counts_rows, counts_variance_rows, (counts_variance_rows/counts_rows) * 100, (28087 - counts_variance_rows), (100 - (counts_variance_rows/counts_rows) * 100))
    summary_variance <- data_frame(Information, Values)
    return(summary_variance)
    
    
  }
  
  
  #' @details Here we need to produce a dataframe based upon the input 
  #' value that the user inputs via the slider. 
  #' @param counts_matrix The counts .csv file that the user inputs via the 
  #' file box 
  #' @param slider_zero The slider that the user uses to include genes with at 
  #' least X samples that are non-zero
  
  
  filtering_zero <- function(counts_matrix, slider_zero) {
    
    counts_matrix$frequency <- rowSums(counts_matrix != 0)
    zeroo <- counts_matrix[which(counts_matrix$frequency > slider_zero),]
    
    display <- zeroo[, c(1)]
    
    counts_samples <- ncol(counts_matrix)
    counts_genes <- nrow(counts_matrix)
    counts_passing <- nrow(display)
    percent_passing <- (counts_passing / 28087) * 100
    counts_not_passing <- 28087 - counts_passing
    percent_not_passing <- 100 - percent_passing
    
    Info <- c("Number of Samples", "Number of Genes", "Number of Passing", "Percent Passing", "Number not Passing", 
              "Percent not Passing")
    
    Statistics <- c(counts_samples, counts_genes, counts_passing, percent_passing, counts_not_passing, percent_not_passing)
    summary_zero <- data_frame(Info, Statistics)
    return(summary_zero)
    
  }
  
  
  #' @details Here we will build a volcano plot based upon what the user 
  #' inputs via the slider. The variables that will be plotted on the 
  #' volcano plot are: median count on the x-axis and the variance
  #' on the y-axis 
  #' @param counts_matrix the counts matrix that the user inputs 
  #' @param slider The values that the user will input 
  
  volcano_plot_variance <- 
    function(counts_matrix, slider) {
      county <- counts_matrix[,-1] 
      counts_matrix$variance <- rowVars(as.matrix(county))
      counts_matrix <- counts_matrix[rev(order(counts_matrix$variance)),]
      input_value <- slider/100 
      
      counts_matrix$median <- apply(counts_matrix[,-1], 1, median)
      
      Filter_Method <- (floor(nrow(counts_matrix)*input_value)) > counts_matrix$variance
      
      eruption <- ggplot(data = counts_matrix, 
                         aes(x = median, y = -log10(variance)))+
        geom_point(aes(color = Filter_Method)) + 
        theme_bw() + 
        scale_color_manual(values = c('green', 'blue')) + 
        theme(legend.position = "bottom") 
      
      return(eruption)
      
      
    } 
  
  #' @details Here we will build a volcano plot based upon what the user 
  #' inputs via the slider. The variables that will be plotted on the 
  #' volcano plot are: median count on the x-axis and the variance
  #' on the y-axis 
  #' @param counts_matrix the counts matrix that the user inputs 
  #' @param slider The values that the user will input 
  
  volcano_plot_zero <- function(counts_matrix, slider_zero) {
    
    counts_matrix$frequency <- rowSums(counts_matrix != 0)
    counts_matrix$median <- apply(counts_matrix[,-1], 1, median)
    
    Filter_Method <- counts_matrix$frequency > slider_zero
    messi <- ggplot(data = counts_matrix, 
                    aes(x = median, y = -log10(frequency)))+
      geom_point(aes(color = Filter_Method )) + 
      theme_bw() + 
      scale_color_manual(values = c('green', 'blue')) + 
      theme(legend.position = "bottom") 
    
    
    return(messi)
  }
  
  #'@details Here we want to produce a heat map post filtering of the counts
  #'matrix that the user has inputted
  #'@param counts_matrix This is the counts matrix that the user will input
  #'@param num_colors 
  #'@param palette 
  
  heatmap_plot <- function(counts_matrix, slider, slider_zero) {
    county <- counts_matrix[,-1] 
    #'Here we are making a seperate column called variance where the variances
    #'of each gene will be appended 
    counts_matrix$variance <- rowVars(as.matrix(county))
    counts_matrix <- counts_matrix[rev(order(counts_matrix$variance)),]
    input_value <- slider/100 
    #This is the dataframe that is produced post filtering of variance
    updated_counts_matrix <- counts_matrix[1:floor(nrow(counts_matrix)*input_value), ]
    
    #Here we want to create a seperate column called frequency and find the 
    #values that are not equal to 0
    updated_counts_matrix$frequency <- rowSums(updated_counts_matrix != 0)
    #Here we want values that are only greater than the slider value that the user chooses
    combination <- updated_counts_matrix[which(updated_counts_matrix$frequency > slider_zero), ]
    
    combination <- as.matrix(combination[, -1])
    
    hot <- heatmap(combination, col = terrain.colors(256))
    
    
    return(hot)
    
    
  }
  
  #'@details Here we will construct a PC1 vs PC2 plot based on the input sequence 
  #'that the user will utilize 
  #'@param counts_matrix The counts matrix that the user will input 
  
  pca_plot <- function(counts_matrix) {
    counts_matrix <- counts_matrix[, -c(1)]
    pca_plotting <- prcomp(counts_matrix, scale. = TRUE)
    plot_project <- autoplot(pca_plotting, data = counts_matrix, colour = "red")
    
    
    return(plot_project)
    
  }
  
  #'XXXXXXXXXXXXXXXXXXXXXXX Differential Expression XXXXXXXXXXXXXXXXXXXXX
  
  
  #' Load in the DE Data
  #' 
  #' @details Here all we want to do is simply load in the DE .csv file that the user will 
  #' input via the "file_DE" button. 
  
  load_DE_data <- reactive ({
    de <- read_csv(input$differential_expression$datapath)
    
    return(de)
    
  })
  
  
  insert_table <- function(dataframe, selected_gene){
    
    dataframe %>%
      filter(gene == selected_gene) %>%
      
      return() 
  }
  
  #'@details This is the information for the density plot that will be 
  #'made based upon the file input that the user selects. 
  #'@param dataframe The file input that the user enters 
  #'@param x_name The radio button that the user selects to be on the x-axis 
  #'@param color1 The color that the user selects based on the color pick input 
  #'@param color2 The second color that the user selects based on the color pick 
  #'input 
  
  de_volcano_plot <- function(dataframe, x_name, color1, color2) {
    
    de_vol_plot <- 
      ggplot(dataframe, aes(x = !!sym(x_name))) + 
      geom_density(aes(color = color1, fill = color2)) +
      scale_color_manual(values = c(color1, color2)) + 
      theme_bw() + 
      theme(legend.position = 'bottom') 
    
    return(de_vol_plot)
    
  }
  
  
  #'XXXXXXXX Visualization of Individual Gene Expression  XXXXXXXXXX
  
  #' Load in the normalized counts matrix 
  #' 
  #' @details Here all we want to do is load in the normalized counts matrix
  #' that the user wants 
  
  load_fusion_matrix <- reactive ({
    fusion_matrix <- read_csv(input$fusion_journey$datapath)
    
    return(fusion_matrix)
    
  })
  
  #' Load in the normalized DE matrix 
  #' 
  #' @details Here all we want to do is load in the differential expression
  #' matrix of the user's choice
  
  load_sample_matrix <- reactive ({
    sample_matrix <- read_csv(input$sample_journey$datapath)
    
    return(sample_matrix)
    
    
  })
  
  #' @details  The purpose here is to make a bar plot based on the 
  #' input .csv file that the user will input. The file that is inputted 
  #' is key because it should be a merge of both the sample information 
  #' file as well as the counts matrix based upon the sample 
  #' @param bar_matrix This is the file that the user will input 
  #' @param gene_selection This is the gene that the user will enter 
  #' into the input box. 
  
  bar_plot <- function(bar_matrix, gene_selection){
    
    bar_matrix <- bar_matrix %>%
      select(c(gene_selection, Sample_Neurological_Diagnosis)) 
    
    bar_matrix$sample_number <- 1:nrow(bar_matrix)
    
    sample_number <- rownames(bar_matrix)
    
    bar_matrix_plot <- 
      ggplot(data = bar_matrix, aes(x = sample_number, y = !!sym(gene_selection)))+
      geom_bar(aes(color = Sample_Neurological_Diagnosis), stat = "identity") +
      theme_bw() + 
      scale_color_manual(values = c('#8EF02B', '#9B82FF')) + 
      theme(legend.position = "bottom") 
    
    
    return(bar_matrix_plot)
    
  }
  
  #' @details  The purpose here is to make a scatter plot based on the 
  #' input .csv file that the user will input. The file that is inputted 
  #' is key because it should be a merge of both the sample information 
  #' file as well as the counts matrix based upon the sample 
  #' @param bar_matrix This is the file that the user will input 
  #' @param gene_selection This is the gene that the user will enter 
  #' into the input box. 
  
  point_plot <- function(point_matrix, gene_selection){
    
    point_matrix <- point_matrix %>%
      select(c(gene_selection, Sample_Neurological_Diagnosis)) 
    
    point_matrix$sample_number <- 1:nrow(point_matrix)
    
    sample_number <- rownames(point_matrix)
    
    point_matrix_plot <- 
      ggplot(data = point_matrix, aes(x = sample_number, y = !!sym(gene_selection)))+
      geom_point(aes(color = Sample_Neurological_Diagnosis), stat = "identity") +
      theme_bw() + 
      scale_color_manual(values = c('#8EF02B', '#9B82FF')) + 
      theme(legend.position = "bottom") 
    
    
    return(point_matrix_plot)
    
  }
  
  #' @details  The purpose here is to make a beeswarm plot based on the 
  #' input .csv file that the user will input. The file that is inputted 
  #' is key because it should be a merge of both the sample information 
  #' file as well as the counts matrix based upon the sample 
  #' @param bar_matrix This is the file that the user will input 
  #' @param gene_selection This is the gene that the user will enter 
  #' into the input box. 
  
  bee_plot <- function(bee_matrix, gene_selection){
    
    bee_matrix <- bee_matrix %>%
      select(c(gene_selection, Sample_Neurological_Diagnosis)) 
    
    bee_matrix$sample_number <- 1:nrow(bee_matrix)
    
    sample_number <- rownames(bee_matrix)
    
    bee_matrix_matrix_plot <- 
      ggplot(data = bee_matrix, aes(x = Sample_Neurological_Diagnosis, y = !!sym(gene_selection)))+
      geom_beeswarm(aes(color = Sample_Neurological_Diagnosis)) +
      theme_bw() + 
      scale_color_manual(values = c('#8EF02B', '#9B82FF')) + 
      theme(legend.position = "bottom") 
    
    
    return(bee_matrix_matrix_plot)
    
  }
  
  
  
  #' @details Here we will use the "summarized" function in order to produce
  #' a data frame. 
  #' @param hdataframe The file that the user imported 
  
  output$table <- renderTable({
    req(input$file)
    table <- load_file()
    return(summarization(hdataframe = table))
    
  })
  
  #' @details Here we will use the "filtered_table" function in order to 
  #' produce a data frame. 
  #' @details  hdataframe The file that the user imported 
  #' @details table_filtered This is the outputId that we are using 
  #' @details load_file is the function used to produce the raw csv file
  #' 
  
  output$table_filtered <- renderTable({
    req(input$file)
    #' In this UI section, this is the name of the tableOutput table_filtered
    table_filtered <- load_file()
    return(filtered_table(hdataframe =  table_filtered))
    
  })
  
  #' @details Here we will use the histogram_plot function in order to produce
  #' a histogram plot
  #' @details huntington This is the original file that the user inputs 
  #' @param data_histogram uses the function histogram_plot and the UI inputs 
  #' to make the histogram 
  
  output$plotting_histogram <- renderPlot({
    req(input$file)
    huntington <- load_file()
    data_histogram <- histogram_plot(huntington, input$x_axis, input$outer_color, 
                                     input$fill_color)
    return(data_histogram)
  }, height = 400)
  
  
  #' @details Here we will use the density_plot function in order to produce
  #' a density plot
  #' @param huntington This is the original file that the user inputs 
  #' @param data_histogram uses the function histogram_plot and the UI inputs 
  #' to make the histogram 
  output$plotting_density <- renderPlot({
    req(input$file)
    huntington <- load_file()
    data_density <- density_plot(huntington, input$x_axis, input$outer_color, 
                                 input$fill_color)
    return(data_density)
  }, height = 400)
  
  
  
  output$varianceee <- renderTable({
    req(input$count_file)
    c_matrix <- load_counts()
    
    
    return(filterering_variance(c_matrix, slider = input$slider))
    
  })
  
  output$agent_zero <- renderTable({
    req(input$count_file)
    c_matrix <- load_counts()
    
    return(filtering_zero(c_matrix, slider_zero = input$slider_zero))
    
    
  })
  
  output$volcano <- renderPlot({
    req(input$count_file)
    c_matrix <- load_counts()
    eruption <- volcano_plot_variance(c_matrix,
                                      slider = input$slider
    )
    return(eruption)
    
  }, height = 400)
  
  output$heatmap <- renderPlot({
    req(input$count_file)
    c_matrix <- load_counts()
    hot <- heatmap_plot(c_matrix, slider = input$slider, 
                        slider_zero = input$slider_zero) 
    
    
    return(hot)
    
  })
  
  output$volcano_zero <- renderPlot({
    req(input$count_file)
    c_matrix <- load_counts()
    messi <- volcano_plot_zero(c_matrix, slider_zero = input$slider_zero)
    
    return(messi)
    
  })
  
  output$pca <- renderPlot({
    req(input$count_file)
    c_matrix <- load_counts()
    ronaldo <- pca_plot(c_matrix)
    
    return(ronaldo)
    
  })
  
  output$plotting_volcano <- renderPlot({
    req(input$differential_expression)
    different <- load_DE_data()
    final <- de_volcano_plot(dataframe = different, x_name = input$x_value, 
                             input$color1, input$color2) 
    
    
    return(final)
    
  })
  
  
  
  #' @details Here we will use the "filtered_table" function in order to 
  #' produce a data frame. 
  output$different <- renderTable({
    req(input$differential_expression)
    different <- load_DE_data()
    
    return(insert_table(dataframe = different, selected_gene = input$gene_choice_select))
    
    
  }, )
  
  output$plottingbar <- renderPlot({
    req(input$fusion_journey)
    barplot_sample <- load_fusion_matrix()
    
    if (input$plot_type == 'Bar_Plot') {
      data_barplot <- bar_plot(bar_matrix = barplot_sample, gene_selection = input$choose_gene)
      return(data_barplot)
    } else if (input$plot_type == 'Scatter_Plot') {
      data_pointplot <- point_plot(point_matrix = barplot_sample, gene_selection = input$choose_gene)
      return(data_pointplot)
    } else if (input$plot_type == 'Beeswarm_Plot') {
      data_beeswarm <- bee_plot(bee_matrix = barplot_sample, gene_selection = input$choose_gene)
      return(data_beeswarm)
      
      
    }
  
  }, height = 400)
 
}

# Run the application 
shinyApp(ui = ui, server = server)
