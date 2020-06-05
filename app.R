library(shiny)
library(ggplot2)
library(scales)

ui <- fluidPage(
  titlePanel("Predicting Binding Curves"),
  sidebarLayout(
    sidebarPanel(
      #sliderInput("concInput","Concentration Range 10^x nM", -10,10,c(-3,3)),
      numericInput("dilInput", "Number of Dilutions in Experiment", 16),
      numericInput("upperInput", "Highest Sample Concentration (nM)", 1000),
      numericInput("factorInput", "Dilution factor for samples (0-1)", 0.5),
      numericInput("kdInput", "Predicted Kd for Interaction (nM)", 10),
      numericInput("ligandInput","Ligand Concetration (constant) (nM)", 50),
      radioButtons("plotInput", "Plot Type", c("Line","Point", "Both"),selected = "Point")
    ),
    mainPanel(
      plotOutput("bindingplots"),
      br(), br(),
      p(strong("The solid black line indicates the specified constant concentration. The red line is predicted Kd.")),
      p("This graph shows two plots. The 'Ideal' is what a binding curve for an interaction with the specficied
      kd will look like if the experiment is set up properly. The 'Predicted' curve is what your binding experiment 
      will look like with a constant 
      concentration as stated."),
        p("If there is discrepancy between the two, then you are conducting your experiment 
      at a concentration too high to accurately get data for the interaction.")
    )
  )
)

server <- function(input, output) {
  
# Functions doing the Heavy Lifting
  hill <- function(b,m,P,Kd) {
    b + (m-b)*(1/(1+Kd/P))
  }
  hillq <- function(b,m,R,P,Kd){
    scaler <- b+(m-b)
    root = sqrt(((R+P+Kd)^2)-(4*R*P))
    scaler*(R+P+Kd-root)/(2*R)
  }
  
  outputdataframe <- reactive({
  
  # setting variables for plots
 
  b <- 0 # This is minimum value (0% bound)
  m <- 1 # This is maximum value (100% bound)
  R <- input$ligandInput # This is concentration of the constant partner / ligand
  
samplenum <- 1:(input$dilInput-1) # List of how many samples there are.
totaldilution <- input$upperInput # Highest starting concentration

totaldilution[2:input$dilInput] <- input$factorInput^(samplenum) * input$upperInput
  P <- totaldilution

  #function for producing dataframe and plot from inputs of hill and hillq
  plotbindhill <- function(input){
    temp1 <- data.frame(P,hill(b,m,P,input))
    temp1$curve <- "Ideal"
    temp2 <- data.frame(P,hillq(b,m,R,P,input))
    temp2$curve <- "Predicted"
    colnames(temp1) <- c("pconc","fbound","curve")
    colnames(temp2) <- c("pconc","fbound","curve")
    temp <- merge(temp1,temp2, all = TRUE)
    temp$kd <- input
    temp
    }

  list <- c(0.1,1,10)
    
  kd.data.group  <- plotbindhill(input$kdInput)

  })
  
  #sets up the graph with required parameters and data.
  mainplot <- reactive({
    R <- input$ligandInput
    type <- c("geom_line()", "geom_point()", "geom_line() + goem_point()")
    names(type) <- c("Line", "Point", "Both")
    plotype <- type[input$plotInput]
    ggplot(outputdataframe(), aes(pconc,fbound, colour=curve)) +
      labs(title = "Predicted Binding Curves", subtitle = paste("Binding curve with kd =", input$kdInput, "nM")) +
      xlab("Variable Concentration (nM)") + 
      ylab("Fraction Bound") +
      geom_linerange(aes(x=input$kdInput, ymax=0.5, ymin=-0.05), colour="gray80") +
      geom_vline(xintercept = R, colour="#FFBF00", alpha=0.5) +
      # geom_hline(yintercept = 0.5, colour="gray80") +
      # geom_vline(xintercept = input$kdInput, color = "red") + 
      theme_classic(base_size = 18) + 
      scale_x_log10(labels = trans_format("log10", function(x) 10^x)) + 
      theme(aspect.ratio = 9/16, legend.position = "right") + coord_cartesian(ylim=c(0,1))
  })
  
  #Creates the 3 types of graphs
  plotpoint <- reactive({
    mainplot() + geom_point(size = 3)
  })
  plotline <- reactive({
    mainplot() + geom_line()
  })
  plotboth <- reactive({
    mainplot() + geom_point(size = 3) + geom_line()
  })
  
  #Selects what type of graph.
  graphInput <- reactive({
    switch(input$plotInput,
           "Point" = plotpoint(),
           "Line" = plotline(),
           "Both" = plotboth()
    )
  })
  
  #renders the actual graph
  output$bindingplots <- renderPlot({ 
    graphInput()
  })

}

shinyApp(ui = ui, server = server)