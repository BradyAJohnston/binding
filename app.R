library(shiny)
library(ggplot2)
library(scales)

bcl <- read.csv("bcl-data.csv", stringsAsFactors = FALSE)

ui <- fluidPage(
  titlePanel("Predicting Binding Curves"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("concInput","Concentration Range 10^x nM", -10,10,c(-3,3)),
      numericInput("dilInput", "Number of Dilutions in Experiment", 16),
      #numericInput("upperInput", "Upper Concentration (nM)", 1000),
      #numericInput("lowerInput", "Lower concentration (")
      numericInput("kdInput", "Predicted Kd for Interaction (nM)", 10),
      numericInput("ligandInput","Ligand Concetration (constant) (nM)", 50),
      radioButtons("plotInput", "Plot Type", c("Line","Point", "Both"),selected = "Point")
    ),
    mainPanel(
      plotOutput("bindingplots"),
      br(), br(),
      p(strong("The solid black line indicates the specified constant concentration.")),
      p("This graph shows two plots. The 'Ideal' is what a binding curve for an interaction with the specficied
      kd will look like. The 'Predicted' curve is what your binding experiment will look like with a constant 
      concentration as stated. If there is discrepancy between the two, then you are conducting your experiment 
      a concentration too high to accurately get data for the interaction.")
    )
  )
)

server <- function(input, output) {
  
# Functions doing the Heavy Lifting
  #####
  hill <- function(b,m,P,Kd) {
    b + (m-b)*(1/(1+Kd/P))
  }
  hillq <- function(b,m,R,P,Kd){
    scaler <- b+(m-b)
    root = sqrt(((R+P+Kd)^2)-(4*R*P))
    scaler*(R+P+Kd-root)/(2*R)
  }
  
  outputdataframe <- reactive({
  
  ## setting variables for plots (will include definitions)
 
  b <- 0
  m <- 1
  R <- input$ligandInput
  logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n)) 
  P <- logspace(input$concInput[1], input$concInput[2], input$dilInput)
  #P <- logspace(-3, 3, 50)
  #  P<-1
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
    assign(paste("kd.data.", input, sep = ""), temp, envir = .GlobalEnv)
  }

  #rm(kd.data.group)
  list <- c(0.1,1,10)
    
  kd.data.group  <- plotbindhill(input$kdInput)
  
  #for(i in list) {
  #  plotbindhill(i)
  #}
  
  
  #kd.data.group <- Reduce(function(x,y){
  #  merge(x=x,y=y, all=TRUE)
  #}, list(kd.data.0.1,kd.data.1,kd.data.10))

  })

##############
  
  mainplot <- reactive({
    R <- input$ligandInput
    type <- c("geom_line()", "geom_point()", "geom_line() + goem_point()")
    names(type) <- c("Line", "Point", "Both")
    plotype <- type[input$plotInput]
    ggplot(outputdataframe(), aes(pconc,fbound,colour=curve)) +
      labs(title = "Predicted Binding Curves", subtitle = paste("Binding curve with kd =", input$kdInput, "nM")) +
      xlab("Variable Concentration (nM)") + 
      ylab("Fraction Bound") +
      geom_vline(xintercept = R) +
      scale_x_log10(labels = comma) + 
      theme(aspect.ratio = 9/16, legend.position = "bottom") #+
    #facet_wrap(vars(kd))
    
  })
  
  plotpoint <- reactive({
    mainplot() + geom_point(size = 3)
  })
  plotline <- reactive({
    mainplot() + geom_line()
  })
  plotboth <- reactive({
    mainplot() + geom_point(size = 3) + geom_line()
  })
  
  graphInput <- reactive({
    switch(input$plotInput,
           "Point" = plotpoint(),
           "Line" = plotline(),
           "Both" = plotboth()
    )
  })
  
  output$bindingplots <- renderPlot({ 
    graphInput()
  })
  
  
  

}

shinyApp(ui = ui, server = server)