library(shiny)
library(shinythemes)
library(dplyr)
library(readr)

dagger <- function(X){
  t(Conj(X))
}

compareMatrix <- function(a,b){
  max(isTRUE(all.equal(a,b)),isTRUE(all.equal(a,1i*b)),
      isTRUE(all.equal(a,-1*b)),isTRUE(all.equal(a,-1i*b)))
}

# Pauli group
I <- matrix(c(1+0i,0,0,1),nrow=2,ncol=2,byrow = TRUE)
X <- matrix(c(0+0i,1,1,0),nrow=2,ncol=2,byrow = TRUE)
Z <- matrix(c(1+0i,0,0,-1),nrow=2,ncol=2,byrow = TRUE)
Y <- matrix(c(0+0i,-1i,1i,0),nrow=2,ncol=2,byrow = TRUE)

# Clifford group
ClifM01 <- matrix(c(1+0i,0,0,1) ,nrow = 2,byrow = TRUE)
ClifM02 <- matrix(exp(-1i*pi/4)*c(1,0,0,1i),nrow = 2,byrow = TRUE)
ClifM03 <- matrix(-1i*c(1,0,0,-1),nrow = 2,byrow = TRUE)
ClifM04 <- matrix(exp(1i*pi/4)*c(1,0,0,-1i),nrow = 2,byrow = TRUE)
ClifM05 <- matrix((-1)*c(0+0i,1,-1,0),nrow = 2,byrow = TRUE)
ClifM06 <- matrix(-exp(1i*pi/4)*c(0,1,1i,0),nrow = 2,byrow = TRUE)
ClifM07 <- matrix((-1i)*c(0,1,1,0),nrow = 2,byrow = TRUE)
ClifM08 <- matrix(exp(-1i*pi/4)*c(0,1,-1i,0),nrow = 2,byrow = TRUE)
ClifM09 <- matrix((-1i/(sqrt(2)))*c(1,1,1,-1),nrow = 2,byrow = TRUE)
ClifM10 <- matrix((1/(sqrt(2)))*c(1+0i,1,-1,1),nrow = 2,byrow = TRUE)
ClifM11 <- matrix((exp(-1i*pi/4)/(sqrt(2)))*c(1,1,-1i,1i),nrow = 2,byrow = TRUE)
ClifM12 <- matrix((-exp(1i*pi/4)/(sqrt(2)))*c(1,1,1i,-1i),nrow = 2,byrow = TRUE)
ClifM13 <- matrix((1i/(sqrt(2)))*c(1,-1,-1,-1),nrow = 2,byrow = TRUE)
ClifM14 <- matrix((exp(-1i*pi/4)/(sqrt(2)))*c(1,-1,1i,1i),nrow = 2,byrow = TRUE)
ClifM15 <- matrix((1/(sqrt(2)))*c(1+0i,-1,1,1),nrow = 2,byrow = TRUE)
ClifM16 <- matrix((exp(1i*pi/4)/(sqrt(2)))*c(1,-1,-1i,-1i),nrow = 2,byrow = TRUE)
ClifM17 <- matrix((exp(-1i*pi/4)/(sqrt(2)))*c(1,1i,-1,1i),nrow = 2,byrow = TRUE)
ClifM18 <- matrix((exp(1i*pi/4)/(sqrt(2)))*c(1,1i,1,-1i),nrow = 2,byrow = TRUE)
ClifM19 <- matrix((1i/(sqrt(2)))*c(1,1i,-1i,-1),nrow = 2,byrow = TRUE)
ClifM20 <- matrix((1/(sqrt(2)))*c(1,1i,1i,1),nrow = 2,byrow = TRUE)
ClifM21 <- matrix((exp(1i*pi/4)/(sqrt(2)))*c(1,-1i,-1,-1i),nrow = 2,byrow = TRUE)
ClifM22 <- matrix((1/(sqrt(2)))*c(1,-1i,-1i,1),nrow = 2,byrow = TRUE)
ClifM23 <- matrix((-1i/(sqrt(2)))*c(1,-1i,1i,-1),nrow = 2,byrow = TRUE)
ClifM24 <- matrix((exp(-1i*pi/4)/(sqrt(2)))*c(1,-1i,1,1i),nrow = 2,byrow = TRUE)

GetGroupstructure <- function(MatrixList){
  m <- length(MatrixList)
  GroupStructure <- matrix(rep(0,m^2),nrow = m)
  for(i in 1:m){ # Current state
    for(k in 1:m){ # Gate
      ki <- MatrixList[[k]] %*% MatrixList[[i]]
      for(j in 1:m){ # Gate
        if(compareMatrix(ki, MatrixList[[j]])){
          GroupStructure[i,j] <- k
          break
        }
      }
    }
  }
  return(GroupStructure)
}

GetErrorStructure <- function(MatrixList){
  m <- length(MatrixList)
  ErrorGroupStructure <- array(rep(0,m^3),rep(m,3))
  for(i in 1:m){ # Current error state
    for(k in 1:m){ # Gate used
      ki <- MatrixList[[k]] %*% MatrixList[[i]]
      for(l in 1:m){ # Error gate used
        lki <- MatrixList[[l]] %*%  ki
        for(j in 1:m){ # Goal error state
          if(compareMatrix(lki, MatrixList[[j]])){
            ErrorGroupStructure[k,i,j] <- l
            break
          }
        }
      }
    }
  }
  return(ErrorGroupStructure)
}

try(load("myRdata.Rdata"))

ui <- fluidPage(theme = shinytheme("lumen"),
   
   titlePanel("Transition Matrices of Markov chains for error accumulation in quantum circuits"),
   
   sidebarLayout(
      sidebarPanel(
        selectInput(inputId = "Group", label=strong("Pauli or Clifford Group"),
                    choices = c("Pauli","Clifford"),
                    selected = "Pauli"
        ),
        selectInput(inputId = "PorQChoice", label=strong("P or Q matrix"),
                    choices = c("P","Q"),
                    selected = "P"
        ),
        conditionalPanel(condition = "input.PorQChoice=='Q'",
          selectInput(inputId = "GateChoice", label=strong("Which gate is used"),
                      choices = c("I","X","Y","Z"),
                      selected = "I"
          )
        ),
        conditionalPanel(condition = "input.PorQChoice=='P'",
          checkboxInput(inputId = "UniformGate", label = strong("Uniform gate probabilites"),TRUE
          ),
          conditionalPanel(condition = "input.UniformGate == false",
                         textInput("GateProb", 'Enter a vector of gate probabilities (comma delimited)', "")
          )
        ),
        checkboxInput(inputId = "UniformError", label = strong("Equal probability for all errors"),TRUE
        ),
        conditionalPanel(condition = "input.UniformError == false",
                         textInput("ErrorProb", 'Enter a vector of error probabilities, including the probability for no error (comma delimited)', "")
        ),
        conditionalPanel(condition = "input.UniformError == true",
                         textInput("TotErrorProb", 'Enter the total probability of an error occuring', "0.1")
        ),
        actionButton("letsgo", "Run code"
        ),
        downloadButton("downloadData", "Download as csv"
        ),
        tags$hr(),
        downloadLink("readme","README"),
        tags$hr()
      ),
      
      mainPanel(
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         h4("The programm is running, please wait ..."),
                         h4("This might take a couple of minutes.")),
        h4(textOutput("caption")),
        conditionalPanel(condition = "input.Group=='Clifford'",
                         h4(textOutput("subcaption"))),
        tableOutput("mytable"),
        h4("Press the \"DOWNLOAD AS CSV\" button to download the full matrix in csv format"),
        h4("The README contains more information about the program and how to use it")
      )
   )
)

server <- function(input, output, session) {
  
  GateChoices <- reactive({
    if(input$Group=="Pauli"){
      c("I","X","Y","Z")
    }else{
      as.character(1:24)
    }
  })
  
  DefaultChoice <- reactive({
    if(input$Group=="Pauli"){
      "I"
    }else{
      "1"
    }
  })
  
  observe({
  updateSelectInput(session, "GateChoice", choices = GateChoices(),selected = DefaultChoice())
  })
  SetupList <- reactive({
    if(!exists("myList")){
    MatrixList <- list(I,X,Y,Z)
    PauliGS <- GetGroupstructure(MatrixList)
    PauliES <- GetErrorStructure(MatrixList)
    
    MatrixList <- list(ClifM01,ClifM02,ClifM03,ClifM04,ClifM05,ClifM06,ClifM07,ClifM08,ClifM09,ClifM10,ClifM11,ClifM12,
                       ClifM13,ClifM14,ClifM15,ClifM16,ClifM17,ClifM18,ClifM19,ClifM20,ClifM21,ClifM22,ClifM23,ClifM24)
    CliffordGS <- GetGroupstructure(MatrixList)
    CliffordES <- GetErrorStructure(MatrixList)
    myList <- list(PauliGS,PauliES,CliffordGS,CliffordES)
    save(myList, file="myRdata.Rdata")
    }
    myList
  })
  
  
  MyPmatrix <- reactive({
    
    input$letsgo
    # Get all variables
    Group <- isolate(input$Group)
    UniformGate <- isolate(input$UniformGate)
    GateProb <- isolate(input$GateProb)
    UniformError <- isolate(input$UniformError)
    ErrorProb <- isolate(input$ErrorProb)
    TotErrorProb <- isolate(input$TotErrorProb)

    if(Group == "Pauli"){
      GroupStructure <- SetupList()[[1]]
      ErrorGroupStructure <- SetupList()[[2]]
      m <- 4
    }else{
      GroupStructure <- SetupList()[[3]]
      ErrorGroupStructure <- SetupList()[[4]]
      m <- 24
    }

    if(UniformGate){
      Pgate <- rep(1/m,m)
    }else{
      Pgate <- unlist(lapply(unlist(strsplit(GateProb,",")),function(x) eval(parse(text=x))))
    }
    
    if(UniformError){
      r <- eval(parse(text=TotErrorProb))
      Perror <- matrix(rep( c(1-r,rep(r/(m-1),m-1)),m),nrow = m, byrow = TRUE)
    }else{
      errorlist <- unlist(lapply(unlist(strsplit(ErrorProb,",")),function(x) eval(parse(text=x))))
      Perror <- matrix(rep(errorlist,m),nrow = m, byrow = TRUE)
    }
    
    Pmatrix <- matrix(rep(0,m^4),nrow = m^2)
    for(i in 1:m){ # Current state
      for(j in 1:m){ # Current error state
        for(k in 1:m){ # Goal state
          for(l in 1:m){ # Goal error state
            Pmatrix[m*i-m+j,m*k-m+l] <- Pgate[GroupStructure[i,k]]*Perror[GroupStructure[i,k],ErrorGroupStructure[GroupStructure[i,k],j,l]]
          }
        }
      }
    }
    
    Pmatrix
    
  })
  
  MyQmatrix <- reactive({
    
    input$letsgo
    # Get all variables
    Group <- isolate(input$Group)
    GateChoice <- isolate(input$GateChoice)
    UniformError <- isolate(input$UniformError)
    ErrorProb <- isolate(input$ErrorProb)
    TotErrorProb <- isolate(input$TotErrorProb)
    
    if(Group == "Pauli"){
      GroupStructure <- SetupList()[[1]]
      ErrorGroupStructure <- SetupList()[[2]]
      m <- 4
    }else{
      GroupStructure <- SetupList()[[3]]
      ErrorGroupStructure <- SetupList()[[4]]
      m <- 24
    }
    
    if(UniformError){
      r <- eval(parse(text=TotErrorProb))
      Perror <- matrix(rep( c(1-r,rep(r/(m-1),m-1)),m),nrow = m, byrow = TRUE)
    }else{
      errorlist <- unlist(lapply(unlist(strsplit(ErrorProb,",")),function(x) eval(parse(text=x))))
      Perror <- matrix(rep(errorlist,m),nrow = m, byrow = TRUE)
    }
    
    Gate <- switch(GateChoice,
                   "I" = 1,
                   "X" = 2,
                   "Y" = 3,
                   "Z" = 4,
                   "1" = 1,
                   "2" = 2,
                   "3" = 3,
                   "4" = 4,
                   "5" = 5,
                   "6" = 6,
                   "7" = 7,
                   "8" = 8,
                   "9" = 9,
                   "10" = 10,
                   "11" = 11,
                   "12" = 12,
                   "13" = 13,
                   "14" = 14,
                   "15" = 15,
                   "16" = 16,
                   "17" = 17,
                   "18" = 18,
                   "19" = 19,
                   "20" = 20,
                   "21" = 21,
                   "22" = 22,
                   "23" = 23,
                   "24" = 24
                   )
    
    Qmatrix <- matrix(rep(0,m^2),nrow=m)
    for(i in 1:m){ # Current error state
      for(j in 1:m){ # Goal error state
        Qmatrix[i,j] <- Perror[Gate,ErrorGroupStructure[Gate,i,j]]
      }
    }
    
    Qmatrix
    
  })
  
  DisplayMatrix <- reactive({
    input$letsgo
    PorQchoice <- isolate(input$PorQChoice)
    Group <- isolate(input$Group)
    if(PorQchoice=="P"){
      if(Group=="Pauli"){
        MyPmatrix()
      }else{
        MyPmatrix()[1:12,1:12]
      }
    }else{
      if(Group=="Pauli"){
        MyQmatrix() 
      }else{
        MyQmatrix()[1:6,1:6]
      }
    }
  })
  
  FullMatrix <- reactive({
    input$letsgo
    PorQchoice <- isolate(input$PorQChoice)
    if(PorQchoice=="P"){
      MyPmatrix()
    }else{
      MyQmatrix()
    }
  })
  
  output$caption <- renderText({
    input$letsgo
    PorQchoice <- isolate(input$PorQChoice)
    Group <- isolate(input$Group)
    GateChoice <- isolate(input$GateChoice)
    if(PorQchoice=="P"){
      if(Group=="Pauli"){
        paste("The transition matrix P for the process Z_t = (X_t,Y_t) for the Pauli group and the given probabilities")
      }else{
        paste("The first 12x12 elements of the transition matrix P for the process Z_t = (X_t,Y_t) for the Clifford group and the given probabilities. WARNING: Not the full matrix")
      }
    }else{
      if(Group=="Pauli"){
        paste("The transition matrix Q_t for the process Y_t when gate",GateChoice,"is applied")
      }else{
        paste("The first 6x6 elements of the transition matrix Q_t for the process Y_t when gate",GateChoice,"is applied. WARNING: Not the full matrix")
      }
    }
  })
  
  output$mytable <- renderTable({
    DisplayMatrix()
  }, digits = 5,colnames = FALSE)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Matrix.csv")
    },
    content = function(file) {
      MyDf <- as.data.frame(FullMatrix())
      names(MyDf) <- c()
      write.csv(MyDf, file, row.names = FALSE)
    }
  )
  
  output$readme <- downloadHandler(
    filename = "README.pdf",
    content = function(file) {
      file.copy("www/readmefile.pdf", file)
    }
  )

}

# Run the application 
shinyApp(ui = ui, server = server)