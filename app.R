#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(dplyr)
library(readr)

dagger <- function(X){
  t(Conj(X))
}

compareMatrix <- function(a,b){
  isTRUE(all.equal(a,b))
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

# Define UI for application
ui <- fluidPage(theme = shinytheme("lumen"),
   
   # Application title
   titlePanel("Transition Matrix for Qbit error accumulation"),
   
   # Sidebar Input
   sidebarLayout(
      sidebarPanel(
        selectInput(inputId = "Group", label=strong("Group"),
                    choices = c("Pauli","Clifford"),
                    selected = "Pauli"
        ),
        
        checkboxInput(inputId = "ManualPsi", label = strong("Manual Psi input")
        ),
        conditionalPanel(condition = "input.ManualPsi == true",
                    textInput("Psi", 'Enter a Psi as vector (comma delimited)', "1.2,2+1i")
        ),
        checkboxInput(inputId = "UniformGate", label = strong("Uniform gate probabilites"),TRUE
        ),
        conditionalPanel(condition = "input.UniformGate == false",
                         textInput("GateProb", 'Enter a vector of gate probabilities (comma delimited)', "")
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
        downloadButton("downloadData", "Download as csv")
        
      ),
      
      # Mainpanel with outcome
      mainPanel(
        h4("Look at this matrix"),
        tableOutput("mytable")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  

  MyTable <- reactive({
    
    input$letsgo # Trigger a update
    
    # Get all variables
    Group <- isolate(input$Group)
    ManualPsi <- isolate(input$ManualPsi)
    Psi <- isolate(input$Psi)
    UniformGate <- isolate(input$UniformGate)
    GateProb <- isolate(input$GateProb)
    UniformError <- isolate(input$UniformError)
    ErrorProb <- isolate(input$ErrorProb)
    TotErrorProb <- isolate(input$TotErrorProb)
    
    
    if(ManualPsi){
      phi <- Psi
    }else{
      phi <-c(complex(real=runif(1),imaginary=runif(1)),complex(real=runif(1),imaginary=runif(1))) # Vector of two complex numbers
    }
    phi <- phi/sqrt(sum(Mod(phi)^2)) # Rescale
    rho0 <- phi %*% dagger(phi) # Create rho0
    
    if(Group == "Pauli"){
      MatrixList <- list(I,X,Y,Z) # Set of gates
    }else{
      MatrixList <- list(ClifM01,ClifM02,ClifM03,ClifM04,ClifM05,ClifM06,ClifM07,ClifM08,ClifM09,ClifM10,ClifM11,ClifM12,
                         ClifM13,ClifM14,ClifM15,ClifM16,ClifM17,ClifM18,ClifM19,ClifM20,ClifM21,ClifM22,ClifM23,ClifM24)
    }
    stateSpace <- lapply(MatrixList, function(x) x%*%rho0%*%dagger(x)) # Compute all states !! check for double states !! breath/depth search?
    n <- length(stateSpace)
    m <- length(MatrixList)
    
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
    
    # Step 3: determine transition matrix for non-errors
    GroupStructure <- matrix(rep(0,n^2),nrow = n)
    for(i in 1:n){ # Current state
      for(j in 1:n){ # Goal state
        for(k in 1:m){ # Gate
          if(compareMatrix(MatrixList[[k]] %*% stateSpace[[i]] %*% dagger(MatrixList[[k]]),stateSpace[[j]])){
            GroupStructure[i,j] <- k
          }
        }
      }
    }
    
    # Step 4: determine corresponding error for gate-couple
    ErrorGroupStructure <- array(rep(0,n^3),rep(n,3))
    for(i in 1:n){ # Current error state
      for(j in 1:n){ # Goal error state
        for(k in 1:m){ # Gate
          for(l in 1:m){ # Error gate
            if(compareMatrix(MatrixList[[l]] %*% MatrixList[[k]] %*% stateSpace[[i]] %*% dagger(MatrixList[[k]]) %*% dagger(MatrixList[[l]]), stateSpace[[j]])){
              ErrorGroupStructure[k,i,j] <- l # gate,current,goal
            }
          }
        }
      }
    }
    
    # Step 5: combine probabilities, compute P matrix
    Pmatrix <- matrix(rep(0,n^4),nrow = n^2)
    for(i in 1:n){ # Current state
      for(j in 1:n){ # Current error state
        for(k in 1:n){ # Goal state
          for(l in 1:n){ # Goal error state
            Pmatrix[m*i-m+j,m*k-m+l] <- Pgate[GroupStructure[i,k]]*Perror[GroupStructure[i,k],ErrorGroupStructure[GroupStructure[i,k],j,l]]
          }
        }
      }
    }
    
    # Step 6: enjoy
    Pmatrix
    
  })
  
  output$mytable <- renderTable({
    MyTable()[1:16,1:16]
  }, digits = 5,colnames = FALSE)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Pmatrix.csv")
    },
    content = function(file) {
      MyDf <- as.data.frame(MyTable())
      names(MyDf) <- c()
      write.csv(MyDf, file, row.names = FALSE)
    }
  )

}

# Run the application 
shinyApp(ui = ui, server = server)
