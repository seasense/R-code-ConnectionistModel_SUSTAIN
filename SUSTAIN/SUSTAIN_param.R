library(GA)

setwd("/Users/lewismarsh/Downloads/SUSTAIN")
setwd("/Users/pseiti/Dropbox/Projekte/CEITER/Internships/Intern_LM_July18/SUSTAIN/genetic_algorithms")
source("SUSTAIN_fx.R")


data_set <- NA
file <- "/Users/lewismarsh/Downloads/SUSTAIN/best_sols.txt"
file <- "/Users/pseiti/Dropbox/Projekte/CEITER/Internships/Intern_LM_July18/SUSTAIN/genetic_algorithms.txt"

shepard_structure_I <- list(c(1,1,1,1),c(1,1,2,1),c(1,2,1,1),c(1,2,2,1),c(2,1,1,2),c(2,1,2,2),c(2,2,1,2),c(2,2,2,2))
nAgents <- 100
stimulus_structure <- shepard_structure_I


set_data_set <- function(data){
  data_set <- data
}

run_SUSTAIN <- function(par){
  labelQuery <- TRUE
  m <- 4
  V <- c(2,2,2,2)
  
  structure_name <- "Type1"
  feature_coding <- c(0,1)
  
  ResponseSet <- c("A","B")
  agents <- paste("A",c(1:nAgents),sep="")
  termCrit <- 4
  nBlocks <- 32
  
  run <- sapply(agents, agent_fx, params=par, nBlocks=nBlocks, termCrit=termCrit,
                labelQuery=labelQuery, ResponseSet=ResponseSet,
                structure_name=structure_name, stimulus_structure=stimulus_structure, m=m, V=V,
                feature_coding=feature_coding)
  
  errorProb_perAgentAndBlock <- run["probOfError_block",]
  errorProb_perAgentAndBlock <- lapply(errorProb_perAgentAndBlock,function(x){c(x,rep(NA,(nBlocks -length(x))))})
  errorProb_perAgentAndBlock <- as.data.frame(errorProb_perAgentAndBlock)
  rownames(errorProb_perAgentAndBlock) <- paste("Block",c(1:nrow(errorProb_perAgentAndBlock)),sep="")
  
  errorProb_perAgentAndBlock[is.na(errorProb_perAgentAndBlock)] <- 0
  averageErrorProb_perBlock <- rowMeans(errorProb_perAgentAndBlock,na.rm=T)
  
  averageErrorProb_perBlock
}

RMSD_fitness <- function(x){
  beta <- x[1]   ###
  r <- x[2]
  d <- x[3]
  eta <- x[4]
  
  par <- c(r, beta, d, eta)
  names(par) <- c("r", "beta", "d", "eta")
  
  data <- run_SUSTAIN(par) #runs SUSTAIN on given model with paramenters
  
  data_points <- min(c(length(data), length(data_set)))
  sum <- 0
  
  for(i in 1:data_points){ # calculates RMSD^2
    sum <- sum + (data_points[[i]] - data[[i]])^2
  }
  
  sum <- -sqrt(sum) # Returns negative as fitness function gets maximised
  sum
}

run_genetic_algorithm <- function(popSize = 10, maxGenerations = 1000, survivalRate = 0.1, saveSols = T){
  if(saveSols){
    write("beta\t r\t\t d\t eta", file = file)
  }
  
  GA <- ga(type = "real-valued", 
           fitness =  RMSD_fitness,
           min = c(0, 0, 0, 0), max = c(14, 8, 30, 1), #Use depricated min and max if lower and upper don't work
           popSize = popSize,
           maxiter = maxGenerations,
           maxFitness = -0.01,
           keepBest = saveSols,
           postFitness = write_to_file,   # Writes best data point to a file after each generation if saveSols = T
           elitism = base::max(1, round(popSize*survivalRate)))
  
  names(GA$solutions) <- c("beta", "r", "d", "eta")
    
  GA  # call GA$solutions to retrieve parameters with best fitness (beta, r, d, eta)
}

write_to_file <- function(ga_object){
  suppressWarnings(
    if(!is.na(ga_object@bestSol)){
      index <- max(which(!is.null(ga_object@bestSol)))
      write(ga_object@bestSol[[index]], file = file, append = T)
    }
  )
  
  ga_object
}

