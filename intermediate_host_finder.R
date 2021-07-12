#' Parameterize and graph a SIR model of a zoonosis with an intermediate host
#' 
#' The wrapper function that simulates the progression of the epidemic through each
#' intermediate host species as listed in an external .xlsx document.
#' 
#' @param disease The pathogen of interest
#' @return species_name The intermediate host species that produced the most severe epidemic
#' @export 

#necessary packages et al.
#rm(list = ls(all.names = TRUE)) #clears environment
#setwd('../katherineroyce/Desktop/new research stuffs/')
# library("readxl") #for loading excel docs
# library("deSolve")

# source("r code/IH_sir.R")
# source("r code/jacobian.R")
# source("r code/nextgen.R")

#valid inputs, based on my file scheme: MERS, SARS, COVID
intermediate_host_finder <- function(disease){
  
  ## Set background birth/death parameters for all compartments
  ## all parameters are per day unless otherwise noted
  b_w <- 0.0625             # wild birth rate
  m_w <- 0.0625             # wild background death rate
  b_d <- 1/17               #birth and mortality rates are at replacement 
  m_d <- 1/17               #individuals/day
  b_h <- 0.0118             # human birth rate
  m_h <- 0.0118             # human background death rate
  
  # naive susceptibles (used in calculating R0 later)
  S_w_naive <- b_w/m_w;
  S_d_naive <- b_d/m_d;
  S_h_naive <- b_h/m_h;

    if (disease == 'MERS'){
      #MERS parameters
      beta_w <- 0.08;            # MERS is NOT highly transmissible in reservoir, but mild illness
      gamma_w <- 0.014;
      gamma_h <- 0.035;          # 1/day (28.4 day recovery period)
    } else {
    #SARS and COVID
      beta_w <- 0.25;            # wild transmission rate
      gamma_w <- 0.14;           # wild recovery rate; SARS R0_w = 1.78
      gamma_h <- 0.035;          # human recovery rate
    }

  ## Read in parameters for intermediate host candidates
  filename = paste("./data/",disease,"_species_params.xlsx",sep='')
  T <- read_excel(filename)
  
  ## Create space to store future info--in R you may not need this lol
  h <- nrow(T); #number of potential species
  newcol <- integer(h) 
  T <- cbind(T,Beta_h = newcol, Final_Ih = newcol, Time_To_Infection = newcol, Max_Ih = newcol, Time_To_Max = newcol, R0 = newcol, Rank = newcol)
 
  ## Seed initial infection in wild compartment
  init_Sw <- 0.9; # seed initial infection in reservoir at equilibrium Sw
  init_Iw <- 1-init_Sw; # (same for all IH simulations)

  ## Test each IH species
  for (i in c(1:h)) {
    
    ## Set IH species parameters using provided information
  
    #biological--measures biological factors that make IH similar to other species
    reservoir_similarity <- T[i, "res_sim"] # unitless
    human_similarity <- T[i, "hum_sim"] # unitless
  
    #contact--measures how frequently IH comes into contact with other species
    reservoir_contact <- T[i, "res_con"] #contact/day
    human_contact <- T[i, "hum_con"] #contact/day
  
    #transmission risk--measures if average contact between IH and other species is likely to transmit virus
    reservoir_risk <- T[i, "res_risk"] # infections/contact
    human_risk <- T[i, "hum_risk"]   # infections/contact
  
    ## Set model parameters as function of IH species parameters
    beta_dw <- reservoir_similarity*beta_w;         # transmission rate of wildtype in IH
    gamma_dw <- gamma_w;                             # recovery rate from wildtype in IH, 1/day
    beta_dm <- beta_w;                               # transmission rate of mutant strain in IH
    gamma_dm <- gamma_w;                             # recovery rate from mutant in IH
    wild_transmission <- reservoir_contact*reservoir_risk;  # rate of spillover wild -> intermediate (infections/day)
    human_transmission <- human_contact*human_risk;         # rate of spillover intermediate -> humans
    mu <- human_similarity*human_contact*human_risk;       # mutation rate (infections/day)
    beta_h <- human_similarity*beta_dw;                     # transmission rate in humans
  
    ## Run simulation
  
    # length of simulation
    time <- c(0:3000)
  
    # initial population values (proportions in each compartment) 
    population_vector <- c(init_Sw,init_Iw,0,    1,0,0,0,     1,0,0) 
    
    # collect all parameters!
    parms <- c(b_w, m_w, beta_w, gamma_w, b_d, m_d, beta_dw, gamma_dw, beta_dm, gamma_dm, wild_transmission, b_h, m_h, beta_h, gamma_h, human_transmission, mu)
    
    # wrapper function for ODE solver
    sirWrapper <- function(t,x, parms){IH_sir(t,x,parms)}
  
    # solves system of ODEs (using ode45)
    #ode(population_vector, time, sirWrapper, parms, method = "ode45");
    x <- ode(population_vector, time, sirWrapper, parms, method = "ode45"); 

    # in Matlab, stored simulation in a cell--maybe feasible here?
  
    ## Start filling in IH species-specific epidemic data based on simulation
    end <- nrow(x)
    S_w <- x[end, 2] #number is INDEX of column, not column name omg
    I_w <- x[end, 3]
    S_d <- x[end, 5]
    I_d <- x[end, 6]
    T_d <- x[end, 7]
    S_h <- x[end, 9]
    I_h <- x[end, 10]
     
    T[i,8] <- beta_h; #record beta_h
    T[i,9] <- I_h;    #record equilibrium I_h
  
    #index where #of humans susceptible drops below 0.99 for first time
    #(assuming spillover into population of at least 100 individuals, this assures one 'real' individual has it)
    init_time <- which(x[,9]<=0.99)[1]
  
    max_time <- which.max(x[,10])
    max_Ih <- max(x[,10])
    T[i,11] <- max_Ih; # max proportion of humans infected
  
    # record time to max
    if (max_Ih > 10^(-2)){ # the epidemic establishes itself in humans
      T[i,10] <- init_time;
      T[i,12] <- max_time; # time to max proportion of humans infected
    } else { # no significant foothold in humans, so indicate with non-integer timestep
      T[i,10] <- 0.5;
      init_time <- 0.5;
      max_time <- 0.5;
      T[i,12] <- 0.5;
    }

    # calculate stability of equilibria
    # disease-free equilibria
    dfe <- jacobian(m_w, beta_w, gamma_w, m_d, beta_dw, gamma_dw, beta_dm, gamma_dm, m_h, beta_h, gamma_h, wild_transmission, human_transmission, mu, 1,0,  1,0,0,  1,0)
    dfeTest <- which(Re(dfe)<0) #stable if all eigenvalues have real part < 0 ie this array is all 1s

    # endemic equilibria 
    ee <- jacobian(m_w, beta_w, gamma_w, m_d, beta_dw, gamma_dw, beta_dm, gamma_dm, m_h, beta_h, gamma_h, wild_transmission, human_transmission, mu, S_w,I_w,    S_d,I_d,T_d,     S_h,I_h)
    eeTest <- which(Re(ee)<0) #stable if all eigenvalues have real part < 0 ie this array is all 1s

    # calculate R_0
    A <- nextgen(beta_w, gamma_w, m_w, beta_dw, gamma_dw, beta_dm, gamma_dm, m_d, mu, beta_h, gamma_h, m_h, wild_transmission, human_transmission, S_w_naive, S_d_naive, S_h_naive)     
    global_r0 <- max(eigen(A)$values)
    T[i,13] <- global_r0; #store R0
    T[i,14] <- (global_r0+max_Ih*100)/2; #rank is average of global r0 and max percent infected humans

    } 
    
    ## After simulating epidemic spreading through all potential IH species, find the most potent
    maximum <- max(T$Rank)
    idx <- which.max(T$Rank)
    species_name <- T$Species[idx]
    
    ## save data on all potential IH species
    saved_file <- paste("./data/",disease,"_species_epidemics.RData",sep='')
    save(T, file = saved_file) #figure out how to rank the data frame in descending order lol
    ## Report on all IH species that produced epidemics above certain threshold
    #rows <- (T$Max_Humans_Infected > .3 | T$R0 > 1.5)
    #cols <- c(1,9:14)
    #T(rows,cols)
    return(species_name)
    }

