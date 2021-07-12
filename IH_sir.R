#' IH mutation model
#' 
#' The 10 differential equations that describe the IH mutation model
#' 
#' @param t Timestep
#' @param x Input vector
#' @param parms List of disease-specific parameters
#' @return species_name The intermediate host species that produced the most severe epidemic
#' @export 

IH_sir <- function (t, x, parms){
  b_w <- parms[1]
  m_w <- parms[2]
  beta_w <- parms[3]
  gamma_w <- parms[4] 
  b_d <- parms[5]
  m_d <- parms[6]
  beta_dw <- parms[7]
  gamma_dw <- parms[8]
  beta_dm <- parms[9]
  gamma_dm <- parms[10]
  p_d <- parms[11]
  b_h <- parms[12]
  m_h<- parms[13]
  beta_h <- parms[14]
  gamma_h <-parms[15]
  p_h <- parms[16]
  mu <- parms[17]
  
  #output vector
  y<-numeric(10)

  # traditional model in wild compartment
  y[1] <- b_w-beta_w*x[1]*x[2]-m_w*x[1]                                #S_r
  y[2] <- beta_w*x[1]*x[2]-gamma_w*x[2]-m_w*x[2]                       #I_r
  y[3] <- gamma_w*x[2]-m_w*x[3]                                        #R_r
  
  # 2 types of infections in domestic animals: caught from wild animals
  # OR transmissible to humans
  y[4]<-b_d-beta_dw*x[4]*x[5]-p_d*x[2]*x[4]-beta_dm*x[4]*x[6]-m_d*x[4]  #S_i
  y[5]<-beta_dw*x[4]*x[5]+p_d*x[2]*x[4]-mu*x[5]-gamma_dw*x[5]-m_d*x[5]  #I_i 
  y[6]<-mu*x[5]+beta_dm*x[4]*x[6]-gamma_dm*x[6]-m_d*x[6]                #T_i
  y[7]<-gamma_dw*x[5]+gamma_dm*x[6]-m_d*x[7]                            #R_i
  
  # traditional model in human compartment seeded by type T_i
  y[8]<-b_h-beta_h*x[8]*x[9]-p_h*x[6]*x[8]-m_h*x[8]                   #S_h
  y[9]<-beta_h*x[8]*x[9]+p_h*x[6]*x[8]-gamma_h*x[9]-m_h*x[9]          #I_h
  y[10]<-gamma_h*x[9]-m_h*x[10]                                       #R_h
  return(list(y,parms)) 
}
