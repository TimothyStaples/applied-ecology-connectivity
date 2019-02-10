# create all possible two-way interactions from a vector of strings (representing covariates)
two.way.interactions<-function(variables){
  interactions.temp<-combn(unique(c(variables, variables)), 2)
  return(paste(interactions.temp[1,], interactions.temp[2,], sep=":"))
}