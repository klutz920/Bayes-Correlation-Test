#######################################################################
####  Dirichlet-Multinomial (DM) Model for Bayes Correlation Test #####
#######################################################################

## Required ##

## 1.  counts:  a matrix of count data (rows are samples, columns are taxa)
## 2.  group: a binary vector indicating the group for each sample.
## 3.  target: a numerical vector for the measurement of interest 

## Note: The group order must be the same for both counts and target. 

## (1) The function will return the posterior mean and lower and upper bounds of the 
## 95% credible interval for the Spearman correlation between counts and target 
## are sent to the global environment (GE) for both groups. 
## These are called "group1_summary" and "group2_summary. 

## (2) The storage containing the parameter, alpha, of the DM will be sent to the GE.
## This is called "alpha_store".

## (3) The console will print which group each one belongs to. 

BayesDMC = function(
  counts,
  group,
  target
){
  
  # reduced data likelihood for Metropolis Hastings ratio
  LDMdens = function(y_i, alpha_i){
    lgamma(sum(alpha_i))-lgamma( sum(y_i+alpha_i) )+
      sum( lgamma(y_i+alpha_i) - lgamma(alpha_i) )
  }
  
  # MCMC settings
  counts = as.matrix(counts)
  tau = 1 
  T = 10000    # mcmc iterations
  B = T/2      # burnins
  n = nrow(counts)
  n_taxa = ncol(counts)
  
  # store results for the DM parameter and sample acceptance rates
  alpha_store <- array(rep(NA, B*n*n_taxa), dim=c(B, n, n_taxa))
  acc_taxa = array(rep(NA, n*B), dim = c(n,B))
  
  # keep a count for process monitoring
  counter = 10
  
  # 1. Run the DM
  # algorithm runs each sample i over t iterations
  for(i in 1:n){
    
    #initial settings for each sample
    y_i = counts[i,]   # sample row vector
    init = rep(1,n_taxa) #same initial alpha values for each sample
    alpha = init
    
    # Monitor the Proces
    if(i*100/n >= counter){
      print(paste0(counter, "% has been completed."))
      counter = counter + 10
    }
    
    # do T iterations of MCMC for sample i
    for(t in 1:T){
      alpha_copy = alpha # copy of alpha[t-1]
      accept_taxa = 0
      
      # taxa-wise sequential update in sample i for alpha_ij at iteration t
      for(j in 1:n_taxa){
        # draw candidate for alpha_ij
        alpha_ij = exp(rnorm(1,mean=log(alpha[j]),sd=tau)) 
        
        #store in copy of alpha
        alpha_copy[j] = alpha_ij
        
        #calculate MH ratio/decision to keep alpha_ij or previous alpha value
        s = LDMdens(y_i=y_i, alpha_i = alpha_copy)-LDMdens(y_i=y_i, alpha_i = alpha)
        v = runif(1,0,1)
        if(s>=log(v)){
          #keep new alpha_ij
          alpha[j] = alpha_ij    # update alpha
          if(t>B){
            accept_taxa = accept_taxa + 1  #acceptance rate of entire sample per iteration
          }
        } else{
          #keep old alpha_ij
          alpha_copy[j] = alpha[j]
        }
        
      } 
      
      # store after burn-in results
      if(t>B){
      alpha_store[t-B,i,] = alpha
      }
      
      # record acceptance rate of alpha_star for sample i, iteration t
      if(t>B){
        acc_taxa[i,t-B] = accept_taxa/n_taxa
      }
    } 
  } 
  
  # 2. Calculate correlations after DM has completed
  # index the groups
  g = sort(unique(group))
  c1 = which(group == g[1])
  c2 = which(group == g[2])
  
  # counts separated by group
  taxa_1 = counts[c1,]
  taxa_2 = counts[c2,]
  n1 = nrow(taxa_1)
  n2 = nrow(taxa_2)
  taxa_name = colnames(counts)
  
  # covariate separated by group
  target_1 = target[c1]
  target_2 = target[c2]
  
  # thin the posterior samples
  thin = seq(from=10,to=B,by=10)
  t_n = length(thin)
  alpha = alpha_store[thin,,]
  
  # storage for correlation estimate at each after burn in iteration
  cor1_store = array(rep(NA, t_n*n_taxa), dim = c(t_n,n_taxa))
  cor2_store = array(rep(NA, t_n*n_taxa), dim = c(t_n,n_taxa))
  
  # calculate Spearman correlation at each after burn in iteration
  for(b in 1:t_n){
    for(c in 1:n_taxa){
    x_1 = alpha[b,c1,c]  #alpha's of iter b, group c1, taxa c
    x_2 = alpha[b,c2,c]  #alpha's of iter b, group c2, taxa c
    cor1_store[b,c] = cor(x_1, target_1, method = "spearman")
    cor2_store[b,c] = cor(x_2, target_2, method = "spearman")
    }
  }
  
  # storage for posterior correlation results
  group1_summary = array(rep(NA, n_taxa*3), dim = c(n_taxa, 3))
  group2_summary = array(rep(NA, n_taxa*3), dim = c(n_taxa, 3))
  rownames(group1_summary) = taxa_name
  rownames(group2_summary) = taxa_name
  colnames(group1_summary) = c("Lower", "Mean", "Upper")
  colnames(group2_summary) = c("Lower", "Mean", "Upper")
  
  #posterior summaries
  group1_summary[,1] = apply(cor1_store, 2, quantile, 0.025)
  group1_summary[,2] = apply(cor1_store, 2, mean)
  group1_summary[,3] = apply(cor1_store, 2, quantile, 0.975)
  
  group2_summary[,1] = apply(cor2_store, 2, quantile, 0.025)
  group2_summary[,2] = apply(cor2_store, 2, mean)
  group2_summary[,3] = apply(cor2_store, 2, quantile, 0.975)
  
  # 3. output results to global environment
  name1 = paste0("group1_summary")
  assign(name1,group1_summary, envir = .GlobalEnv) 
  name2 = paste0("group2_summary")
  assign(name2,group2_summary, envir = .GlobalEnv) 
  name3 = paste0("alpha")
  assign(name3, alpha_store, envir = .GlobalEnv)
  
  print(paste0("group1_summary corresponds to the group labeled as: ", g[1]))
  print(paste0("group2_summary corresponds to the group labeled as: ", g[2]))
} 
