# remove unsignificant edges from plotting

# function to remove unsignificant edges for plotting
clean.semplot <- function(sem_fit_object, pvalue_cutoff=0.05){
  obj <- semPlot:::semPlotModel(sem_fit_object)
  original_Pars <- obj@Pars
  check_Pars <- obj@Pars %>% dplyr::filter(!(edge %in% c("int","<->") | lhs == rhs))
  keep_Pars <- obj@Pars %>% dplyr::filter(edge %in% c("int","<->") | lhs == rhs) # list of parameters to keep untouched
  test_against <- lavaan::standardizedSolution(sem_fit_object) %>% dplyr::filter(pvalue < pvalue_cutoff, rhs != lhs)
  test_against_rev <- test_against %>% 
    rename(rhs2 = lhs,  lhs = rhs) %>%  # rhs and lhs are reversed in the standardizedSolution() output
    rename(rhs = rhs2)
  checked_Pars <- check_Pars %>% 
    semi_join(test_against, by = c("lhs", "rhs")) %>% 
    bind_rows(
      check_Pars %>% semi_join(test_against_rev, by = c("lhs", "rhs"))
    )
  obj@Pars <- keep_Pars %>% bind_rows(checked_Pars)
  
  #verify by looking at the list of the edges removed from the object
  print(anti_join(original_Pars,obj@Pars))
  return(obj)
}



# function to relabel edges to include thresholds of significance
sem_sig_labels <- function(sem_fit_object){
  table2<-parameterEstimates(sem_fit_object,standardized=TRUE)[!is.na(parameterEstimates(sem_fit_object)$pvalue) & 
                                                            parameterEstimates(sem_fit_object)$op!=":=",]
  sig <- rep(" ",nrow(table2))
  sig[table2$pvalue<=0.05] <- "*"
  sig[table2$pvalue<=0.01] <- "**"
  sig[table2$pvalue<=0.001] <- "***"
  b <- paste0(round(table2$std.all,2), sig)
  return(b)
}
