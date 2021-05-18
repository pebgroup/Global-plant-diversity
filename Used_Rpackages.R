# used R packages in code
library(NCmisc)

scripts <- dir("publish_scripts/")

funcs <- c()
for(i in 1:length(scripts)){
  print(scripts[i])
  print(list.functions.in.file(paste0("publish_scripts/", scripts[i])))
  func <- names(list.functions.in.file(paste0("publish_scripts/", scripts[i])))
  temp <- unlist(str_split(func, ","))
  temp <- str_extract(temp, ":[a-z|A-Z]*")
  funcs <- c(funcs, temp)
  #funcs <- c(funcs, names(list.functions.in.file(paste0("publish_scripts/", scripts[i])))) 
  
}
funcs <- gsub(":", "", funcs)
sort(unique(funcs))
