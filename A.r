

data <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/v6.csv", h = T)

mods <- c("cf.type","genre","Age","ed.level")

A <- setNames(lapply(seq_along(mods), function(i) table(data[[mods[i]]], dnn = NULL)), mods)

target <- sapply(seq_along(A), function(i) any(A[[i]] <= 4))
A <- A[target]
                 
low <- setNames(lapply(seq_along(A), function(i) A[[i]][which(A[[i]] <= 4)]), names(A))

lst <- Filter(length, lapply(split(data[names(low)], data$study.name), 
                              function(dat) Filter(nrow, Map(function(x, y) 
                                merge(x, y[setdiff(names(y), "values")], by = "ind"), lapply(dat, 
                                function(x) stack(table(x))), lapply(low, stack)))))

## Desired Output:                                                                                             
do.call(rbind, c(Map(cbind, study.name = names(lst), lapply(lst,                              
                 function(x) do.call(rbind, c(Map(cbind, x, mod.name = names(x)),
                 make.row.names = FALSE)))), make.row.names = FALSE))
