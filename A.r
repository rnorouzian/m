

data <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/v6.csv", h = T)

mods <- c("cf.type","genre","Age","ed.level")

A <- setNames(lapply(seq_along(mods), function(i) table(data[[mods[i]]], dnn = NULL)), mods)

target <- sapply(seq_along(A), function(i) any(A[[i]] <= 4))
A <- A[target]
low <- setNames(lapply(seq_along(A), function(i) A[[i]][which(A[[i]] <= 4)]), names(A))

# Notice `low` (above) has an entry for "cf.type" with `ind` == 15 and `values` == 2.

lst1 <- Filter(length, lapply(split(data[mods], data$study.name), 
                              function(dat) Filter(nrow, Map(merge, lapply(dat, 
                               function(x) stack(table(x))), lapply(low, stack)))))

# It doesn't show "cf.type" with `ind` == 15 and `values` == 2, which is in the `low` above?
do.call(rbind, c(Map(cbind, study.name = names(lst1), lapply(lst1, 
                 function(x) do.call(rbind, c(Map(cbind, x, mod.name = names(x)),
                 make.row.names = FALSE)))), make.row.names = FALSE))
