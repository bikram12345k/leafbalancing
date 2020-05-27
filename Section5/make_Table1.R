#  This code computes the numbers in Table 1 
#	based on run results from 
#	run_Table1_normal and run_Table1_laplace
## 	
res_tab <- read.csv('sim_10k_laplace_6.csv', header=FALSE)
res_tab <- res_tab[,-1]	

res_tab1 <- res_tab[res_tab[,74]>= n*.00,,drop=FALSE]
taus = c(0, .5, 0.75, 1, 1.25, 1.5)
between <- function(x, y)
 (y[1] <= x) & (x<=y[2])
 
table <- list()

for(j in 1:length(taus)){
	table[[j]] = rbind(
		sapply(1:4, function(k) mean(res_tab1[, 12*(j-1)+(k-1)*3 + 1], na.rm=TRUE) - taus[j]),
		sapply(1:4, function(k) sqrt(mean((res_tab1[, 12*(j-1)+(k-1)*3 + 1] - taus[j])^2, na.rm=TRUE))),
		sapply(1:4, function(k) mean(apply(res_tab1[, 12*(j-1)+(k-1)*3 + 2:3], 1, function(y) between(taus[j],y) ), na.rm=TRUE)),
		sapply(1:4, function(k) 1-mean(apply(res_tab1[, 12*(j-1)+(k-1)*3 + 2:3], 1, function(y) between(0,y) ), na.rm=TRUE)))
	rownames(table[[j]]) = c('bias', 'rmse', 'coverage', 'power')
	colnames(table[[j]]) = c("Prop", "Leaf", "PairCaliper", "LeafPair")
}

names(table) = paste('tau = ', taus)
table

## Columns "Prop", "Leaf", "PairCaliper", "LeafPair"
## rows bias, rmse, coverage, power