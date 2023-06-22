# DEG analysis of de novo transcriptome assembly from Trinity

(Example)
```R
install.packages("knitr")
library(TCC)
library(ROC)
library(plotly)
library(magrittr)
in_f <- "ctrl_fi_fltrexpr_contig_cnt.csv"    # summary of contig count from Trinity
out_f1 <- "DEG_profile.csv"                  # output DEG profile
out_f2 <- "DEG_profile.png"                  # output DEG M-A plot
param_G1 <- 4                          # number of sample in group 1
param_G2 <- 4                          # number of sample in group 2
param_FDR <- 0.1                      # false discovery rate (FDR)
param_fig <- c(1200, 1140)               # plot size

data <- read.table(in_f, 
                   header=TRUE, 
                   row.names=1, 
                   sep=",", 
                   quote="")

knitr::kable(head(data, n=10))# check input file
### ctrl_fi_fltrexpr_contig_cnt.csv
|                       | ctrl_rep1| ctrl_rep2| ctrl_rep3| ctrl_rep4| fi_rep1| fi_rep2| fi_rep3| fi_rep4|
|:----------------------|---------:|---------:|---------:|---------:|-------:|-------:|-------:|-------:|
|TRINITY_DN0_c0_g1_i1   |   1299.73|   1253.00|   1638.00|   1549.00| 1291.00| 1260.00| 1100.00| 1243.00|
|TRINITY_DN0_c0_g2_i1   |    320.00|    328.00|    675.00|    666.61|  368.00|  398.00|  402.00|  415.00|
|TRINITY_DN0_c12_g1_i2  |      4.00|     15.00|     10.00|      5.00|   15.00|   20.00|    9.00|   10.00|
|TRINITY_DN0_c12_g2_i1  |      8.87|      3.32|      5.23|      8.29|    2.06|    4.99|    6.03|    0.00|
|TRINITY_DN0_c14_g1_i1  |    293.93|    296.93|    240.96|    258.04|  250.97|  228.15|  174.09|  186.54|
|TRINITY_DN0_c163_g1_i1 |      0.00|      1.04|      0.00|      3.39|    2.00|    1.00|    1.00|    0.00|
|TRINITY_DN0_c174_g1_i1 |     67.00|     62.00|     92.00|     84.00|   61.00|   71.00|   44.00|   51.00|
|TRINITY_DN0_c19_g1_i1  |     37.36|     48.53|     62.29|     71.06|   51.34|   55.01|   51.37|   76.69|
|TRINITY_DN0_c19_g2_i3  |      9.21|      6.00|     10.00|      6.00|    8.00|   13.00|    5.00|   10.94|
|TRINITY_DN0_c19_g3_i1  |     51.00|     61.00|    150.70|    128.00|   84.67|   88.00|   70.00|   87.00|

###
data.cl <- c(rep(1, param_G1), rep(2, param_G2))
data.cl
tcc <- new("TCC", data, data.cl)     
tcc <- calcNormFactors(tcc, 
                       norm.method="tmm",
                       test.method="edger",
                       iteration=30, 
                       FDR=0.1, 
                       floorPDEG=0.05)
tcc
tcc <- estimateDE(tcc, 
                  test.method="edger",
                  FDR=param_FDR)
tcc
result <- getResult(tcc, sort=FALSE)
knitr::kable(head(result, n=10)) 
###
tmp <- cbind(rownames(tcc$count), tcc$count, result)
knitr::kable(head(tmp, n=10))
write.csv(tmp, out_f1)
# (Example)
# rownames(tcc$count)	ctrl_rep1	ctrl_rep2	ctrl_rep3	ctrl_rep4	fi_rep1	fi_rep2	fi_rep3	fi_rep4	gene_id	a.value	m.value	p.value	q.value	rank	estimatedDEG
# TRINITY_DN7038_c0_g1_i5	TRINITY_DN7038_c0_g1_i5	757	708	573	639	646	818	594	546	TRINITY_DN7038_c0_g1_i5	9.41934178077124	0.370386035284893	0.245918609070331	0.435771919592758	147160	0
###

sum(tcc$stat$q.value < param_FDR)      # number of DEG

# Method 1
png(out_f2, 
	pointsize=3, 
	width=param_fig[1], 
	height=param_fig[2],
	res=350)

# M-A plot
plot.TCC(tcc, FDR=param_FDR)
legend("topright", 
       c(paste("DEG(FDR<", param_FDR, ")", sep=""), 
         "non-DEG"),           
       col=c("magenta", "black"), 
       pch=20)
dev.off()
quit()
```