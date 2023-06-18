#setwd("/Users/liux/Documents/Ken_analysis/data/")
setwd("/home/katagirf/liux4215/Documents/Ken_analysis/data")
load("col_count_data.Rdata")
#source("/Users/liux/Documents/fit_function_fumi.R")
source("/home/katagirf/liux4215/Documents/fit_function_fumi.R")
all.genes <- rownames(col_count_data)
treatment <- c()
time <- c()
for (mylib in colnames(col_count_data)){
  treatment <- c(treatment, strsplit(mylib,split = "_")[[1]][4])
  time <- c(time, strsplit(mylib,split = "_")[[1]][5])
}
time <- as.numeric(substr(time,start=1,stop=2))
after_3_indexes <- which(time >= 3)
time_after3 <- time[after_3_indexes]
treatment_after3 <- treatment[after_3_indexes]
col_offset_after3 <- col_offset[after_3_indexes]
offset_mock_after3 <- col_offset_after3[treatment_after3 == "mock"]
offset_Rpt2_after3 <- col_offset_after3[treatment_after3 == "AvrRpt2"]
offset_EV_after3 <- col_offset_after3[treatment_after3 == "EV"]
myoffset_list <- list(Rpt2 = offset_Rpt2_after3,mock = offset_mock_after3)

load("mock_start.vals.Rdata")
load("AvrRpt2_genes.Rdata")
rownames(mock_coef) <- all.genes
k_vec = c(1.8, 0.5487, 1.1733, 1.1253, 0.6492, 0.6853, 0.6957, 0.6987, 0.4949,
          0.5065, 0.513, 0.5174, 0.2133)
#k_vec = c(1.8, 2.4573, 1.0775, 1.3624, 0.6624, 0.7903, 0.646, 0.6757, 0.5669, 0.5814,
#          0.5003, 0.5094)
k_upper_bound = 10
a_upper_bound = k_upper_bound * 50
a_lower_bound = -k_upper_bound * 10
#compart_vec <- c(1,3,5,7,9,11)
args = commandArgs(trailingOnly=TRUE)
gene_index_start <- as.numeric(args[1])
gene_index_end <- as.numeric(args[2])
input_lognormal <- args[3]
result_name <- args[4]
load(input_lognormal)
start_mat <- myresult
rm(myresult)

cores <- 8
start_time <- date()
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
myresult <- foreach(m = 1:cores, .combine = "rbind", .packages = c("functional","Rmpfr","pracma")) %dopar%
  {
    results <- c()
    mycompartment <- m + 10
    compart1<- compartment_function.add(mycompartment)[1]
    compart2 <- compartment_function.add(mycompartment)[2]
    func1_str <- func_string(compart1)
    func1 <- eval(parse(text =func1_str))
    func2_str <- func_string(compart2)
    func2 <- eval(parse(text =func2_str))
    
    myindex <- 1
    for (mygene in AvrRpt2_mock_positive_genes[gene_index_start:gene_index_end]){
      y_after3 <- as.numeric(col_count_data[mygene,after_3_indexes])
      y_Rpt2 <- y_after3[treatment_after3 == "AvrRpt2"]
      y_mock <- y_after3[treatment_after3 == "mock"]
      mycount_list <- list(Rpt2 = y_Rpt2,mock = y_mock)
      
      x <- as.numeric(time_after3[treatment_after3 == "mock"])
      time_list <- list(Rpt2 = x,mock = x)
      #x is the same across three treatments
      
      poly_start = mock_coef[mygene,1:6]
      
      compart_index <- intersect(which(start_mat[,15] == mycompartment), which(start_mat[,16] == myindex))
      compart_start <- start_mat[compart_index,c(7,8,9)]
      
      mytheta = mock_coef[mygene,7]
      likelihood_function_curry <- Curry(likelihood_function_Rpt2, y_Rpt2 = y_Rpt2,
                                         time = x, offset_Rpt2 = offset_Rpt2_after3, 
                                         theta = mytheta, poly_para = poly_start,
                                         k_vec = k_vec)
      
      fit <- tryCatch(optim(likelihood_function_curry,par =  compart_start, lower = c(a_lower_bound,0,0.1),
                            method = "L-BFGS-B", upper = c(a_upper_bound,a_upper_bound,k_upper_bound)),
                      error = function(err){return(1)})
      
      if (inherits(fit,"list")){
        results <- rbind(results, c(fit$par,poly_start,fit$convergence,fit$value,mycompartment,myindex))
      }else{
        results <- rbind(results, c(rep(NA,3),poly_start,NA,NA,mycompartment,myindex))
      }
      myindex <- myindex + 1
    }
    results
  }
end_time <- date()
stopImplicitCluster()
stopCluster(cl)
save(myresult,start_time,end_time,file = result_name)


