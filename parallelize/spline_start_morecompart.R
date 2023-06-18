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
myoffset_list <- list(Rpt2 = offset_Rpt2_after3,EV = offset_EV_after3,mock = offset_mock_after3)

load("mock_start.vals.Rdata")
load("AvrRpt2_genes.Rdata")
rownames(mock_coef) <- all.genes
k_vec = c(1.8, 0.5487, 1.1733, 1.1253, 0.6492, 0.6853, 0.6957, 0.6987, 0.4949,
          0.5065, 0.513, 0.5174, 0.2133)
#k_vec = c(1.8, 2.4573, 1.0775, 1.3624, 0.6624, 0.7903, 0.646, 0.6757, 0.5669, 0.5814,
#          0.5003, 0.5094)
k_upperbound = 1.5
#compart_vec <- c(1,3,5,7,9,11)
args = commandArgs(trailingOnly=TRUE)
gene_index_start <- as.numeric(args[1])
gene_index_end <- as.numeric(args[2])
result_name <- args[3]

cores <- 8
start_time <- date()
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
myresult <- foreach(m = 1:cores, .combine = "rbind", .packages = c("functional","Rmpfr","pracma")) %dopar%
  {
    results <- c()
    mycompartment <- m + 10
    compart1 <- compartment_function.add(mycompartment)[1]
    compart2 <- compartment_function.add(mycompartment)[2]
    func1_str <- func_string(compart1)
    func1 <- eval(parse(text =func1_str))
    func2_str <- func_string(compart2)
    func2 <- eval(parse(text =func2_str))
    
    myindex <- 1
    for (mygene in AvrRpt2_mock_positive_genes[gene_index_start:gene_index_end]){
      y_after3_log <- log(as.numeric(col_count_data[mygene,after_3_indexes]))
      y_Rpt2_log <- y_after3_log[treatment_after3 == "AvrRpt2"] - offset_Rpt2_after3
      y_EV_log <- y_after3_log[treatment_after3 == "EV"] - offset_EV_after3
      y_mock_log <- y_after3_log[treatment_after3 == "mock"] - offset_mock_after3
      
      x <- as.numeric(time_after3[treatment_after3 == "mock"])
      
      poly_start = mock_coef[mygene,1:6]
      mock_subtract <- poly_function(p=poly_start,x = x)
      
      Rpt2_subtract <- y_Rpt2_log - mock_subtract
      EV_subtract <- y_EV_log - mock_subtract
      
      Rpt2.sq.sp.fun = splinefun(sqrt(x), Rpt2_subtract)
      EV.sq.sp.fun = splinefun(sqrt(x), EV_subtract)

      Rpt2.sp.val = Rpt2.sq.sp.fun(sqrt(seq(3,24,0.8)))
      EV.sp.val = EV.sq.sp.fun(sqrt(seq(3,24,0.8)))
      
      Rpt2.spline_func_curry <- Curry(LS_function,x = seq(3,24,0.8),y = Rpt2.sp.val,k_vec = k_vec)
      EV.spline_func_curry <- Curry(LS_function,x = seq(3,24,0.8),y = EV.sp.val,k_vec = k_vec)
      
      EV.fit_spline <- tryCatch(optim(EV.spline_func_curry,par = c(3,3,0.3), lower = c(-10 * k_upperbound,0,0.1),
                                        method = "L-BFGS-B", upper = c(50 * k_upperbound ,50 * k_upperbound ,k_upperbound )),error = function(err){return(1)})
      
      Rpt2.fit_spline <- tryCatch(optim(Rpt2.spline_func_curry,par = c(3,3,0.3), lower = c(-10 * k_upperbound,0,0.1),
                                        method = "L-BFGS-B", upper = c(50 * k_upperbound,50 * k_upperbound,k_upperbound)),error = function(err){return(1)})
      
      if (inherits(Rpt2.fit_spline,"list")){
        Rpt2.start_spline <- Rpt2.fit_spline$par
      }else{
        Rpt2.start_spline <- c(3,3,0.3)
      } 
      if (inherits(EV.fit_spline,"list")){
        EV.start_spline <- EV.fit_spline$par
      }else{
        EV.start_spline <- c(3,3,0.3)
      }
      
      
      Rpt2.func_curry <- Curry(LS_function,x = x, y = Rpt2_subtract, k_vec = k_vec)
      EV.func_curry <- Curry(LS_function,x = x, y = EV_subtract, k_vec = k_vec)
      
      EV.fit <- tryCatch(optim(EV.func_curry,par = EV.start_spline, lower = c(-10 * k_upperbound,0,0.1),
                               method = "L-BFGS-B", upper = c(50 * k_upperbound,50 * k_upperbound,k_upperbound)),error = function(err){return(1)})
      
      Rpt2.fit <- tryCatch(optim(Rpt2.func_curry,par = Rpt2.start_spline, lower = c(-10 * k_upperbound,0,0.1),
                                 method = "L-BFGS-B", upper = c(50 * k_upperbound,50 * k_upperbound,k_upperbound)),error = function(err){return(1)})
      
      if (inherits(Rpt2.fit,"list")){
        Rpt2.fit_para <- Rpt2.fit$par
        Rpt2.fit_value <- Rpt2.fit$value
      }else{
        Rpt2.fit_para <- c(3,3,0.3)
        Rpt2.fit_value <- NA
      }
      
      if (inherits(EV.fit,"list")){
        EV.fit_para <- EV.fit$par
        EV.fit_value <- EV.fit$value
      }else{
        EV.fit_para <- c(3,3,0.3)
        EV.fit_value <- NA
      }
    
      results <- rbind(results, c(Rpt2.start_spline,EV.start_spline,Rpt2.fit_para,EV.fit_para,Rpt2.fit_value,EV.fit_value,mycompartment,myindex))
      myindex <- myindex + 1
    }
    results
  }
end_time <- date()
stopImplicitCluster()
stopCluster(cl)
save(myresult,start_time,end_time,file = result_name)

