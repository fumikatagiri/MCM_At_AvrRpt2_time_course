rm(list=ls())
load("./data/col_count_data.Rdata")   # RNA-seq read count data
load("./data/mock_start.vals.Rdata")  # polynomial regression model for mock time course

source("MCM_function_copy.r")    # collection of MCM-related functions
all.genes = rownames(col_count_data)
treatment = c()
time = c()
for (mylib in colnames(col_count_data)){
  treatment = c(treatment, strsplit(mylib,split = "_")[[1]][4])
  time = c(time, strsplit(mylib,split = "_")[[1]][5])
}
time = as.numeric(substr(time,start=1,stop=2))
after_3_indexes = which(time >= 3)
time_after3 = time[after_3_indexes]
treatment_after3 = treatment[after_3_indexes]
col_offset_after3 = col_offset[after_3_indexes]
offset_mock_after3 = col_offset_after3[treatment_after3 == "mock"]
offset_Rpt2_after3 = col_offset_after3[treatment_after3 == "AvrRpt2"]
offset_EV_after3 = col_offset_after3[treatment_after3 == "EV"]
myoffset_list = list(Rpt2 = offset_Rpt2_after3,EV = offset_EV_after3,mock = offset_mock_after3)

rownames(mock_coef) = all.genes
k_vec = c(1.8, 0.5487, 1.1733, 1.1253, 0.6492, 0.6853, 0.6957, 0.6987, 0.4949,
          0.5065, 0.513, 0.5174, 0.2133)  # Set1 decay rates
#k_vec = c(1.8, 2.4573, 1.0775, 1.3624, 0.6624, 0.7903, 0.646, 0.6757, 0.5669, 0.5814,
#          0.5003, 0.5094)  # Set2 decay rates
k_upperbound = 1.5

####
demo.gene.set = c("AT1G10865", "AT1G22650")  # two genes for a demo
result_name = "MCM.fit.continous.demo.RData" # output .RData file name
result_name = paste0('./data/', result_name)
comp.conv = 1:6*2-1

start_time = date()
results = c()
for (mygene in demo.gene.set) {
  for (mycompartment in 1:10) {
    compart1 = comp.conv[compartment_function(mycompartment)[1]]
    compart2 = comp.conv[compartment_function(mycompartment)[2]]
    func1_str = func_string(compart1)
    func1 = eval(parse(text =func1_str))  # MCM for compartment 1
    func2_str = func_string(compart2)
    func2 = eval(parse(text =func2_str))  # MCM for compartment 2
    
    y_after3_log = log(as.numeric(col_count_data[mygene,after_3_indexes]))  # consinuous approximation of read count using log
    y_Rpt2_log = y_after3_log[treatment_after3 == "AvrRpt2"] - offset_Rpt2_after3  # between-libraries offsets
    y_EV_log = y_after3_log[treatment_after3 == "EV"] - offset_EV_after3
    y_mock_log = y_after3_log[treatment_after3 == "mock"] - offset_mock_after3
    
    x = as.numeric(time_after3[treatment_after3 == "mock"])  # time points in hpi
    
    poly_start = mock_coef[mygene,1:6]
    mock_subtract = poly_function(p=poly_start,x = x)
    
    Rpt2_subtract = y_Rpt2_log - mock_subtract
    EV_subtract = y_EV_log - mock_subtract
    
    Rpt2.sq.sp.fun = splinefun(sqrt(x), Rpt2_subtract)
    EV.sq.sp.fun = splinefun(sqrt(x), EV_subtract)
    
    Rpt2.sp.val = Rpt2.sq.sp.fun(sqrt(seq(3,24,0.8)))
    EV.sp.val = EV.sq.sp.fun(sqrt(seq(3,24,0.8)))
    
    Rpt2.spline_func_curry = Curry(LS_function,x = seq(3,24,0.8),y = Rpt2.sp.val,k_vec = k_vec)
    EV.spline_func_curry = Curry(LS_function,x = seq(3,24,0.8),y = EV.sp.val,k_vec = k_vec)
    
    EV.fit_spline = tryCatch(optim(EV.spline_func_curry,par = c(3,3,0.3), lower = c(-10 * k_upperbound,0,0.1),
                                   method = "L-BFGS-B", upper = c(50 * k_upperbound ,50 * k_upperbound ,k_upperbound )),error = function(err){return(1)})
    
    Rpt2.fit_spline = tryCatch(optim(Rpt2.spline_func_curry,par = c(3,3,0.3), lower = c(-10 * k_upperbound,0,0.1),
                                     method = "L-BFGS-B", upper = c(50 * k_upperbound,50 * k_upperbound,k_upperbound)),error = function(err){return(1)})
    
    if (inherits(Rpt2.fit_spline,"list")){
      Rpt2.start_spline = Rpt2.fit_spline$par
    }else{
      Rpt2.start_spline = c(3,3,0.3)
    } 
    if (inherits(EV.fit_spline,"list")){
      EV.start_spline = EV.fit_spline$par
    }else{
      EV.start_spline = c(3,3,0.3)
    }
    
    Rpt2.func_curry = Curry(LS_function,x = x, y = Rpt2_subtract, k_vec = k_vec)
    EV.func_curry = Curry(LS_function,x = x, y = EV_subtract, k_vec = k_vec)
    
    EV.fit = tryCatch(optim(EV.func_curry,par = EV.start_spline, lower = c(-10 * k_upperbound,0,0.1),
                            method = "L-BFGS-B", upper = c(50 * k_upperbound,50 * k_upperbound,k_upperbound)),error = function(err){return(1)})
    
    Rpt2.fit = tryCatch(optim(Rpt2.func_curry,par = Rpt2.start_spline, lower = c(-10 * k_upperbound,0,0.1),
                              method = "L-BFGS-B", upper = c(50 * k_upperbound,50 * k_upperbound,k_upperbound)),error = function(err){return(1)})
    
    if (inherits(Rpt2.fit,"list")){
      Rpt2.fit_para = Rpt2.fit$par
      Rpt2.fit_value = Rpt2.fit$value
    }else{
      Rpt2.fit_para = c(3,3,0.3)
      Rpt2.fit_value = NA
    }
    
    if (inherits(EV.fit,"list")){
      EV.fit_para = EV.fit$par
      EV.fit_value = EV.fit$value
    }else{
      EV.fit_para = c(3,3,0.3)
      EV.fit_value = NA
    }
    
    results = rbind(results, c(mygene, mycompartment, Rpt2.start_spline,EV.start_spline,Rpt2.fit_para,EV.fit_para,Rpt2.fit_value,EV.fit_value))
  }
}  # 9.5 min
end_time = date()

save(results,start_time,end_time,file = result_name)

