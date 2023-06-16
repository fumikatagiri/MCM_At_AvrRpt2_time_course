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
k_upper_bound = 10
a_upper_bound = k_upper_bound * 50
a_lower_bound = -k_upper_bound * 10

input_lognormal = "MCM.fit.continous.demo.RData"   # preliminary MCM fit results to continuous data
input_lognormal = paste0('./data/', input_lognormal)
result_name = 'MCM.fit.read.count.demo.RData'      # output file name - results of MCM fit
result_name = paste0('./data/', result_name)

load(input_lognormal)
demo.gene.set = as.character(unique(results[,1]))
comp.conv = 1:6*2-1
#####

start_time = date()
myresults = c()
for (mygene in demo.gene.set) {
  sel.mat = results[results[,1] == mygene,]
  for (sel.row in 1:nrow(sel.mat)) {
    app.fit = sel.mat[sel.row,]
    mycompartment = as.integer(app.fit[2])
    compart1 = comp.conv[compartment_function(mycompartment)[1]]
    compart2 = comp.conv[compartment_function(mycompartment)[2]]
    func1_str = func_string(compart1)
    func1 = eval(parse(text =func1_str))  # MCM for compartment 1
    func2_str = func_string(compart2)
    func2 = eval(parse(text =func2_str))  # MCM for compartment 2
    
    y_after3 = as.numeric(col_count_data[mygene,after_3_indexes])  # read count data
    y_Rpt2 = y_after3[treatment_after3 == "AvrRpt2"]
    y_mock = y_after3[treatment_after3 == "mock"]
    mycount_list = list(Rpt2 = y_Rpt2,mock = y_mock)
    
    x = as.numeric(time_after3[treatment_after3 == "mock"])
    time_list = list(Rpt2 = x,mock = x)

    poly_start = mock_coef[mygene,1:6]
    compart_start = app.fit[9:11]  # Rpt2.fit_para
    
    mytheta = mock_coef[mygene,7]
    likelihood_function_curry = Curry(likelihood_function_Rpt2, y_Rpt2 = y_Rpt2,
                                      time = x, offset_Rpt2 = offset_Rpt2_after3, 
                                      theta = mytheta, poly_para = poly_start,
                                      k_vec = k_vec)
    
    fit = tryCatch(optim(likelihood_function_curry,par =  compart_start, lower = c(a_lower_bound,0,0.1),
                         method = "L-BFGS-B", upper = c(a_upper_bound,a_upper_bound,k_upper_bound)),
                   error = function(err){return(1)})
    
    if (inherits(fit,"list")){
      myresults = rbind(myresults, c(fit$par,poly_start,fit$convergence,fit$value,mycompartment, mygene))
    }else{
      myresults = rbind(myresults, c(rep(NA,3),poly_start,NA,NA,mycompartment, mygene))
    }
  }
}  # 2 min
end_time = date()

save(myresults,start_time,end_time,file = result_name)

############
#### selection of best MCM based on the likelihood (column 11 is -likelihood, so mininum is the best)
parameter_best = c()
for (mygene in demo.gene.set) {
  sel.mat = myresults[myresults[,13] == mygene,]
  best.m = sel.mat[which.min(sel.mat[,11]),]
  parameter_best = rbind(parameter_best, best.m)
}
colnames(parameter_best) = c('a1','a2','k','p5','p4','p3','p2','p1','p0',
                             'convergence', 'likelihood','compart_index','gene')
gene.n = parameter_best[,'gene']
parameter_best = parameter_best[,-13]
parameter_best = apply(parameter_best, 2, as.numeric)
rownames(parameter_best) = gene.n

save(parameter_best, file = './data/merged_result_demo.RData')
