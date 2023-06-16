library(functional)
library(Rmpfr)
library(pracma)
library(doParallel)
library(foreach)
library(numDeriv)

op.f2 = function(n=3) {
  n.func = c()
  for (i in 0:n) {
    ai.inv = c()
    for (j in 0:n) {
      if (i == j) next
      ai.inv = c(ai.inv,paste0('(k[',j+1,']-k[',i+1,'])'))
    }
    ai.inv = paste(ai.inv, collapse = '*')
    ai.inv = paste0('(',ai.inv,')')
    ai.t = paste0('exp(-k[',i+1,']*t)/', ai.inv)
    nf = paste0('1/(k[',i+1,']*', ai.inv,')')
    n.func = c(n.func, ai.t)
  }
  n.func = paste(n.func, collapse = '+')
  norm.fac = "k[1]";
  for (q in 2:(n+1)){
    norm.fac = paste(norm.fac,"*k[",q,"]",sep = "")
  }
  return(list(func=n.func, norm=norm.fac))
}  

poly_function <- function(p,x){
  x.sqrt <- sqrt(x)
  range.x = range(x.sqrt)
  x.sqrt.sd = (x.sqrt - mean(range.x)) / (range.x[2] - range.x[1]) * 2
  mymean <- polyval(p,x.sqrt.sd)
  return(mymean)
}
mypoly_function <- function(p,x){
  x.sqrt <- sqrt(x)
  x.sqrt.sd = (x.sqrt - mean(c(sqrt(24),sqrt(3)))) / (sqrt(24) - sqrt(3)) * 2
  mymean <- polyval(p,x.sqrt.sd)
  return(mymean)
}
likelihood_poly <- function(p,y,x,theta,offset){
  mymean <- poly_function(p,x)
  mymu <- exp(mymean + offset)
  likelihood_sum <- lapply(1:length(mymu),function(x){
    dnbinom(y[x],size = theta, mu = mymu[x],log = T)})
  -sum(unlist(likelihood_sum))
}
func_string <- function(compart){
  compart1 = compart[1]
  op.f2.p1.func = op.f2(compart1 + 1)$func
  op.f2.p1.norm = op.f2(compart1)$norm
  all.f.p1.str=paste0('(', op.f2.p1.func, ')*(', op.f2.p1.norm, ')')
  last.k = paste0('k\\[',compart1+2,'\\]')
  all.f.p1.str = gsub(last.k, 'kx', all.f.p1.str)
  if (length(compart) == 1){
    all.f.str = paste0('a*(',all.f.p1.str,')')
    arg = 't, k, kx, a'
    if (compart < 8){
      func.str = paste0('multi.c.m = function(',arg,') {', all.f.str,'}')
    }else{
      func.str = paste0('multi.c.m = function(',arg,') {k = mpfr(k,100);kx = mpfr(kx, 100);result=', 
                        all.f.str,';return(asNumeric(result)) }')
    }
  }else if(length(compart) == 2){
    compart2 = compart[2]
    op.f2.p2.func = op.f2(compart2 + 1)$func
    op.f2.p2.norm = op.f2(compart2)$norm
    all.f.p2.str=paste0('(', op.f2.p2.func, ')*(', op.f2.p2.norm, ')')
    last.k = paste0('k\\[',compart2+2,'\\]')
    all.f.p2.str = gsub(last.k, 'kx', all.f.p2.str)
    all.f.str = paste0('a1*(',all.f.p1.str,')+a2*(',all.f.p2.str,')')
    arg = 't, k, kx, a1, a2'
    func.str = paste0('multi.c.m = function(',arg,') {k = mpfr(k,100);kx = mpfr(kx, 100);result= ', 
                      all.f.str,';return(asNumeric(result)) }')
  }
  return(func.str)
}
compartment_output <- function(t,p,k_vec){
  if (length(p) == 3){
    a1_fit <- p[1]
    a2_fit <- p[2]
    k_fit <- p[3]
    r1 = func1(t - 3, k = k_vec, kx = k_fit, a = a1_fit)
    r2 = func2(t - 3, k = k_vec, kx = k_fit, a = a2_fit) 
    return(r1 + r2)
  }else if(length(p) == 2){
    a1_fit <- p[1]
    a2_fit <- p[2]
    k_fit <- 0.8
    r1 = func1(t - 3, k = k_vec, kx = k_fit, a = a1_fit)
    r2 = func2(t - 3, k = k_vec, kx = k_fit, a = a2_fit) 
    return(r1 + r2)
  }else{
    return("error")
  }
}
peak_time_function <- function(p,k_vec){
  t = seq(0.5,30,0.01) + 3
  outputs = compartment_output(t=t,p=p,k_vec=k_vec)
  first_series = outputs[1:length(outputs) - 1]
  second_series = outputs[2:length(outputs)]
  diff = second_series - first_series
  diff_point = diff[1:length(diff) - 1] * diff[2:length(diff)] < 0
  maxima = c()
  minima = c()
  for (checkpoint in which(diff_point == T)){
    if (first_series[checkpoint + 1] > first_series[checkpoint]){
      if (first_series[checkpoint + 1] > first_series[checkpoint + 2]){
        maxima = c(maxima,checkpoint/length(t) * 29.5 + 3.5)
      }
    }else{
      if (first_series[checkpoint + 1] < first_series[checkpoint + 2]){
        minima = c(minima,checkpoint/length(t) * 29.5 + 3.5)
      }
    }
  }
  return(list(maxima = maxima, minima = minima))
}
model_output <- function(t, p , k_vec,mock){
  parameter <- p[4:length(p)]
  pf = mypoly_function(parameter,t)
  if (mock == 1){
    return(pf)
  }else if(mock == 3){
    return(p[1] + pf)
  }else{
    a1_fit <- p[1]
    a2_fit <- p[2]
    k_fit <- p[3]
    #k_vec <- mpfr(k_vec,80)
    r1 = func1(t - 3, k = k_vec, kx = k_fit, a = a1_fit)
    r2 = func2(t - 3, k = k_vec, kx = k_fit, a = a2_fit) 
    return(r1 + r2 + pf)
  }
}
model_output.after <- function(t, p , k_vec,mock){
  parameter <- p[4:length(p)]
  pf = mypoly_function(parameter,t)
  if (mock == 1){
    return(pf)
  }else if(mock == 3){
    return(p[1] + pf)
  }else{
    a1_fit <- p[1]
    a2_fit <- p[2]
    k_fit <- p[3]
    #k_vec <- mpfr(k_vec,80)
    r1 = func1.after(t - 3, k = k_vec, kx = k_fit, a = a1_fit)
    r2 = func2.after(t - 3, k = k_vec, kx = k_fit, a = a2_fit) 
    return(r1 + r2 + pf)
  }
}
model_output_sp <- function(t, p , k_vec){
  parameter <- p[3:length(p)]
  pf = mypoly_function(parameter,t)
  a_fit <- p[1]
  k_fit <- p[2]
  #k_vec <- mpfr(k_vec,80)
  r = func(t - 3, k = k_vec, kx = k_fit, a = a_fit)
  return(r + pf)
}
likelihood_function <- function(p, y_list, t_list, offset_list,
                                k_vec, theta, mutant_no,compart,poly_order){
  likelihood <- 0
  poly_para <- p[(2 * mutant_no + 2):length(p)]
  mypoly_para <- rep(0,6)
  mypoly_para[poly_order] <- poly_para
  for (thecondi in 1:(mutant_no + 1)){
    if (thecondi == (mutant_no + 1)){
      mock_option = 1
      parameters <- c(10,10,5,mypoly_para)
    }else{
      mock_option = 0
      a1 <- p[2 * thecondi - 1]
      a2 <- p[2 * thecondi]
      k <- p[2 * mutant_no + 1]
      parameters <- c(a1,a2,k,mypoly_para)
    }
    t_condi <- t_list[[thecondi]]
    y_condi <- y_list[[thecondi]]
    offset_condi <- offset_list[[thecondi]]
    
    mymean <- model_output(t=t_condi,p=parameters,k_vec = k_vec,mock= mock_option)
    mymu <- exp(mymean + offset_condi)
    likelihood_sum <- lapply(1:length(mymu),function(x){
      dnbinom(y_condi[x],size = theta, mu = mymu[x],log = T)})
    likelihood <- likelihood + sum(unlist(likelihood_sum))
  }
  -likelihood
}
likelihood_function_a2k_same <- function(p, y_list, t_list, offset_list,
                                         k_vec, theta, mutant_no, compart, poly_order){
  likelihood <- 0
  poly_para <- p[(2 * mutant_no + 1):length(p)]
  mypoly_para <- rep(0,6)
  mypoly_para[poly_order] <- poly_para
  for (thecondi in 1:(mutant_no + 1)){
    if (thecondi == (mutant_no + 1)){
      mock_option = 1
      parameters <- c(10,10,5,mypoly_para)
    }else{
      mock_option = 0
      a1 <- p[thecondi]
      a2 <- p[3]
      k <- p[4]
      parameters <- c(a1,a2,k,mypoly_para)
    }
    t_condi <- t_list[[thecondi]]
    y_condi <- y_list[[thecondi]]
    offset_condi <- offset_list[[thecondi]]
    
    mymean <- model_output(t=t_condi,p=parameters,k_vec = k_vec,mock= mock_option)
    mymu <- exp(mymean + offset_condi)
    likelihood_sum <- lapply(1:length(mymu),function(x){
      dnbinom(y_condi[x],size = theta, mu = mymu[x],log = T)})
    likelihood <- likelihood + sum(unlist(likelihood_sum))
  }
  -likelihood
}
likelihood_function_a2_same <- function(p, y_list, t_list, offset_list,
                                        k_vec, theta, mutant_no, compart, poly_order){
  likelihood <- 0
  poly_para <- p[(2 * mutant_no + 2):length(p)]
  mypoly_para <- rep(0,6)
  mypoly_para[poly_order] <- poly_para
  for (thecondi in 1:(mutant_no + 1)){
    if (thecondi == (mutant_no + 1)){
      mock_option = 1
      parameters <- c(10,10,5,mypoly_para)
    }else{
      mock_option = 0
      a1 <- p[thecondi]
      a2 <- p[3]
      k <- p[thecondi + 3]
      parameters <- c(a1,a2,k,mypoly_para)
    }
    t_condi <- t_list[[thecondi]]
    y_condi <- y_list[[thecondi]]
    offset_condi <- offset_list[[thecondi]]
    
    mymean <- model_output(t=t_condi,p=parameters,k_vec = k_vec,mock= mock_option)
    mymu <- exp(mymean + offset_condi)
    likelihood_sum <- lapply(1:length(mymu),function(x){
      dnbinom(y_condi[x],size = theta, mu = mymu[x],log = T)})
    likelihood <- likelihood + sum(unlist(likelihood_sum))
  }
  -likelihood
}
likelihood_function_k_flex <- function(p, y_list, t_list, offset_list,
                                       k_vec, theta, mutant_no, compart, poly_order){
  likelihood <- 0
  poly_para <- p[(2 * mutant_no + 3):length(p)]
  mypoly_para <- rep(0,6)
  mypoly_para[poly_order] <- poly_para
  for (thecondi in 1:(mutant_no + 1)){
    if (thecondi == (mutant_no + 1)){
      mock_option = 1
      parameters <- c(10,10,5,mypoly_para)
    }else{
      mock_option = 0
      a1 <- p[3 * thecondi - 2]
      a2 <- p[3 * thecondi - 1]
      k <- p[3 * thecondi]
      parameters <- c(a1,a2,k,mypoly_para)
    }
    t_condi <- t_list[[thecondi]]
    y_condi <- y_list[[thecondi]]
    offset_condi <- offset_list[[thecondi]]
    
    mymean <- model_output(t=t_condi,p=parameters,k_vec = k_vec,mock= mock_option)
    mymu <- exp(mymean + offset_condi)
    likelihood_sum <- lapply(1:length(mymu),function(x){
      dnbinom(y_condi[x],size = theta, mu = mymu[x],log = T)})
    likelihood <- likelihood + sum(unlist(likelihood_sum))
  }
  -likelihood
}
likelihood_function_Rpt2 <- function(p, y_Rpt2, time, offset_Rpt2, poly_para,
                                     k_vec, theta, weight = NA, sp = F){
  mypoly_para <- poly_para
  if (sp == T){
    mymean <- model_output_sp(t=time,p=c(p,poly_para),k_vec = k_vec)
  }else{
    mymean <- model_output(t=time,p=c(p,poly_para),k_vec = k_vec,mock= 0)
  }
  mymu <- exp(mymean + offset_Rpt2)
  likelihood_sum <- lapply(1:length(mymu),function(x){
    dnbinom(y_Rpt2[x],size = theta, mu = mymu[x],log = T)})
  if (is.na(weight)){
    likelihood <- sum(unlist(likelihood_sum))
  }else{
    likelihood <- sum(unlist(likelihood_sum) * weight)
  }
  -likelihood
}

likelihood_function_intercept <- function(p, y_list, t_list, offset_list,
                                          k_vec, theta, mutant_no,compart,poly_order){
  likelihood <- 0
  poly_para <- p[(2 * mutant_no + 1):length(p)]
  mypoly_para <- rep(0,6)
  mypoly_para[poly_order] <- poly_para
  for (thecondi in 1:(mutant_no + 1)){
    if (thecondi == (mutant_no + 1)){
      mock_option = 1
      parameters <- c(10,10,5,mypoly_para)
    }else if(thecondi == 1){
      mock_option = 0
      a1 <- p[1]
      a2 <- p[2]
      k <- p[3]
      parameters <- c(a1,a2,k,mypoly_para)
    }else if(thecondi == 2){
      mock_option = 3
      intercept = p[4]
      parameters <- c(rep(intercept,3),mypoly_para)
    }
    t_condi <- t_list[[thecondi]]
    y_condi <- y_list[[thecondi]]
    offset_condi <- offset_list[[thecondi]]
    
    mymean <- model_output(t=t_condi,p=parameters,k_vec = k_vec,mock= mock_option)
    mymu <- exp(mymean + offset_condi)
    likelihood_sum <- lapply(1:length(mymu),function(x){
      dnbinom(y_condi[x],size = theta, mu = mymu[x],log = T)})
    likelihood <- likelihood + sum(unlist(likelihood_sum))
  }
  -likelihood
}

LS_function <- function(p,x,y,k_vec){
  mymean <- compartment_output(t = x, p = p, k_vec = k_vec)
  sse <- sum((mymean - y) ^ 2)
  return(sse)
}
mps = function(coefs,order = 6) {
  cum.term = rep(0,order)
  n.cf = names(coefs)
  for (i in 1:length(coefs)) {
    if (grepl('^\\(Inter', n.cf[i])) {
      cum.term[1] <- coefs[i]
    }else if(grepl('^I\\(', n.cf[i])){
      ep = as.integer(substr(n.cf[i],13,13)) + 1
      cum.term[ep] = coefs[i]
    }else{
      cum.term[2] = coefs[i]
    }
  }
  cum.term <- rev(cum.term)
  return(cum.term)
}
compartment_function <- function(mycompartment){
  if (mycompartment %in% c(1,2,3,4)){
    compart1 <- 1
    compart2 <- 1 + mycompartment + 1
  }else if(mycompartment %in% c(5,6,7)){
    compart1 <- 2
    compart2 <- 2 + mycompartment - 3
  }else if(mycompartment %in% c(8,9)){
    compart1 <- 3
    compart2 <- 3 + mycompartment - 6
  }else{
    compart1 <- 4
    compart2 <- 6
  }
  return(c(compart1,compart2))
}
#no longer used function
compartment_function.add <- function(mycompartment){
  compart1_vec = c(2,2,2,2,2,4,4,4)
  compart2_vec = c(4,5,7,9,11,7,9,11)
  compart1 = compart1_vec[mycompartment - 10]
  compart2 = compart2_vec[mycompartment - 10]
  return(c(compart1, compart2))
}
#It looks like I missed a compartment combination here, compart1 = 1, compart2 = 4
#But since we decided not to use the full time course data, I will not change this to keep the raw record.
#Note: 2021/06/02
compartment_function.add.16h <- function(mycompartment){
  compart1_vec = c(1,1,2,2,2,3,4)
  compart2_vec = c(3,4,4,5,7,5,7)
  compart1 = compart1_vec[mycompartment]
  compart2 = compart2_vec[mycompartment]
  return(c(compart1, compart2))
}




#compart_transform <- function(mycompartment){
#  if (mycompartment %in% c(1,2,3,4,5)){
#    compart1 <- 1
#    compart2 <- 1 + mycompartment + 1
#  }else if(mycompartment %in% c(6,7,8,9)){
#    compart1 <- 2
#    compart2 <- 2 + mycompartment - 4
#  }else if(mycompartment %in% c(10,11,12)){
#    compart1 <- 3
#    compart2 <- 3 + mycompartment - 8
#  }else if(mycompartment %in% c(13,14)){
#    compart1 <- 4
#    compart2 <- 4 + mycompartment - 11
#  }else{
#    compart1 <- 5
#    compart2 <- 7
#  }
#  return(c(compart1,compart2))
#}

model_output_spa_hete <- function(t, p , k_vec,mock){
  parameter <- p[3:length(p)]
  pf = poly_function(parameter,t)
  if (mock == 1){
    return(pf)
  }else if(mock == 3){
    return(p[1] + pf)
  }else{
    a1_fit <- p[1]
    a2_fit <- p[2]
    k_fit <- 0.8
    #k_vec <- mpfr(k_vec,80)
    r1 = func1(t - 3, k = k_vec, kx = k_fit, a = a1_fit)
    r2 = func2(t - 3, k = k_vec, kx = k_fit, a = a2_fit) 
    return(r1 + r2 + pf)
  }
}

likelihood_function_spa_hete <- function(p, y_list, t_list, offset_list,
                                         k_vec, theta, mutant_no,compart,poly_order){
  likelihood <- 0
  poly_para <- p[(2 * mutant_no + 1):length(p)]
  mypoly_para <- rep(0,6)
  mypoly_para[poly_order] <- poly_para
  for (thecondi in 1:(mutant_no + 1)){
    if (thecondi == (mutant_no + 1)){
      mock_option = 1
      parameters <- c(10,10,mypoly_para)
    }else{
      mock_option = 0
      a1 <- p[2 * thecondi - 1]
      a2 <- p[2 * thecondi]
      parameters <- c(a1,a2,mypoly_para)
    }
    t_condi <- t_list[[thecondi]]
    y_condi <- y_list[[thecondi]]
    offset_condi <- offset_list[[thecondi]]
    
    mymean <- model_output_spa_hete(t=t_condi,p=parameters,k_vec = k_vec,mock= mock_option)
    mymu <- exp(mymean + offset_condi)
    likelihood_sum <- lapply(1:length(mymu),function(x){
      dnbinom(y_condi[x],size = theta, mu = mymu[x],log = T)})
    likelihood <- likelihood + sum(unlist(likelihood_sum))
  }
  -likelihood
}

library(evd)
single_peak_gumbel = function(x,p){
  a1 = p[1]
  mu1 = p[2]
  w = p[3]
  a1 * dgumbel(x - 3, loc = mu1 - 3, scale = w * (mu1 - 3))
}
double_peak_gumbel = function(x,p){
  a1 = p[1]
  a2 = p[2]
  mu1 = p[3]
  mu2 = p[4]
  w = p[5]
  a1 * dgumbel(x - 3, loc = mu1 -3, scale = w * (mu1 - 3)) + a2 * dgumbel(x - 3, loc = mu2 - 3 , scale = w * (mu2 - 3))
}
gumbel_single_fit_LS = function(x,y,p){
  output = single_peak_gumbel(x=x,p=p)
  sse <- sum((output - y) ^ 2)
  return(sse)
}
gumbel_double_fit_LS = function(x,y,p){
  output = double_peak_gumbel(x=x,p=p)
  sse <- sum((output - y) ^ 2)
  return(sse)
}

first_derive = function(t, y){
  if (length(t) != length(y)){
    print("error")
    return(0)
  }else{
    time = t[1:(length(y) - 1)]
    derive = y[2:length(y)] - y[1:(length(y) - 1)]
    return(rbind(time,derive))
  }
}
second_derive = function(t,y){
  time = first_derive(t = t, y = y)[1,]
  first_derivative = first_derive(t=t,y=y)[2,]
  thetime = time[1:(length(time)-1)]
  thederive = first_derivative[2:length(first_derivative)] - first_derivative[1:(length(first_derivative)-1)]
  return(rbind(thetime, thederive))
}

compartment_output_sp <- function(t,p,k_vec){
  a <- p[1]
  k <- p[2]
  r = func(t - 3, k = k_vec, kx = k, a = a)
  return(r)
}
LS_function_sp <- function(p,x,y,k_vec){
  mymean <- compartment_output_sp(t = x, p = p, k_vec = k_vec)
  sse <- sum((mymean - y) ^ 2)
  return(sse)
}
