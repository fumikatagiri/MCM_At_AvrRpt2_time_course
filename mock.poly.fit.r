library(MASS)

load("./data/col_count_data.Rdata")
#source("MCM_function_copy.r")
all.genes = rownames(col_count_data)
treatment = c()
time = c()
for (mylib in colnames(col_count_data)){
  treatment = c(treatment, strsplit(mylib,split = "_")[[1]][4])
  time = c(time, strsplit(mylib,split = "_")[[1]][5])
}
time = as.numeric(substr(time,start=1,stop=2))
after_3_indexes = which(time >= 3)
mock_indices = which(treatment == 'mock')
mock_3l_indices = intersect(after_3_indexes, mock_indices)

mock.read.c = col_count_data[, mock_3l_indices]
time_3l = time[mock_3l_indices]
x.sqrt <- sqrt(time_3l)
range.x = range(x.sqrt)
tss = (x.sqrt - mean(range.x)) / (range.x[2] - range.x[1]) * 2

offset_mock = col_offset[mock_3l_indices]

all.zero = rep(0, 7)
names(all.zero) = c('I(tss^5)', 'I(tss^4)', 'I(tss^3)', 'I(tss^2)', 'I(tss^1)', 
                    '(Intercept)', 'theta')

## for unknown reason, the step function doesn't work with row 2983 "AT5G64170"
date()
mock_coef = t(apply(mock.read.c[-2983,], 1, function(mrc) {
  gnb = glm.nb(mrc ~ I(tss^5) + I(tss^4) +I(tss^3) +I(tss^2) +I(tss^1) + 
               offset(offset_mock))
  gnb.sel = step(gnb, trace = F)
  cef = gnb.sel$coefficients
  coef = all.zero
  coef[names(cef)] = cef
  coef['theta'] = gnb.sel$theta
  return(coef)
}))  # 1.3 min
date()

# run for row 2983 "AT5G64170"
mrc=as.numeric(mock.read.c[2983,])
gnb = glm.nb(mrc ~ I(tss^5) + I(tss^4) +I(tss^3) +I(tss^2) +I(tss^1) + 
               offset(offset_mock))
cef = gnb$coefficients
coef = all.zero
coef[names(cef)] = cef
coef['theta'] = gnb$theta

mock_coef_fin = rbind(mock_coef[1:2982,], "AT5G64170" = coef, mock_coef[2983:nrow(mock_coef), ])
colnames(mock_coef_fin)[1:6] = paste0('p', 5:0)

mock_coef = mock_coef_fin
save(mock_coef, file = './data/mock_start.vals.Rdata')
