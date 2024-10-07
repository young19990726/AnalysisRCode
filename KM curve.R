#### Libraries ####
library(cowplot)
library(ggplot2)
library(ggpubr)
library(glmnet)
library(gridExtra) 
library(magrittr)
library(mice, lib.loc = "/usr/local/lib/R/site-library")
library(pROC)
library(scales) 
library(survival)
library(xgboost)

#### Pathways ####
setwd('/home/young19990726/RCODE/AI-CXR for aortic aneurysm')
data_path = './patient_data.RData'
plot_path = './KM curve.pdf'

#### Functions ####
KM_PLOT = function (time, status, show_risk_table = TRUE, x, cont_x = NULL, age = NULL, sex = NULL, MAIN = 'KM curve', km_curve = TRUE, x_lim = c(0, 2, 4, 6), y_lim = 0.5, xlab_name = 'months', ylab_name = 'Cumulative incidence', col_list = c('#F79952', '#EA7A81', '#7C2451'), x_name = expression(paste('Low risk'), paste('Median risk'), paste('High risk'))) {
  
  if (!km_curve & (is.null(age) | is.null(sex))) {stop('Need to add age and sex')}
  
  sub_data = data.frame(time, status, x)
  sub_data[, 'x'] = factor(sub_data[, 'x'])
  
  if (!is.null(cont_x)) {sub_data[, 'cont_x'] = cont_x}
  if (!is.null(age)) {sub_data[, 'age'] = age}
  if (!is.null(sex)) {sub_data[, 'sex'] = factor(sex)}
  
  sub_data = sub_data[apply(is.na(sub_data), 1, sum) %in% 0,]
  sub_data[sub_data[, 'time'] >= max(x_lim), 'status'] = 0L
  sub_data[sub_data[, 'time'] >= max(x_lim), 'time'] = max(x_lim)
  
  if (sum(sub_data[, 'status'] == 1) > 2) {
    
    if (!is.null(cont_x)) {
      
      C_idx_model = coxph(Surv(time, status) ~ cont_x, data = sub_data[, colnames(sub_data)[!colnames(sub_data) %in% 'x']])
      concordance = concordance(C_idx_model)
      C_idx_txt = paste0('C-index = ', formatC(concordance[['concordance']], format = 'f', 3), ' (', formatC(concordance[['concordance']] + qnorm(0.025) * sqrt(concordance[['var']]), format = 'f', 3), '-', formatC(concordance[['concordance']] + qnorm(0.975) * sqrt(concordance[['var']]), format = 'f', 3), ')')
      
    }
    
    survival_model = coxph(Surv(time, status) ~ x + ., data = sub_data[, colnames(sub_data)[!colnames(sub_data) %in% 'cont_x']])
    new_data = data.frame(x = levels(sub_data[, 'x']), age = mean(sub_data[, 'age']), sex = levels(sub_data[, 'sex'])[1])
    
    if (km_curve) {
      
      Predict.Surv = survfit(as.formula(paste0('Surv(time, status) ~ x')), data = sub_data)
      group_surv_data = data.frame(group = rep(letters[1:length(x_name)], Predict.Surv$strata), time = Predict.Surv$time, risk = (1 - Predict.Surv$surv) * 100) 
      
    } else {
      
      Predict.Surv = survfit(survival_model, newdata = new_data)
      group_surv_data = data.frame(group = rep(letters[1:length(x_name)], each = length(Predict.Surv$time)), time = rep(Predict.Surv$time, 2), risk = (1 - as.numeric(Predict.Surv$surv)) * 100)
      
    }
    
    group_surv_data = rbind(group_surv_data, data.frame(group = letters[1:length(x_name)], time = 0,  risk = 0))
    
    for (group in letters[1:length(x_name)]) {
      
      group_surv_data = rbind(group_surv_data, data.frame(group = group, time = max(x_lim), risk = max(group_surv_data[group_surv_data[,'group'] %in% group,'risk'])))
      
    }
    
    group_surv_data = group_surv_data[order(group_surv_data[, 'time']),]
    group_surv_data = group_surv_data[order(group_surv_data[, 'group']),]
    group_surv_data[group_surv_data[, 'time'] > max(x_lim), 'time'] = max(x_lim)
    group_surv_data = group_surv_data[!duplicated(group_surv_data[, c('group', 'time')]),]
    
    max_time = tapply(group_surv_data[, 'time'], group_surv_data[, 'group'], max)
    
    y.HR = summary(survival_model)[['coefficients']][1:(length(x_name) - 1), 1]
    se.HR = summary(survival_model)[['coefficients']][1:(length(x_name) - 1), 3]
    y.HR[y.HR > 6] = 1e5
    se.HR[y.HR > 6] = 1e6
    HR_txt = paste0(formatC(exp(y.HR), 2, format = 'f'), ' (', formatC(exp(y.HR - qnorm(0.975) * se.HR), 2, format = 'f'), ', ', formatC(exp(y.HR + qnorm(0.975) * se.HR), 2, format = 'f'), ')')
    HR_txt = c('Reference', HR_txt)
    
    KM_surv = survfit(as.formula(paste0('Surv(time, status) ~ x')), data = sub_data)
    
    AT_RISK_TABLE = data.frame(group = rep(letters[1:length(x_name)], KM_surv$strata), time = KM_surv$time, n.risk = KM_surv$n.risk, surv = KM_surv$surv, stringsAsFactors = FALSE)
    AT_RISK_TABLE = rbind(AT_RISK_TABLE, data.frame(group = letters[1:length(x_name)], time = 0, n.risk = as.numeric(table(sub_data[,'x'])), surv = 1))
    AT_RISK_TABLE = AT_RISK_TABLE[order(AT_RISK_TABLE[, 'time']),]
    AT_RISK_TABLE = AT_RISK_TABLE[order(AT_RISK_TABLE[, 'group']),]  
    AT_RISK_SUMMARY = data.frame(x = rep(1:length(x_name), each = length(x_lim)),
                                 color = rep(col_list[1:length(x_name)], each = length(x_lim)),
                                 group = rep(letters[1:length(x_name)], each = length(x_lim)),
                                 time = x_lim,
                                 n.risk = NA,
                                 surv = NA,
                                 stringsAsFactors = FALSE)
    
    for (m in 1:length(x_name)) {
      max_time = max(AT_RISK_TABLE[AT_RISK_TABLE[, 'group'] %in% letters[m], 'time'])
      for (l in x_lim) {
        summary_pos = which(AT_RISK_SUMMARY[, 'group'] %in% letters[m] & AT_RISK_SUMMARY[, 'time'] %in% l)
        if ((l - max_time) > 3) {
          AT_RISK_SUMMARY[summary_pos, 'n.risk'] = 0
          AT_RISK_SUMMARY[summary_pos, 'surv'] = NA
        } else {
          table_pos = which(AT_RISK_TABLE[, 'group'] %in% letters[m] & AT_RISK_TABLE[, 'time'] >= min(l, max_time))[1]
          AT_RISK_SUMMARY[summary_pos, 'n.risk'] = AT_RISK_TABLE[table_pos, 'n.risk']
          AT_RISK_SUMMARY[summary_pos, 'surv'] = 1 - AT_RISK_TABLE[table_pos, 'surv']
        }
      }
    }
    
    AT_RISK_SUMMARY[, 'txt'] = paste0(AT_RISK_SUMMARY[, 'n.risk'], '\n(', formatC(AT_RISK_SUMMARY[, 'surv'] * 100, 1, format = 'f'), '%)')
    AT_RISK_SUMMARY[AT_RISK_SUMMARY[,'n.risk'] %in% 0, 'txt'] = ''
    
    time_diff = max(x_lim) - min(x_lim)
    
    gg_p = ggplot(AT_RISK_SUMMARY, aes(x = time, y = x, group = group))
    gg_p = gg_p + geom_text(label = AT_RISK_SUMMARY[, 'txt'], color = AT_RISK_SUMMARY[, 'color'], size = 3, fontface = 2)
    gg_p = gg_p + ylim(c(0.5, length(x_name) + 0.5))
    gg_p = gg_p + xlim(c(min(x_lim) - time_diff * 0.05, max(x_lim) + time_diff * 0.05))
    gg_p = gg_p + ggtitle('Number at risk/event rate (%)')
    gg_p = gg_p + theme(plot.title = element_text(color = "#000000", size = 5), legend.position = "none")
    table_p = gg_p + theme_void()
    
    gg_p = ggplot(group_surv_data, aes(x = time, y = risk, group = group))
    gg_p = gg_p + geom_step(aes(color = group), size = 1)
    gg_p = gg_p + theme_bw()
    gg_p = gg_p + scale_color_manual(values = paste0(col_list[1:length(x_name)], 'A0'))
    gg_p = gg_p + xlab(xlab_name)
    gg_p = gg_p + ylab(paste0(ylab_name, " (%)"))
    gg_p = gg_p + scale_x_continuous(limits = c(0, max(x_lim)), breaks = x_lim, labels = x_lim)
    gg_p = gg_p + scale_y_continuous(limits = c(0, y_lim * 100), breaks = c(seq(0, y_lim * 100, length.out = 6)))
    gg_p = gg_p + ggtitle(MAIN)
    
    for (m in 1:length(x_name)) {
      gg_p = gg_p + annotate(geom = "text", x = 0, y = y_lim * (55 + 8.5 * m), label = x_name[m], size = 3.5, color = substr(col_list[m], 1, 7), hjust = 0, fontface = 2)
      gg_p = gg_p + annotate(geom = "text", x = max(x_lim) * 0.45, y = y_lim * (55 + 8.5 * m), label = HR_txt[m], size = 3.5, color = 'black', hjust = 0)
    }
    
    gg_p = gg_p + annotate(geom = "text", x = max(x_lim) * 0.45, y = y_lim * (55 + 8.5 * (length(x_name) + 1)), label = 'Adjusted HR:', size = 3.5, color = 'black', hjust = 0, fontface = 2)
    if (!is.null(cont_x)) {gg_p = gg_p + annotate(geom = "text", x = 0, y = y_lim * (55 + 8.5 * (length(x_name) + 2)), label = C_idx_txt, size = 3.5, color = 'black', hjust = 0, fontface = 2)}
    gg_p = gg_p + theme(plot.title = element_text(color = "#000000", size = 14), legend.position = "none")
    
    survival_p = ggdraw()
    survival_p = survival_p + draw_plot(gg_p, x = 0, y = 0.30, width = 0.93, height = 0.70)
    survival_p = survival_p + draw_plot(table_p, x = 0.12, y = 0, width = 0.83, height = 0.30)
    
    if (show_risk_table) {return(survival_p)} else {return(gg_p)}
    
  }
  
}

#### Settings ####
set.seed(0)

DATASET_NAME = c('internal-test' = 'Internal validation set', 'community-test' = 'External validation set')

x_list = c('x1', 'x2')
x_name_list = c("predition [[CCT] Thoracic aortic aneurysm]", "predition [[ACT] Abdominal aortic aneurysm]")

#### Processes ####
load(data_path)

## Randomly shuffle data and remove duplicates
patient_data = patient_data[sample(1:nrow(patient_data)),]
patient_data = patient_data[!duplicated(patient_data[, 'CNO']),]

## Risk 
patient_data[, 'x1'] = 'Low risk'
patient_data[patient_data[, 'group [[CCT] Thoracic aortic aneurysm]'] %in% c('Median risk', 'High risk'), 'x1'] = 'High/Median risk'
patient_data[, 'x1'] = factor(patient_data[, 'x1'], levels = c('Low risk', 'High/Median risk'))

patient_data[, 'x2'] = 'Low risk'
patient_data[patient_data[, 'group [[ACT] Abdominal aortic aneurysm]'] %in% c('Median risk', 'High risk'), 'x2'] = 'High/Median risk'
patient_data[, 'x2'] = factor(patient_data[, 'x2'], levels = c('Low risk', 'High/Median risk'))

## Group 
# patient_data = patient_data[patient_data[, '[CCT] Thoracic aortic aneurysm'] %in% 0,]
# patient_data = patient_data[patient_data[, '[ACT] Abdominal aortic aneurysm] %in% 0,]

## Survival analysis:30-day mortality
critical = 30

patient_data[, 'event'] = patient_data[, 'event[all-cause mortality]']
patient_data[, 'time'] = patient_data[, 'time[all-cause mortality]']
patient_data = patient_data[!is.na(patient_data[, 'event']) & !is.na(patient_data[, 'time']),]
patient_data[, 'time'] = (patient_data[, 'time'] + 0.5) 
patient_data = patient_data[patient_data[, 'time'] >= 0,] 
patient_data[patient_data[, 'time'] >= critical, 'event'] = 0L
patient_data[patient_data[, 'time'] >= critical, 'time'] = critical

dat_list = list()
x_lim_list = list()

dat_list[['All-cause mortality (Thoracic aortic aneurysm)']] = impute_data
x_lim_list[[1]] = c(0, 10, 20, 30)

dat_list[['All-cause mortality (Abdominal aortic aneurysm)']] = impute_data
x_lim_list[[2]] = c(0, 10, 20, 30)

#### plotting ####
km_p_list = list()

for (j in 1:length(DATASET_NAME)) {
  km_p_list[[DATASET_NAME[j]]] = list()
  
  for (i in 1:length(dat_list)) {
    sub_data = dat_list[[i]]
    sub_data = sub_data[sub_data[, 'DATASET'] %in% names(DATASET_NAME)[j],]
    
    km_p_list[[DATASET_NAME[j]]][[names(dat_list)[i]]] = KM_PLOT(time = sub_data[, 'time'], status = sub_data[, 'event'], age = sub_data[, 'AGE'], sex = sub_data[, 'GENDER'],
                                                                 x = sub_data[, x_list[i]], cont_x = sub_data[, x_name_list[i]], 
                                                                 MAIN = names(dat_list)[i], km_curve = TRUE, 
                                                                 x_lim = x_lim_list[[i]], y_lim = 0.2,
                                                                 xlab_name = 'days', ylab_name = 'Cumulative incidence',
                                                                 # col_list =c('#F79952', '#EA7A81', '#7C2451'),
                                                                 col_list =c('#F79952', '#7C2451'),
                                                                 x_name = levels(sub_data[, x_list[i]]))
    
  }
  
}

#### Merge plots ####
merge_km_p.1 = arrangeGrob(km_p_list[[1]][[1]], km_p_list[[1]][[2]], nrow = 2)
merge_km_p.2 = arrangeGrob(km_p_list[[2]][[1]], km_p_list[[2]][[2]], nrow = 2)

final_p = ggdraw()
final_p = final_p + draw_plot(merge_km_p.1, x = 0, y = 0, width = 0.5, height = 0.95)
final_p = final_p + draw_plot(merge_km_p.2, x = 0.5, y = 0, width = 0.5, height = 0.95)
final_p = final_p + draw_plot_label(c("Internal validation", "External validation"), c(0.005, 0.495), c(1, 1), size = 14, colour = 'black', hjust = 0, fontface = 2)

pdf(plot_path, width = 12, height = 8)
print(final_p)
dev.off()