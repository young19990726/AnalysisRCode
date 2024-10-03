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
setwd('/home/young19990726/RCODE/010 AI-CXR for Future compression fracture/New')
data_path = './patient_data_after_mice.RData'
plot_path_1 = './c_index_bar.pdf'

#### Functions ####
c_index_test = function (concordance.1, concordance.2, r) {
  if (concordance.1$n != concordance.2$n) {stop("the concordance indices are computed from different number of samples!")}
  if (is.na(concordance.1$var) || is.na(concordance.2$var)) {stop("the concordance indices must be computed using method noether!")}
  if (is.na(r)) {r = 0}
  
  n = concordance.1$n
  total_var = concordance.1$var + concordance.2$var - 2 * r * sqrt(concordance.2$var * concordance.2$var)
  
  if (total_var > 1e-15) {
    t.stat = (concordance.1$concordance - concordance.2$concordance)/sqrt(concordance.1$var + concordance.2$var - 2 * r * sqrt(concordance.2$var * concordance.2$var))
    diff.ci.p = pt(q = t.stat, df = n - 1, lower.tail = TRUE)
    
    return(diff.ci.p)
  } else {return(1)}
}

#### Settings ####
set.seed(0)

DATASET_NAME = c('internal-test' = 'Internal validation set', 'community-test' = 'External validation set')

col_list.1 = c('#9B0AFAFF', '#9B0AFADF', '#9B0AFABF', '#9B0AFA9F')
col_list.2 = c('#D42300', '#003D9E')
col_list.3 = c('#404040')
col_list.4 = c('#15ED32', '#D6BD09', '#FA9057')
col_list.5 = c('#6893F0', '#003D9E')

y_name_list = c('Y1')
x_name_list = c('predition [[BASIC] Future compression fracture risk]')

critical_cut = 5

intrested_var = c('V0', 'SEX_G', 'AGE_G')

#### Processes ####
load(data_path)

patient_data = current_data

## Select criteria
patient_data = patient_data[!patient_data[, '[DXA] T-score'] %in% NA,]
OP_CNOs = patient_data[patient_data[, '[DXA] T-score'] <= -2.5, 'CNO'] %>% unique()
patient_data[, 'OP_HIS'] = factor((patient_data[, 'CNO'] %in% OP_CNOs) + 0L, levels = 0:1)

patient_data = patient_data[patient_data[, 'CXR[POSITION]'] %in% c('OPD', 'PEC'),]
patient_data = patient_data[patient_data[, 'BMD-drug'] %in% 0,]
patient_data = patient_data[patient_data[, 'Fracture'] %in% 0,]
patient_data = patient_data[!patient_data[, 'event[Fracture_outcome]'] %in% NA & !patient_data[, 'time[Fracture_outcome]'] %in% NA,]

## Define censor
drug_before_end_pos = which(patient_data[, 'event[BMD-drug]'] %in% 1 & patient_data[, 'time[BMD-drug]'] < patient_data[, 'time[Fracture_outcome]'])
patient_data[drug_before_end_pos, 'time[Fracture_outcome]'] = patient_data[drug_before_end_pos, 'time[BMD-drug]']
patient_data[drug_before_end_pos, 'event[Fracture_outcome]'] = 0L

## Randomly shuffle data, sort by time and remove duplicates
patient_data = patient_data[sample(1:nrow(patient_data)),]
patient_data = patient_data[order(patient_data[, 'CXR[TIME]']),]
patient_data = patient_data[!duplicated(patient_data[, 'CNO']) | CXR_data[, 'DATASET'] %in% 'train',]

## Survival analysis
valid_data = patient_data[patient_data[, 'DATASET'] %in% 'valid',]

# new_pred.1：GENDER, AGE 
cox_model.1 = coxph(Surv(valid_data[, 'time[all-cause mortality]'], valid_data[, 'event[all-cause mortality]']) ~ ., data = valid_data[, c('GENDER', 'AGE')])
patient_data[, 'new_pred.1'] = predict(cox_model.1, patient_data)

# new_pred.2：GENDER, AGE, CXR_pred.BMD
cox_model.2 = coxph(Surv(valid_data[, 'time[all-cause mortality]'], valid_data[, 'event[all-cause mortality]']) ~ ., data = valid_data[, c('GENDER', 'AGE', 'predition [[BASIC] Future compression fracture risk]')])
patient_data[, 'new_pred.2'] = predict(cox_model.2, patient_data) 

## Format variables
patient_data = patient_data[!patient_data[, 'RheumatoidArthritis'] %in% NA,]

patient_data[, 'GENDER'] = factor(patient_data[, 'GENDER'])
patient_data[, 'CXR[MODALITY]'] = factor(patient_data[, 'CXR[MODALITY]'])
patient_data[, 'CXR[VIEW]'] = factor(patient_data[, 'CXR[VIEW]'])

# patient_data[, 'Fracture'] = factor(patient_data[, 'Fracture'])
patient_data[, 'RheumatoidArthritis'] = factor(patient_data[, 'RheumatoidArthritis'])
patient_data[, 'Hyperparathyroidism'] = factor(patient_data[, 'Hyperparathyroidism'])
patient_data[, 'Hyperthyroidism'] = factor(patient_data[, 'Hyperthyroidism'])
patient_data[, 'CushingDisease'] = factor(patient_data[, 'CushingDisease'])
patient_data[, 'CeliacDisease'] = factor(patient_data[, 'CeliacDisease'])
patient_data[, 'SecondaryOsteoporosis'] = factor((patient_data[, 'Hyperparathyroidism'] %in% 1 | patient_data[, 'Hyperthyroidism'] %in% 1 | patient_data[, 'CushingDisease'] %in% 1 | patient_data[, 'CeliacDisease'] %in% 1) + 0L, levels = c(0:1))

final_C_table = NULL
SUB_DATASET_NAME = gsub(' .*', '', DATASET_NAME)

for (j in 1:length(DATASET_NAME)) {
  
  test_data = patient_data[patient_data[, 'DATASET'] %in% names(DATASET_NAME)[j],]
  
  SUB_DATASET_NAME[j] = paste0(SUB_DATASET_NAME[j], '\nn = ', nrow(test_data))
  
  C_table = data.frame(x_name = c('Age, sex', 'CXR-risk', 'Age, sex + CXR-risk'),
                       x = 1:3 + (j - 1) * 4.5,
                       y = 0,
                       y_low = 0,
                       y_up = 0,
                       txt = '',
                       p_val = '',
                       bar_pos = NA,
                       group = factor(1:3, levels = 1:3), stringsAsFactors = FALSE)
  
  sub_model.0 = coxph(Surv(test_data[, 'time[all-cause mortality]'], test_data[, 'event[all-cause mortality]']) ~ ., data = test_data[, c('new_pred.1'), drop = FALSE])
  sub_model.1 = coxph(Surv(test_data[, 'time[all-cause mortality]'], test_data[, 'event[all-cause mortality]']) ~ ., data = test_data[, c('predition [[BASIC] Future compression fracture risk]'), drop = FALSE])
  sub_model.2 = coxph(Surv(test_data[, 'time[all-cause mortality]'], test_data[, 'event[all-cause mortality]']) ~ ., data = test_data[, c('new_pred.2'), drop = FALSE])
  
  concordance.0 = concordance(sub_model.0)
  concordance.1 = concordance(sub_model.1)
  concordance.2 = concordance(sub_model.2)
  
  C_table[1, 'y'] = concordance.0[['concordance']]
  C_table[1, 'y_low'] = concordance.0[['concordance']] + qnorm(0.025) * sqrt(concordance.0[['var']])
  C_table[1, 'y_up'] = concordance.0[['concordance']] + qnorm(0.975) * sqrt(concordance.0[['var']])
  
  C_table[2, 'y'] = concordance.1[['concordance']]
  C_table[2, 'y_low'] = concordance.1[['concordance']] + qnorm(0.025) * sqrt(concordance.1[['var']])
  C_table[2, 'y_up'] = concordance.1[['concordance']] + qnorm(0.975) * sqrt(concordance.1[['var']])
  
  p_compare = c_index_test(concordance.1 = concordance.0, concordance.2 = concordance.1, r = cor(sub_model.0[['linear.predictors']], sub_model.1[['linear.predictors']], use = 'pairwise.complete.obs'))
  if (p_compare < 0.001) {C_table[2, 'p_val'] = '***'} else if (p_compare < 0.01) {C_table[2, 'p_val'] = '**'} else if (p_compare < 0.05) {C_table[2, 'p_val'] = '*'}
  
  C_table[3, 'y'] = concordance.2[['concordance']]
  C_table[3, 'y_low'] = concordance.2[['concordance']] + qnorm(0.025) * sqrt(concordance.2[['var']])
  C_table[3, 'y_up'] = concordance.2[['concordance']] + qnorm(0.975) * sqrt(concordance.2[['var']])
  
  p_compare = c_index_test(concordance.1 = concordance.2, concordance.2 = concordance.1, r = cor(sub_model.1[['linear.predictors']], sub_model.2[['linear.predictors']], use = 'pairwise.complete.obs'))
  if (p_compare < 0.001) {C_table[3, 'p_val'] = '***'} else if (p_compare < 0.01) {C_table[3, 'p_val'] = '**'} else if (p_compare < 0.05) {C_table[3, 'p_val'] = '*'}
  
  bar_pos.y = max(C_table[, 'y_up']) + .03
  
  for (k in 2:3) {
    if (C_table[k, 'p_val'] != '') {
      C_table[k, 'bar_pos'] = bar_pos.y
    }
    if (k %in% 3) {bar_pos.y = bar_pos.y + .06}
  }
  
  final_C_table = rbind(final_C_table, C_table)
  
}

final_C_table[, 'txt'] = paste0(formatC(final_C_table[, 'y'], 3, format = 'f'), ' (', formatC(final_C_table[, 'y_low'], 3, format = 'f'), '-', formatC(final_C_table[, 'y_up'], 3, format = 'f'), ')')

#### Plotting ####
gg_p = ggplot(final_C_table, aes(x = x, y = y, fill = group))
gg_p = gg_p + geom_bar(position = "dodge", stat = "identity")
gg_p = gg_p + geom_errorbar(aes(ymin = y_low, ymax= y_up), width = .3, position = position_dodge(.9))
gg_p = gg_p + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.22), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
gg_p = gg_p + scale_x_continuous(limits = c(0.5, max(final_C_table[, 'x']) + 0.5), breaks = final_C_table[, 'x'] - 0.18, labels = final_C_table[, 'x_name'])
gg_p = gg_p + ggtitle('')
gg_p = gg_p + xlab('')
gg_p = gg_p + ylab('C index')
gg_p = gg_p + theme_classic()
gg_p = gg_p + scale_fill_manual(values = col_list.4)
gg_p = gg_p + theme(legend.position = "none",
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_text(size = 12, face = "bold"),
                    axis.text.x = element_text(color = col_list.1[final_C_table[, 'group']], angle = 45, hjust = 1, size = 10, face = "bold"),
                    axis.text.y = element_text(angle = 0, hjust = 1, size = 8, face = "bold"))
gg_p = gg_p + annotate(geom = "text", x = final_C_table[, 'x'], y = 0.02, label = final_C_table[, 'txt'], size = 4, color = "white", angle = 90, fontface = 2, hjust = 0)

for (k in 1:nrow(final_C_table)) {
  if (!is.na(final_C_table[k, 'bar_pos'])) {
    if (k %% 4 == 0) {
      x_pos <- c(final_C_table[k-3, 'x'] + 0.05, final_C_table[k, 'x'] - 0.05)
    } else {
      x_pos <- c(final_C_table[k-1, 'x'] + 0.05, final_C_table[k, 'x'] - 0.05)
    }
    
    gg_p = gg_p + annotate(geom = "line", x = c(x_pos[1], x_pos[1], x_pos[2], x_pos[2]), y = final_C_table[k, 'bar_pos'] + c(0, .02, .02, 0), color = 'black', size = 0.3)
    gg_p = gg_p + annotate(geom = "text", x = mean(x_pos), y = final_C_table[k, 'bar_pos'] + .04, label = final_C_table[k, 'p_val'], color = 'black', size = 3, fontface = 2)
    
  }
}

gg_p = gg_p + annotate(geom = "text", x = 2.5, y = 1.13, label = SUB_DATASET_NAME[1], color = 'black', size = 4, fontface = 2)
bar_p = gg_p + annotate(geom = "text", x = 7, y = 1.13, label = SUB_DATASET_NAME[2], color = 'black', size = 4, fontface = 2)
bar_p

#### Merge plots ####
