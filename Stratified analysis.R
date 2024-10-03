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
plot_path_1 = './Figure_5_Stratified analysis_c_index.pdf'
plot_path_2 = './Figure_5_Stratified analysis_auc.pdf'
plot_path_3 = './Figure_5_Stratified analysis_prauc.pdf'

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
patient_data = patient_data[!duplicated(patient_data[, 'CNO']),]

## Survival analysis
patient_data[, 'Y1'] = 0L
patient_data[patient_data[, 'event[Fracture_outcome]'] %in% 1 & patient_data[, 'time[Fracture_outcome]'] <= 365.25 * critical_cut, 'Y1'] = 1L
patient_data[patient_data[, 'event[Fracture_outcome]'] %in% 0 & patient_data[, 'time[Fracture_outcome]'] < 365.25 * critical_cut, 'Y1'] = NA

patient_data[, 'event'] = patient_data[, 'event[Fracture_outcome]']
patient_data[, 'time'] = patient_data[, 'time[Fracture_outcome]']
patient_data[, 'time'] = (patient_data[, 'time'] + 0.5) / 365.25
patient_data = patient_data[patient_data[, 'time'] >= 0,] 
patient_data[patient_data[, 'time'] >= critical_cut, 'event'] = 0L
patient_data[patient_data[, 'time'] >= critical_cut, 'time'] = critical_cut

## Format variables
patient_data = patient_data[!patient_data[, 'RheumatoidArthritis'] %in% NA,]

patient_data[, 'V0'] = factor('all samples', levels = 'all samples')
patient_data[, 'AGE_G'] = cut(patient_data[, 'AGE'], breaks = c(-Inf, 59.999, 69.999, 79.999, Inf), labels = c('<60 y/o', '60-69 y/o', '70-79 y/o', '>79 y/o'))
patient_data[, 'SEX_G'] = factor(patient_data[, 'GENDER'])
patient_data[, 'CXR[MODALITY]'] = factor(patient_data[, 'CXR[MODALITY]'], labels = c('computed radiography', 'digital radiography'), levels = c('CR', 'DX'))
patient_data[, 'CXR[VIEW]'] = factor(patient_data[, 'CXR[VIEW]'])
# patient_data[, 'SEXAGE_G'] = (as.integer(patient_data[, 'AGE_G']) * 2 + as.integer(patient_data[, 'SEX_G']) - 2) %>% factor()
# levels(patient_data[, 'SEXAGE_G']) = c('female <60 y/o', 'male <60 y/o', 'female 60-69 y/o', 'male 60-69 y/o', 'female 70-79 y/o', 'male 70-79 y/o', 'female >79 y/o', 'male >79 y/o')

# patient_data[, 'Fracture'] = factor(patient_data[, 'Fracture'])
patient_data[, 'RheumatoidArthritis'] = factor(patient_data[, 'RheumatoidArthritis'])
patient_data[, 'Hyperparathyroidism'] = factor(patient_data[, 'Hyperparathyroidism'])
patient_data[, 'Hyperthyroidism'] = factor(patient_data[, 'Hyperthyroidism'])
patient_data[, 'CushingDisease'] = factor(patient_data[, 'CushingDisease'])
patient_data[, 'CeliacDisease'] = factor(patient_data[, 'CeliacDisease'])
patient_data[, 'SecondaryOsteoporosis'] = factor((patient_data[, 'Hyperparathyroidism'] %in% 1 | patient_data[, 'Hyperthyroidism'] %in% 1 | patient_data[, 'CushingDisease'] %in% 1 | patient_data[, 'CeliacDisease'] %in% 1) + 0L, levels = c(0:1))

#### Plotting ####
# C-index bar plot

p_list_1 = list()

for (h in 1:length(DATASET_NAME)) {
  
  sub_data = patient_data[patient_data[, 'DATASET'] %in% names(DATASET_NAME)[h],]
  
  p_list_1[[names(DATASET_NAME)[h]]] = list()
  
  for (i in 1:length(y_name_list)) {
    
    x_name = x_name_list
    y_name = y_name_list[i]
    
    x_pos = 1
    summary_list = list()
    
    for (l in 1:length(intrested_var)) {
      
      lvl.var = patient_data[, intrested_var[l]] %>% factor() %>% levels() 
      
      if (length(lvl.var) %in% 4) {
        col_list = col_list.1
      } else if (length(lvl.var) %in% 1) {
        col_list = col_list.3
      } else if (length(lvl.var) %in% 3) {
        col_list = col_list.4
      } else if (!grepl('SEX', intrested_var[l])) {
        col_list = col_list.5
      } else {
        col_list = paste0(rep(col_list.2, 4), rep(c('FF', 'DF', 'BF', '9F'), each = 2))
      }
      
      for (m in 1:length(lvl.var)) {
        
        sub_sub_test_data = sub_data[sub_data[, intrested_var[l]] %in% lvl.var[m],]
        
        C_idx_model = coxph(Surv(time, event) ~ `predition [[BASIC] Future compression fracture risk]`, data = sub_sub_test_data)
        concordance = concordance(C_idx_model)
        
        if (!'try-error' %in% class(concordance)) {
          
          summary_list[[length(summary_list) + 1]] = data.frame(var1 = factor(intrested_var[l], levels = intrested_var),
                                                                var2 = lvl.var[m],
                                                                col = col_list[m],
                                                                x = x_pos,
                                                                val = concordance[['concordance']],
                                                                ci.l = concordance[['concordance']] + qnorm(0.025) * sqrt(concordance[['var']]),
                                                                ci.u = concordance[['concordance']] + qnorm(0.975) * sqrt(concordance[['var']]),
                                                                stringsAsFactors = FALSE)
          
        }
        x_pos = x_pos + 1
      }
      x_pos = x_pos + 0.5
    }
    
    my_summary = do.call('rbind', summary_list)
    my_summary = my_summary[order(my_summary[, 'x']),]
    
    for (c in 1:nrow(my_summary)) {if (!my_summary[c, 'ci.u'] %in% NaN & my_summary[c, 'ci.u'] > 1) {my_summary[c, 'ci.u'] = 1}}
    
    my_summary[, 'x_name'] = my_summary[, 'var2']
    my_summary[, 'txt'] = paste0(formatC(my_summary[, 'val'], 3, format = 'f'), ' (', formatC(my_summary[, 'ci.l'], 3, format = 'f'), '-', formatC(my_summary[, 'ci.u'], 3, format = 'f'), ')')
    
    main_txt = main_name_list[i]
    bar_p = ggplot(my_summary, aes(x = x, y = val, fill = col))
    bar_p = bar_p + geom_bar(position = "dodge", stat = "identity")
    bar_p = bar_p + geom_errorbar(aes(ymin = ci.l, ymax = ci.u), width = .4, position = position_dodge(.9))
    bar_p = bar_p + theme_classic()
    bar_p = bar_p + labs(title = '', x = '', y = 'C-index (95% CI)')
    bar_p = bar_p + theme(plot.title = element_text(color = "#000000", size = 12),
                          legend.position = "none",
                          axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, size = 8)) + labs(fill = '')
    bar_p = bar_p + scale_fill_manual(values = factor(my_summary[,"col"]) %>% levels())
    bar_p = bar_p + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
    bar_p = bar_p + scale_x_continuous(limits = c(0.5, max(my_summary[,'x']) + 0.5),
                                       breaks = my_summary[, 'x'],
                                       labels = my_summary[, 'x_name'])
    bar_p = bar_p + annotate(geom = "text", x = my_summary[, 'x'], y = 0.02, label = my_summary[, 'txt'], size = 2.5, color = "white", angle = 90, fontface = 2, hjust = 0)
    p_list_1[[names(DATASET_NAME)[h]]][[y_name]] = bar_p
  }
}

c_idx_list = p_list_1

# AUC bar plot

# p_list_2 = list()
# 
# for (h in 1:length(DATASET_NAME)) {
#   
#   sub_data = patient_data[patient_data[, 'DATASET'] %in% names(DATASET_NAME)[h],]
#   
#   p_list_2[[names(DATASET_NAME)[h]]] = list()
#   
#   for (i in 1:length(y_name_list)) {
#     
#     x_name = x_name_list
#     y_name = y_name_list[i]
#     
#     x_pos = 1
#     summary_list = list()
#     
#     for (l in 1:length(intrested_var)) {
#       
#       lvl.var = patient_data[, intrested_var[l]] %>% factor() %>% levels()
#       
#       if (length(lvl.var) %in% 4) {
#         col_list = col_list.1
#       } else if (length(lvl.var) %in% 1) {
#         col_list = col_list.3
#       } else if (length(lvl.var) %in% 3) {
#         col_list = col_list.4
#       } else if (!grepl('SEX', intrested_var[l])) {
#         col_list = col_list.5
#       } else {
#         col_list = paste0(rep(col_list.2, 4), rep(c('FF', 'DF', 'BF', '9F'), each = 2))
#       }
#       
#       for (m in 1:length(lvl.var)) {
#         
#         sub_sub_test_data = sub_data[sub_data[, intrested_var[l]] %in% lvl.var[m],]
#         roc_curve.test = try(roc(response = sub_sub_test_data[, y_name], predictor = sub_sub_test_data[, x_name]), silent = TRUE)
#         
#         if (!'try-error' %in% class(roc_curve.test)) {
#           
#           summary_list[[length(summary_list) + 1]] = data.frame(var1 = factor(intrested_var[l], levels = intrested_var),
#                                                                 var2 = lvl.var[m],
#                                                                 col = col_list[m],
#                                                                 x = x_pos,
#                                                                 val = ci.auc(roc_curve.test)[2],
#                                                                 ci.l = ci.auc(roc_curve.test)[1],
#                                                                 ci.u = ci.auc(roc_curve.test)[3],
#                                                                 stringsAsFactors = FALSE)
#           
#         }
#         x_pos = x_pos + 1
#       }
#       x_pos = x_pos + 0.5
#     }
#     
#     my_summary = do.call('rbind', summary_list)
#     my_summary = my_summary[order(my_summary[, 'x']),]
#     
#     my_summary[, 'x_name'] = my_summary[, 'var2']
#     my_summary[, 'txt'] = paste0(formatC(my_summary[, 'val'], 3, format = 'f'), ' (', formatC(my_summary[, 'ci.l'], 3, format = 'f'), '-', formatC(my_summary[, 'ci.u'], 3, format = 'f'), ')')
#     
#     main_txt = main_name_list[i]
#     bar_p = ggplot(my_summary, aes(x = x, y = val, fill = col))
#     bar_p = bar_p + geom_bar(position = "dodge", stat = "identity")
#     bar_p = bar_p + geom_errorbar(aes(ymin = ci.l, ymax = ci.u), width = .4, position = position_dodge(.9))
#     bar_p = bar_p + theme_classic()
#     bar_p = bar_p + labs(title = '', x = '', y = 'AUC (95% CI)')
#     bar_p = bar_p + theme(plot.title = element_text(color = "#000000", size = 12),
#                           legend.position = "none",
#                           axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, size = 8)) + labs(fill = '')
#     bar_p = bar_p + scale_fill_manual(values = factor(my_summary[,"col"]) %>% levels())
#     bar_p = bar_p + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
#     bar_p = bar_p + scale_x_continuous(limits = c(0.5, max(my_summary[,'x']) + 0.5),
#                                        breaks = my_summary[, 'x'],
#                                        labels = my_summary[, 'x_name'])
#     bar_p = bar_p + annotate(geom = "text", x = my_summary[, 'x'], y = 0.02, label = my_summary[, 'txt'], size = 2.5, color = "white", angle = 90, fontface = 2, hjust = 0)
#     p_list_2[[names(DATASET_NAME)[h]]][[y_name]] = bar_p
#   }
# }
# 
# auc_p_list = p_list_2

# PRAUC bar plot

# p_list_3 = list()
# cut_list = rep(NA, length(x_name_list))
# 
# for (h in 1:length(DATASET_NAME)) {
#   
#   sub_data = patient_data[patient_data[, 'DATASET'] %in% names(DATASET_NAME)[h],]
#   sub_valid_data = patient_data[patient_data[, 'DATASET'] %in% 'valid',]
#   
#   p_list_3[[names(DATASET_NAME)[h]]] = list()
#   
#   for (i in 1:length(y_name_list)) {
#     
#     x_name = x_name_list
#     y_name = y_name_list[i]
#     
#     x_pos = 1
#     summary_list = list()
#     
#     if (cut_list[h] %in% NA) {
#       
#       roc_curve.valid = roc(response = sub_valid_data[, y_name], predictor = sub_valid_data[, x_name])
#       
#       tab = table(factor(sub_valid_data[, y_name], levels = 0:1))
#       
#       roc_curve.valid[['ppv']] = tab[2] * roc_curve.valid[['sensitivities']] / (tab[1] * (1 - roc_curve.valid[['specificities']]) + tab[2] * roc_curve.valid[['sensitivities']] + 1e-15)
#       roc_curve.valid[['f1']] = 2 * roc_curve.valid[['ppv']] * roc_curve.valid[['sensitivities']] / (roc_curve.valid[['ppv']] + roc_curve.valid[['sensitivities']] + 1e-15)
#       
#       pb = txtProgressBar(max = length(roc_curve.valid[['thresholds']]), style = 3)
#       for (l in 2:length(roc_curve.valid[['thresholds']])) {
#         if (roc_curve.valid[['ppv']][l] < max(roc_curve.valid[['ppv']][-l:-length(roc_curve.valid[['thresholds']])])) {roc_curve.valid[['ppv']][l] = max(roc_curve.valid[['ppv']][-l:-length(roc_curve.valid[['thresholds']])])}
#         setTxtProgressBar(pb, l)
#       }
#       close(pb)
#       
#       best_pos = which.max(roc_curve.valid[['f1']])
#       cut_list[h] = roc_curve.valid[['thresholds']][best_pos]
#       
#     }
#     
#     for (l in 1:length(intrested_var)) {
#       
#       lvl.var = patient_data[, intrested_var[l]] %>% factor() %>% levels()
#       
#       if (length(lvl.var) %in% 4) {
#         col_list = col_list.1
#       } else if (length(lvl.var) %in% 1) {
#         col_list = col_list.3
#       } else if (length(lvl.var) %in% 3) {
#         col_list = col_list.4
#       } else if (!grepl('SEX', intrested_var[l])) {
#         col_list = col_list.5
#       } else {
#         col_list = paste0(rep(col_list.2, 4), rep(c('FF', 'DF', 'BF', '9F'), each = 2))
#       }
#       
#       for (m in 1:length(lvl.var)) {
#         
#         sub_sub_test_data = sub_data[sub_data[, intrested_var[l]] %in% lvl.var[m],]
#         current_table = table(sub_sub_test_data[, x_name] >= cut_list[h], factor(sub_sub_test_data[, y_name]))
#         prev = prop.table(current_table)[, 2] %>% sum()
#         ppv = prop.table(current_table[2,])[2]
#         n = current_table[2,] %>% sum()
#         
#         summary_list[[length(summary_list) + 1]] = data.frame(var1 = factor(intrested_var[l], levels = intrested_var),
#                                                               var2 = lvl.var[m],
#                                                               col = col_list[m],
#                                                               x = x_pos,
#                                                               val = ppv,
#                                                               ci.l = max(ppv + qnorm(0.025) * sqrt(ppv * (1 - ppv) / n), 0),
#                                                               ci.u = min(ppv + qnorm(0.975) * sqrt(ppv * (1 - ppv) / n), 1),
#                                                               prev = prev,
#                                                               n = n,
#                                                               stringsAsFactors = FALSE)
#         
#         x_pos = x_pos + 1
#       }
#       x_pos = x_pos + 0.5
#     }
#     
#     my_summary = do.call('rbind', summary_list)
#     my_summary = my_summary[order(my_summary[, 'x']),]
#     
#     my_summary[, 'x_name'] = my_summary[, 'var2']
#     my_summary[, 'txt'] = paste0(formatC(my_summary[, 'val'], 3, format = 'f'), ' (', formatC(my_summary[, 'ci.l'], 3, format = 'f'), '-', formatC(my_summary[, 'ci.u'], 3, format = 'f'), ')')
#     my_summary[, 'txt2'] = paste0('Prev = ', formatC(my_summary[, 'prev'] * 100, 1, format = 'f'), '% (n = ', my_summary[, 'n'], ')')
#     
#     main_txt = main_name_list[i]
#     bar_p = ggplot(my_summary, aes(x = x, y = val, fill = col))
#     bar_p = bar_p + geom_bar(position = 'dodge', stat = 'identity')
#     bar_p = bar_p + geom_errorbar(aes(ymin = ci.l, ymax = ci.u), width = .4, position = position_dodge(.9))
#     bar_p = bar_p + theme_classic()
#     bar_p = bar_p + labs(title = main_txt, x = '', y = 'Positive predictive value')
#     bar_p = bar_p + theme(plot.title = element_text(color = "#000000", size = 12),
#                           legend.position = 'none',
#                           axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, size = 8)) + labs(fill = '')
#     
#     bar_p = bar_p + scale_fill_manual(values = factor(my_summary[, 'col']) %>% levels())
#     bar_p = bar_p + scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
#     bar_p = bar_p + scale_x_continuous(limits = c(0.5, max(my_summary[, 'x']) + 0.5),
#                                        breaks = my_summary[, 'x'],
#                                        labels = my_summary[, 'x_name'])
#     
#     pos = my_summary[, 'val'] >= 0.4
#     
#     bar_p = bar_p + annotate(geom = "text", x = my_summary[pos, 'x'], y = 0.02, label = my_summary[pos, 'txt'], size = 2.5, color = "white", angle = 90, fontface = 2, hjust = 0)
#     bar_p = bar_p + annotate(geom = "text", x = my_summary[!pos, 'x'], y = my_summary[!pos, 'ci.u'] + 0.02, label = my_summary[!pos, 'txt'], size = 2.5, color = "black", angle = 90, fontface = 2, hjust = 0)
#     bar_p = bar_p + annotate(geom = "text", x = my_summary[,'x'], y = my_summary[,'ci.u'] + 0.02, label = my_summary[,'txt2'], size = 2.5, color = "black", angle = 90, fontface = 2, hjust = 0)
#     
#     p_list_3[[names(DATASET_NAME)[h]]][[y_name]] = bar_p
#   }
# }
# 
# ppv_p_list = p_list_3

#### Merge plots ####
# # pp.1 = arrangeGrob(c_idx_list[[1]][[1]], c_idx_list[[2]][[1]], ncol = 2)
# # 
# # final_p = ggdraw()
# # final_p = final_p + draw_plot(pp.1, x = 0, y = 0, width = 1, height = 0.97)
# # final_p = final_p + draw_plot_label(DATASET_NAME, c(0.005, 0.505), c(1, 1), size = 15, hjust = 0, fontface = 2)
# 
# # pdf(plot_path_1, width = 8, height = 6)
# # print(final_p)
# # dev.off()

# pp.2 = arrangeGrob(auc_p_list[[1]][[1]], auc_p_list[[2]][[1]], nrow = 2)
# 
# final_p = ggdraw()
# final_p = final_p + draw_plot(pp.2, x = 0, y = 0, width = 1, height = 0.97)
# final_p = final_p + draw_plot_label(DATASET_NAME, c(0.005, 0.005), c(1, 0.505), size = 15, hjust = 0, fontface = 2)

# pdf(plot_path_2, width = 8, height = 12)
# print(final_p)
# dev.off()

# pp.3 = arrangeGrob(ppv_p_list[[1]][[1]], ppv_p_list[[2]][[1]], nrow = 2)
# 
# final_p = ggdraw()
# final_p = final_p + draw_plot(pp.3, x = 0, y = 0, width = 1, height = 0.97)
# final_p = final_p + draw_plot_label(DATASET_NAME, c(0.005, 0.005), c(1, 0.505), size = 15, hjust = 0, fontface = 2)

# pdf(plot_path_3, width = 8, height = 12)
# print(final_p)
# dev.off()