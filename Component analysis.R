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
data_path = '/home/chinlin/Multi-center-database/002. data/002. project data/048. project-048/analysis/study 014/patient_data.RData'
plot_path_1 = './Component_analysis_Thoracic.pdf'

#### Settings ####
set.seed(0)

n_importance = 8

DATASET_NAME = c('internal-test' = 'Internal validation set', 'community-test' = 'External validation set')

x_name = 'predition [[CCT] Thoracic aortic aneurysm]'

intrested_var_list = list()

intrested_var_list[[1]] = list(var = c(
                                       'CXR[POSITION]', 'AGE', 'GENDER', 
                                       'DM', 'HTN', 'CKD', 'HLP', 'HF', 'CAD', 'COPD'
                                       ),
                               name = c(
                                        'Data source', 'Age', 'Sex', 
                                        'History of diabetes mellitus', 'History of hypertension', 'History of chronic kidney disease', 'History of hyperlipidemia', 'History of heart failure', 'History of coronary artery disease', 'History of chronic obstructive pulmonary disease'
                                        ))

intrested_var_list[[2]] = list(var = c(
                                       'CXR[MODALITY]', 'CXR[VIEW]', 
                                       '[RADIO] Consolidation change', '[RADIO] Pneumonia', '[RADIO] Emphysematous change', '[RADIO] Pneumothorax', 
                                       '[RADIO] Atelectasis', '[RADIO] Scalloping of the diaphragm', '[RADIO] Costophrenic angle blunting', 
                                       '[RADIO] Pleural effusion', '[RADIO] Atherosclerosis', '[RADIO] Cardiomegaly', '[RADIO] Prominence of hilar shadow',
                                       '[RADIO] Pulmonary edema', '[RADIO] Aneurysm', '[RADIO] Degenerative joint disease', '[RADIO] Fracture', 
                                       '[RADIO] Spondylosis', '[RADIO] Osteophyte formation', '[RADIO] Osteoporosis', '[RADIO] Osteoarthritis',
                                       '[RADIO] Widening of the mediastinum', '[RADIO] Malignancy', '[RADIO] Inflammatory', '[RADIO] Pigtail or drainage',
                                       '[RADIO] Sternotomy', '[RADIO] Port a implantation', '[RADIO] Perm catheter insertion', '[RADIO] Pacemaker',
                                       '[RADIO] Tracheostomy', '[RADIO] Vertebroplasty', '[RADIO] Endotracheal tube', '[RADIO] Nasogastric tube'
                                       ),
                               name = c(
                                        'Modality (DR or CR)', 'PA or AP view', 
                                        'Consolidation change', 'Pneumonia', 'Emphysematous change', 'Pneumothorax', 
                                        'Atelectasis', 'Scalloping of the diaphragm', 'Costophrenic angle blunting',
                                        'Pleural effusion', 'Atherosclerosis', 'Cardiomegaly', 'Prominence of hilar shadow',
                                        'Pulmonary edema', 'Aneurysm', 'Degenerative joint disease', 'Fracture',
                                        'Spondylosis', 'Osteophyte formation', 'Osteoporosis', 'Osteoarthritis',
                                        'Widening of the mediastinum', 'Malignancy', 'Inflammatory', 'Pigtail or drainage',
                                        'Sternotomy', 'Port a implantation', 'Perm catheter insertion', 'Pacemaker',
                                        'Tracheostomy', 'Vertebroplasty', 'Endotracheal tube', 'Nasogastric tube'
                                        ))
 
intrested_var_list[[3]] = list(var = c(
                                       'CXR[POSITION]', 'AGE', 'GENDER', 
                                       'DM', 'HTN', 'CKD', 'HLP', 'HF', 'CAD', 'COPD',
                                       'CXR[MODALITY]', 'CXR[VIEW]', 
                                       '[RADIO] Consolidation change', '[RADIO] Pneumonia', '[RADIO] Emphysematous change', '[RADIO] Pneumothorax', 
                                       '[RADIO] Atelectasis', '[RADIO] Scalloping of the diaphragm', '[RADIO] Costophrenic angle blunting', 
                                       '[RADIO] Pleural effusion', '[RADIO] Atherosclerosis', '[RADIO] Cardiomegaly', '[RADIO] Prominence of hilar shadow',
                                       '[RADIO] Pulmonary edema', '[RADIO] Aneurysm', '[RADIO] Degenerative joint disease', '[RADIO] Fracture', 
                                       '[RADIO] Spondylosis', '[RADIO] Osteophyte formation', '[RADIO] Osteoporosis', '[RADIO] Osteoarthritis',
                                       '[RADIO] Widening of the mediastinum', '[RADIO] Malignancy', '[RADIO] Inflammatory', '[RADIO] Pigtail or drainage',
                                       '[RADIO] Sternotomy', '[RADIO] Port a implantation', '[RADIO] Perm catheter insertion', '[RADIO] Pacemaker',
                                       '[RADIO] Tracheostomy', '[RADIO] Vertebroplasty', '[RADIO] Endotracheal tube', '[RADIO] Nasogastric tube'
                                       ),
                               name = c(
                                        'Data source', 'Age', 'Sex', 
                                        'History of diabetes mellitus', 'History of hypertension', 'History of chronic kidney disease', 'History of hyperlipidemia', 'History of heart failure', 'History of coronary artery disease', 'History of chronic obstructive pulmonary disease',
                                        'Modality (DR or CR)', 'PA or AP view',
                                        'Consolidation change', 'Pneumonia', 'Emphysematous change', 'Pneumothorax', 
                                        'Atelectasis', 'Scalloping of the diaphragm', 'Costophrenic angle blunting',
                                        'Pleural effusion', 'Atherosclerosis', 'Cardiomegaly', 'Prominence of hilar shadow',
                                        'Pulmonary edema', 'Aneurysm', 'Degenerative joint disease', 'Fracture',
                                        'Spondylosis', 'Osteophyte formation', 'Osteoporosis', 'Osteoarthritis',
                                        'Widening of the mediastinum', 'Malignancy', 'Inflammatory', 'Pigtail or drainage',
                                        'Sternotomy', 'Port a implantation', 'Perm catheter insertion', 'Pacemaker',
                                        'Tracheostomy', 'Vertebroplasty', 'Endotracheal tube', 'Nasogastric tube'
                                        ))

intrested_var_list[[4]] = list(var = c(
                                       'CXR[POSITION]', 'AGE', 'GENDER', 
                                       'DM', 'HTN', 'CKD', 'HLP', 'HF', 'CAD', 'COPD',
                                       'predition [[CCT] Thoracic aortic aneurysm]'
                                       ),
                               name = c(
                                        'Data source', 'Age', 'Sex', 
                                        'History of diabetes mellitus', 'History of hypertension', 'History of chronic kidney disease', 'History of hyperlipidemia', 'History of heart failure', 'History of coronary artery disease', 'History of chronic obstructive pulmonary disease',
                                        'DLM prediction'
                                        ))

names(intrested_var_list) = c('Patient characteristics', 'CXR features', 'Combination', 'DLM + patient data')

#### Processes ####
load(data_path)

## Randomly shuffle data and remove duplicates
patient_data = patient_data[sample(1:nrow(patient_data)),]
patient_data = patient_data[!duplicated(patient_data[, 'CNO']),]

## Survival analysis
patient_data[ ,'Y1'] = 0L
patient_data[patient_data[, '[CCT] Thoracic aortic aneurysm'] %in% 1, 'Y1'] = 1L
patient_data[patient_data[, '[CCT] Thoracic aortic aneurysm'] %in% NA, 'Y1'] = NA

patient_data = patient_data[!patient_data[, 'Y1'] %in% NA,]

#### Plotting ####
XGB_list = list()

for (j in 1:length(intrested_var_list)) {
  
  intrested_var = intrested_var_list[[j]][['var']]
  intrested_var_name = intrested_var_list[[j]][['name']]
  names(intrested_var_name) = intrested_var
  
  for (var in intrested_var) {
    if (is.factor(patient_data[[var]])) {
      if (length(levels(patient_data[[var]])) < 2) {
        stop(paste("Factor variable", var, "has less than 2 levels"))
      }
    }
  }
  
  train_data = patient_data[patient_data[, 'DATASET'] %in% 'train', intrested_var] %>% model.matrix(~ ., data = .)
  train_label = patient_data[patient_data[, 'DATASET'] %in% 'train', 'Y1']
  
  test_data = patient_data[patient_data[, 'DATASET'] %in% 'valid', intrested_var] %>% model.matrix(~ ., data = .)
  test_label = patient_data[patient_data[, 'DATASET'] %in% 'valid', 'Y1']
  
  all_data = patient_data[, intrested_var] %>% model.matrix(~ ., data = .)
  all_label = patient_data[, 'Y1']
  
  dfold.1 = xgb.DMatrix(data = as.matrix(train_data), label = train_label)
  dfold.2 = xgb.DMatrix(data = as.matrix(test_data), label = test_label)
  dfold.3 = xgb.DMatrix(data = as.matrix(all_data), label = all_label)
  
  XGB_model.1 = xgb.train(data = dfold.1,
                          watchlist = list(train = dfold.1, test = dfold.2),
                          eta = 0.05,  nrounds = 300, max.depth = 5, early_stopping_rounds = 10,
                          nthread = 5, verbose = FALSE, eval_metric = "auc", objective = "binary:logistic")
  
  patient_data[, paste0('combine_', j)] = predict(XGB_model.1, dfold.3)
  
  XGB_importance = xgb.importance(model = XGB_model.1)
  XGB_data.1 = data.frame(name = XGB_importance[['Feature']], importance = XGB_importance[['Gain']], sign = sign(XGB_importance[['Gain']]), stringsAsFactors = FALSE)
  XGB_data.1[, 'importance'] = XGB_data.1[, 'importance'] / max(XGB_data.1[, 'importance'])
  
  XGB_data.2 = XGB_data.1
  
  XGB_data = merge(XGB_data.1, XGB_data.2, by = 'name', all = TRUE)
  XGB_data[XGB_data[, 'importance.x'] %in% NA, 'importance.x'] = 0
  XGB_data[XGB_data[, 'importance.y'] %in% NA, 'importance.y'] = 0
  XGB_data[, 'importance'] = (XGB_data[, 'importance.x'] + XGB_data[, 'importance.y']) / 2
  XGB_data[, 'revise_name'] = XGB_data[, 'name'] %>% gsub('^`', '', .) %>% gsub('`.*$', '', .)
  XGB_data[, 'revise_name'] = intrested_var_name[XGB_data[, 'revise_name']]
  
  for (i in 1:nrow(XGB_data)) {
    if (XGB_data[i, 'revise_name'] %in% NA) {
      for (q in 1:length(intrested_var_name)) {
        if (grepl(names(intrested_var_name)[q], XGB_data[i, 'name'])) {
          XGB_data[i, 'revise_name'] = intrested_var_name[q]
        }
      }
    }
  }
  
  SUMMARY_XGB = tapply(XGB_data[,'importance'], XGB_data[,'revise_name'], sum)
  
  FIANL_XGB_DATA = data.frame(name = names(SUMMARY_XGB), importance = as.numeric(SUMMARY_XGB), stringsAsFactors = FALSE)
  FIANL_XGB_DATA = rbind(data.frame(name = rep('', 20), importance = rep(0, 20), stringsAsFactors = FALSE), FIANL_XGB_DATA)
  FIANL_XGB_DATA[, 'importance'] = FIANL_XGB_DATA[, 'importance'] / max(FIANL_XGB_DATA[, 'importance'])
  FIANL_XGB_DATA = FIANL_XGB_DATA[order(FIANL_XGB_DATA[, 'importance'], decreasing = FALSE),]
  FIANL_XGB_DATA = tail(FIANL_XGB_DATA, n_importance)
  
  FIANL_XGB_DATA[, 'x'] = 1:nrow(FIANL_XGB_DATA)
  FIANL_XGB_DATA[, 'txt'] = formatC(FIANL_XGB_DATA[, 'importance'] * 100, 1, format = 'f') %>% paste0(., '%')
  FIANL_XGB_DATA[FIANL_XGB_DATA[, 'name'] %in% '', 'txt'] = ''
  
  if (j %in% 1) {
    
    FIANL_XGB_DATA[, 'group'] = '1'
    input_col = c("#4DBBD592")
    
  } else if (j %in% 2) {
    
    FIANL_XGB_DATA[, 'group'] = '1'
    input_col = c("#E64B3592")
    
  } else if (j %in% 3) {
    
    FIANL_XGB_DATA[, 'group'] = '1'
    FIANL_XGB_DATA[FIANL_XGB_DATA[, 'name'] %in% intrested_var_list[[2]][['name']], 'group'] = '2'
    input_col = c("#4DBBD592", "#E64B3592")
    
  }
  
  p_XGB = ggplot(FIANL_XGB_DATA, aes(x = as.factor(x), y = importance, fill = as.factor(group), group = as.factor(group))) 
  p_XGB = p_XGB + geom_bar(stat = "identity", position = "dodge")
  p_XGB = p_XGB + ylab("Related importance") + xlab("") + ggtitle(names(intrested_var_list)[j]) + labs(fill = '')
  p_XGB = p_XGB + theme_classic()
  p_XGB = p_XGB + theme(plot.title = element_text(color = "#000000", size = 15, hjust = 0, margin = margin(), debug = FALSE),
                        legend.position = "none",
                        axis.text.x = element_text(size = 10, face = 1, margin = margin(), debug = FALSE),
                        axis.text.y = element_text(size = 10, face = 1, margin = margin(), debug = FALSE))
  p_XGB = p_XGB + scale_fill_manual(values = c(input_col))
  
  p_XGB = p_XGB + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.32), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
  p_XGB = p_XGB + scale_x_discrete(labels = FIANL_XGB_DATA[, 'name'])
  p_XGB = p_XGB + annotate('text', FIANL_XGB_DATA[, 'x'], FIANL_XGB_DATA[, 'importance'] + 0.02, label = FIANL_XGB_DATA[, 'txt'], hjust = 0)
  p_XGB = p_XGB + coord_flip()
  
  XGB_list[[j]] = p_XGB
  
}

thoracic_p_auc_list = list()

auc_data = data.frame(name = c(paste0('combine_', 1:4), x_name, intrested_var_list[[3]][['var']]),
                      show_name = c(paste0('Xgboost[', names(intrested_var_list), ']'), 'DLM prediction', intrested_var_list[[3]][['name']]),
                      colour = c('1', '1', '1', '2', '2', rep('3', intrested_var_list[[1]][['var']] %>% length()), rep('4', intrested_var_list[[2]][['var']] %>% length())),
                      auc = NA,
                      auc.l = NA,
                      auc.u = NA,
                      stringsAsFactors = FALSE)

for (SET in names(DATASET_NAME)) {
  
  for (q in 1:nrow(auc_data)) {
    
    if (auc_data[q, 'colour'] %in% c('1', '2')) {
      
      patient_data[patient_data[, 'DATASET'] %in% SET, 'x'] = patient_data[patient_data[, 'DATASET'] %in% SET, auc_data[q, 'name']]
      
    } else {
      
      train_data = patient_data[patient_data[, 'DATASET'] %in% 'train', auc_data[q, 'name'], drop = FALSE] %>% model.matrix(~ ., data = .)
      train_label = patient_data[patient_data[, 'DATASET'] %in% 'train', 'Y1']
      
      valid_data = patient_data[patient_data[, 'DATASET'] %in% 'valid', auc_data[q, 'name'], drop = FALSE] %>% model.matrix(~ ., data = .)
      valid_label = patient_data[patient_data[, 'DATASET'] %in% 'valid', 'Y1']
      
      test_data = patient_data[patient_data[, 'DATASET'] %in% SET, auc_data[q, 'name'], drop = FALSE] %>% model.matrix(~ ., data = .)
      test_label = patient_data[patient_data[, 'DATASET'] %in% SET, 'Y1']
      
      dfold.1 = xgb.DMatrix(data = as.matrix(train_data), label = train_label)
      dfold.2 = xgb.DMatrix(data = as.matrix(valid_data), label = valid_label)
      dfold.3 = xgb.DMatrix(data = as.matrix(test_data), label = test_label)
      
      XGB_model.1 = xgb.train(data = dfold.1,
                              watchlist = list(train = dfold.1, test = dfold.2),
                              eta = 0.05,  nrounds = 300, max.depth = 5, early_stopping_rounds = 10,
                              nthread = 5, verbose = FALSE, eval_metric = "auc", objective = "binary:logistic")
      
      patient_data[patient_data[, 'DATASET'] %in% SET, 'x'] = predict(XGB_model.1, dfold.3)
      
    }
    
    
    
    sub_data = patient_data[patient_data[, 'DATASET'] %in% SET,]
    
    # AUC
    roc_curve = roc(response = sub_data[, 'Y1'], predictor = as.numeric(sub_data[, 'x']))

    auc_data[q, 'auc'] = roc_curve[['auc']]
    auc_data[q, 'auc.l'] = ci.auc(roc_curve)[1]
    auc_data[q, 'auc.u'] = ci.auc(roc_curve)[3]
    
    # C-index
    # C_idx_model = coxph(Surv(time, event) ~ x, data = sub_data)
    # concordance = concordance(C_idx_model)
    # 
    # auc_data[q, 'auc'] = concordance[['concordance']]
    # auc_data[q, 'auc.l'] = concordance[['concordance']] + qnorm(0.025) * sqrt(concordance[['var']])
    # auc_data[q, 'auc.u'] = concordance[['concordance']] + qnorm(0.975) * sqrt(concordance[['var']])
    
  }
  
  auc_data[auc_data[, 'auc.u'] > 1, 'auc.u'] = 1
  auc_data = auc_data[order(auc_data[, 'auc'], decreasing = TRUE),]
  auc_data[, 'txt'] = paste0(formatC(auc_data[, 'auc'], 3, format = 'f'), ' (', formatC(auc_data[, 'auc.l'], 3, format = 'f'), '-', formatC(auc_data[, 'auc.u'], 3, format = 'f'), ')')
  
  sub_auc_data = auc_data[!auc_data[, 'auc.l'] %in% NA,]
  sub_auc_data = sub_auc_data[!duplicated(sub_auc_data[, 'show_name']),]
  sub_auc_data = head(sub_auc_data, 10)
  sub_auc_data[, 'x'] = 1:nrow(sub_auc_data)
  
  p_auc = ggplot(sub_auc_data, aes(x = as.factor(x), y = auc, fill = as.factor(colour), group = as.factor(colour)))
  p_auc = p_auc + geom_bar(stat = "identity", position = "dodge", color = "black")
  p_auc = p_auc + geom_errorbar(aes(ymin = auc.l, ymax= auc.u), width = .3, position = position_dodge(.9), color = '#303030A0')
  p_auc = p_auc + theme_classic()
  p_auc = p_auc + ylab("AUC (95% CI)") + xlab("") + ggtitle(DATASET_NAME[SET]) + labs(fill = '')
  p_auc = p_auc + scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = paste0(seq(0, 100, by = 20), '%'))
  p_auc = p_auc + scale_x_discrete(breaks = sub_auc_data[, 'x'], labels = sub_auc_data[, 'show_name']) 
  p_auc = p_auc + annotate('text', sub_auc_data[, 'x'], 0.02, label = sub_auc_data[, 'txt'], size = 3.5, hjust = 0, colour = 'black', angle = 90, fontface = 1)
  p_auc = p_auc + scale_fill_manual(values =  c("#00A08792", '#7E614892', "#4DBBD592", "#E64B3592"))
  p_auc = p_auc + theme(plot.title = element_text(size = 15, face = 2),
                        legend.position = "none",
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(size = 14, face = 1),
                        axis.text.x = element_text(hjust = 1, angle = 45, size = 12, face = 1),
                        axis.text.y = element_text(angle = 0, size = 12, face = 1)) 
  
  thoracic_p_auc_list[[SET]] = p_auc
  
}

#### Merge plots ####
auc_p.1 = arrangeGrob(thoracic_p_auc_list[[1]], thoracic_p_auc_list[[2]], ncol = 2)

final_p = ggdraw()
final_p = final_p + draw_plot(XGB_list[[1]], x = -0.01, y = 0.5, width = 0.31, height = 0.48)
final_p = final_p + draw_plot(XGB_list[[2]], x = 0.29, y = 0.5, width = 0.36, height = 0.48)
final_p = final_p + draw_plot(XGB_list[[3]], x = 0.65, y = 0.5, width = 0.33, height = 0.48)
final_p = final_p + draw_plot(auc_p, x = 0.02, y = 0, width = 0.98, height = 0.48)
final_p = final_p + draw_plot_label(LETTERS[1:2], c(0.005, 0.005), c(1, 0.5), size = 20, hjust = 0)

pdf(plot_path_1, width = 18, height = 12)
print(final_p)
dev.off()