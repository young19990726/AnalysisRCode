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
plot_path_1 = './Figure_4_ROC_curve_CXR.pdf'
plot_path_2 = './Figure_4_ROC_curve_DXA.pdf'
plot_path_3 = './Figure_4_PRAUC_curve_CXR.pdf'
plot_path_4 = './Figure_4_PRAUC_curve_DXA.pdf'

#### Functions ####
PLOT_ROC_CURVE = function (x, y, best_cut, w = NULL, title = 'ROC curve', col = '#F8766D', show_PRAUC = FALSE, show_cut = FALSE) {
  
  if (is.null(w)) {
    
    X = x[!y %in% NA & !x %in% NA]
    Y = y[!y %in% NA & !x %in% NA]
    
  } else {
    
    X = x[!y %in% NA & !x %in% NA & !w %in% NA]
    Y = y[!y %in% NA & !x %in% NA & !w %in% NA]
    W = w[!y %in% NA & !x %in% NA & !w %in% NA]
    
  }
  
  if (!0 %in% Y | !1 %in% Y) {stop('y need to include 0 and 1')} else {
    
    if (!is.null(w)) {
      
      unique_sort_x = sort(x) %>% unique() # %>% sample(., size = 1000) %>% sort()
      
      roc_curve = list()
      roc_curve[['thresholds']] = c(-Inf, unique_sort_x[-1] - diff(unique_sort_x), Inf)
      roc_curve[['sensitivities']] = rep(NA, length(roc_curve[['thresholds']]))
      roc_curve[['specificities']] = rep(NA, length(roc_curve[['thresholds']]))
      
      pb = txtProgressBar(max = length(roc_curve[['thresholds']]), style = 3)
      
      for (q in 1:length(roc_curve[['thresholds']])) {
        
        group = factor(as.integer(factor(X >= roc_curve[['thresholds']][q], levels = c(FALSE, TRUE))) + Y * 10, c(1, 2, 11, 12))
        w_ori = tapply(W, group, sum, na.rm = TRUE)
        w_ori[w_ori %in% NA] = 0L
        
        w_tab = w_ori
        dim(w_tab) = c(2, 2)
        
        roc_curve[['sensitivities']][q] = w_tab[2,2] / sum(w_tab[,2])
        roc_curve[['specificities']][q] = w_tab[1,1] / sum(w_tab[,1])
        setTxtProgressBar(pb, q)
        
      }
      
      close(pb)
      
      roc_curve[['auc']] = sum(diff(roc_curve[['specificities']]) * (roc_curve[['sensitivities']][-1] + roc_curve[['sensitivities']][-length(roc_curve[['sensitivities']])]) / 2)
      
    } else {
      
      roc_curve = roc(response = Y, predictor = X)
      
    }
    
    if (show_PRAUC) {
      
      tab = table(factor(Y, levels = 0:1))
      
      roc_curve[['ppv']] = tab[2] * roc_curve[['sensitivities']] / (tab[1] * (1 - roc_curve[['specificities']]) + tab[2] * roc_curve[['sensitivities']] + 1e-15)
      roc_curve[['f1']] = 2 * roc_curve[['ppv']] * roc_curve[['sensitivities']] / (roc_curve[['ppv']] + roc_curve[['sensitivities']] + 1e-15)
      
      pb = txtProgressBar(max = length(roc_curve[['thresholds']]), style = 3)
      for (l in 2:length(roc_curve[['thresholds']])) {
        if (roc_curve[['ppv']][l] < max(roc_curve[['ppv']][-l:-length(roc_curve[['thresholds']])])) {roc_curve[['ppv']][l] = max(roc_curve[['ppv']][-l:-length(roc_curve[['thresholds']])])}
        setTxtProgressBar(pb, l)
      }
      close(pb)
      
    }
    
    if (!is.null(w)) {
      
      group = factor(as.integer(factor(X >= best_cut, levels = c(FALSE, TRUE))) + Y * 10, c(1, 2, 11, 12))
      w_ori = tapply(W, group, sum, na.rm = TRUE)
      w_ori[w_ori %in% NA] = 0L
      
      tab_test = w_ori
      dim(tab_test) = c(2, 2)
      
      sens_test = tab_test[2,2] / sum(tab_test[,2])
      spec_test = tab_test[1,1] / sum(tab_test[,1])
      ppv_test = tab_test[2,2] / sum(tab_test[2,])
      npv_test = tab_test[1,1] / sum(tab_test[1,])
      F1_test = (ppv_test^(-1) + sens_test^(-1) + 1e-15)^(-1) * 2
      
    } else {
      
      tab_test = table(factor(X > best_cut, levels = c(TRUE, FALSE)), factor(Y, levels = c(1, 0)))
      
      sens_test = tab_test[1,1] / sum(tab_test[,1])
      spec_test = tab_test[2,2] / sum(tab_test[,2])
      ppv_test = tab_test[1,1] / sum(tab_test[1,])
      npv_test = tab_test[2,2] / sum(tab_test[2,])
      F1_test = (ppv_test^(-1) + sens_test^(-1) + 1e-15)^(-1) * 2
      
    }
    
    # Show text
    
    if (show_PRAUC) {
      
      roc_data = data.frame(ppv = roc_curve[['ppv']], sens = roc_curve[['sensitivities']])
      roc_data = rbind(data.frame(ppv = 0, sens = 1), roc_data, data.frame(ppv = 1, sens = 0))
      roc_data = roc_data[order(roc_data[, 'ppv'], decreasing = TRUE),]
      roc_data = roc_data[order(roc_data[, 'sens']),]
      rownames(roc_data) = 1:nrow(roc_data)
      
      ppv_test_plot = max(roc_data[which(roc_data[, 'sens'] %in% sens_test), 'ppv'])
      
      if (show_cut) {
        
        roc_txt.1 = paste0('Cut-off', 
                           '\nPRAUC ', 
                           '\nF score', 
                           '\nSens. ', 
                           '\nSpec. ', 
                           '\nPPV ', 
                           '\nNPV ', 
                           '')
        
        roc_txt.2 = paste0('', formatC(abs(best_cut), 1 + (best_cut >= -1) * 2, format = 'f'),
                           '\n', formatC(diff(roc_data[, 'sens']) %*% roc_data[, 'ppv'][-nrow(roc_data)], 4, format = 'f'),
                           '\n', formatC(F1_test, 4, format = 'f'),
                           '\n', formatC(sens_test * 100, 1, format = 'f'),
                           '%\n', formatC(spec_test * 100, 1, format = 'f'),
                           '%\n', formatC(ppv_test * 100, 1, format = 'f'),
                           '%\n', formatC(npv_test * 100, 1, format = 'f'),
                           '%')
        
      } else {
        
        roc_txt.1 = paste0('PRAUC ', 
                           '\nF score', 
                           '\nSens. ', 
                           '\nSpec. ', 
                           '\nPPV ', 
                           '\nNPV ', 
                           '')
        
        roc_txt.2 = paste0(formatC(diff(roc_data[, 'sens']) %*% roc_data[, 'ppv'][-nrow(roc_data)], 4, format = 'f'),
                           '\n', formatC(F1_test, 4, format = 'f'),
                           '\n', formatC(sens_test * 100, 1, format = 'f'),
                           '%\n', formatC(spec_test * 100, 1, format = 'f'),
                           '%\n', formatC(ppv_test * 100, 1, format = 'f'),
                           '%\n', formatC(npv_test * 100, 1, format = 'f'),
                           '%')
        
      }
      
    } else {
      
      roc_data = data.frame(spec = roc_curve[['specificities']], sens = roc_curve[['sensitivities']])
      roc_data = rbind(data.frame(spec = 0, sens = 1), roc_data, data.frame(spec = 1, sens = 0))
      roc_data = roc_data[order(roc_data[, 'spec'], decreasing = TRUE),]
      roc_data = roc_data[order(roc_data[, 'sens']),]
      rownames(roc_data) = 1:nrow(roc_data)
      
      if (show_cut) {
        
        roc_txt.1 = paste0('Cut-off', 
                           '\nAUC ', 
                           '\nSens. ', 
                           '\nSpec. ', 
                           '\nPPV ', 
                           '\nNPV ', 
                           '')
        
        roc_txt.2 = paste0('', formatC(abs(best_cut), 1 + (best_cut >= -1) * 2, format = 'f'),
                           '\n', formatC(roc_curve[['auc']], 4, format = 'f'),
                           '\n', formatC(sens_test * 100, 1, format = 'f'),
                           '%\n', formatC(spec_test * 100, 1, format = 'f'),
                           '%\n', formatC(ppv_test * 100, 1, format = 'f'),
                           '%\n', formatC(npv_test * 100, 1, format = 'f'),
                           '%')
        
      } else {
        
        roc_txt.1 = paste0('AUC ', 
                           '\nSens. ', 
                           '\nSpec. ', 
                           '\nPPV ', 
                           '\nNPV ', 
                           '')
        
        roc_txt.2 = paste0(formatC(roc_curve[['auc']], 3, format = 'f'),
                           '\n', formatC(sens_test * 100, 1, format = 'f'),
                           '%\n', formatC(spec_test * 100, 1, format = 'f'),
                           '%\n', formatC(ppv_test * 100, 1, format = 'f'),
                           '%\n', formatC(npv_test * 100, 1, format = 'f'),
                           '%')
        
      }
      
    }
    
    # ROC curve
    
    if (show_PRAUC) {
      
      roc_p = ggplot(data = roc_data, aes(x = sens, y = ppv))
      roc_p = roc_p + geom_line(colour = col, size = 1.5)
      roc_p = roc_p + theme_bw()
      roc_p = roc_p + coord_equal()
      roc_p = roc_p + ggtitle(title)
      roc_p = roc_p + ggtitle(title) + xlab('Sensitivity') + ylab('Positive predictive value')
      roc_p = roc_p + annotate(geom = "point", x = sens_test, y = ppv_test_plot, shape = 21, size = 5, fill = paste0(col, 'A0'), color = '#000000')
      roc_p = roc_p + annotate(geom = "text", x = 0.1, y = 0.05, label = roc_txt.1, size = 5, fontface = 2, colour = '#00000080', hjust = 0, vjust = 0)
      roc_p = roc_p + annotate(geom = "text", x = 0.9, y = 0.05, label = roc_txt.2, size = 5, fontface = 2, colour = '#00000080', hjust = 1, vjust = 0)
      
    } else {
      
      roc_p = ggplot(data = roc_data, aes(x = spec, y = sens))
      roc_p = roc_p + geom_line(colour = col, size = 1)
      roc_p = roc_p + theme_bw()
      roc_p = roc_p + coord_equal()
      roc_p = roc_p + ggtitle(title)
      roc_p = roc_p + ggtitle(title) + xlab('Specificity') + ylab('Sensitivity')
      roc_p = roc_p + annotate(geom = "point", x = spec_test, y = sens_test, shape = 21, size = 5, fill = paste0(col, 'A0'), color = '#000000')
      roc_p = roc_p + annotate(geom = "text", x = 0.05, y = 0.00, label = roc_txt.1, size = 4, fontface = 2, colour = col, hjust = 0, vjust = 0)
      roc_p = roc_p + annotate(geom = "text", x = 0.60, y = 0.00, label = roc_txt.2, size = 4, fontface = 2, colour = col, hjust = 1, vjust = 0)
      
    }
    
    roc_p = roc_p + theme(plot.title = element_text(color = "#000000", size = 14),
                          axis.title = element_text(color = "#000000", size = 12),
                          legend.position = "none")
    
    return(roc_p)
    
  }
  
}

#### Settings ####
set.seed(0)

DATASET_NAME = c('internal-test' = 'Internal validation set', 'community-test' = 'External validation set')

col_list = c('internal-test' = '#F8766D', 'community-test' = '#619CFF')

y_name_list = c('Y1', 'Y2', 'Y3')
x_name_list = c('predition [[BASIC] Future compression fracture risk]', 'T-score')
main_name_list = c('1 year', '3 year', '5 year')

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
patient_data[patient_data[, 'event[Fracture_outcome]'] %in% 1 & patient_data[, 'time[Fracture_outcome]'] <= 365.25 * 1, 'Y1'] = 1L
patient_data[patient_data[, 'event[Fracture_outcome]'] %in% 0 & patient_data[, 'time[Fracture_outcome]'] < 365.25 * 1, 'Y1'] = NA

patient_data[, 'Y2'] = 0L
patient_data[patient_data[, 'event[Fracture_outcome]'] %in% 1 & patient_data[, 'time[Fracture_outcome]'] <= 365.25 * 3, 'Y2'] = 1L
patient_data[patient_data[, 'event[Fracture_outcome]'] %in% 0 & patient_data[, 'time[Fracture_outcome]'] < 365.25 * 3, 'Y2'] = NA

patient_data[, 'Y3'] = 0L
patient_data[patient_data[, 'event[Fracture_outcome]'] %in% 1 & patient_data[, 'time[Fracture_outcome]'] <= 365.25 * 5, 'Y3'] = 1L
patient_data[patient_data[, 'event[Fracture_outcome]'] %in% 0 & patient_data[, 'time[Fracture_outcome]'] < 365.25 * 5, 'Y3'] = NA

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

patient_data[, 'T-score'] = 1 - rescale(patient_data[, '[DXA] T-score'], to = c(0, 1))

#### Plotting ####
roc_list = list()
cut_list = rep(NA, length(x_name_list))

for (SET in names(DATASET_NAME)) {
  
  sub_valid_data = patient_data[patient_data[, 'DATASET'] %in% 'valid',]
  sub_test_data = patient_data[patient_data[, 'DATASET'] %in% SET,]
  
  roc_list[[SET]] = list()
  
  for (i in 1:length(y_name_list)) {
    
    y_name = y_name_list[i]
    main_txt = DATASET_NAME[SET]
    
    roc_list[[SET]][[y_name]] = list()
    
    for (j in 1:length(x_name_list)) {
      
      x_name = x_name_list[j]
      
      if (cut_list[j] %in% NA) {
        
        roc_curve.valid = roc(response = sub_valid_data[, y_name], predictor = sub_valid_data[, x_name])
        best_pos = which.max(roc_curve.valid[['sensitivities']] + roc_curve.valid[['specificities']])
        cut_list[j] = roc_curve.valid[['thresholds']][best_pos]
        
      }
      
      roc_list[[SET]][[y_name]][[x_name]] = PLOT_ROC_CURVE(x = sub_test_data[, x_name],
                                                           y = sub_test_data[, y_name],
                                                           best_cut = cut_list[j],
                                                           title = main_name_list[i],
                                                           col = col_list[SET],
                                                           show_PRAUC = FALSE,
                                                           show_cut = TRUE)
    }
  }
}

prroc_list = list()
cut_list = rep(NA, length(x_name_list))

for (SET in names(DATASET_NAME)) {
  
  sub_valid_data = patient_data[patient_data[, 'DATASET'] %in% 'valid',]
  sub_test_data = patient_data[patient_data[, 'DATASET'] %in% SET,]
  
  prroc_list[[SET]] = list()
  
  for (i in 1:length(y_name_list)) {
    
    y_name = y_name_list[i]
    main_txt = DATASET_NAME[SET]
    
    prroc_list[[SET]][[y_name]] = list()
    
    for (j in 1:length(x_name_list)) {
      
      x_name = x_name_list[j]
      
      # cut_list 要從 sub_model_info <- model_info[model_info[, 'label'] %in% var_of_interested,] 找 best_cut.1
      if (cut_list[j] %in% NA) {
        
        roc_curve.valid = roc(response = sub_valid_data[, y_name], predictor = sub_valid_data[, x_name])
        
        tab = table(factor(sub_valid_data[, y_name], levels = 0:1))
        
        roc_curve.valid[['ppv']] = tab[2] * roc_curve.valid[['sensitivities']] / (tab[1] * (1 - roc_curve.valid[['specificities']]) + tab[2] * roc_curve.valid[['sensitivities']] + 1e-15)
        roc_curve.valid[['f1']] = 2 * roc_curve.valid[['ppv']] * roc_curve.valid[['sensitivities']] / (roc_curve.valid[['ppv']] + roc_curve.valid[['sensitivities']] + 1e-15)
        
        pb = txtProgressBar(max = length(roc_curve.valid[['thresholds']]), style = 3)
        for (l in 2:length(roc_curve.valid[['thresholds']])) {
          if (roc_curve.valid[['ppv']][l] < max(roc_curve.valid[['ppv']][-l:-length(roc_curve.valid[['thresholds']])])) {roc_curve.valid[['ppv']][l] = max(roc_curve.valid[['ppv']][-l:-length(roc_curve.valid[['thresholds']])])}
          setTxtProgressBar(pb, l)
        }
        close(pb)
        
        best_pos = which.max(roc_curve.valid[['f1']])
        cut_list[j] = roc_curve.valid[['thresholds']][best_pos]
        
      }
      
      prroc_list[[SET]][[y_name]][[x_name]] = PLOT_ROC_CURVE(x = sub_test_data[, x_name],
                                                             y = sub_test_data[, y_name],
                                                             best_cut = cut_list[j],
                                                             title = main_txt,
                                                             col = col_list[SET],
                                                             show_PRAUC = TRUE,
                                                             show_cut = TRUE)
      
      
    }
  }
}

#### Merge plots ####
pp.1 = arrangeGrob(roc_list[[1]][[1]][[1]], roc_list[[1]][[2]][[1]], roc_list[[1]][[3]][[1]], 
                   roc_list[[2]][[1]][[1]], roc_list[[2]][[2]][[1]], roc_list[[2]][[3]][[1]], ncol = 3)

pp.2 = arrangeGrob(roc_list[[1]][[1]][[2]], roc_list[[1]][[2]][[2]], roc_list[[1]][[3]][[2]], 
                   roc_list[[2]][[1]][[2]], roc_list[[2]][[2]][[2]], roc_list[[2]][[3]][[2]], ncol = 3)

final_p = ggdraw()
final_p = final_p + draw_plot(pp.1, x = 0, y = 0, width = 1, height = 0.97)
final_p = final_p + draw_plot_label(c("Internal validation set", "External validation set"), c(0.005, 0.005), c(1, 0.505), size = 18, hjust = 0, fontface = 2)

# pdf(plot_path_1, width = 10, height = 8)
# print(final_p)
# dev.off()

final_p = ggdraw()
final_p = final_p + draw_plot(pp.2, x = 0, y = 0, width = 1, height = 0.97)
final_p = final_p + draw_plot_label(c("Internal validation set", "External validation set"), c(0.005, 0.005), c(1, 0.505), size = 18, hjust = 0, fontface = 2)

# pdf(plot_path_2, width = 10, height = 8)
# print(final_p)
# dev.off()

pp.1 = arrangeGrob(prroc_list[[1]][[1]][[1]], prroc_list[[1]][[2]][[1]], prroc_list[[1]][[3]][[1]], 
                   prroc_list[[2]][[1]][[1]], prroc_list[[2]][[2]][[1]], prroc_list[[2]][[3]][[1]], ncol = 3)

pp.2 = arrangeGrob(prroc_list[[1]][[1]][[2]], prroc_list[[1]][[2]][[2]], prroc_list[[1]][[3]][[2]], 
                   prroc_list[[2]][[1]][[2]], prroc_list[[2]][[2]][[2]], prroc_list[[2]][[3]][[2]], ncol = 3)

final_p = ggdraw()
final_p = final_p + draw_plot(pp.1, x = 0, y = 0, width = 1, height = 0.97)
final_p = final_p + draw_plot_label(c("Internal validation set", "External validation set"), c(0.005, 0.005), c(1, 0.505), size = 18, hjust = 0, fontface = 2)

# pdf(plot_path_3, width = 10, height = 8)
# print(final_p)
# dev.off()

final_p = ggdraw()
final_p = final_p + draw_plot(pp.2, x = 0, y = 0, width = 1, height = 0.97)
final_p = final_p + draw_plot_label(c("Internal validation set", "External validation set"), c(0.005, 0.005), c(1, 0.505), size = 18, hjust = 0, fontface = 2)

# pdf(plot_path_4, width = 10, height = 8)
# print(final_p)
#  dev.off()