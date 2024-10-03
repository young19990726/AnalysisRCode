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
plot_path_1 = './Scatter.pdf'

#### Functions ####
get_density = function(x, y, ...) {
  dens = MASS::kde2d(x, y, ...)
  ix = findInterval(x, dens$x)
  iy = findInterval(y, dens$y)
  ii = cbind(ix, iy)
  return(dens$z[ii])
}

#### Settings ####
set.seed(0)

DATASET_NAME = c('internal-test' = 'Internal validation set', 'community-test' = 'External validation set')

y_name_list = c('Y1')
x_name_list = c('predition [[BASIC] Future compression fracture risk]')

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
scattor_p_list = list()

for (SET in names(DATASET_NAME)) {
  
  sub_test_data = patient_data[patient_data[, 'DATASET'] %in% SET,]
  
  x1 = sub_test_data[, 'predition [[BASIC] Future compression fracture risk]']
  x2 = sub_test_data[, '[DXA] T-score']
  
  current_cor = cor(x1, x2)
  current_MAE = abs(x1 - x2) %>% mean()
  current_mean = (x1 - x2) %>% mean()
  current_sd = (x1 - x2) %>% sd()
  
  current_txt = paste0('Diff = ', formatC(current_mean, 2, format = 'f'), '(+/-)', formatC(current_sd, 2, format = 'f'), '\nr = ', formatC(current_cor, 2, format = 'f'), '\nMAE = ', formatC(current_MAE, 2, format = 'f'))

  myColor = rev(RColorBrewer::brewer.pal(11, "Spectral"))
  
  scattor_dat = data.frame(x = x1, y = x2)
  scattor_dat = scattor_dat[apply(is.na(scattor_dat), 1, sum) == 0,]
  scattor_dat$density = get_density(scattor_dat$x, scattor_dat$y, n = 100)
  scattor_dat$density = scattor_dat$density / max(scattor_dat$density)
  
  scattor_p = ggplot(data = scattor_dat, aes(x, y))
  scattor_p = scattor_p + geom_point(aes(x, y, col = density), size = 1)
  scattor_p = scattor_p + scale_colour_gradientn(colours = myColor)
  scattor_p = scattor_p + theme_bw()
  scattor_p = scattor_p + theme(legend.position = "none")
  scattor_p = scattor_p + xlim(c(-5, 5)) + ylim(c(-5, 5))
  scattor_p = scattor_p + xlab("CXR-T score") + ylab('DXA-T score') + ggtitle(DATASET_NAME[SET]) 
  scattor_p = scattor_p + coord_equal()
  scattor_p = scattor_p + geom_smooth(method = lm, se = TRUE, linetype = "dashed", color = "#303030A0",  fullrange = TRUE)
  scattor_p = scattor_p + annotate(geom = "text",
                                   x = -5 + 10 * 0.01,
                                   y = -5 + 10 * 0.9,
                                   label = current_txt, size = 4, fontface = 2, colour = '#000000B0', hjust = 0)
  
  scattor_p = scattor_p + theme(plot.title = element_text(color = "#000000", size = 14, face = 1, hjust = 0),
                                legend.position = "none",
                                axis.title.x = element_text(color = "#000000", size = 14),
                                axis.title.y = element_text(color = "#000000", size = 14),
                                axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, size = 8, face = "bold"),
                                axis.text.y = element_text(color = "#000000", angle = 45, hjust = 1, size = 8, face = "bold"))
  
  scattor_p_list[[SET]] = scattor_p
}

#### Merge plots ####
final_p = ggdraw()
final_p = final_p + draw_plot(scattor_p_list[[1]], x = 0, y = 0, width = 0.5, height = 1)
final_p = final_p + draw_plot(scattor_p_list[[2]], x = 0.5, y = 0, width = 0.5, height = 1)

# pdf(plot_path_1, width = 8, height = 4)
print(final_p)
# dev.off()