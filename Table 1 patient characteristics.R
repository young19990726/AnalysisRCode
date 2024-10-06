#### Library ####
library(magrittr)

#### Pathways ####
setwd('/home/young19990726/RCODE/AI-CXR for aortic aneurysm')
data_path = '/home/chinlin/Multi-center-database/002. data/002. project data/048. project-048/analysis/study 014/patient_data.RData'
table_path = './Table 01.doc'

#### Function ####
source('/home/young19990726/RCODE/000 Table function.R')

#### Settings ####
set.seed(0)

DATASET_NAME = c('internal-test' = 'Internal validation set', 'community-test' = 'External validation set')

group_name.1 = 'DATASET'
group_name.2 = 'Thoracic aortic aneurysm'
group_name.3 = 'Abdominal aortic aneurysm'

intrested_var.1 = c(
  'AGE', 'GENDER', 'CXR[MODALITY]', 'CXR[POSITION]', 'CXR[VIEW]', 
  'DM', 'HTN', 'CKD', 'HLP', 'HF', 'CAD', 'COPD'
)

intrested_var.2 = c("[RADIO] Consolidation change", "[RADIO] Pneumonia", "[RADIO] Emphysematous change", "[RADIO] Pneumothorax", "[RADIO] Atelectasis", "[RADIO] Scalloping of the diaphragm", "[RADIO] Costophrenic angle blunting", "[RADIO] Pleural effusion", "[RADIO] Atherosclerosis", "[RADIO] Cardiomegaly", "[RADIO] Prominence of hilar shadow", "[RADIO] Pulmonary edema", "[RADIO] Aneurysm", "[RADIO] Degenerative joint disease", "[RADIO] Fracture", "[RADIO] Spondylosis", "[RADIO] Osteophyte formation", "[RADIO] Osteoporosis", "[RADIO] Osteoarthritis", "[RADIO] Widening of the mediastinum", "[RADIO] Malignancy", "[RADIO] Inflammatory", "[RADIO] Pigtail or drainage", "[RADIO] Sternotomy", "[RADIO] Port a implantation", "[RADIO] Perm catheter insertion", "[RADIO] Pacemaker", "[RADIO] Tracheostomy", "[RADIO] Vertebroplasty", "[RADIO] Endotracheal tube", "[RADIO] Nasogastric tube")

intrested_var = c(intrested_var.1, intrested_var.2)

#### Processes ####
load(data_path)

## Randomly shuffle data and remove duplicates
patient_data = patient_data[sample(1:nrow(patient_data)),]
patient_data = patient_data[!duplicated(patient_data[, 'CNO']),]

patient_data[patient_data[, 'CXR[POSITION]'] %in% 'Unknown', 'CXR[POSITION]'] = 'OPD'
patient_data[patient_data[, 'CXR[POSITION]'] %in% c('OPD', 'PEC'), 'CXR[POSITION]'] = 'OPD'
patient_data[patient_data[, 'CXR[POSITION]'] %in% 'ER', 'CXR[POSITION]'] = 'ED'

## Preprocessing
patient_data[, group_name.1] = factor(patient_data[, group_name.1], levels = c('train', 'valid', names(DATASET_NAME)))
patient_data[, group_name.2] = factor(patient_data[, 'group [[CCT] Thoracic aortic aneurysm]'], levels = c('Low risk', 'Median risk', 'High risk'), labels = c('Low risk', 'Median risk', 'High risk'))
patient_data[, group_name.3] = factor(patient_data[, 'group [[ACT] Abdominal aortic aneurysm]'], levels = c('Low risk', 'Median risk', 'High risk'), labels = c('Low risk', 'Median risk', 'High risk'))

for (m in 1:length(intrested_var)) {
  
  len_lvl = factor(patient_data[, intrested_var[m]]) %>% levels %>% length()
  
  if (len_lvl > 5) {
    
    patient_data[, intrested_var[m]] = as.numeric(patient_data[, intrested_var[m]])
    
  } else {patient_data[, intrested_var[m]] = factor(patient_data[, intrested_var[m]])}
}

## Build table
table_list = list()

table_list[[1]] = Table1(X = patient_data[, group_name.1], Y.matrix = patient_data[, c(group_name.2, group_name.3, intrested_var)], x.name = "Dataset")

sub_data = patient_data[patient_data[,group_name.1] %in% names(DATASET_NAME)[1],]
table_list[[2]] = Table1(X = sub_data[, group_name.2], Y.matrix = sub_data[, c(group_name.2, group_name.3, intrested_var)], x.name = "Thoracic aortic aneurysm")
table_list[[3]] = Table1(X = sub_data[, group_name.3], Y.matrix = sub_data[, c(group_name.2, group_name.3, intrested_var)], x.name = "Abdominal aortic aneurysm")

sub_data = patient_data[patient_data[,group_name.1] %in% names(DATASET_NAME)[2],]
table_list[[4]] = Table1(X = sub_data[, group_name.2], Y.matrix = sub_data[, c(group_name.2, group_name.3, intrested_var)], x.name = "Thoracic aortic aneurysm")
table_list[[5]] = Table1(X = sub_data[, group_name.3], Y.matrix = sub_data[, c(group_name.2, group_name.3, intrested_var)], x.name = "Abdominal aortic aneurysm")

for (i in 1:length(table_list)) {
  total_table = table_list[[i]]
  for (j in 1:nrow(total_table)) {
    if (total_table[j, 'p-value'] == '') {
      total_table[j, 'p-value'] = total_table[j-1, 'p-value']
    }
  }
  save_pos = total_table[nrow(total_table):1, 'Variable'] %>% gsub('\\:.*', '', .) %>% duplicated()
  save_pos[length(save_pos):(length(save_pos)-10)] = FALSE
  total_table = total_table[!save_pos[length(save_pos):1],]
  table_list[[i]] = total_table
}

Table2doc(table_list = table_list, filename = table_path)