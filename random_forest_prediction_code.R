# load libraries, загрузить библиотеки
library("tidyverse")
library("randomForest")

#load samples taxonomic profiles, загрузить таблицу с таксономическими профилями образцов
read.table("path/to/table/data_cleaned.csv", sep = "\t", header = T)

# load random forest model, загрузить модель случайного леса 
rf <- readRDS("path to rds")

# make random forest predictions, предсказать уровень загрязнения моделью
predicted_pollution <- predict(my_rf, data_cleaned) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% as_tibble()  %>%
  mutate(across(2, ~ ., .names = "predicted pollution level")) %>%
  select(`predicted pollution level`, sample) 

# save predicted values to a new table, сохранить предсказанное в новую таблицу 
write.table(predicted_pollution,
            "path/to/predicted_pollution.csv", sep = '\t', row.names = F)