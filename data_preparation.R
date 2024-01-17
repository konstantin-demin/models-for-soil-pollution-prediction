# load libraries, загрузить библиотеки
library("tidyverse")
library("randomForest")

# set wd, установить рабочую директорию
setwd("~/")

# set dir with kraken2 outputs, задать директорию с профилями, созданными kraken2
mydir = "path/to/profiles"

# load predicting_labels dataset, загрузить таблицу predicting_labels 
labels <- read.table("path/to/predicting_labels.csv", header = T, sep = "\t")

# list taxonomic profiles in "myfiles"
myfiles = list.files(path=mydir, pattern="*", full.names=TRUE)

# create empty tibble
z <- tibble()

# load taxonomic profiles and put them into the tibble
for (i in 1:length(myfiles)) {
  z<-rbind(z, (read.csv(myfiles[i], header = F) %>%
                 as_tibble() %>% 
                 mutate(sample = myfiles[i])))
}

# clean sample names, убрать пути из имен файлов
z1 <- z %>% mutate(sample = str_replace_all(sample, "path/to/profiles/", ""))

# check the sample list, проверить список загруженных образцов
z1$sample %>% unique()

# introduce kingdom to prokaryotic rows, создать новую переменную kingdom в датасет
z1 <- z1 %>%
  mutate(V1 = str_replace_all(V1, "d__Bacteria", "d__Bacteria|k__Bacteria"),
         V1 = str_replace_all(V1, "d__Archaea", "d__Archaea|k__Archaea"))

# split by taxonomic level, разделить датасет по таксономическим категориям
z2 <- z1 %>% 
  mutate(V1 = str_replace_all(V1, "\t", "#")) %>% 
  mutate(V1 = str_replace_all(V1, "\\|", "#")) %>% 
  mutate(delimcount = pmap_int(., ~ str_count(c(...), "#") %>% 
                                 sum)) %>% 
  filter(grepl("k__Bacteria", V1) |
           grepl("k__Archaea", V1) |
           grepl("k__Fungi", V1)) %>% 
  mutate(V1 = case_when(delimcount == 1 ~ gsub("#([0-9]*)$", "########\\1", V1),
                        delimcount == 2 ~ gsub("#([0-9]*)$", "#######\\1", V1),
                        delimcount == 3 ~ gsub("#([0-9]*)$", "######\\1", V1),
                        delimcount == 4 ~ gsub("#([0-9]*)$", "#####\\1", V1),
                        delimcount == 5 ~ gsub("#([0-9]*)$", "####\\1", V1),
                        delimcount == 6 ~ gsub("#([0-9]*)$", "###\\1", V1),
                        delimcount == 7 ~ gsub("#([0-9]*)$", "##\\1", V1),
                        delimcount == 8 ~ V1)) %>% 
  separate_wider_delim(cols = V1,
                       names = c("d", "k", "p", "c", "o", "f", "g", "s", "n"),
                       delim = "#", too_few = "align_start") %>% 
  mutate(n = as.numeric(n))

# вывести сумму классифицированных прокариот и грибов в отдельный объект
allmicrobes <- z2 %>%
  group_by(sample) %>%
  filter(grepl("k__", k) & !grepl(".", p)) %>% 
  reframe(allmicrobes = sum(n))

# посчитать число родов в датасете
g_in_g <-z2 %>% 
  filter(grepl("g__", g) & !grepl(".", s)) %>%
  select(sample, g, n)

g_in_f <-z2 %>% 
  filter(grepl("g__", f) & !grepl(".", g)) %>%
  mutate(g = f) %>% select(sample, g, n)

g_in_o <-z2 %>% 
  filter(grepl("g__", o) & !grepl(".", f)) %>%
  mutate(g = o) %>% select(sample, g, n)

g_in_c <-z2 %>% 
  filter(grepl("g__", c) & !grepl(".", o)) %>%
  mutate(g = c) %>% select(sample, g, n)

all_g <- rbind(g_in_g, g_in_f, g_in_o, g_in_c)

all_g_share <- all_g %>% 
  inner_join(allmicrobes, by = "sample") %>% 
  mutate(share = (n/allmicrobes)*100) %>% 
  select(sample, g, share)%>% 
  mutate(g = str_replace_all(g, "g__Candidatus ", ""),
         g = word(g, 1),
         g = str_replace_all(g, "g__", ""),
         tax = g) %>% 
  select(sample, tax, share)

# посчитать число семейств в датасете 
f_in_f <-z2 %>% 
  filter(grepl("f__", f) & !grepl(".", g)) %>%
  select(sample, f, n)

f_in_o <-z2 %>% 
  filter(grepl("f__", o) & !grepl(".", f)) %>%
  mutate(f = o) %>% select(sample, f, n)

f_in_c <-z2 %>% 
  filter(grepl("f__", c) & !grepl(".", o)) %>%
  mutate(f = c) %>% select(sample, f, n)

f_in_p <-z2 %>% 
  filter(grepl("f__", p) & !grepl(".", c)) %>%
  mutate(f = p) %>% select(sample, f, n)

all_f <- rbind(f_in_p, f_in_c, f_in_o, f_in_f)

all_f_share <- all_f %>% 
  inner_join(allmicrobes, by = "sample") %>% 
  mutate(share = (n/allmicrobes)*100) %>% 
  select(sample, f, share) %>% 
  mutate(f = str_replace_all(f, "f__Candidatus ", ""),
         f = word(f, 1),
         f = str_replace_all(f, "f__", ""),
         tax = f) %>% 
  select(sample, tax, share)

# объединить число родов и семейств в единый объект
f_and_g <- rbind(all_f_share, all_g_share)

# join genus-family dataset and predicting labels dataset saving mutual taxa,
# отфильтровать итоговый датасет, сохраняя таксоны из файла predicting_labels
data_cleaned <- f_and_g %>%
  inner_join(tax_labels, by = "tax") %>%
  select(sample, label, share) %>% 
  group_by(sample, label) %>% 
  reframe(share = sum(share)) %>% 
  pivot_wider(names_from = label, values_from = share, values_fill = 0) %>%
  as.data.frame() %>% 
  column_to_rownames(var = "sample")

# save samples taxonomic profiles, сохранить таксономические профили образцов
write.table(data_cleaned, "path/to/save/data_cleaned.csv", sep = "\t", row.names = F)