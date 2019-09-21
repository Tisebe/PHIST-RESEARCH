# PHIST-RESEARCH
R scripts for my analysis on Antibody responses
Stepwise R scripts for all my analysis measuring antibody responses for expressed and purified recombinant PHISTb antigens.
# Measuring antibody responses for children aged 3 to 5 years in two different geographical locations (Siaya and Takaungu)
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)

mydata1 <- read.csv("C:/Users/Tony Isebe/Desktop/Backup/PHIST RESEARCH/csv raw/Response_age.csv")
view(mydata1)
mydata1 %>%
  filter(AGE.YEARS. >=3,AGE.YEARS. <=5)%>%
  ggplot(aes(location, SCHIZONT.EXTRACT.OD.VALUES))+
  geom_boxplot()+
  geom_quasirandom()+
 stat_compare_means( method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme_bw()
  
# reading in phist antigens and ELISA data
sukuta = read.csv("C:/Users/Tony Isebe/Desktop/Backup/PHIST RESEARCH/csv raw/SUKUTA_D.csv", stringsAsFactors = F, header = T)
elisa_data = read.csv("C:/Users/Tony Isebe/Desktop/Backup/PHIST RESEARCH/csv raw/Copy of ELISA_withmetadata.csv", stringsAsFactors = F, header = T)

# merging data; phist antigens with elisa data
elisa_data$PF3D7_0532400 = as.character(elisa_data$PF3D7_0532400)
sukuta$PF3D7_0532400 = as.character(suguta$PF3D7_0532400)
phist_ags_elisa = merge(suguta, elisa_data, by.x = "PF3D7_0532400", by.y = "PF3D7_0532400")

# measuring responses across different age groups
mydata <- read.csv("C:\\Users/John/Desktop/tony/Response_age.csv")

mydata$age_categ <- cut(mydata$AGE.YEARS., c(0,3,6,10))

mydata %>% 
            ggplot(aes(age_categ, PF3D7_0532400))+
            geom_boxplot()
            
# adding statistical comparisons
library(tidyverse)
library(ggpubr)

mydata <- read.csv("C:\\Users/John/Desktop/tony/Response_age.csv")

mydata$age_categ <- cut(mydata$AGE.YEARS., c(0,3,6,10))

comparisons <- list(c('(0,3]', '(3,6]'), c('(3,6]', '(6,10]'), c('(0,3]', '(6,10]'))

mydata %>% 
            ggplot(aes(age_categ, PF3D7_0532400))+
            geom_boxplot()+
            stat_compare_means(comparisons = comparisons)
# comparing responses in locations

library(tidyverse)
library(ggpubr)
library(ggbeeswarm)

mydata <- read.csv("C:\\Users/John/Desktop/tony/Response_age.csv")

mydata$age_categ <- cut(mydata$AGE.YEARS., c(0,3,6,10))

comparisons <- list(c('(0,3]', '(3,6]'), c('(3,6]', '(6,10]'), c('(0,3]', '(6,10]'))

mydata %>% 
            ggplot(aes(age_categ, PF3D7_0532400))+
            geom_boxplot()+
            geom_quasirandom()+
            stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
            stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
            facet_wrap(~location, scales = 'free')+
            theme_bw()
# comparing responses based on transmission intensity
my_comparisons <- list(c('JUNJU','SIAYA'),c('JUNJU','SUKUTA'),c('JUNJU','TAKAUNGU'),c('SIAYA','SUKUTA'),c('SIAYA','TAKAUNGU'),c('SUKUTA','TAKAUNGU'))
nytoy %>% 
  ggplot(aes(location, SCHIZONT.EXTRACT.OD.VALUES))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme_light()
  # generating beeswarm plots
  boxplot(PF3D7_0532400 ~ PARA, data = nytoy, 
          outline = FALSE,     ## avoid double-plotting outliers, if any
          main = 'Comparing Antibody Levels in Children with  Parasitemia
')
beeswarm(PF3D7_0532400 ~ location, data = nytoy, 
         col = 4, pch = 16, add = TRUE)+
  stat_compare_means( method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme_bw()
  # Add geom lines to box plots
  tony<-tony %>% gather(Variables,Values,c("SCHIZONT.EXTRACT.OD.VALUES","PF3D7_0532400","PF3D7_1401600","PF3D7_1102500")) 

tony$Variables<-factor(tony$Variables,levels = unique(tony$Variables))
tony %>% 
  ggplot(aes(as.factor(AGECODE),Values))+geom_boxplot()+facet_wrap(~Variables)+
  xlab('AGE(YEARS)') + ylab('OD@492nm')+
  stat_summary(fun.y = mean, geom = 'line', aes(group=1))+
  stat_summary(fun.y = mean, geom = 'point')
  
