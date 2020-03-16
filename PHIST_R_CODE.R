##FOR ANYTHING DATA RELATED

# https://www.r-graph-gallery.com/269-ggplot2-boxplot-with-average-value.html

library(tidyverse) 
library(ggpubr)
library(ggbeeswarm)
library(ggplot2)

##DATA
Phist_All_sites <- read.csv("C:/Users/Tony Isebe/Desktop/Backup/PHIST RESEARCH/csv raw/PHISTb_ELISA_Metadata_All_Locations.csv")
View(Phist_All_sites)


PHIST_Data_without_Sukuta <- read.csv("C:/Users/Tony Isebe/Desktop/Backup/PHIST RESEARCH/csv raw/PHIST_Data_without_Sukuta.csv")
view(PHIST_Data_without_Sukuta)

xy <- read.csv("C:/Users/Tony Isebe/Desktop/Backup/PHIST RESEARCH/csv raw/xy.csv")
view(xy)

xy %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")

ggplot(data = xy)+
  geom_boxplot(mapping = aes())

##PLOTS
#Antibody responses vs. transmission intensity

# comparing responses based on transmission intensity
#for schizont extract
my_comparisons <- list(c('SIAYA','TAKAUNGU'),c('SIAYA','SUKUTA'),c('TAKAUNGU', 'SUKUTA'))
Phist_All_sites %>% 
  ggplot(aes(location, SCHIZONT.EXTRACT.OD.VALUES))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=median, geom="point", shape=20, size=5, color="red", fill="red") +
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "grey")
  )+
  xlab('Location') + ylab('Crude_schizont_extract (OD 492nm')

#for PF3D7_0532400
my_comparisons <- list(c('SIAYA','TAKAUNGU'),c('SIAYA','SUKUTA'),c('TAKAUNGU', 'SUKUTA'))
Phist_All_sites %>% 
  ggplot(aes(location, PF3D7_0532400))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "grey")
  )+
  xlab('Location') + ylab('Pf3D7_0532400 (OD 492nm)')


#for PF3D7_1102500
my_comparisons <- list(c('SIAYA','TAKAUNGU'),c('SIAYA','SUKUTA'),c('TAKAUNGU', 'SUKUTA'))
Phist_All_sites %>% 
  ggplot(aes(location, PF3D7_1102500))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "grey")
  )+
  xlab('Location') + ylab('Pf3D7_1102500 (OD 492nm)')


#for PF3D7_1401600
my_comparisons <- list(c('SIAYA','TAKAUNGU'),c('SIAYA','SUKUTA'),c('TAKAUNGU', 'SUKUTA'))
Phist_All_sites %>% 
  ggplot(aes(location, PF3D7_1401600))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "grey")
  )+
  xlab('Location') + ylab('Pf3D7_1401600 (OD 492nm)')


#Antibody responses vs. Age Groups

#for crude schizont extract
PHIST_Data_without_Sukuta$age_categ <- cut(PHIST_Data_without_Sukuta$AGE.YEARS., c(0,2,4,6,8,10))

comparisons <- list(c('(0,2]', '(2,4]'), c('(0,2]', '(4,6]'), c('(0,2]', '(6,8]'), c('(0,2]', '(8,10]'), c('(2,4]', '(4,6]'), c('(2,4]', '(6,8]'), c('(2,4]', '(8,10]'), c('(4,6]', '(6,8]'), c('(4,6]', '(8,10]'),c('(6,8]', '(8,10]'))          

PHIST_Data_without_Sukuta %>% 
  ggplot(aes(age_categ, SCHIZONT.EXTRACT.OD.VALUES))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Crude_schizont extract (OD 492nm)') + xlab('Age_Group (Years)')

# for PF3D7_1401600
PHIST_Data_without_Sukuta$age_categ <- cut(PHIST_Data_without_Sukuta$AGE.YEARS., c(0,2,4,6,8,10))

comparisons <- list(c('(0,2]', '(2,4]'), c('(0,2]', '(4,6]'), c('(0,2]', '(6,8]'), c('(0,2]', '(8,10]'), c('(2,4]', '(4,6]'), c('(2,4]', '(6,8]'), c('(2,4]', '(8,10]'), c('(4,6]', '(6,8]'), c('(4,6]', '(8,10]'),c('(6,8]', '(8,10]'))          

PHIST_Data_without_Sukuta %>% 
  ggplot(aes(age_categ, PF3D7_1401600))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Pf3D7_1401600 (OD 492nm)') + xlab('Age_Group (Years)')

#For PF3D7_1102500

PHIST_Data_without_Sukuta$age_categ <- cut(PHIST_Data_without_Sukuta$AGE.YEARS., c(0,2,4,6,8,10))

comparisons <- list(c('(0,2]', '(2,4]'), c('(0,2]', '(4,6]'), c('(0,2]', '(6,8]'), c('(0,2]', '(8,10]'), c('(2,4]', '(4,6]'), c('(2,4]', '(6,8]'), c('(2,4]', '(8,10]'), c('(4,6]', '(6,8]'), c('(4,6]', '(8,10]'),c('(6,8]', '(8,10]'))          

PHIST_Data_without_Sukuta %>% 
  ggplot(aes(age_categ, PF3D7_1102500))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Pf3D7_1102500 (OD 492nm)') + xlab('Age_Group (Years)')

#For PF3D7_0532400

PHIST_Data_without_Sukuta$age_categ <- cut(PHIST_Data_without_Sukuta$AGE.YEARS., c(0,2,4,6,8,10))

comparisons <- list(c('(0,2]', '(2,4]'), c('(0,2]', '(4,6]'), c('(0,2]', '(6,8]'), c('(0,2]', '(8,10]'), c('(2,4]', '(4,6]'), c('(2,4]', '(6,8]'), c('(2,4]', '(8,10]'), c('(4,6]', '(6,8]'), c('(4,6]', '(8,10]'),c('(6,8]', '(8,10]'))          

PHIST_Data_without_Sukuta %>% 
  ggplot(aes(age_categ, PF3D7_0532400))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Pf3D7_0532400 (OD 492nm)') + xlab('Age_Group (Years)')

#Scaling down
#For schizont

PHIST_Data_without_Sukuta$age_categ <- cut(PHIST_Data_without_Sukuta$AGE.YEARS., c(0,2,4,6,8,10))

comparisons <- list(c('(0,2]', '(2,4]'), c('(0,2]', '(4,6]'), c('(0,2]', '(6,8]'), c('(0,2]', '(8,10]'), c('(2,4]', '(4,6]'), c('(2,4]', '(6,8]'), c('(2,4]', '(8,10]'), c('(4,6]', '(6,8]'), c('(4,6]', '(8,10]'),c('(6,8]', '(8,10]'))          

PHIST_Data_without_Sukuta %>% 
  ggplot(aes(age_categ, SCHIZONT.EXTRACT.OD.VALUES))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Crude_schizont_extract (OD 492nm)') + xlab('Age_Group (Years)')+ ylim(0, 3)

PHIST_Data_without_Sukuta$age_categ <- cut(PHIST_Data_without_Sukuta$AGE.YEARS., c(0,2,4,6,8,10))

#For 0532400

comparisons <- list(c('(0,2]', '(2,4]'), c('(0,2]', '(4,6]'), c('(0,2]', '(6,8]'), c('(0,2]', '(8,10]'), c('(2,4]', '(4,6]'), c('(2,4]', '(6,8]'), c('(2,4]', '(8,10]'), c('(4,6]', '(6,8]'), c('(4,6]', '(8,10]'),c('(6,8]', '(8,10]'))          

PHIST_Data_without_Sukuta %>% 
  ggplot(aes(age_categ, PF3D7_0532400))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Pf3D7_0532400 (OD 492nm)') + xlab('Age_Group (Years)')+ ylim(0, 3)

#for 1102500

PHIST_Data_without_Sukuta$age_categ <- cut(PHIST_Data_without_Sukuta$AGE.YEARS., c(0,2,4,6,8,10))

comparisons <- list(c('(0,2]', '(2,4]'), c('(0,2]', '(4,6]'), c('(0,2]', '(6,8]'), c('(0,2]', '(8,10]'), c('(2,4]', '(4,6]'), c('(2,4]', '(6,8]'), c('(2,4]', '(8,10]'), c('(4,6]', '(6,8]'), c('(4,6]', '(8,10]'),c('(6,8]', '(8,10]'))          

PHIST_Data_without_Sukuta %>% 
  ggplot(aes(age_categ, PF3D7_1102500))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Pf3D7_1102500 (OD 492nm)') + xlab('Age_Group (Years)')+ ylim(0, 1.5)

#for 1401600

PHIST_Data_without_Sukuta$age_categ <- cut(PHIST_Data_without_Sukuta$AGE.YEARS., c(0,2,4,6,8,10))

comparisons <- list(c('(0,2]', '(2,4]'), c('(0,2]', '(4,6]'), c('(0,2]', '(6,8]'), c('(0,2]', '(8,10]'), c('(2,4]', '(4,6]'), c('(2,4]', '(6,8]'), c('(2,4]', '(8,10]'), c('(4,6]', '(6,8]'), c('(4,6]', '(8,10]'),c('(6,8]', '(8,10]'))          

PHIST_Data_without_Sukuta %>% 
  ggplot(aes(age_categ, PF3D7_1401600))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Pf3D7_1401600 (OD 492nm)') + xlab('Age_Group (Years)')+ ylim(0, 1)

#RECATEGORIZING BASED ON THE BOB SNOW DATA

View(PHIST_Data_without_Sukuta)
library(ggplot2)

#for 1401600

PHIST_Data_without_Sukuta$age_categ <- cut(PHIST_Data_without_Sukuta$AGE.YEARS., c(1,2,3,5,6,8,9,11))

comparisons <- list(c('(1,2]', '(3,5]'), c('(1,2]', '(6,8]'), c('(1,2]', '(9,11]'), c('(3,5]', '(6,8]'), c('(3,5]', '(9,11]'), c('(6,8]', '(9,11]'))
PHIST_Data_without_Sukuta %>% 
  ggplot(aes(age_categ, PF3D7_1401600))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red", na.rm = TRUE) +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Pf3D7_1401600 (OD 492nm)') + xlab('Age_Group (Years)')+ ylim(0, 1)

library(ggpubr)
library("tidyverse")
data<-read.csv("PHIST_Data_without_Sukuta.csv",header=TRUE)
data

#REFLECT BOB SNOW DATA
#for 1401600
data$age_categ <- cut(data$AGE.YEARS., c(1,2,3,5,6,8,9,11))
data2=data[!is.na(data$age_categ),]
data3=data2[data2$age_categ!='(2,3]',]
data4=data3[data3$age_categ!='(5,6]',]
data5=data4[data4$age_categ!='(8,9]',]
comparisons <- list(c('(1,2]', '(3,5]'), c('(1,2]', '(6,8]'), c('(1,2]', '(9,11]'), c('(3,5]', '(6,8]'), c('(3,5]', '(9,11]'), c('(6,8]', '(9,11]'))
data5 %>%
  ggplot(aes(age_categ, PF3D7_1401600))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red", na.rm = TRUE) +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Pf3D7_1401600 (OD 492nm)') + xlab('Age_Group (Years)')+ ylim(0, 1)

#1102500


data$age_categ <- cut(data$AGE.YEARS., c(1,2,3,5,6,8,9,11))
data2=data[!is.na(data$age_categ),]
data3=data2[data2$age_categ!='(2,3]',]
data4=data3[data3$age_categ!='(5,6]',]
data5=data4[data4$age_categ!='(8,9]',]
comparisons <- list(c('(1,2]', '(3,5]'), c('(1,2]', '(6,8]'), c('(1,2]', '(9,11]'), c('(3,5]', '(6,8]'), c('(3,5]', '(9,11]'), c('(6,8]', '(9,11]'))
data5 %>%
  ggplot(aes(age_categ, PF3D7_1102500))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red", na.rm = TRUE) +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Pf3D7_1102500 (OD 492nm)') + xlab('Age_Group (Years)')+ ylim(0, 1)

#for 0532400
data$age_categ <- cut(data$AGE.YEARS., c(1,2,3,5,6,8,9,11))
data2=data[!is.na(data$age_categ),]
data3=data2[data2$age_categ!='(2,3]',]
data4=data3[data3$age_categ!='(5,6]',]
data5=data4[data4$age_categ!='(8,9]',]
comparisons <- list(c('(1,2]', '(3,5]'), c('(1,2]', '(6,8]'), c('(1,2]', '(9,11]'), c('(3,5]', '(6,8]'), c('(3,5]', '(9,11]'), c('(6,8]', '(9,11]'))
data5 %>%
  ggplot(aes(age_categ, PF3D7_0532400))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red", na.rm = TRUE) +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Pf3D7_0532400 (OD 492nm)') + xlab('Age_Group (Years)')+ ylim(0, 2)

#for schizont extract

data$age_categ <- cut(data$AGE.YEARS., c(1,2,3,5,6,8,9,11))
data2=data[!is.na(data$age_categ),]
data3=data2[data2$age_categ!='(2,3]',]
data4=data3[data3$age_categ!='(5,6]',]
data5=data4[data4$age_categ!='(8,9]',]
comparisons <- list(c('(1,2]', '(3,5]'), c('(1,2]', '(6,8]'), c('(1,2]', '(9,11]'), c('(3,5]', '(6,8]'), c('(3,5]', '(9,11]'), c('(6,8]', '(9,11]'))
data5 %>%
  ggplot(aes(age_categ, SCHIZONT.EXTRACT.OD.VALUES))+
  geom_boxplot()+
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red", na.rm = TRUE) +
  facet_wrap(~location, scales = 'free')+
  theme_bw() +
  ylab('Crude_schizont_extract (OD 492nm)') + xlab('Age_Group (Years)')+ ylim(0, 3)






##ANTIBODY RESPONSES VS. PARASITE POSITIVE AND NEGATIVE

# generating beeswarm plots
boxplot(PF3D7_0532400 ~ PARA, data = PHIST_Data_without_Sukuta, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'Comparing Antibody Levels in Children with  Parasitemia
')+
beeswarm(PF3D7_0532400 ~ location, data = PHIST_Data_without_Sukuta, 
         col = 4, pch = 16, add = TRUE)+
  stat_compare_means( method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme_bw()

boxplot(PF3D7_1102500 ~ PARA, data = PHIST_Data_without_Sukuta, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'Comparing Antibody Levels in Children with  Parasitemia
')+
  beeswarm(PF3D7_1102500 ~ PARA, data = PHIST_Data_without_Sukuta, 
           col = 4, pch = 16, add = TRUE)+
  stat_compare_means( method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme_bw()

boxplot(SCHIZONT.EXTRACT.OD.VALUES ~ PARA, data = PHIST_Data_without_Sukuta, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'Comparing Antibody Levels in Children with  Parasitemia
')+
  beeswarm(SCHIZONT.EXTRACT.OD.VALUES ~ PARA, data = PHIST_Data_without_Sukuta, 
           col = 4, pch = 16, add = TRUE)+
  stat_compare_means( method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme_bw()

boxplot(PF3D7_0532400 ~ PARA, data = PHIST_Data_without_Sukuta, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'Comparing Antibody Levels in Children with  Parasitemia
')+
  beeswarm(PF3D7_0532400 ~ PARA, data = PHIST_Data_without_Sukuta, 
           col = 4, pch = 16, add = TRUE)+
  stat_compare_means( method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme_bw()

boxplot(PF3D7_0532400 ~ PARA, data = PHIST_Data_without_Sukuta, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        main = 'Comparing Antibody Levels in Children with  Parasitemia
')+
  beeswarm(PF3D7_0532400 ~ PARA, data = PHIST_Data_without_Sukuta, 
           col = 4, pch = 16, add = TRUE)+
  stat_compare_means( method = 'wilcox.test')+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme_bw()





##REGRESSION MODELS
mydata <- read.csv("C:/Users/Tony Isebe/Desktop/Backup/PHIST RESEARCH/csv raw/Junju_Metadata_PHISTb_Antigens.csv")

library(ggplot2)

#Univariate model with elisa od only 0532400
Pf053 <- glm (manytm6 ~ PF3D7_0532400, data = mydata, family = poisson)
str(Pf053)
Pf053$deviance
#converting the co-efficients to odd ratios
exp(coef(Pf053))
confint(Pf053)

exp(cbind(OR = coef(Pf053), confint(Pf053)))

#Univariate model with elisa od only 1102500
Pf110 <- glm (manytm6 ~ PF3D7_1102500, data = mydata, family = poisson)
str(Pf110)
Pf110$deviance
#converting the co-efficients to odd ratios
exp(coef(Pf110))
confint(Pf110)

exp(cbind(OR = coef(Pf110), confint(Pf110)))

#Univariate model with elisa od only 1401600
Pf140 <- glm (manytm6 ~ PF3D7_1401600, data = mydata, family = poisson)
str(Pf140)
Pf140$deviance
#converting the co-efficients to odd ratios
exp(coef(Pf140))
confint(Pf140)

exp(cbind(OR = coef(Pf140), confint(Pf140)))


#m

#pf3d7_0532400

mydata.maypos<-subset(mydata,maypfpos=="0")

Pf053_1 <- glm (tm6 ~ PF3D7_0532400, data = mydata.maypos, family = binomial)
summary(Pf053_1)

#converting the co-efficients to odd ratios
exp(coef(Pf053_1))
confint(Pf053_1)

exp(cbind(OR = coef(Pf053_1), confint(Pf053_1)))

#pf3d7_1102500

mydata.maypos<-subset(mydata,maypfpos=="0")

Pf110_1 <- glm (tm6 ~ PF3D7_1102500, data = mydata.maypos, family = binomial)
summary(Pf110_1)

#converting the co-efficients to odd ratios
exp(coef(Pf110_1))
confint(Pf110_1)

exp(cbind(OR = coef(Pf110_1), confint(Pf110_1)))

#pf3d7_1401600

mydata.maypos<-subset(mydata,maypfpos=="0")

Pf140_1 <- glm (tm6 ~ PF3D7_1401600, data = mydata.maypos, family = binomial)
summary(Pf140_1)

#converting the co-efficients to odd ratios
exp(coef(Pf140_1))
confint(Pf140_1)

exp(cbind(OR = coef(Pf140_1), confint(Pf140_1)))

#0532400
cutN <- function(X , n = 4){
  cut(
    X ,
    include.lowest = TRUE ,
    breaks = quantile(
      X , 
      probs = (0:n)/n ,
      na.rm = TRUE ))}

mydata$tert<-cutN(mydata$PF3D7_0532400 , n = 3 )
mydata$tert

mydata$tert <- ordered(mydata$tert,
                       labels = c("Low", "Medium", "High"))





hilow<-subset(mydata, tert!="Medium")

Pf053_2 <- glm (tm6 ~ tert, data = hilow, family = binomial)
summary(Pf053_2)

#converting the co-efficients to odd ratios
exp(coef(Pf053_2))
confint(Pf053_2)

exp(cbind(OR = coef(Pf053_2), confint(Pf053_2)))

s<-ggplot(mydata,aes(tert,PF3D7_0532400))+theme_bw(base_size = 25)
s+geom_dotplot(binaxis="y",stackdir = "center",aes(fill=as.factor(tm6)))+
  xlab("Antibody level")+ylab("Pf3D7_0532400")+ scale_fill_discrete(name = "Malaria",
                                                                    labels = c("Neg", "Pos"))
#1102500

cutN <- function(X , n = 4){
  cut(
    X ,
    include.lowest = TRUE ,
    breaks = quantile(
      X , 
      probs = (0:n)/n ,
      na.rm = TRUE ))}

mydata$tert<-cutN(mydata$PF3D7_1102500 , n = 3 )
mydata$tert

mydata$tert <- ordered(mydata$tert,
                       labels = c("Low", "Medium", "High"))





hilow<-subset(mydata, tert!="Medium")

Pf110_2 <- glm (tm6 ~ tert, data = hilow, family = binomial)
summary(Pf053_2)

#converting the co-efficients to odd ratios
exp(coef(Pf110_2))
confint(Pf110_2)

exp(cbind(OR = coef(Pf110_2), confint(Pf110_2)))

s<-ggplot(mydata,aes(tert,PF3D7_1102500))+theme_bw(base_size = 25)
s+geom_dotplot(binaxis="y",stackdir = "center",aes(fill=as.factor(tm6)))+
  xlab("Antibody level")+ylab("Pf3D7_1102500")+ scale_fill_discrete(name = "Malaria",
                                                                    labels = c("Neg", "Pos"))


#1401600

cutN <- function(X , n = 4){
  cut(
    X ,
    include.lowest = TRUE ,
    breaks = quantile(
      X , 
      probs = (0:n)/n ,
      na.rm = TRUE ))}

mydata$tert<-cutN(mydata$PF3D7_1401600 , n = 3 )
mydata$tert

mydata$tert <- ordered(mydata$tert,
                       labels = c("Low", "Medium", "High"))





hilow<-subset(mydata, tert!="Medium")

Pf140_2 <- glm (tm6 ~ tert, data = hilow, family = binomial)
summary(Pf053_2)

#converting the co-efficients to odd ratios
exp(coef(Pf140_2))
confint(Pf140_2)

exp(cbind(OR = coef(Pf140_2), confint(Pf140_2)))

s<-ggplot(mydata,aes(tert,PF3D7_1401600))+theme_bw(base_size = 25)
s+geom_dotplot(binaxis="y",stackdir = "center",aes(fill=as.factor(tm6)))+
  xlab("Antibody level")+ylab("Pf3D7_1401600")+ scale_fill_discrete(name = "Malaria",
                                                                    labels = c("Neg", "Pos"))











