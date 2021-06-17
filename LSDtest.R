
rm(list = ls(all = TRUE))
graphics.off()
shell("cls")

data <- read_excel("marker.xlsx")
head(data)
str(data)
data$status <- as.factor(x = data$status)
data$drug <- as.factor(x = data$drug)
data$sex <- as.factor(x = data$sex)
str(data)
attach(data)
library(stats)
#Fit an Analysis of Variance Model
model.aov <- aov(biomarker ~ 
                 status +
                 drug +
                 sex +
                 status:drug +
                 status:sex +
                 drug:sex +
                 status:drug:sex)

anova(model.aov)

model.lm <- lm(biomarker ~ 
                   status +
                   drug +
                   sex +
                   status:drug +
                   status:sex +
                   drug:sex +
                   status:drug:sex)

anova(model.lm)

library(agricolae)

# First factor variable (status)
LSDstatus = LSD.test(y = biomarker,
                 trt = status,
                 DFerror = model.aov$df.residual,
                 MSerror = deviance(model.aov)/model.aov$df.residual,
                 alpha = 0.05,
                 p.adj = "bonferroni",
                 group = TRUE,
                 console = TRUE)

# second factor variable (drug)
LSDdrug = LSD.test(y = biomarker,
                     trt = drug,
                     DFerror = model.aov$df.residual,
                     MSerror = deviance(model.aov)/model.aov$df.residual,
                     alpha = 0.05,
                     p.adj = "bonferroni",
                     group = TRUE,
                     console = TRUE)

# 3rd factor variable (sex)
LSDsex = LSD.test(y = biomarker,
                     trt = sex,
                     DFerror = model.aov$df.residual,
                     MSerror = deviance(model.aov)/model.aov$df.residual,
                     alpha = 0.05,
                     p.adj = "bonferroni",
                     group = TRUE,
                     console = TRUE)

# 4th factor variable (status:drug)
LSD4 = LSD.test(y = biomarker,
                trt = status:drug,
                DFerror = model.aov$df.residual,
                MSerror = deviance(model.aov)/model.aov$df.residual,
                alpha = 0.05,
                p.adj = "bonferroni",
                group = TRUE,
                console = TRUE)

# 5th factor variable (status:sex)
LSD5 = LSD.test(y = biomarker,
                trt = status:sex,
                DFerror = model.aov$df.residual,
                MSerror = deviance(model.aov)/model.aov$df.residual,
                alpha = 0.05,
                p.adj = "bonferroni",
                group = TRUE,
                console = TRUE)
# 6th factor variable (drug:sex)
LSD6 = LSD.test(y = biomarker,
                trt = drug:sex,
                DFerror = model.aov$df.residual,
                MSerror = deviance(model.aov)/model.aov$df.residual,
                alpha = 0.05,
                p.adj = "bonferroni",
                group = TRUE,
                console = TRUE)

# 7th factor variable (status:drug:sex)
LSD7 = LSD.test(y = biomarker,
                trt = status:drug:sex,
                DFerror = model.aov$df.residual,
                MSerror = deviance(model.aov)/model.aov$df.residual,
                alpha = 0.05,
                p.adj = "bonferroni",
                group = TRUE,
                console = TRUE)

library(dplyr) 

MeanSE1 = data %>%
  group_by(status) %>%
  summarise(avg_1 = mean(biomarker),
            se = sd(biomarker)/sqrt(length(biomarker)))

print(MeanSE1)
attach(MeanSE1)

MeanSE2 = data %>%
  group_by(drug) %>%
  summarise(avg_2 = mean(biomarker),
            se = sd(biomarker)/sqrt(length(biomarker)))

print(MeanSE2)
attach(MeanSE2)

MeanSE3 = data %>%
  group_by(sex) %>%
  summarise(avg_3 = mean(biomarker),
            se = sd(biomarker)/sqrt(length(biomarker)))

print(MeanSE3)
attach(MeanSE3)

MeanSE4 = data %>%
  group_by(status, drug) %>%
  summarise(avg_4 = mean(biomarker), .groups = 'drop',
            se = sd(biomarker)/sqrt(length(biomarker)))

print(MeanSE4)
attach(MeanSE4)

MeanSE5 = data %>%
  group_by(status, sex) %>%
  summarise(avg_5 = mean(biomarker), .groups = 'drop',
            se = sd(biomarker)/sqrt(length(biomarker)))

print(MeanSE5)
attach(MeanSE5)

MeanSE6 = data %>%
  group_by(drug, sex) %>%
  summarise(avg_6 = mean(biomarker), .groups = 'drop',
            se = sd(biomarker)/sqrt(length(biomarker)))

print(MeanSE6)
attach(MeanSE6)

MeanSE7 = data %>%
  group_by(status, drug, sex) %>%
  summarise(avg_7 = mean(biomarker), .groups = 'drop',
            se = sd(biomarker)/sqrt(length(biomarker)))

print(MeanSE7)
attach(MeanSE7)

## Original order of groups (for lettering)

ascend1 = LSDstatus$groups %>%
  group_by(rownames(LSDstatus$groups)) %>%
  arrange(rownames(LSDstatus$groups))
print(ascend1)

ascend2 = LSDdrug$groups %>%
  group_by(rownames(LSDdrug$groups)) %>%
  arrange(rownames(LSDdrug$groups))
print(ascend2)

ascend3 = LSDsex$groups %>%
  group_by(rownames(LSDsex$groups)) %>%
  arrange(rownames(LSDsex$groups))
print(ascend3)

ascend4 = LSD4$groups %>%
  group_by(rownames(LSD4$groups)) %>%
  arrange(rownames(LSD4$groups))
print(ascend4)

ascend5 = LSD5$groups %>%
  group_by(rownames(LSD5$groups)) %>%
  arrange(rownames(LSD5$groups))
print(ascend5)

ascend6 = LSD6$groups %>%
  group_by(rownames(LSD6$groups)) %>%
  arrange(rownames(LSD6$groups))
print(ascend6)

ascend7 = LSD7$groups %>%
  group_by(rownames(LSD7$groups)) %>%
  arrange(rownames(LSD7$groups))
print(ascend7)

## Create plotting object

library(ggplot2)

# plot for status
p1 = ggplot(MeanSE1, aes(x = status, y = avg_1))
print(p1)

plot1 = p1 + 
  geom_bar(stat = "identity",
           color = "black",
           position = position_dodge(width=0.9))
# Adding error bars       
plot1a = plot1 + 
  geom_errorbar(aes(ymax = avg_1 + se,
                    ymin = avg_1 - se), 
                position = position_dodge(width=0.9), 
                width = 0.25)
# Changing main title, X & Y labels
plot1b = plot1a +
  labs(title = "",
       x = "status",
       y = "serBilir Biomarker")

# Adding lettering from the test applied (LSD$group)

plot1c = plot1b +
  geom_text(aes(x = status,
                y = avg_1 + se,
                label = as.matrix(ascend1$groups)),
            position = position_dodge(width = 0.9),
            vjust = -(0.5))
plot1c

..........................................................
# plot for drug factor variable
p2 = ggplot(MeanSE2, aes(x = drug, y = avg_2))
print(p2)

plot2 = p2 + 
  geom_bar(stat = "identity",
           color = "black",
           position = position_dodge(width=0.9),
           width = 0.8)
# Adding error bars       
plot2a = plot2 + 
  geom_errorbar(aes(ymax = avg_2 + se,
                    ymin = avg_2 - se), 
                position = position_dodge(width=0.9), 
                width = 0.25)
# Changing main title, X & Y labels
plot2b = plot2a +
  labs(title = "",
       x = "drug",
       y = "serBilir Biomarker")
# Adding lettering from the test applied (LSD$group)

plot2c = plot2b +
  geom_text(aes(x = drug,
                y = avg_2 + se,
                label = as.matrix(ascend2$groups)),
            position = position_dodge(width = 0.9),
            vjust = -(0.5))
plot2c

....................................................................
### plot for sex factor variable
p3 = ggplot(MeanSE3, aes(x = sex, y = avg_3))
print(p3)

plot3 = p3 + 
  geom_bar(stat = "identity",
           color = "black",
           position = position_dodge(width=0.9),
           width = 0.8)
# Adding error bars       
plot3a = plot3 + 
  geom_errorbar(aes(ymax = avg_3 + se,
                    ymin = avg_3 - se), 
                position = position_dodge(width=0.9), 
                width = 0.25)
# Changing main title, X & Y labels
plot3b = plot3a +
  labs(title = "",
       x = "sex",
       y = "serBilir Biomarker")
# Adding lettering from the test applied (LSD$group)

plot3c = plot3b +
  geom_text(aes(x = sex,
                y = avg_3 + se,
                label = as.matrix(ascend3$groups)),
            position = position_dodge(width = 0.9),
            vjust = -(0.5))
plot3c

............................................................................
### plot for status and drug interaction
p4 = ggplot(MeanSE4, aes(x = drug, y = avg_4, fill = factor(status)))
print(p4)

plot4 = p4 + 
  geom_bar(stat = "identity",
           color = "black",
           position = position_dodge(width=0.9))
# Adding color for fill and changing legend text
plot4a = plot4 + 
  scale_fill_manual(values = gray(1:3/3),
                    labels = c("Dead", "Alive", "transplanted"))

# Adding error bars       

plot4b = plot4a + 
  geom_errorbar(aes(ymax = avg_4 + se,
                    ymin = avg_4 - se), 
                position = position_dodge(width=0.9), 
                width = 0.25)
# Changing main title, X & Y labels
plot4c = plot4b +
  labs(title = "",
       x = "drug",
       y = "serBilir Biomarker",
       fill = "status")
# Adding lettering from the test applied (LSD$group)

plot4d = plot4c +
  geom_text(aes(x = drug,
                y = avg_4 + se,
                label = as.matrix(ascend4$groups)),
            position = position_dodge(width = 0.9),
            vjust = -(0.5))
plot4d

.............................................................................
### plot for status and sex interaction
# Adding lettering from the test applied (LSD$group)

p5 = ggplot(MeanSE5, aes(x = sex, y = avg_5, fill = factor(status)))
print(p5)

plot5 = p5 + 
  geom_bar(stat = "identity",
           color = "black",
           position = position_dodge(width=0.9))
# Adding color for fill and changing legend text
plot5a = plot5 + 
  scale_fill_manual(values = gray(1:3/3),
                    labels = c("Dead", "Alive", "Transplanted"))
# Adding error bars       
plot5b = plot5a + 
  geom_errorbar(aes(ymax = avg_5 + se,
                    ymin = avg_5 - se), 
                position = position_dodge(width=0.9), 
                width = 0.25)
# Changing main title, X & Y labels
plot5c = plot5b +
  labs(title = "",
       x = "sex",
       y = "serBilir Biomarker",
       fill = "status")

plot5d = plot5c +
  geom_text(aes(x = sex,
                y = avg_5 + se,
                label = as.matrix(ascend5$groups)),
            position = position_dodge(width = 0.9),
            vjust = -(0.5))
plot5d

...............................................................................
### plot for drug and sex interaction
# Adding lettering from the test applied (LSD$group)

p6 = ggplot(MeanSE6, aes(x = sex, y = avg_6, fill = factor(drug)))
print(p6)

plot6 = p6 + 
  geom_bar(stat = "identity",
           color = "black",
           position = position_dodge(width=0.9))
# Adding color for fill and changing legend text
plot6a = plot6 + 
  scale_fill_manual(values = gray(1:4/4),
                    labels = c("D-penicil", "placebo"))
# Adding error bars       
plot6b = plot6a + 
  geom_errorbar(aes(ymax = avg_6 + se,
                    ymin = avg_6 - se), 
                position = position_dodge(width=0.9), 
                width = 0.25)
# Changing main title, X & Y labels
plot6c = plot6b +
  labs(title = "",
       x = "sex",
       y = "serBilir Biomarker",
       fill = "drug")

plot6d = plot6c +
  geom_text(aes(x = sex,
                y = avg_6 + se,
                label = as.matrix(ascend6$groups)),
            position = position_dodge(width = 0.9),
            vjust = -(0.5))
plot6d

............................................................................................
### plot for drug and sex interaction
# Adding lettering from the test applied (LSD$group)

p7 = ggplot(MeanSE7, aes(x = sex, y = avg_7, fill = factor(status:drug)))
print(p7)

plot7 = p7 + 
  geom_bar(stat = "identity",
           color = "black",
           position = position_dodge(width=0.9))
# Adding color for fill and changing legend text
plot7a = plot7 + 
  scale_fill_manual(values = gray(1:6/6),
                    labels = c("alive:D-penicil", "alive:placebo", "dead:D-penicil",
                    "dead:placebo", "transplanted:D-penicil", "transplanted:placebo"))
# Adding error bars       
plot7b = plot7a + 
  geom_errorbar(aes(ymax = avg_7 + se,
                    ymin = avg_7 - se), 
                position = position_dodge(width=0.9), 
                width = 0.25)
# Changing main title, X & Y labels
plot7c = plot7b +
  labs(title = "",
       x = "sex",
       y = "serBilir",
       fill = "status:drug")

plot7d = plot7c +
  geom_text(aes(x = sex,
                y = avg_7 + se,
                label = as.matrix(ascend7$groups)),
            position = position_dodge(width = 0.9),
            vjust = -(0.5))
plot7d

-------------------------------------------------------------------
p7 = ggplot(MeanSE7, aes(x = status, y = avg_7, fill = factor(sex:drug)))
print(p7)

plot7 = p7 + 
  geom_bar(stat = "identity",
           color = "black",
           position = position_dodge(width=0.9))
# Adding color for fill and changing legend text
plot7a = plot7 + 
  scale_fill_manual(values = gray(1:6/6),
                    labels = c("male:D-penicil",  "female:D-penicil", "male:placebo",
                               "female:placebo"))
# Adding error bars       
plot7b = plot7a + 
  geom_errorbar(aes(ymax = avg_7 + se,
                    ymin = avg_7 - se), 
                position = position_dodge(width=0.9), 
                width = 0.25)
# Changing main title, X & Y labels
plot7c = plot7b +
  labs(title = "",
       x = "status",
       y = "serBilir biomarker",
       fill = "sex:drug")

plot7d = plot7c +
  geom_text(aes(x = status,
                y = avg_7 + se,
                label = as.matrix(ascend7$groups)),
            position = position_dodge(width = 0.9),
            vjust = -(0.5))
plot7d
