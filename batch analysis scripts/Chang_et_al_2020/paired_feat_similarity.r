#load them data
dat = read.csv("/home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/R/paired.csv") # load data
dat$coupled <- factor(dat$coupled, levels = c(0, 1), labels = c("non", "coup"))
dat$cue <- factor(dat$cue, levels = c(1, 2, 3), labels = c("3", "2", "1"))

#atanh transform
dat$r <- atanh(dat$r)
head(dat)

# plot data
library("ggpubr")
ggboxplot(dat, x = "cue", y = "r", color = "coupled")

#run two-way anova
# note to future self: I know this shit looks as though as it should be a mixed anova, but these are actually two between factors
md.aov <- aov(r ~ cue * coupled, data = dat)
summary(md.aov)

#check homogeneity - satisfied
plot(md.aov, 1)
library(car)
leveneTest(r ~ cue * coupled, data = dat)

#check normality - satisfied
plot(md.aov, 2)
md.residuals <- residuals(object = md.aov)
shapiro.test(x = md.residuals)

#plot diagnostics
par(mfrow=c(2,2))
plot(md.aov)
par(mfrow=c(1,1))

#type II anova for unequal samples, since no interaction
Anova(md.aov, type = "II")

#post hocs
library(multcomp)
md.mult <- glht(md.aov, linfct = mcp(coupled = "Tukey"))
summary(md.mult)
TukeyHSD(md.aov)