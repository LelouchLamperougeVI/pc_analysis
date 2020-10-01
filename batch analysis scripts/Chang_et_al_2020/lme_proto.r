# enough screwing around, this is the real thing
dat = read.csv("/home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/inkscape/MATLAB/test.csv") # load data
anovaModelRM <- aov(response ~ type*time + Error(subject/time), dat)
summary(anovaModelRM)

dat$interac <- interaction(dat$type, dat$time)
require(nlme)
#mixed_model = lme(response ~ interac, random = ~1 | subject/time, dat)
mixed_model = lme(response ~ interac, random = ~1 | subject, dat)
#mixed_model = lme(response ~ interac, random = ~ 0 | subject, dat)
anova(mixed_model)

require(multcomp)
summary(glht(mixed_model, linfct=mcp(interac = "Tukey")), test = adjusted(type = "bonferroni"))


library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 18))

ggplot(dat, aes(sample=response)) +
  stat_qq()

# sphericity test
neo <- matrix(dat$response, nrow = 321, ncol = 4)
mauchly.test(lm(neo ~ 1), X = ~ 1)

d_bycond = na.omit(dat) %>%
  group_by(type, time) %>%
  summarise(mean_response = mean(response))

ggplot(d_bycond, aes(x=time, y=mean_response, 
                     colour=type, group=type)) +
  geom_line(size=2) + geom_point(size=5, shape=21, fill="white")