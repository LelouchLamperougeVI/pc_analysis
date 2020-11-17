require(nlme)
require(multcomp)
require(dplyr)
library(ggplot2)

# enough screwing around, this is the real thing
dat = read.csv("/home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/inkscape/MATLAB/test.csv") # load data
dat$time <- factor(dat$time)
dat$type <- factor(dat$type)
dat$subject <- factor(dat$subject)
levels(dat$time) <- c('before', 'early', 'during', 'after')
levels(dat$type) <- c('none', 'traj', 'cue')

md.aov <- aov(response ~ type*time + Error(subject/(type*time)), dat)
summary(md.aov)

dat$interac <- interaction(dat$type, dat$time)
md.lme = lme(response ~ interac, random = ~ 1 | subject/time, dat)
summary(md.lme)
anova(md.lme)

test_results <- summary(glht(md.lme, linfct=mcp(interac = "Tukey")), test = adjusted(type = "bonferroni"))

se <- function(x) {
  if(!(is.double(x))) { return(NA) }
  return( sd(x)/sqrt(length(x)) )
}
dat$prediction <- predict(md.lme)
dat_summary <- dat %>% group_by(time, type) %>% summarise(dat_mu = mean(response), dat_se = se(response), p_mu = mean(prediction), p_se = se(prediction))

ggplot(dat_summary, aes(x=time, y=dat_mu, ymin=dat_mu-dat_se, ymax=dat_mu+dat_se, colour=type)) + geom_pointrange() + geom_line(mapping=aes(x=time, y=p_mu, group=type, color=type))

# data exploration
ggplot(aes(time, response, color=time), data=dat) + geom_boxplot() + facet_wrap(~ type)

library(lme4)
md.lmer <- lmer(response ~ type*time + (1 | time), data=dat)
summary(md.lmer)
plot(md.lmer)

qqnorm(resid(md.lmer))
qqline(resid(md.lmer))

ggplot(aes(time, response, color=time), data=dat) + geom_point() + facet_wrap(~ type) + geom_boxplot(data = cbind(dat, pred = predict(md.lmer)), aes(y=pred))

# todo: a better model would be response ~ type*time + (1 | session/ensemble/time) accounting for the recording session
# better even, response ~ type*time + (1 | session/ensemble/neuron/time) to perform the analysis at the single neuron level

### new table
dat = read.csv("/home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/inkscape/MATLAB/r1.csv") # load data
dat <- dat %>% tidyr::gather(key = time, value = response, t0:t3) %>% arrange(neuron)
dat$type <- factor(dat$type)
dat$session <- factor(dat$session)
dat$ensemble <- factor(dat$ensemble)
levels(dat$type) <- c('none', 'cue', 'traj')

library(lmerTest)
library(lme4)
md.lmer <- lmer(response ~ type*time + (1 | session/ensemble) + (1 | time), data=dat)
summary(md.lmer)
plot(md.lmer)

as.numeric(VarCorr(md.lmer))

qqnorm(resid(md.lmer))
qqline(resid(md.lmer))

ggplot(aes(time, response, color=time), data=dat) + geom_point() + facet_wrap(~ type) + geom_boxplot(data = cbind(dat, pred = predict(md.lmer)), aes(y=pred))

dat$interac <- interaction(dat$type, dat$time)
md.lme = lme(response ~ interac, random = ~ 1 | session/ensemble, dat)
summary(glht(md.lme, linfct=mcp(interac = "Tukey")), test = adjusted(type = "bonferroni"))