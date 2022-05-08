library(tidyverse)
library(ggpubr)
library(rstatix)
library(MASS)

#load them data
dat = read.csv("/home/loulou/Documents/my_docs/Manuscripts/Chang_et_al_2020/R/err.csv")

# arrange them data
dat$treatment <- factor(dat$lesion, levels = c(0, 1), labels = c("sham", "lesion"))
dat$session <- c(1:nrow(dat))
dat <- dat %>%
  gather(key = "position", value = "error", x1:x50) %>%
  convert_as_factor(session, treatment)
dat$position <- as.numeric(str_extract(dat$position, "\\d+"))
dat$position <- (dat$position - 1) / 49 * 150

# boxcox transform data for normality
out <- boxcox(dat$error ~ dat$treatment, lambda = seq(-2, 2, .001))
lambda <- out$x[which.max(out$y)]
transformed <- (dat$error^lambda - 1) / lambda

dat$error <- transformed

# plot them data
ggline(dat, x = "position", y = "error", color = "treatment", add = c("mean_se", "jitter"))

# check for outliers
md.outliers <- dat %>% group_by(treatment, position) %>% identify_outliers(error)
any(md.outliers$is.extreme)

# shapiro test for normality
md.shapiro <- dat %>% group_by(treatment, position) %>% shapiro_test(error)
any(md.shapiro$p < .05)

# ANOVA them data (two-way repeated measures)
md.aov <- anova_test(data = dat, dv = error, wid = session, within = position, between = treatment)
md.aov_tab <- get_anova_table(md.aov)
md.aov_tab

md.residuals <- residuals(attr(md.aov, "args")$model)
qqnorm(md.residuals)
qqline(md.residuals)

md.mult <- dat %>%
  group_by(position) %>%
  anova_test(dv = error, wid = session, between = treatment) %>%
  get_anova_table() %>%  
  adjust_pvalue(method = "fdr") %>%
  #adjust_pvalue(method = "bonferroni") %>%
  add_significance(p.col = "p.adj", cutpoints = c(0, .05, 1), symbols = c("*", ""))
md.mult

# plot them anova
md.summary <- dat %>%
  group_by(treatment, position) %>%
  summarize(mean = mean(error), se = sd(error)/sqrt(length(error)))

ggplot(data = md.summary, mapping = aes(x = position, y = mean, group = treatment, color = treatment)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.1) +
  geom_text(data = md.mult, aes(x = position, y = 7, label = p.adj.signif, group = 1, color = NULL, size = 100)) +
  theme_classic2()

