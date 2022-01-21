#!/usr/bin/env Rscript

## version history -------------------------------------------------------------
# v1 by J.T. 6 Jan 2022, v2 by B.P. 7 Jan 2022, fig1 by B.P. 9 Jan 2022
# updated input data, analyzed logistic fits
# [JT] v3 by J. T. 11 Jan 2022
# [JT] more input and output formatting; code testing and commenting
# B.P. v4, updated CIs to include dates b/w 11-30 and 12-31
# B.P. v5, added GISAID data and MA curves
# B.P. v6, added MA data post-university dates
# B.P. v7, added MA and NE new cases per day
# B.P. v8, added inverse predictions

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/Volumes/GoogleDrive/My Drive/Omicron/omicron_repo")
suppressPackageStartupMessages(library(tidyverse)) # [JT] load libs in script
options(stringsAsFactors=FALSE)
theme_set(theme_classic())
library('msm')

# global variables
cols <- c(Delta="#bdbdbd", Omicron="#ca0020")
casecols <- c(colleges="#ca0020", MA="black", NE='darkgrey')

## inputs ----------------------------------------------------------------------
# institutions: every case
inst <- read.csv("input-data-freeze/combined-inst.csv",
                 na.strings="",
                 colClasses=c(Date="Date")) %>% # [JT] set col classes while reading
        filter(Variant %in% c("Delta", "Omicron"), # [JT] only keep delta and omicron
               Date >= "2021-12-02",
               Date <= "2021-12-21")

# MA data (GISAID)
ma <- read.csv("input-data-freeze/gisaid-MA.csv",
                 na.strings="",
                 colClasses=c(date="Date")) %>% # [JT] set col classes while reading
  filter(variant %in% c("Delta", "Omicron"), # [JT] only keep delta and omicron
         date >= "2021-12-01",
         date <= "2022-01-01")

# NE data (GISAID)
ne <- read.csv("input-data-freeze/gisaid-NE.csv",
               na.strings="",
               colClasses=c(date="Date")) %>% # [JT] set col classes while reading
  filter(variant %in% c("Delta", "Omicron"), # [JT] only keep delta and omicron
         date >= "2021-12-01",
         date <= "2022-01-01")

# add MA data to inst
inst <- ma %>%
  mutate(Institution="MA",
         Date=date,
         Variant=variant) %>%
  select(Institution, Date, Variant) %>%
  full_join(inst, by=c("Institution", "Date", "Variant"))

# add NE data to inst
inst <- ne %>%
  mutate(Institution="NE",
         Date=date,
         Variant=variant) %>%
  select(Institution, Date, Variant) %>%
  full_join(inst, by=c("Institution", "Date", "Variant"))
inst$bin<-ifelse(inst$Variant == "Omicron",1,0)

# format as daily percentages
prct <- inst %>%
        group_by(Institution, Date) %>%
        summarise(Omicron=sum(Variant=="Omicron"),
                  Delta=sum(Variant=="Delta"),
                  Cases=n(),
                  .groups="drop") %>%
        mutate(Fraction=Omicron/Cases,
               Percent=100*Fraction) %>%
        mutate(Institution=factor(Institution, levels=c("BU", "HU", "NEU", "MA", "NE")))

#add total case counts for MA and state from CDC
states <- read.csv("input-data-freeze/cdc-cases-state.csv",
                   na.strings="",
                   colClasses=c(date="Date")) %>% # [JT] set col classes while reading
  filter(date >= "2021-12-01",
         date <= "2022-01-15")
ma_case = subset(states, states$state == 'MA')

# add cases across NE states
ne_case <- states %>%
  group_by(date) %>%
  summarise(cases=sum(tot_cases))
rm(states)

## logistic fits ---------------------------------------------------------------
# figure 1C
p3 <- inst %>%
      ggplot(aes(Date, bin)) +
      stat_smooth(data=inst, aes(col=Institution), 
                  method="glm", se=FALSE, fullrange=TRUE, 
                  method.args=list(family=binomial)) +
      geom_point(data = prct, aes(Date, Fraction, fill=Institution), pch=21, size=3)  +
      scale_color_manual(values=c("BU"="#1b9e77","HU"="#d95f02", "NEU"="#7570b3", "MA"="black", "NE"="darkgray")) +
      scale_fill_manual(values=c("BU"="#1b9e77","HU"="#d95f02", "NEU"="#7570b3", "MA"="black", "NE"="darkgray")) +
      labs(x="date", 
           y="Omicron fraction",
           fill="institution",
           col="institution") +
      scale_x_continuous(breaks=as.Date(c("2021-11-30", "2021-12-05", "2021-12-10", 
                                          "2021-12-15", "2021-12-20", "2021-12-25", "2021-12-30", "2022-01-4")), 
                         limits=as.Date(c("2021-11-30", "2022-01-06")),
                         labels=c("Nov 30","Dec 5", "Dec 10","Dec 15","Dec 20",
                                  "Dec 25","Dec 30", "Jan 4")) +
      guides(col="none") # [JT] removed extra legend
p3
ggsave("outputs/logfit.png", units="cm", width=20, height=10)
#ggsave("outputs/logfit.svg", units="cm", width=15, height=10)

# logistic fits
fits <- unique(inst$Institution) %>%
  lapply(function(i) {
    inst %>%
      filter(Institution==i) %>%
      glm(bin ~ Date, data=., family="binomial")
  }) 
names(fits) <- unique(inst$Institution)

# null model to determine McFadden's R^2
null <- unique(inst$Institution) %>%
  lapply(function(i) {
    inst %>%
      filter(Institution==i) %>%
      glm(bin ~ 1, data=., family="binomial")
  }) 
names(null) <- unique(inst$Institution)

# summary statistics for each logistic model
# BU
summary(fits$BU)
confint(fits$BU, '(Intercept)', 0.95)
confint(fits$BU, 'Date', 0.95)
r2 <- 1-logLik(fits$BU)/logLik(null$BU)

# HU
summary(fits$HU)
confint(fits$HU, '(Intercept)', 0.95)
confint(fits$HU, 'Date', 0.95)
r2 <- 1-logLik(fits$HU)/logLik(null$HU)

# NEU
summary(fits$NEU)
confint(fits$NEU, '(Intercept)', 0.95)
confint(fits$NEU, 'Date', 0.95)
r2 <- 1-logLik(fits$NEU)/logLik(null$NEU)

# MA
summary(fits$MA)
confint(fits$MA, '(Intercept)', 0.95)
confint(fits$MA, 'Date', 0.95)
r2 <- 1-logLik(fits$MA)/logLik(null$MA)

# NE
summary(fits$NE)
confint(fits$NE, '(Intercept)', 0.95)
confint(fits$NE, 'Date', 0.95)
r2 <- 1-logLik(fits$NE)/logLik(null$NE)

rm(fits, null, pred, r2)

## inference--------------------------------------------------------------------
# functions
logit <- function(p){log(p/(1-p))}

# generate point estimates and 95% CIs for O10, O50, O90 using inverse pred
# adapted from https://www.talkstats.com/threads/inverse-prediction-from-binary-logistic-regression.52121/
omi_cis = data.frame()
# loop over each institution x Omicron fraction
for (ins in unique(inst$Institution)){
  for (p in c(0.1, 0.5, 0.9)){
    mod <- glm(bin ~ Date, data=subset(inst, inst$Institution == ins), family="binomial")
    est <- unname((logit(p) - coef(mod)[1])/coef(mod)[2])
    se <- deltamethod(~ (log(p/(1-p)) - x1)/x2, coef(mod), vcov(mod))
    ci <- est + c(-1, 1) * qt(0.975, sum(prct$Institution == ins) - 1) * se
    cis <- data.frame(inst = ins, prct = p, est = as.Date(est), e = est, stderr = se, 
                     ci_low = as.Date(ci[1]), ci_high = as.Date(ci[2]), n = sum(prct$Institution == ins))
    omi_cis = rbind(omi_cis, cis)}}
write.csv(omi_cis[c("inst", "prct", "est", "ci_low", "ci_high")], 
    'outputs/omicron_cis_inv_pred.csv', row.names = FALSE)
rm (ci, cis, est, ins, mod, p, se)


# generate p-vals and time differences for O10, O50, O90 by inst using inverse pred
ttest = data.frame(inst1 = combn(unique(inst$Institution), 2)[1,], 
        inst2 = combn(unique(inst$Institution), 2)[2,]) %>%
        slice(rep(1:n(), each = 3))
ttest$prct = rep(c(0.1, 0.5, 0.9), times = length(ttest$inst1)/3)

t = data.frame()
for (r in 1:length(ttest$inst1)){
  x = subset(omi_cis, inst == ttest$inst1[r] & prct == ttest$prct[r])
  y = subset(omi_cis, inst == ttest$inst2[r] & prct == ttest$prct[r])
  tt = tsum.test(x$e, x$stderr, x$n, y$e, y$stderr, y$n, conf.level = 0.95)
  t = rbind(t, data.frame(est = tt[["estimate"]][1] - tt[["estimate"]][2],
                          ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], 
                          pval = tt$p.value))}
ttest = cbind(ttest, t)
ttest$padj = p.adjust(ttest$pval, method = "BH")
write.csv(ttest, 'outputs/omicron_ttest.csv', row.names = FALSE)
rm(r, t, tt, x, y)


# generate CIs and time differences for deltaO(90-10) by inst using inverse pred
omi_delta = data.frame()
# loop over each institution
for (ins in unique(inst$Institution)){
  x = subset(omi_cis, inst == ins & prct == 0.9)
  y = subset(omi_cis, inst == ins & prct == 0.1)
  tt = tsum.test(x$e, x$stderr, x$n, y$e, y$stderr, y$n, conf.level = 0.95)
  omi_delta = rbind(omi_delta, data.frame(inst = ins, 
                            est = tt[["estimate"]][1] - tt[["estimate"]][2],
                            ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], 
                            pval = tt$p.value))}
rm(ins, tt, x, y)
write.csv(omi_delta, 'outputs/omicron_delta_ttest.csv', row.names = FALSE)

#clean up
rm(logit, omi_cis, omi_delta, ttest)

## observational figs-----------------------------------------------------------
# remove CDC and get total cases per day
tots <- prct %>%
        filter(!Institution %in% c("NE", "MA")) %>%
        group_by(Date) %>%
        summarise(Omicron=sum(Omicron),
                  Delta=sum(Delta),
                  Total=sum(Cases),
                  .groups="drop") %>%
        reshape2::melt(id.vars=c("Date", "Total"),
                       variable.name="Variant",
                       value.name="Cases") %>%
        mutate(Fraction=Cases/Total)

# college cases by variant
# [JT] fixed axis
p1 <- tots %>%
      filter((Date >= as.Date("2021-12-02")) & (Date <= as.Date("2021-12-21"))) %>%
      ggplot(aes(Date, Cases, fill=Variant)) +
      geom_col(col="black", position="stack") +
      scale_fill_manual(values=cols) +
      labs(x="date", 
           y="SARS-CoV-2-positive individuals",
           fill="variant") +
      scale_x_continuous(breaks=as.Date(c("2021-12-02", "2021-12-06", "2021-12-11", 
                                          "2021-12-16", "2021-12-21")), 
                         limits=c(as.Date("2021-12-01"), 
                                  as.Date("2021-12-22")), 
                         labels=c("Dec 02", "Dec 06", "Dec 11", 
                                  "Dec 16", "Dec 21"))
p1
ggsave("outputs/cases.png", units="cm", width=15, height=10)
#ggsave("outputs/cases.svg", units="cm", width=15, height=10)


# frequency
# [JT] fixed axis
p2 <- tots %>%
  mutate(Variant=factor(Variant, levels=c("Delta", "Omicron"))) %>%
  filter((Date >= as.Date("2021-12-02")) & (Date <= as.Date("2021-12-21"))) %>%
  ggplot(aes(Date, Fraction, fill=Variant)) +
  geom_col(col="black", position="stack") +
  scale_fill_manual(values=cols) +
  labs(x="date", 
       y="Omicron fraction",
       fill="variant")  +
  scale_x_continuous(breaks=as.Date(c("2021-12-02", "2021-12-06", "2021-12-11", 
                                      "2021-12-16", "2021-12-21")), 
                     limits=c(as.Date("2021-12-01"), 
                              as.Date("2021-12-22")), 
                     labels=c("Dec 02", "Dec 06", "Dec 11", 
                              "Dec 16", "Dec 21"))
p2
ggsave("outputs/frequency.png", units="cm", width=15, height=10)
#ggsave("outputs/frequency.svg", units="cm", width=15, height=10)


# combine college, MA, NE cases
college = data.frame(group = "colleges", cases = tots$Total, date = tots$Date,
                     roll = rollmean(tots$Total, 3, na.pad = TRUE))
mass = data.frame(group = "MA", cases = ma_case$tot_cases, date = ma_case$date,
                  roll = rollmean(ma_case$tot_cases, 7, na.pad = TRUE))
region = data.frame(group = "NE", cases = ne_case$cases, date = ne_case$date,
                    roll = rollmean(ne_case$cases, 7, na.pad = TRUE))
total = rbind(college, mass, region)
rm(college, mass, region)
total$merge <- ifelse(is.na(total$roll),total$cases,total$roll)
total$plot <-ifelse(total$group %in% c("MA", "NE"),total$merge/100, total$merge)

 
## plot total case counts over time
p4 = total %>%
    filter((date >= as.Date("2021-12-02")) & (date <= as.Date("2022-01-14"))) %>%
    ggplot(aes(date, plot, fill = group)) +
    geom_line(aes(date, plot, col = group), show.legend = F, size = 2) +
    geom_point(aes(date, plot), pch=21, size=3)  +
    scale_color_manual(values=casecols) +
    scale_fill_manual(values=casecols) +
    scale_y_continuous("Cases (colleges) ",
    sec.axis = sec_axis(~ . * 100, name = "Cases (MA, NE)")) +
    labs(x="date",
        y="case count",
        fill="Group")  +
    theme(legend.position = c(0.1, 0.8)) +
    scale_x_continuous(breaks=as.Date(c("2021-12-02", "2021-12-06", "2021-12-11",
                                        "2021-12-16", "2021-12-21", "2021-12-26",
                                        "2021-12-31", "2022-01-05", "2022-01-10", 
                                        "2022-01-15")),
                      limits=c(as.Date("2021-12-01"),
                                as.Date("2022-1-15")),
                      labels=c("Dec 02", "Dec 06", "Dec 11",
                                "Dec 16", "Dec 21", "Dec 26",
                               "Dec 31", "Jan 5", "Jan 10", "Jan 15"))
p4
ggsave("outputs/tot_cases.png", units="cm", width=15, height=10)
#ggsave("outputs/tot_cases.svg", units="cm", width=15, height=10)

# frequency
# [JT] fixed axis
p2 <- tots %>%
      mutate(Variant=factor(Variant, levels=c("Delta", "Omicron"))) %>%
      filter((Date >= as.Date("2021-12-02")) & (Date <= as.Date("2021-12-21"))) %>%
      ggplot(aes(Date, Fraction, fill=Variant)) +
      geom_col(col="black", position="stack") +
      scale_fill_manual(values=cols) +
      labs(x="date", 
           y="Omicron fraction",
           fill="variant")  +
      scale_x_continuous(breaks=as.Date(c("2021-12-02", "2021-12-06", "2021-12-11", 
                                          "2021-12-16", "2021-12-21")), 
                         limits=c(as.Date("2021-12-01"), 
                                  as.Date("2021-12-22")), 
                         labels=c("Dec 02", "Dec 06", "Dec 11", 
                                  "Dec 16", "Dec 21"))
p2
ggsave("outputs/frequency.png", units="cm", width=15, height=10)
#ggsave("outputs/frequency.svg", units="cm", width=15, height=10)

# concatenate
p <- cowplot::plot_grid(p1+theme(legend.position="none"), p2,
                        ncol=2, rel_widths=c(1, 1.2),
                        labels=c("B", "C"))
p <- cowplot::plot_grid(p4, p, nrow=2,
                        labels=c("A", ""))
p <- cowplot::plot_grid(p, p3, nrow=2,
                        labels=c("", "D"), rel_heights=c(2, 1))
p
ggsave("outputs/figure1.png", units="cm", width=25, height=35)
#ggsave("outputs/figure1.svg", units="cm", width=20, height=20)

# clean up
rm(p1, p2, p3, p)

## fin! ------------------------------------------------------------------------
sessionInfo()

