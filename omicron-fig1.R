#!/usr/bin/env Rscript

## version history -------------------------------------------------------------
# v1 by J.T. 6 Jan 2022, v2 by B.P. 7 Jan 2022, v3 by B.P. 9 Jan 2022
# updated input data, analyzed logistic fits
# [JT] v4 by J. T. 11 Jan 2022
# [JT] v5 more input and output formatting; code testing and commenting
# B.P. v6, updated CIs to include dates b/w 11-30 and 12-31
# B.P. v7, added GISAID data and MA curves
# B.P. v8, added MA data post-university dates
# B.P. v9, added MA and NE new cases per day
# B.P. v10, added inverse predictions
# B.P. v11, added cases per day panel
# B.P. v12, input and output formatting; code testing and commenting

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/Volumes/GoogleDrive/My Drive/Omicron/omicron_repo")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(zoo))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())
library('BSDA') #for t tests
library('msm') #for delta method
library('zoo') #for working with dates

# global variables
instcols <- c(IHEs="#ca0020", MA="black", NE='darkgrey')
instvarcols <- c("Delta (IHEs)" = "#984ea3", "Omicron (IHEs)" = "#ca0020", 
                 "Delta (MA)" = "black", "Omicron (MA)" = "darkgrey")

## inputs ----------------------------------------------------------------------
# institutions: every case
inst <- read.csv("input-data-freeze/combined-inst.csv",
                 na.strings="",
                 colClasses=c(Date="Date")) %>% # [JT] set col classes while reading
        filter(Variant %in% c("Delta", "Omicron"), # [JT] only keep delta and omicron
               Date >= "2021-12-02", # only keep dates w/ data for all IHEs
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
         division != "Massachusetts", # NE minus MA (e.g., CT, ME, NH, RI, VT)
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
  mutate(Institution="NE (without MA)",
         Date=date,
         Variant=variant) %>%
  select(Institution, Date, Variant) %>%
  full_join(inst, by=c("Institution", "Date", "Variant")) %>%
  mutate(bin = ifelse(Variant == "Omicron",1,0))
rm(ma, ne)

# format as daily percentages
prct <- inst %>%
        group_by(Institution, Date) %>%
        summarise(Omicron=sum(Variant=="Omicron"),
                  Delta=sum(Variant=="Delta"),
                  Cases=n(),
                  .groups="drop") %>%
        mutate(Fraction=Omicron/Cases,
               Percent=100*Fraction) %>%
        mutate(Institution=factor(Institution, levels=c("BU", "HU", "NU", "MA", "NE (without MA)")))

#daily case counts for MA (CDC)
ma <- read.csv("input-data-freeze/cdc-cases-state.csv",
                   na.strings="",
                   colClasses=c(date="Date")) %>% # [JT] set col classes while reading
  filter(date >= "2021-12-01",
         date <= "2022-01-15",
         state == "MA") %>%
  mutate(day = weekdays(date)) %>% # note weekday vs. weekend variation
  arrange(date)

#daily case counts for NE (CDC)
ne <- read.csv("input-data-freeze/cdc-cases-state.csv",
               na.strings="",
               colClasses=c(date="Date")) %>% # [JT] set col classes while reading
  filter(date >= "2021-12-01",
         date <= "2022-01-15") %>%
  mutate(day = weekdays(date)) %>% # note weekday vs. weekend variation
  arrange(date) %>%
  group_by(date) %>%
  summarise(cases=sum(tot_cases))

# population sizes
ma_pop = 7033469 #from Wikipedia
ne_pop = 15116205 #from Wikipedia
BU = 43904 #from BU
HU = 38434 #from HUCL
NU = 30602 #from NU
col_pop = NU + BU + HU #total college testing program size
rm(BU, HU, NU)

## logistic fits ---------------------------------------------------------------
# figure 1B
p2 <- inst %>%
      ggplot(aes(Date, bin)) +
      stat_smooth(data=inst, aes(col=Institution), 
                  method="glm", se=FALSE, fullrange=TRUE, 
                  method.args=list(family=binomial)) +
      geom_point(data = prct, aes(Date, Fraction, fill=Institution), pch=21, size=3)  +
      scale_color_manual(values=c("BU"="#1b9e77","HU"="#d95f02", "NU"="#7570b3", "MA"="black", "NE (without MA)"="darkgray")) +
      scale_fill_manual(values=c("BU"="#1b9e77","HU"="#d95f02", "NU"="#7570b3", "MA"="black", "NE (without MA)"="darkgray")) +
      labs(x="Date", 
           y="Omicron fraction",
           fill="Institution",
           col="Institution") +
      theme(legend.text=element_text(size=12), legend.title=element_text(size=14)) +
      scale_x_continuous(breaks=as.Date(c("2021-12-01", "2021-12-06", "2021-12-12", 
                                          "2021-12-18", "2021-12-24", "2022-01-01")),
                         labels=c("Dec 1","Dec 6", "Dec 12","Dec 18","Dec 24", "Jan 1")) + guides(col="none") # [JT] removed extra legend
p2
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

# NU
summary(fits$NU)
confint(fits$NU, '(Intercept)', 0.95)
confint(fits$NU, 'Date', 0.95)
r2 <- 1-logLik(fits$NU)/logLik(null$NU)

# MA
summary(fits$MA)
confint(fits$MA, '(Intercept)', 0.95)
confint(fits$MA, 'Date', 0.95)
r2 <- 1-logLik(fits$MA)/logLik(null$MA)

# NE
summary(fits$'NE (without MA)')
confint(fits$'NE (without MA)', '(Intercept)', 0.95)
confint(fits$'NE (without MA)', 'Date', 0.95)
r2 <- 1-logLik(fits$'NE (without MA)')/logLik(null$'NE (without MA)')

# get MA Omicron fraction per day to estimate variant-spec case counts
ma_prct <- data.frame(Date = seq(as.Date('2021-12-01'), as.Date('2022-01-15'), by = 1))
ma_prct$prct <- predict(fits$MA, ma_prct, se.fit = F, type = 'response')
rm(fits, null, pred, r2)

## inference--------------------------------------------------------------------
# functions
logit <- function(p){log(p/(1-p))}

# point estimates and 95% CIs for O10, O50, O90 using delta method
# adapted from www.talkstats.com/threads/inverse-prediction-from-binary-logistic-regression.52121/
omi_cis = data.frame()
# loop over each institution x Omicron fraction
for (ins in unique(inst$Institution[inst$Institution != "NE (without MA)"])){
  for (p in c(0.1, 0.5, 0.9)){
    mod <- glm(bin ~ Date, data=subset(inst, inst$Institution == ins), family="binomial")
    est <- unname((logit(p) - coef(mod)[1])/coef(mod)[2])
    se <- deltamethod(~ (log(p/(1-p)) - x1)/x2, coef(mod), vcov(mod)) #delta method to estimate standard error
    ci <- est + c(-1, 1) * qt(0.975, sum(prct$Institution == ins) - 1) * se #95% CI via student's t distribution
    cis <- data.frame(inst = ins, prct = p, est = as.Date(est), e = est, stderr = se, 
                     ci_low = as.Date(ci[1]), ci_high = as.Date(ci[2]), n = sum(prct$Institution == ins))
    omi_cis = rbind(omi_cis, cis)}}
write.csv(omi_cis[c("inst", "prct", "est", "ci_low", "ci_high")], 
    'outputs/omicron_cis_inv_pred.csv', row.names = FALSE)
rm (ci, cis, est, ins, mod, p, se)


# generate delta t and pval for O10, O50, O90 by inst
ttest = data.frame(inst1 = combn(unique(inst$Institution[inst$Institution != "NE (without MA)"]), 2)[1,], 
        inst2 = combn(unique(inst$Institution[inst$Institution != "NE (without MA)"]), 2)[2,]) %>%
        slice(rep(1:n(), each = 3))
ttest$prct = rep(c(0.1, 0.5, 0.9), times = length(ttest$inst1)/3)

t = data.frame()
for (r in 1:length(ttest$inst1)){
  x = subset(omi_cis, inst == ttest$inst1[r] & prct == ttest$prct[r])
  y = subset(omi_cis, inst == ttest$inst2[r] & prct == ttest$prct[r])
  tt = tsum.test(x$e, x$stderr, x$n, y$e, y$stderr, y$n, conf.level = 0.95) #two-sample student's t test
  t = rbind(t, data.frame(est = tt[["estimate"]][1] - tt[["estimate"]][2],
                          ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], 
                          pval = tt$p.value))}
ttest = cbind(ttest, t)
ttest$padj = p.adjust(ttest$pval, method = "BH")
write.csv(ttest, 'outputs/omicron_ttest.csv', row.names = FALSE)
rm(r, t, tt, x, y)


# generate point estimates and CIs for deltaO(90-10) by inst
omi_delta = data.frame()
# loop over each institution
for (ins in unique(inst$Institution[inst$Institution != "NE (without MA)"])){
  x = subset(omi_cis, inst == ins & prct == 0.9)
  y = subset(omi_cis, inst == ins & prct == 0.1)
  tt = tsum.test(x$e, x$stderr, x$n, y$e, y$stderr, y$n, conf.level = 0.95) #95% CI via student's t distribution
  omi_delta = rbind(omi_delta, data.frame(inst = ins, 
                            est = tt[["estimate"]][1] - tt[["estimate"]][2],
                            ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], 
                            pval = tt$p.value))}
rm(ins, tt, x, y)
write.csv(omi_delta, 'outputs/omicron_delta_ttest.csv', row.names = FALSE)

#clean up
rm(logit, omi_cis, omi_delta, ttest)

## daily cases------------------------------------------------------------------
# total college cases per day
tots <- prct %>%
  filter(!Institution %in% c("NE (without MA)", "MA")) %>%
  group_by(Date) %>%
  summarise(Omicron=sum(Omicron),
            Delta=sum(Delta),
            Total=sum(Cases),
            .groups="drop") %>%
  reshape2::melt(id.vars=c("Date", "Total"),
                 variable.name="Variant",
                 value.name="Cases") %>%
  mutate(Group="IHEs")

# add MA data to tots
ma$Omicron <- ma$tot_cases*ma_prct$prct
ma$Delta <- ma$tot_cases*(1-ma_prct$prct)

tots <- ma %>%
  group_by(date) %>%
  summarise(Omicron=sum(Omicron),
            Delta=sum(Delta),
            Total=sum(tot_cases),
            .groups="drop") %>%
  reshape2::melt(id.vars=c("date", "Total"),
                 variable.name="Variant",
                 value.name="Cases") %>%
  mutate(Group="MA", Date = date) %>%
  select(Cases, Date, Group, Total, Variant) %>%
  full_join(tots, by=c("Cases", "Date", "Group", "Total", "Variant"))

# add NE data to tots
tots <- ne %>%
  mutate(Group = "NE", Date = date, Total = cases) %>%
  select(Group, Date, Total) %>%
  full_join(tots, by = c("Total", "Date", "Group"))
rm(ma, ma_prct, ne)

# reported case counts per 100K
tots$Rep = 100000*c(tots$Total[tots$Group == "NE"]/ne_pop, tots$Total[tots$Group == "MA"]/ma_pop,
                    tots$Total[tots$Group == "IHEs"]/col_pop)

# rolling case counts (due to weekly variation) per 100K
tots$Roll = 100000*c(rollmean(tots$Total[tots$Group == "NE"]/ne_pop, 7, na.pad = T), 
  rollmean(tots$Total[tots$Group == "MA"]/ma_pop, 7, na.pad = T),
  rollmean(tots$Total[tots$Group == "IHEs"]/col_pop, 3, na.pad = T))

## case plots-------------------------------------------------------------------

## plot case counts per 100K over time
ps = tots %>%
    filter((Date >= as.Date("2021-12-01")) & (Date <= as.Date("2022-01-15"))) %>%
    ggplot(aes(Date, Rep, fill = Group)) +
    geom_point(aes(Date, Rep), pch=21, size=3)  +
    scale_color_manual(values=instcols) +
    scale_fill_manual(values=instcols) +
    scale_y_continuous("Reported cases per 100K") + 
    theme(legend.text=element_text(size=14), legend.title=element_text(size=14)) +
    labs(x="Date",
        fill="Institution")  +
    scale_x_continuous(breaks=as.Date(c("2021-12-01", "2021-12-10", "2021-12-20",
                                        "2021-12-30", "2022-01-7", "2022-01-15")),
                      limits=c(as.Date("2021-12-01"),
                                as.Date("2022-1-15")),
                      labels=c("Dec 01", "Dec 10", "Dec 20",
                                "Dec 30", "Jan 07", "Jan 15"))
ps
ggsave("outputs/cases_per_100k.png", units="cm", width=20, height=10)
#ggsave("outputs/cases_per_100k.svg", units="cm", width=15, height=10)
rm(instcols, ps)

# subset variant-specific college and MA cases
tots = tots %>%
  filter(Group != "NE") %>%
  mutate(Label = paste(Variant, ' (' ,Group, ')', sep = ""))

# rolling variant-specific case counts per 100K
tots$Roll = 100000*c(rollmean(tots$Cases[tots$Label == "Omicron (MA)"]/ma_pop, 7, na.pad = T),
                     rollmean(tots$Cases[tots$Label == "Delta (MA)"]/ma_pop, 7, na.pad = T),
                     rollmean(tots$Cases[tots$Label == "Omicron (IHEs)"]/col_pop, 3, na.pad = T),
                     rollmean(tots$Cases[tots$Label == "Delta (IHEs)"]/col_pop, 3, na.pad = T))
rm(col_pop, ma_pop, ne_pop)

## plot variant-specific case counts per 100K
p1 = tots %>%
  filter((Date >= as.Date("2021-12-03")) & (Date <= as.Date("2022-01-12"))) %>%
  ggplot(aes(Date, Roll, fill = Label)) +
  geom_smooth(aes(Date, Roll, col = Label), se = F, show.legend = F, span = 0.75) +
  geom_point(aes(Date, Roll), pch = 21, size=3)  +
  scale_color_manual(values=instvarcols) +
  scale_fill_manual(values=instvarcols) +
  scale_y_continuous("Cases per 100K",limits = c(0, 400)) + 
  theme(legend.text=element_text(size=12), legend.title=element_text(size=14)) +
  labs(x="Date",
       fill="Institution")  +
  scale_x_continuous(breaks=as.Date(c("2021-12-03", "2021-12-13", "2021-12-23",
                                      "2022-01-02", "2022-01-12")),
                     limits=c(as.Date("2021-12-03"),
                              as.Date("2022-01-12")),
                     labels=c("Dec 03", "Dec 13", "Dec 23", "Jan 01", "Jan 12"))
p1
ggsave("outputs/cases_per_100k.png", units="cm", width=15, height=10)
#ggsave("outputs/cases_per_100k.svg", units="cm", width=15, height=10)

# concatenate panels
p <- cowplot::plot_grid(p1, p2, nrow=2, labels=c("A", "B"))
p
ggsave("outputs/figure1_new.png", units="cm", width=25, height=25)
#ggsave("outputs/figure1_new.svg", units="cm", width=25, height=25)

# clean up
rm(p, p1, p2)

## fin! ------------------------------------------------------------------------
sessionInfo()

