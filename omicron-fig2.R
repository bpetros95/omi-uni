#!/usr/bin/env Rscript

## version history -------------------------------------------------------------
# v1 by B.P. 9 Jan 2022
# v2 by J.T. 13 Jan 2022
# v3 by B.P. 14 Jan 2022 [added case counts over time by affiliation]
# v4 by B.P. 19 Jan 2022 [added bootstrapping]
# v5 by B.P. 20 Jan 2022 [updated log reg, added inverse predictions]

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/Volumes/GoogleDrive/My Drive/Omicron/omicron_repo")
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())
library('BSDA') #for t tests
library('msm') #for delta method
library('zoo') #for working with dates

# global variables
cols <- c(Delta="#bdbdbd", Omicron="#ca0020")
affl <- c(Student="#1b9e77", Employee="#99d8c9", MA="black")

## inputs ----------------------------------------------------------------------
# BU student vs. employee data
bu <- read.csv("input-data-freeze/student-v-employee.csv",
               na.strings="",
               colClasses=c(Date="Date")) %>%
      filter(Variant %in% c("Delta", "Omicron"),
             Date >= "2021-12-02",
             Date <= "2021-12-21") %>%
      select(-"Ct_RP")

# MA data (GISAID)
ma <- read.csv("input-data-freeze/gisaid-MA.csv",
               na.strings="",
               colClasses=c(date="Date")) %>% # [JT] set col classes while reading
  filter(variant %in% c("Delta", "Omicron"), # [JT] only keep delta and omicron
         date >= "2021-12-01",
         date <= "2022-01-01")

# add MA data to BU
bu <- ma %>%
  mutate(Affiliation="MA",
         Date=date,
         Variant=variant) %>%
  select(Affiliation, Date, Variant) %>%
  full_join(bu, by=c("Affiliation", "Date", "Variant"))
bu$bin<-ifelse(bu$Variant == "Omicron",1,0)
rm(ma)

# format as daily percentages
prct <- bu %>%
        group_by(Affiliation, Date) %>%
        summarise(Omicron=sum(Variant=="Omicron"),
                  Delta=sum(Variant=="Delta"),
                  Cases=n(),
                  .groups="drop") %>%
        mutate(Fraction=Omicron/Cases,
               Percent=100*Fraction) %>%
        mutate(Affiliation=factor(Affiliation, 
               levels=c("Student", "Employee", "MA")))

## observational figs  ---------------------------------------------------------

# remove MA, NE and get BU total cases per day
tots <- prct %>%
  filter(Affiliation != "MA") %>%
  group_by(Date) %>%
  summarise(Omicron=sum(Omicron),
            Delta=sum(Delta),
            Total=sum(Cases),
            .groups="drop") %>%
  reshape2::melt(id.vars=c("Date", "Total"),
                 variable.name="Variant",
                 value.name="Cases") %>%
  mutate(Fraction=Cases/Total)

# BU cases total
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
ggsave("outputs/cases_bu.png", units="cm", width=15, height=10)
#ggsave("outputs/cases_bu.svg", units="cm", width=15, height=10)


# remove MA and get BU cases per day by affiliation
affil <- prct %>%
  filter(Affiliation != "MA") %>%
  group_by(Date) %>%
  summarise(Student=sum(Cases[Affiliation == 'Student']),
            Employee=sum(Cases[Affiliation == 'Employee']),
            Total=sum(Cases),
            .groups="drop") %>%
  reshape2::melt(id.vars=c("Date", "Total"),
                 variable.name="Affiliation",
                 value.name="Cases") %>%
  mutate(Fraction=Cases/Total ) %>%
  mutate(Affiliation=factor(Affiliation, 
                            levels=c("Employee", "Student")))
  
# BU cases by affiliation
p2 <- affil %>%
  filter((Date >= as.Date("2021-12-02")) & (Date <= as.Date("2021-12-21"))) %>%
  ggplot(aes(Date, Cases, fill=Affiliation)) +
  geom_col(col="black", position="stack") +
  scale_fill_manual(values=affl[1:2]) +
  labs(x="date", 
       y="SARS-CoV-2-positive individuals",
       fill="affiliation") +
  scale_x_continuous(breaks=as.Date(c("2021-12-02", "2021-12-06", "2021-12-11", 
                                      "2021-12-16", "2021-12-21")), 
                     limits=c(as.Date("2021-12-01"), 
                              as.Date("2021-12-22")), 
                     labels=c("Dec 02", "Dec 06", "Dec 11", 
                              "Dec 16", "Dec 21"))
p2
ggsave("outputs/affil_bu.png", units="cm", width=15, height=10)
#ggsave("outputs/affil_bu.svg", units="cm", width=15, height=10)

fish = fisher.test(cbind(c(sum(bu$Affiliation == 'Student' & bu$Variant == 'Delta'), 
                           sum(bu$Affiliation == 'Student' & bu$Variant == 'Omicron')),
                         c(sum(bu$Affiliation == 'Employee' & bu$Variant == 'Delta'), 
                           sum(bu$Affiliation == 'Employee' & bu$Variant == 'Omicron'))))

rm(affil, fish, tots)

## logistic growth -------------------------------------------------------------
# students vs. employees vs. CDC
p3 <- bu %>%
  ggplot(aes(Date, bin)) +
  stat_smooth(data=bu, aes(col=Affiliation), 
              method="glm", se=FALSE, fullrange=TRUE, 
              method.args=list(family=binomial)) +
  geom_point(data = prct, aes(Date, Fraction, fill=Affiliation), pch=21, size=3)  +
  scale_color_manual(values=affl) +
  scale_fill_manual(values=affl) +
  labs(x="date", 
       y="Omicron fraction",
       fill="Affiliation") +
  scale_x_continuous(breaks=as.Date(c("2021-11-30", "2021-12-05", "2021-12-10", 
                                      "2021-12-15", "2021-12-20", "2021-12-25", "2021-12-30", "2022-01-4")), 
                     limits=as.Date(c("2021-11-30", "2022-01-06")),
                     labels=c("Nov 30","Dec 5", "Dec 10","Dec 15","Dec 20",
                              "Dec 25","Dec 30", "Jan 4")) +
  guides(col="none") # [JT] removed extra legend
p3
ggsave("outputs/logfit-bu.png", units="cm", width=20, height=10)
#ggsave("outputs/logfit-bu.svg", units="cm", width=15, height=10)

# concatenate
p <- cowplot::plot_grid(p1, p2, 
                        ncol=2,
                        labels=c("A", "B"))
p <- cowplot::plot_grid(p, p3, 
                        nrow=2,
                        labels=c("", "C"))
p
ggsave("outputs/figure2.png", units="cm", width=20, height=20)
#ggsave("outputs/figure2.svg", units="cm", width=20, height=20)
rm(p, p1, p2, p3)

# logistic fits
fits <- unique(bu$Affiliation) %>%
  lapply(function(i) {
    bu %>%
      filter(Affiliation==i) %>%
      glm(bin ~ Date, data=., family="binomial")
  }) 
names(fits) <- unique(bu$Affiliation)

# null model to determine McFadden's R^2
null <- unique(bu$Affiliation) %>%
  lapply(function(i) {
    bu %>%
      filter(Affiliation==i) %>%
      glm(bin ~ 1, data=., family="binomial")
  }) 
names(null) <- unique(bu$Affiliation)


# summary statistics for each logistic model
# Student
summary(fits$Student)
confint(fits$Student, '(Intercept)', 0.95)
confint(fits$Student, 'Date', 0.95)
r2 <- 1-logLik(fits$Student)/logLik(null$Student)

# Employee
summary(fits$Employee)
confint(fits$Employee, '(Intercept)', 0.95)
confint(fits$Employee, 'Date', 0.95)
r2 <- 1-logLik(fits$Employee)/logLik(null$Employee)

# MA
summary(fits$MA)
confint(fits$MA, '(Intercept)', 0.95)
confint(fits$MA, 'Date', 0.95)
r2 <- 1-logLik(fits$MA)/logLik(null$MA)

# clean up
rm(fits, null, r2)

## inference--------------------------------------------------------------------
# functions
logit <- function(p){log(p/(1-p))}

# point estimates and 95% CIs for O10, O50, O90 using delta method
# adapted from www.talkstats.com/threads/inverse-prediction-from-binary-logistic-regression.52121/
omi_cis = data.frame()
# loop over each affiliation x Omicron fraction
for (ins in unique(bu$Affiliation)){
  for (p in c(0.1, 0.5, 0.9)){
    mod <- glm(bin ~ Date, data=subset(bu, bu$Affiliation == ins), family="binomial")
    est <- unname((logit(p) - coef(mod)[1])/coef(mod)[2])
    se <- deltamethod(~ (log(p/(1-p)) - x1)/x2, coef(mod), vcov(mod)) #delta method to estimate standard error
    ci <- est + c(-1, 1) * qt(0.975, sum(prct$Affiliation == ins) - 1) * se #95% CI via student's t distribution
    cis <- data.frame(affil = ins, prct = p, est = as.Date(est), e = est, stderr = se, 
                      ci_low = as.Date(ci[1]), ci_high = as.Date(ci[2]), n = sum(prct$Affiliation == ins))
    omi_cis = rbind(omi_cis, cis)}}
write.csv(omi_cis[c("affil", "prct", "est", "ci_low", "ci_high")], 
          'outputs/affil_cis_inv_pred.csv', row.names = FALSE)
rm (ci, cis, est, ins, mod, p, se)


# generate delta t and pval for O10, O50, O90 by inst
ttest = data.frame(affil1 = combn(unique(bu$Affiliation), 2)[1,], 
                   affil2 = combn(unique(bu$Affiliation), 2)[2,]) %>%
  slice(rep(1:n(), each = 3))
ttest$prct = rep(c(0.1, 0.5, 0.9), times = length(ttest$affil1)/3)

t = data.frame()
for (r in 1:length(ttest$affil1)){
  x = subset(omi_cis, affil == ttest$affil1[r] & prct == ttest$prct[r])
  y = subset(omi_cis, affil == ttest$affil2[r] & prct == ttest$prct[r])
  tt = tsum.test(x$e, x$stderr, x$n, y$e, y$stderr, y$n, conf.level = 0.95) #two-sample student's t test
  t = rbind(t, data.frame(est = tt[["estimate"]][1] - tt[["estimate"]][2],
                          ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], 
                          pval = tt$p.value))}
ttest = cbind(ttest, t)
ttest$padj = p.adjust(ttest$pval, method = "BH")
write.csv(ttest, 'outputs/affil_ttest.csv', row.names = FALSE)
rm(r, t, tt, x, y)


# generate point estimates and CIs for deltaO(90-10) by inst
omi_delta = data.frame()
# loop over each affiliation
for (ins in unique(bu$Affiliation)){
  x = subset(omi_cis, affil == ins & prct == 0.9)
  y = subset(omi_cis, affil == ins & prct == 0.1)
  tt = tsum.test(x$e, x$stderr, x$n, y$e, y$stderr, y$n, conf.level = 0.95) #95% CI via student's t distribution
  omi_delta = rbind(omi_delta, data.frame(affil = ins, 
                                          est = tt[["estimate"]][1] - tt[["estimate"]][2],
                                          ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], 
                                          pval = tt$p.value))}
rm(ins, tt, x, y)
write.csv(omi_delta, 'outputs/affil_delta_ttest.csv', row.names = FALSE)

#clean up
rm(affl, logit, omi_cis, omi_delta, ttest)

## Ct comparison ---------------------------------------------------------------
# compare across-affiliation Cts
pval = wilcox.test(bu$Ct_N1[bu$Affiliation == 'Student'], bu$Ct_N1[bu$Affiliation == 'Employee'])$p.value
pval = wilcox.test(bu$Ct_N2[bu$Affiliation == 'Student'], bu$Ct_N2[bu$Affiliation == 'Employee'])$p.value

# compare within-affiliation variant vs. primer
ct <- bu %>%
      select(Affiliation, Variant, starts_with("Ct")) %>%
      reshape2::melt(id.vars=c("Affiliation", "Variant"),
                     value.name="Ct",
                     variable.name="Primer") %>%
      mutate(Primer=str_extract(Primer, "(?<=Ct_)[A-z0-9]+")) %>%
      # remove primer dropouts
      # remove crazy Cts
      filter(!is.na(Ct),
             Ct >= 5,
             Ct <= 40) %>%
      # groups are primer, affiliation, and variant
      mutate(Group=paste0(Affiliation, ".", Primer, "-", Variant),
             Label=paste(Affiliation, Primer))

# calculate FDR-adjusted p-value
pval <- pairwise.wilcox.test(ct$Ct, ct$Group, p.adjust.method="BH")

# format adjusted p-values for plotting
pval <- pval$p.value %>%
        as.data.frame() %>%
        rownames_to_column("Rows") %>%
        reshape2::melt(id.vars="Rows",
                       variable.name="Columns",
                       value.name="p.adj") %>%
        mutate(Affiliation=str_extract(Rows, "^[A-z ]+(?=\\.)"),
               Affiliation2=str_extract(Columns, "^[A-z ]+(?=\\.)"),
               Primer=str_extract(Rows, "(?<=\\.)[A-z0-9]+(?=\\-)"),
               Primer2=str_extract(Columns, "(?<=\\.)[A-z0-9]+(?=\\-)"),
               Variant1=str_extract(Rows, "[A-z]+$"),
               Variant2=str_extract(Columns, "[A-z]+$")) %>%
        filter(Affiliation==Affiliation2,
               Primer==Primer2,
               !is.na(p.adj)) %>%
        select(Affiliation, Primer, Variant1, Variant2, p.adj) %>%
        mutate(p.adj=scales::pvalue(p.adj, accuracy=0.0001),
               Label=paste(Affiliation, Primer))


plot.ct <- function(label, d=ct, p=pval, ylim=c(0, 45), pheight=42) {
  # format p-value
  p <- p %>%
       filter(Label==label) %>%
       rename(group1=Variant1,
              group2=Variant2) %>%
       select(group1, group2, p.adj)
  # plot it
  pt <- d %>%
        filter(Label==label) %>%
        ggplot(aes(Variant, Ct)) +
        geom_boxplot(aes(fill=Variant), col="black") +
        scale_fill_manual(values=cols) +
        ggpubr::stat_pvalue_manual(p, y.position=pheight) +
        facet_wrap(~Label) +
        labs(x="variant of concern",
             y=paste(str_extract(label, "[A-z0-9]+$"), "cycle threshold")) +
        ylim(ylim[1], ylim[2]) +
        theme(legend.position="none",
              axis.title.x=element_blank())
  # save it
  f <- tolower(label) %>%
       str_replace(" ", "-") %>%
       paste0("outputs/ct-", .)
  ggsave(paste0(f, ".png"), units="cm", width=10, height=10)
  #ggsave(paste0(f, ".svg"), units="cm", width=10, height=10)
  # return it
  return(pt)
}

# BU N1, BU N2, HU N1, NEU N2
p1 <- plot.ct("Student N1")
p2 <- plot.ct("Employee N1")
p3 <- plot.ct("Student N2")
p4 <- plot.ct("Employee N2")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels="AUTO")
ggsave("outputs/supp_figure3.png", width=12, height=12, units="cm")
#ggsave("outputs/supp_figure3.svg", width=12, height=12, units="cm")

# table
ct %>%
  group_by(Affiliation, Variant, Primer) %>%
  summarise(Median=median(Ct, na.rm = TRUE),
            Mean=mean(Ct, na.rm = TRUE),
            .groups="drop") %>%
  arrange(Affiliation, Primer) %>%
  write.csv("outputs/ct-breakdown-bu.csv", row.names=FALSE)

rm(p1, p2, p3, p4, pval, plot.ct)

## fin! ------------------------------------------------------------------------
sessionInfo()
