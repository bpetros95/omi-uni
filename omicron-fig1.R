#!/usr/bin/env Rscript

## version history -------------------------------------------------------------
# v1 by J.T. 6 Jan 2022, v2 by B.P. 7 Jan 2022, fig1 by B.P. 9 Jan 2022
# updated input data, analyzed logistic fits
# [JT] v3 by J. T. 11 Jan 2022
# [JT] more input and output formatting; code testing and commenting
# B.P. v4, updated CIs to include dates b/w 11-30 and 12-31
# B.P. v5, added GISAID data and MA curves
# B.P. v6, added MA data post-university dates

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/Volumes/GoogleDrive/My Drive/Omicron/omicron_repo")
suppressPackageStartupMessages(library(tidyverse)) # [JT] load libs in script
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# global variables
cols <- c(Delta="#bdbdbd", Omicron="#ca0020")

## inputs ----------------------------------------------------------------------
# institutions: every case
inst <- read.csv("input-data-freeze/combined-inst.csv",
                 na.strings="",
                 colClasses=c(Date="Date")) %>% # [JT] set col classes while reading
        filter(Variant %in% c("Delta", "Omicron"), # [JT] only keep delta and omicron
               Date >= "2021-12-02",
               Date <= "2021-12-21")

# GISAID: MA data
ma <- read.csv("input-data-freeze/gisaid-MA.csv",
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

# format as daily percentages
prct <- inst %>%
        group_by(Institution, Date) %>%
        summarise(Omicron=sum(Variant=="Omicron"),
                  Delta=sum(Variant=="Delta"),
                  Cases=n(),
                  .groups="drop") %>%
        mutate(Fraction=Omicron/Cases,
               Percent=100*Fraction)

# CDC: weekly omicron percent
# [JT] input and formatting
cdc <- read.csv("input-data-freeze/cdc-northeast.csv",
                na.strings="",
                colClasses=c(Week.mid="Date")) %>%
       mutate(CI95.low=CI95.low/100,
              CI95.high=CI95.high/100)

# add CDC percents to prct
prct <- cdc %>%
        mutate(Institution="NE",
               Date=Week.mid,
               Fraction=Percent/100) %>%
        select(Institution, Date, Fraction, Percent, CI95.low, CI95.high) %>%
        full_join(prct, by=c("Institution", "Date", "Fraction", "Percent")) %>%
        mutate(Institution=factor(Institution, levels=c("BU", "HU", "NEU", "MA", "NE")))

## logistic fits ---------------------------------------------------------------
# figure 1C
p3 <- prct %>%
      ggplot(aes(Date, Fraction, ymin = CI95.low, ymax = CI95.high)) +
      stat_smooth(data=subset(prct, Institution!="NE"), aes(col=Institution), 
                  method="glm", se=FALSE, fullrange=TRUE, 
                  method.args=list(family=binomial)) +
      geom_errorbar(data=subset(prct, Institution=="NE"), size=1, width = 1, color = "darkgray") +
      geom_point(aes(fill=Institution), pch=21, size=3) +
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
fits <- unique(subset(prct, Institution!="NE")$Institution) %>%
        lapply(function(i) {
          prct %>%
            filter(Institution==i) %>%
            glm(Fraction ~ Date, data=., family="binomial")
        }) 
names(fits) <- unique(subset(prct, Institution!="NE")$Institution)

# null model to determine McFadden's R^2
null <- unique(subset(prct, Institution!="NE")$Institution) %>%
        lapply(function(i) {
          prct %>%
            filter(Institution==i) %>%
            glm(Fraction ~ 1, data=., family="binomial")
        }) 
names(null) <- unique(subset(prct, Institution!="NE")$Institution)

# summary statistics for each logistic model
# BU
confint(fits$BU, '(Intercept)', 0.95)
confint(fits$BU, 'Date', 0.95)
r2 <- 1-logLik(fits$BU)/logLik(null$BU)
pred <- data.frame(Date=seq(as.Date("2021-11-30"), as.Date("2022-01-05"), by="days"))
pred <- ciTools::add_ci(pred, fits$BU, type="boot", names=c("lwr", "upr"), 
                        alpha=0.05, nSims=500) %>%
        mutate(type="bootstrap")
pred$Institution <- "BU"
ci <- pred

# HU
confint(fits$HU, '(Intercept)', 0.95)
confint(fits$HU, 'Date', 0.95)
r2 <- 1-logLik(fits$HU)/logLik(null$HU)
pred <- data.frame(Date=seq(as.Date("2021-11-30"), as.Date("2022-01-05"), by="days"))
pred <- ciTools::add_ci(pred, fits$HU, type="boot", names=c("lwr", "upr"), 
                        alpha=0.05, nSims=500) %>%
        mutate(type="bootstrap")
pred$Institution <- "HU"
ci <- rbind(ci, pred)

# NEU
confint(fits$NEU, '(Intercept)', 0.95)
confint(fits$NEU, 'Date', 0.95)
r2 <- 1-logLik(fits$NEU)/logLik(null$NEU)
pred <- data.frame(Date=seq(as.Date("2021-11-30"), as.Date("2022-01-05"), by="days"))
pred <- ciTools::add_ci(pred, fits$NEU, type="boot", names=c("lwr", "upr"), 
                        alpha=0.05, nSims=500) %>%
        mutate(type="bootstrap")
pred$Institution <- "NEU"
ci <- rbind(ci, pred)

# MA
confint(fits$MA, '(Intercept)', 0.95)
confint(fits$MA, 'Date', 0.95)
r2 <- 1-logLik(fits$MA)/logLik(null$MA)
pred <- data.frame(Date=seq(as.Date("2021-11-30"), as.Date("2022-01-05"), by="days"))
pred <- ciTools::add_ci(pred, fits$MA, type="boot", names=c("lwr", "upr"), 
                        alpha=0.05, nSims=500) %>%
  mutate(type="bootstrap")
pred$Institution <- "MA"
ci <- rbind(ci, pred)

# write 95% CIs to .csv
ci <- ci %>%
      rename("Point_est"="pred", "CI95_upper"="upr", "CI95_lower"="lwr")
ci <- ci[,!(names(ci) %in% "type")]
write.csv(ci, "outputs/bootstrap_ci_logfit.csv", row.names=FALSE)

# clean up
rm(fits, null, pred, ci)

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

# cases total
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
ggsave("outputs/frequency.svg", units="cm", width=15, height=10)

# concatenate
p <- cowplot::plot_grid(p1+theme(legend.position="none"), p2, 
                        ncol=2, rel_widths=c(1, 1.2),
                        labels=c("A", "B"))
p <- cowplot::plot_grid(p, p3, 
                        nrow=2,
                        labels=c("", "C"))
p
ggsave("outputs/figure1.png", units="cm", width=20, height=20)
#ggsave("outputs/figure1.svg", units="cm", width=20, height=20)

# clean up
rm(p1, p2, p3, p, r2)

## fin! ------------------------------------------------------------------------
sessionInfo()

