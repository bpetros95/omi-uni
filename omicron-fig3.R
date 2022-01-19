#!/usr/bin/env Rscript

## version history -------------------------------------------------------------
# v1 by B.P. 9 Jan 2022
# v2 by J.T. 13 Jan 2022
# v3 by B.P. 14 Jan 2022 [added case counts over time by affiliation]

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/Volumes/GoogleDrive/My Drive/Omicron/omicron_repo")
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# global variables
cols <- c(Delta="#bdbdbd", Omicron="#ca0020")
affl <- c(Student="#1b9e77", Employee="#99d8c9", MA="black", NE="darkgrey")

## inputs ----------------------------------------------------------------------
# BU student vs. employee data
bu <- read.csv("input-data-freeze/student-v-employee.csv",
               na.strings="",
               colClasses=c(Date="Date")) %>%
      filter(Variant %in% c("Delta", "Omicron"),
             Date >= "2021-12-02",
             Date <= "2021-12-21")

# format as daily percentages
prct <- bu %>%
        group_by(Affiliation, Date) %>%
        summarise(Omicron=sum(Variant=="Omicron"),
                  Delta=sum(Variant=="Delta"),
                  Cases=n(),
                  .groups="drop") %>%
        mutate(Fraction=Omicron/Cases,
               Percent=100*Fraction)

# CDC: weekly omicron percent
cdc <- read.csv("input-data-freeze/cdc-northeast.csv",
                na.strings="",
                colClasses=c("Week.start"="Date",
                             "Week.end"="Date",
                             "Week.mid"="Date")) %>%
        mutate(CI95.low=CI95.low/100,
               CI95.high=CI95.high/100,
               Date=Week.mid,
               Fraction=Percent/100,
               Affiliation="NE") %>%
        select(Affiliation, Date, Fraction, Percent, CI95.low, CI95.high)

# GISAID data
gisa <- read.csv("input-data-freeze/gisaid-MA.csv",
                 colClasses=c(date="Date")) %>%
        rename(Date=date, Variant=variant) %>%
        filter(Date >= "2021-12-01",
               Date <= "2022-01-01") %>%  
        select(Date, Variant) %>%
        group_by(Date) %>%
        summarise(Omicron=sum(Variant=="Omicron"),
                  Delta=sum(Variant=="Delta"),
                  Cases=n(),
                  .groups="drop") %>%
        mutate(Fraction=Omicron/Cases,
               Percent=100*Fraction,
               Affiliation="MA") %>%
        select(Affiliation, Date, Fraction, Percent)


# add CDC and GISAID to prct
prct <- prct %>%
        full_join(gisa, by=c("Affiliation", "Date", "Fraction", "Percent")) %>%
        full_join(cdc, by=c("Affiliation", "Date", "Fraction", "Percent")) %>%
        mutate(Affiliation=factor(Affiliation, 
                                  levels=c("Student", "Employee", "MA", "NE")))

## logistic growth -------------------------------------------------------------

# remove MA, CDC and get BU total cases per day
tots <- prct %>%
  filter(!Affiliation %in% c("NE", "MA")) %>%
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
ggsave("outputs/cases_bu.svg", units="cm", width=15, height=10)


# remove MA, CDC and get BU cases per day by affiliation
affil <- prct %>%
  filter(!Affiliation %in% c("NE", "MA")) %>%
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
ggsave("outputs/affil_bu.svg", units="cm", width=15, height=10)

fish = fisher.test(cbind(c(sum(bu$Affiliation == 'Student' & bu$Variant == 'Delta'), 
                           sum(bu$Affiliation == 'Student' & bu$Variant == 'Omicron')),
                         c(sum(bu$Affiliation == 'Employee' & bu$Variant == 'Delta'), 
                           sum(bu$Affiliation == 'Employee' & bu$Variant == 'Omicron'))))

# students vs. employees vs. CDC
p3 <- prct %>%
      ggplot(aes(Date, Fraction, ymin=CI95.low, ymax=CI95.high)) +
      stat_smooth(data=subset(prct, Affiliation!="NE"), aes(col=Affiliation), 
                  method="glm", se = FALSE, fullrange=TRUE, 
                  method.args = list(family=binomial)) +
      geom_errorbar(data=subset(prct, Affiliation=="NE"), size=1, width=1, col="darkgrey") +
      geom_point(aes(fill=Affiliation), pch=21, size=3) +
      scale_color_manual(values=affl) +
      scale_fill_manual(values=affl) +
      labs(x="date", 
           y="Omicron fraction",
           fill="affiliation") +
      scale_x_continuous(breaks=as.Date(c("2021-11-30", "2021-12-05", "2021-12-10", 
                                          "2021-12-15", "2021-12-20", "2021-12-25", "2021-12-30", "2022-01-4")), 
                         limits=as.Date(c("2021-11-30", "2022-01-05")),
                         labels=c("Nov 30","Dec 5", "Dec 10","Dec 15","Dec 20",
                                  "Dec 25","Dec 30", "Jan 4")) +
      guides(col=FALSE) # [JT] removed extra legend
p3
ggsave("outputs/logfit-bu.png", units="cm", width=20, height=10)
ggsave("outputs/logfit-bu.svg", units="cm", width=15, height=10)

# concatenate
p <- cowplot::plot_grid(p1, p2, 
                        ncol=2,
                        labels=c("A", "B"))
p <- cowplot::plot_grid(p, p3, 
                        nrow=2,
                        labels=c("", "C"))
p
ggsave("outputs/figure3.png", units="cm", width=20, height=20)
#ggsave("outputs/figure3.svg", units="cm", width=20, height=20)

# logistic fits
fits <- unique(subset(prct, Affiliation!="NE")$Affiliation) %>%
        lapply(function(i) {
          prct %>%
            filter(Affiliation==i) %>%
            glm(Fraction ~ Date, data=., family="binomial")
        }) 
names(fits) <- unique(subset(prct, Affiliation!="NE")$Affiliation)

# null model to determine McFadden's R^2
null <- unique(subset(prct, Affiliation!="NE")$Affiliation) %>%
        lapply(function(i) {
          prct %>%
            filter(Affiliation==i) %>%
            glm(Fraction ~ 1, data=., family="binomial")
        }) 
names(null) <- unique(subset(prct, Affiliation!="NE")$Affiliation)

# summary statistics for each logistic model
# Student
confint(fits$Student, '(Intercept)', 0.95)
confint(fits$Student, 'Date', 0.95)
r2 <- 1-logLik(fits$Student)/logLik(null$Student)
pred <- data.frame(Date=seq(as.Date("2021-11-30"), as.Date("2021-12-31"), by="days"))
pred <- ciTools::add_ci(pred, fits$Student, type="boot", names=c("lwr", "upr"), 
                        alpha=0.05, nSims=500) %>%
        mutate(type="bootstrap")
pred$Affiliation <- "Student"
ci <- pred
# Employee
confint(fits$Employee, '(Intercept)', 0.95)
confint(fits$Employee, 'Date', 0.95)
r2 <- 1-logLik(fits$Employee)/logLik(null$Employee)
pred <- data.frame(Date=seq(as.Date("2021-11-30"), as.Date("2021-12-31"), by="days"))
pred <- ciTools::add_ci(pred, fits$Employee, type="boot", names=c("lwr", "upr"), 
                        alpha=0.05, nSims=500) %>%
        mutate(type="bootstrap")
pred$Affiliation <- "Employee"
ci <- rbind(ci, pred)
# MA
confint(fits$MA, '(Intercept)', 0.95)
confint(fits$MA, 'Date', 0.95)
r2 <- 1-logLik(fits$MA)/logLik(null$MA)
pred <- data.frame(Date=seq(as.Date("2021-11-30"), as.Date("2021-12-31"), by="days"))
pred <- ciTools::add_ci(pred, fits$MA, type="boot", names=c("lwr", "upr"), 
                        alpha=0.05, nSims=500) %>%
        mutate(type="bootstrap")
pred$Affiliation <- "MA"
ci <- rbind(ci, pred)

# write 95% CIs to .csv
ci <- ci %>%
      rename("Point_est"="pred", "CI95_upper"="upr", "CI95_lower"="lwr")
ci <- ci[,!(names(ci) %in% "type")]
write.csv(ci, "outputs/bootstrap_ci_logfit_affiliation.csv", row.names=FALSE)

# clean up
rm(affil, ci, fish, fits, null, p1, p2, p3, p, pred, r2, tots)

## Ct comparison ---------------------------------------------------------------

# compare across-affiliation Cts
pval = wilcox.test(bu$Ct_N1[bu$Affiliation == 'Student'], bu$Ct_N1[bu$Affiliation == 'Employee'])$p.value
pval = wilcox.test(bu$Ct_N2[bu$Affiliation == 'Student'], bu$Ct_N2[bu$Affiliation == 'Employee'])$p.value
pval = wilcox.test(bu$Ct_RP[bu$Affiliation == 'Student'], bu$Ct_RP[bu$Affiliation == 'Employee'])$p.value

# compare within-affiliation variant vs. primer
ct <- bu %>%
      select(Affiliation, Variant, starts_with("Ct")) %>%
      reshape2::melt(id.vars=c("Affiliation", "Variant"),
                     value.name="Ct",
                     variable.name="Primer") %>%
      mutate(Primer=str_extract(Primer, "(?<=Ct_)[A-z0-9]+")) %>%
      # remove primer dropouts
      # remove spike (SGTF makes it meaningless)
      # remove crazy Cts
      filter(!is.na(Ct),
             Ct >= 5,
             Ct <= 40) %>%
      # groups are primer, affiliation, and variant
      mutate(Group=paste0(Affiliation, ".", Primer, "-", Variant),
             Label=paste(Affiliation, Primer))

# calculate FDR-adjusted p-value
pval <- pairwise.wilcox.test(ct$Ct, ct$Group, p.adjust.method="BH")
pval

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
ggsave("outputs/supp_figure3.svg", width=12, height=12, units="cm")

# supplemental: NEU ORF
p1 = plot.ct("Student RP")
p2 = plot.ct("Employee RP")
cowplot::plot_grid(p1, p2, ncol=2, labels="AUTO")
ggsave("outputs/supp_figure4.png", width=12, height=12, units="cm")
ggsave("outputs/supp_figure4.svg", width=12, height=12, units="cm")

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
