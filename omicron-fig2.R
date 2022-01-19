#!/usr/bin/env Rscript

## version history ----------------------------------------------------------
# v1 by J.T. 6 Jan 2022, v2 by B.P. 7 Jan 2022, fig1 by B.P. 9 Jan 2022
# [JT] v3 by J. T. 11 Jan 2022
# [JT] more input and output formatting; code testing and commenting

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/Volumes/GoogleDrive/My Drive/Omicron/omicron_repo")
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# global variables
cols <- c(Delta="#bdbdbd", Omicron="#ca0020")

## inputs ----------------------------------------------------------------------
# institutions: every case
# [JT] updated read data
inst <- read.csv("input-data-freeze/combined-inst.csv",
                 na.strings="#N/A",
                 colClasses=c(Date="Date")) %>% # set col classes while reading
        filter(Variant %in% c("Delta", "Omicron"), # only keep delta and omicron
               Date >= "2021-12-02",
               Date <= "2021-12-21")

## calculate and format adjusted p-values --------------------------------------
# pairwise comparison for all institutions by primer
# combined w/ FDR-corrected pvals
ct <- inst %>%
      select(-Date) %>%
      reshape2::melt(id.vars=c("Institution", "Variant"),
                     variable.name="Primer",
                     value.name="Ct") %>%
      # get primer names
      mutate(Primer=str_extract(Primer,"(?<=_)[A-z0-9]+(?=$)")) %>%
      # remove primer dropouts
      # remove spike (SGTF makes it meaningless)
      # remove crazy Cts
       filter(!is.na(Ct),
             Primer != "S",
             Ct >= 5,
             Ct <= 40) %>%
      # groups are primer, institution, and variant
      mutate(Group=paste0(Institution, ".", Primer, "-", Variant),
             Label=paste(Institution, Primer))

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
        mutate(Institution=str_extract(Rows, "^[A-z ]+(?=\\.)"),
               Institution2=str_extract(Columns, "^[A-z ]+(?=\\.)"),
               Primer=str_extract(Rows, "(?<=\\.)[A-z0-9]+(?=\\-)"),
               Primer2=str_extract(Columns, "(?<=\\.)[A-z0-9]+(?=\\-)"),
               Variant1=str_extract(Rows, "[A-z]+$"),
               Variant2=str_extract(Columns, "[A-z]+$")) %>%
        filter(Institution==Institution2,
               Primer==Primer2,
               !is.na(p.adj)) %>%
        select(Institution, Primer, Variant1, Variant2, p.adj) %>%
        mutate(p.adj=scales::pvalue(p.adj, accuracy=0.0001),
               Label=paste(Institution, Primer))

## figure ----------------------------------------------------------------------
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
p1 <- plot.ct("BU N1")
p2 <- plot.ct("HU N1")
p3 <- plot.ct("BU N2")
p4 <- plot.ct("NEU N2")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels="AUTO")
ggsave("outputs/figure2.png", width=12, height=12, units="cm")
#ggsave("outputs/figure2.svg", width=12, height=12, units="cm")

# supplemental: NEU ORF
plot.ct("NEU ORF")
plot.ct("BU RP")

# table
ct %>%
  group_by(Institution, Variant, Primer) %>%
  summarise(Median=median(Ct, na.rm = TRUE),
            Mean=mean(Ct, na.rm = TRUE),
            .groups="drop") %>%
  arrange(Institution, Primer) %>%
  write.csv("outputs/ct-breakdown.csv", row.names=FALSE)

## fin! ------------------------------------------------------------------------
sessionInfo()
