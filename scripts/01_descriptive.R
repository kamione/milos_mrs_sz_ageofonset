# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(haven)
library(gt)
library(gtsummary)
library(psych)
library(ggpubr)

# Data I/O ---------------------------------------------------------------------
data <- here("data", "raw", "MILOS_20210621.csv") %>% 
    read_csv(col_types = cols())


# Demographics -----------------------------------------------------------------
basic_demog_table1 <- data %>% 
    select(group2, Gender, Age, Years_of_education, Occupation_recode, dx_sz,
           DUP_months, PANSS_P_Total:PANSS_G_Total, CDSS_TOTAL, MRS_TOTAL,
           SOFAS, CPZ_before_scan) %>%
    filter(group2 %in% c("EOS", "LOS")) %>% 
    mutate(group2 = factor(group2)) %>%
    tbl_summary(
        by = group2,
        missing = "no",
        digits = list(
            where(is.numeric) ~ c(1, 2),
            all_categorical() ~ c(0, 1),
            c("CPZ_before_scan") ~ c(0, 0)
        ),
        statistic = list(
            where(is.numeric) ~ "{mean} ({sd})",
            c("DUP_months", "CPZ_before_scan") ~ "{median} ({IQR})",
            all_categorical() ~ "{n} ({p}%)"
        ),
        value = list(
            Gender ~ "Male",
            Occupation_recode ~ "employed",
            dx_sz ~ "Yes"),
        label = list(Gender ~ "Male",
                     Age ~ "Age",
                     Years_of_education ~ "Education Years",
                     Occupation_recode ~ "Employement",
                     dx_sz ~ "SZ Diagnosis",
                     DUP_months ~ "DUP in months, median (IQR)",
                     PANSS_P_Total ~ "PANSS Positive",
                     PANSS_N_Total ~ "PANSS Negative",
                     PANSS_G_Total ~ "PANSS General",
                     CDSS_TOTAL ~ "CDSS",
                     MRS_TOTAL ~ "YMRS",
                     SOFAS ~ "SOFAS",
                     CPZ_before_scan ~ "CPZ, median (IQR)")
    ) %>% 
    add_overall() %>% 
    modify_footnote(update = everything() ~ NA) %>% 
    add_p(test = list(where(is.numeric) ~ "t.test",
                      all_categorical() ~ "chisq.test")) %>% 
    add_q() %>% 
    bold_p(q = TRUE) %>% 
    modify_column_hide("p.value")

basic_demog_table2 <- data %>% 
    select(group2, Gender, Age, Years_of_education, Occupation_recode) %>% 
    filter(group2 %in% c("EOS-C", "LOS-C")) %>% 
    mutate(group2 = factor(group2)) %>% 
    tbl_summary(
        by = group2,
        missing = "no",
        digits = list(all_categorical() ~ c(0, 1),
                      all_continuous() ~ c(1, 2)),
        statistic = list(all_continuous() ~ "{mean} ({sd})",
                         all_categorical() ~ "{n} ({p}%)"),
        value = list(Gender ~ "Male",
                     Occupation_recode ~ "employed"),
        label = list(Gender ~ "Male",
                     Age ~ "Age",
                     Years_of_education ~ "Education Years",
                     Occupation_recode ~ "Employement")
    ) %>% 
    add_overall() %>% 
    modify_footnote(update = everything() ~ NA) %>% 
    add_p(test = list(all_continuous() ~ "t.test",
                      all_categorical() ~ "chisq.test")) %>% 
    add_q() %>% 
    bold_p(q = TRUE) %>% 
    modify_column_hide("p.value")
         
basic_demog_table3 <- data %>% 
    select(group2, Gender, Age, Years_of_education, Occupation_recode) %>% 
    filter(group2 %in% c("EOS", "EOS-C")) %>% 
    mutate(group2 = factor(group2)) %>% 
    tbl_summary(
        by = group2,
        missing = "no",
        digits = list(all_categorical() ~ c(0, 1),
                      all_continuous() ~ c(1, 2)),
        statistic = list(all_continuous() ~ "{mean} ({sd})",
                         all_categorical() ~ "{n} ({p}%)"),
        value = list(Gender ~ "Male",
                     Occupation_recode ~ "employed"),
        label = list(Gender ~ "Male",
                     Age ~ "Age",
                     Years_of_education ~ "Education Years",
                     Occupation_recode ~ "Employement")
    ) %>% 
    modify_footnote(update = everything() ~ NA) %>% 
    add_p(test = list(all_continuous() ~ "t.test",
                      all_categorical() ~ "chisq.test")) %>% 
    add_q() %>% 
    bold_p(q = TRUE) %>% 
    modify_column_hide(c("p.value", "stat_1", "stat_2"))

basic_demog_table4 <- data %>% 
    select(group2, Gender, Age, Years_of_education, Occupation_recode) %>% 
    filter(group2 %in% c("LOS", "LOS-C")) %>% 
    mutate(group2 = factor(group2)) %>% 
    tbl_summary(
        by = group2,
        missing = "no",
        digits = list(all_categorical() ~ c(0, 1),
                      all_continuous() ~ c(1, 2)),
        statistic = list(all_continuous() ~ "{mean} ({sd})",
                         all_categorical() ~ "{n} ({p}%)"),
        value = list(Gender ~ "Male",
                     Occupation_recode ~ "employed"),
        label = list(Gender ~ "Male",
                     Age ~ "Age",
                     Years_of_education ~ "Education Years",
                     Occupation_recode ~ "Employement")
    ) %>% 
    modify_footnote(update = everything() ~ NA) %>% 
    add_p(test = list(all_continuous() ~ "t.test",
                      all_categorical() ~ "chisq.test")) %>% 
    add_q() %>% 
    bold_p(q = TRUE) %>% 
    modify_column_hide(c("p.value", "stat_1", "stat_2"))
 
basic_demog_table <- tbl_merge(
    tbls = list(basic_demog_table1, 
                basic_demog_table2, 
                basic_demog_table3, 
                basic_demog_table4),
    tab_spanner = c("**First Episode Schizophrenia**", 
                    "**Heathy Controls**",
                    "**EOS vs. EOS-C**",
                    "**LOS vs. LOS-C**"))

basic_demog_table %>% 
    as_gt() %>%
    tab_footnote( # and can modify/add footnotes this way
        footnote = "SZ = Schizophrenia;
            DUP = Duration of untreated psychosis;
            IQR = Interquartile range;
            PANSS = Positive and Negative Syndrome Scale; 
            CDSS = Calgary Depression Scale for Schizophrenia; 
            MRS = Young Mania Rating Scale; 
            SOFAS = Social and Occupational Functioning Scale;
            CPZ = Chlorpromazine equivalent doses (before brain scan);
            EOS = Early-onset schizophrenia;
            LOS = Late-onset schizophrenia; 
            EOS-C = Healthy controls for EOS;
            LOS-C = Healthy controls for LOS",
        locations = cells_column_labels(columns = label)
    ) %>% 
    tab_footnote( # and can modify/add footnotes this way
        footnote = "Pearson's Chi-squared test; Welch Two Sample t-test",
        locations = cells_column_labels(columns = c("q.value_1", "q.value_2", "q.value_3", "q.value_4"))
    ) %>% 
    tab_style(
        style = list(cell_fill(color = "lightgrey")),
        locations = cells_body(columns = c(stat_0_1, stat_0_2))
    ) %>% 
    gtsave(filename = here("outputs", "tables", "basic_demog_table.html"))

# get site differences
site_compare_table <- data %>% 
    mutate(site = if_else(is.na(AP), "HKU", "HKSH")) %>% 
    filter(group == "FEP") %>% 
    select(site, Gender, Age, Years_of_education, Occupation_recode, dx_sz,
           DUP_months, PANSS_P_Total:PANSS_G_Total, CDSS_TOTAL, MRS_TOTAL,
           SOFAS, CPZ_before_scan) %>%
    tbl_summary(
        by = site,
        missing = "no",
        digits = list(all_categorical() ~ c(0, 1),
                      all_continuous() ~ c(1, 2),
                      c("CPZ_before_scan") ~ c(0, 0)),
        statistic = list(all_continuous() ~ "{mean} ({sd})",
                         c("DUP_months", "CPZ_before_scan") ~ "{median} ({IQR})",
                         all_categorical() ~ "{n} ({p}%)"),
        value = list(Gender ~ "Male",
                     Occupation_recode ~ "employed",
                     dx_sz ~ "Yes"),
        label = list(Gender ~ "Male",
                     Age ~ "Age",
                     Years_of_education ~ "Education Years",
                     Occupation_recode ~ "Employement",
                     dx_sz ~ "SZ Diagnosis",
                     DUP_months ~ "DUP in months, median (IQR)",
                     PANSS_P_Total ~ "PANSS Positive",
                     PANSS_N_Total ~ "PANSS Negative",
                     PANSS_G_Total ~ "PANSS General",
                     CDSS_TOTAL ~ "CDSS",
                     MRS_TOTAL ~ "YMRS",
                     SOFAS ~ "SOFAS",
                     CPZ_before_scan ~ "CPZ, median (IQR)")
    ) %>% 
    add_overall() %>% 
    modify_footnote(update = everything() ~ NA) %>% 
    add_p(test = list(all_continuous() ~ "wilcox.test",
                      all_categorical() ~ "chisq.test"))

site_compare_table %>% 
    as_gt() %>%
    tab_footnote( # and can modify/add footnotes this way
        footnote = "SZ = Schizophrenia;
            DUP = Duration of untreated psychosis;
            IQR = Interquartile range;
            PANSS = Positive and Negative Syndrome Scale; 
            CDSS = Calgary Depression Scale for Schizophrenia; 
            MRS = Young Mania Rating Scale; 
            SOFAS = Social and Occupational Functioning Scale;
            CPZ = Chlorpromazine equivalent doses (before brain scan);
            EOS = Early-onset schizophrenia;
            LOS = Late-onset schizophrenia; 
            EOS-C = Healthy controls for EOS;
            LOS-C = Healthy controls for LOS",
        locations = cells_column_labels(columns = label)
    ) %>%
    tab_style(
        style = list(cell_fill(color = "lightgrey")),
        locations = cells_body(columns = c(stat_0))
    ) %>% 
    gtsave(filename = here("outputs", "tables", "site_compare_sz_table.html"))

