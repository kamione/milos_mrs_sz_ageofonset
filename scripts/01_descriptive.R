# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(haven)
library(gt)
library(gtsummary)
library(psych)
library(ggpubr)
library(lubridate)
library(flextable)
library(officer)

sect_properties <- prop_section(
    page_size = page_size(orient = "landscape",
                          width = 8.3, height = 11.7),
    type = "continuous",
    page_margins = page_mar()
)


# Data I/O ---------------------------------------------------------------------
aoo_data <- here("data", "raw", "2022 01 14 MILOS IRAOS.xlsx") %>% 
    readxl::read_xlsx(skip = 1) %>% 
    select(`Case ID`, DOB, `First symptom_date`) %>% 
    mutate(DOB = as.numeric(DOB)) %>% 
    mutate(DOB = lubridate::as_date(DOB, origin = "1900-01-01")) %>% 
    mutate(first_symptom_date = lubridate::as_date(`First symptom_date`, origin = "1900-01-01")) %>% 
    select(-`First symptom_date`) %>% 
    mutate(aoo = (first_symptom_date - DOB) / ddays(365.25)) %>% 
    rename("subjectcode" = "Case ID")
    
data <- here("data", "raw", "MILOS_20210621.csv") %>% 
    read_csv(col_types = cols()) %>% 
    left_join(aoo_data, by = "subjectcode") %>% 
    mutate(group = factor(group, levels = c("FEP", "HC"))) %>% 
    mutate(group2 = factor(
        group2,
        levels = c("EOS", "EOS-C", "LOS", "LOS-C"),
        labels = c("AOS", "AOS-C", "LOS", "LOS-C"))
    ) %>% 
    filter(!(aoo < 40 & group2 == "LOS"))

write_rds(data, here("data", "processed", "data4analysis.rds"))
write_csv(data, here("data", "processed", "data4analysis.csv"))
# data <- read_rds(here("data", "processed", "data4analysis.rds"))

qa_data <-  here("data", "raw", "mrs_qa.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    left_join(data %>% select(subjectcode, group, group2), by = "subjectcode")


# Demographics -----------------------------------------------------------------
basic_demog_table1 <- data %>% 
    select(group2, Gender, Age, Years_of_education, Occupation_recode, dx_sz,
           DUP_months, aoo, PANSS_P_Total:PANSS_G_Total, CDSS_TOTAL, MRS_TOTAL,
           SOFAS, CPZ_before_scan) %>%
    filter(group2 %in% c("AOS", "LOS")) %>% 
    mutate(group2 = droplevels(group2)) %>%
    tbl_summary(
        by = group2,
        missing = "no",
        digits = list(
            where(is.numeric) ~ c(1, 2),
            all_categorical() ~ c(0, 1),
            c("CPZ_before_scan") ~ c(0, 0)
        ),
        type = list(MRS_TOTAL ~ 'continuous'),
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
                     Occupation_recode ~ "Employment",
                     dx_sz ~ "SZ Diagnosis",
                     DUP_months ~ "DUP in months, median (IQR)",
                     aoo ~ "Age of Onset",
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
    add_p(
        test = list(
            where(is.numeric) ~ "t.test", all_categorical() ~ "chisq.test"
        ),
        test.args = all_tests("t.test") ~ list(var.equal = TRUE)
    ) %>% 
    add_q() %>% 
    bold_p(q = TRUE) %>% 
    modify_column_hide("p.value")

basic_demog_table2 <- data %>% 
    select(group2, Gender, Age, Years_of_education, Occupation_recode) %>% 
    filter(group2 %in% c("AOS-C", "LOS-C")) %>% 
    mutate(group2 = droplevels(group2)) %>% 
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
                     Occupation_recode ~ "Employment")
    ) %>% 
    add_overall() %>% 
    modify_footnote(update = everything() ~ NA) %>% 
    add_p(
        test = list(
            where(is.numeric) ~ "t.test", all_categorical() ~ "chisq.test"
        ),
        test.args = all_tests("t.test") ~ list(var.equal = TRUE)
    ) %>% 
    add_q() %>% 
    bold_p(q = TRUE) %>% 
    modify_column_hide("p.value")
         
basic_demog_table3 <- data %>% 
    select(group2, Gender, Age, Years_of_education, Occupation_recode) %>% 
    filter(group2 %in% c("AOS", "AOS-C")) %>% 
    mutate(group2 = droplevels(group2)) %>% 
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
                     Occupation_recode ~ "Employment")
    ) %>% 
    modify_footnote(update = everything() ~ NA) %>% 
    add_p(
        test = list(
            where(is.numeric) ~ "t.test", all_categorical() ~ "chisq.test"
        ),
        test.args = all_tests("t.test") ~ list(var.equal = TRUE)
    ) %>% 
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
                     Occupation_recode ~ "Employment")
    ) %>% 
    modify_footnote(update = everything() ~ NA) %>% 
    add_p(
        test = list(
            where(is.numeric) ~ "t.test", all_categorical() ~ "chisq.test"
        ),
        test.args = all_tests("t.test") ~ list(var.equal = TRUE)
    ) %>% 
    add_q() %>% 
    bold_p(q = TRUE) %>% 
    modify_column_hide(c("p.value", "stat_1", "stat_2"))
 
basic_demog_table <- tbl_merge(
    tbls = list(basic_demog_table1, 
                basic_demog_table2, 
                basic_demog_table3, 
                basic_demog_table4),
    tab_spanner = c("**First Episode Schizophrenia (FES)**", 
                    "**Heathy Control (HC)**",
                    "**AOS vs. AOS-C**",
                    "**LOS vs. LOS-C**")) %>% 
    modify_footnote(c("q.value_1", "q.value_2", "q.value_3", "q.value_4") ~ NA)

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
            AOS = Average-onset schizophrenia;
            LOS = Late-onset schizophrenia; 
            AOS-C = Healthy controls for AOS;
            LOS-C = Healthy controls for LOS",
        locations = cells_column_labels(columns = label)
    ) %>% 
    tab_footnote( # and can modify/add footnotes this way
        footnote = "n (%); Mean (SD) unless specified otherwise",
        locations = cells_column_labels(
            columns = c("stat_0_1", "stat_1_1", "stat_2_1", 
                        "stat_0_2", "stat_1_2", "stat_2_2")
        )
    ) %>% 
    tab_footnote( # and can modify/add footnotes this way
        footnote = "Pearson's Chi-squared test; Two Sample t-test",
        locations = cells_column_labels(
            columns = c("q.value_1", "q.value_2", "q.value_3", "q.value_4")
        )
    ) %>% 
    tab_footnote( # and can modify/add footnotes this way
        footnote = "False discovery rate correction for multiple testing",
        locations = cells_column_labels(
            columns = c("q.value_1", "q.value_2", "q.value_3", "q.value_4")
        )
    ) %>% 
    tab_style(
        style = list(cell_fill(color = "lightgrey")),
        locations = cells_body(columns = c(stat_0_1, stat_0_2))
    ) %>% 
    gtsave(filename = here("outputs", "tables", "basic_demog_table.html"))

basic_demog_table %>% 
    as_flex_table() %>% 
    bold(part = "header") %>% 
    footnote(
        i = 2, j = 1,
        value = as_paragraph(paste(
            "SZ = Schizophrenia;",
            "DUP = Duration of untreated psychosis;",
            "IQR = Interquartile range;",
            "PANSS = Positive and Negative Syndrome Scale;",
            "CDSS = Calgary Depression Scale for Schizophrenia;",
            "MRS = Young Mania Rating Scale;", 
            "SOFAS = Social and Occupational Functioning Scale;",
            "CPZ = Chlorpromazine equivalent doses (before brain scan);",
            "AOS = Average-onset schizophrenia;",
            "LOS = Late-onset schizophrenia;", 
            "AOS-C = Healthy controls for AOS;",
            "LOS-C = Healthy controls for LOS"
        )),
        ref_symbols = "1",
        part = "header"
    ) %>% 
    footnote(
        i = 2, j = c(2, 3, 4, 6, 7, 8),
        value = as_paragraph(c(
            "n (%); Mean (SD) unless specified otherwise"
        )),
        ref_symbols = c("2"),
        part = "header"
    ) %>% 
    footnote(
        i = 2, j = c(5, 9, 10, 11),
        value = as_paragraph(c(
            "Pearson's Chi-squared test; Two Sample t-test"
        )),
        ref_symbols = c("3"),
        part = "header"
    ) %>% 
    footnote(
        i = 2, j = c(5, 9, 10, 11),
        value = as_paragraph(c(
            "False discovery rate correction for multiple testing"
        )),
        ref_symbols = c("4"),
        part = "header"
    ) %>% 
    bg(j = c(2, 6), bg = "grey85", part = "body") %>% 
    save_as_docx(
        path = here("outputs", "tables", "basic_demog_table.docx"),
        pr_section = sect_properties)
    

# get site differences
site_compare_table <- data %>% 
    mutate(site = if_else(is.na(AP), "HKU", "HKSH")) %>% 
    filter(group == "FEP") %>% 
    select(site, Gender, Age, Years_of_education, Occupation_recode, dx_sz,
           DUP_months, aoo, PANSS_P_Total:PANSS_G_Total, CDSS_TOTAL, MRS_TOTAL,
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
                     Occupation_recode ~ "Employment",
                     dx_sz ~ "SZ Diagnosis",
                     DUP_months ~ "DUP in months, median (IQR)",
                     aoo ~ "Age of Onset",
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
    add_p(
        test = list(
            where(is.numeric) ~ "t.test", all_categorical() ~ "chisq.test"
        ),
        test.args = all_tests("t.test") ~ list(var.equal = TRUE)
    ) %>% 
    bold_p()

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
            AOS = Average-onset schizophrenia;
            LOS = Late-onset schizophrenia; 
            AOS-C = Healthy controls for AOS;
            LOS-C = Healthy controls for LOS",
        locations = cells_column_labels(columns = label)
    ) %>%
    tab_style(
        style = list(cell_fill(color = "lightgrey")),
        locations = cells_body(columns = c(stat_0))
    ) %>% 
    gtsave(filename = here("outputs", "tables", "site_compare_sz_table.html"))


# QA table for supplementary
acc_qa_table <- qa_data %>% 
    filter(voi == "acc") %>% 
    select(-c(voi, subjectcode, group2)) %>% 
    tbl_summary(
        by = group,
        type = list(glx_crlb ~ "continuous", naa_crlb ~ "continuous"),
        statistic = list(all_continuous() ~ "{mean} ({sd})"),
        missing = "no",
        label = list(
            fwhm ~ "FWHM",
            sn_ratio ~ "SN Ratio",
            glx_crlb ~ "GLX CRLB",
            naa_crlb ~ "NAA CRLB"
        )
    ) %>% 
    add_p() %>% 
    add_q()

bg_qa_table <- qa_data %>% 
    filter(voi == "bg") %>% 
    select(-c(voi, subjectcode, group2)) %>% 
    tbl_summary(
        by = group,
        type = list(glx_crlb ~ "continuous", naa_crlb ~ "continuous"),
        statistic = list(all_continuous() ~ "{mean} ({sd})"),
        missing = "no",
        label = list(
            fwhm ~ "FWHM",
            sn_ratio ~ "SN Ratio",
            glx_crlb ~ "GLX CRLB",
            naa_crlb ~ "NAA CRLB"
        )
    ) %>% 
    add_p() %>% 
    add_q()
    
list(acc_qa_table, bg_qa_table) %>% 
    tbl_stack(group_header = c("Anterior Cingular Cortex", "Basal Ganglia")) %>% 
    modify_footnote(label ~ paste(
        "FWHM = Full-width at half-maximum;", 
        "SN = Signal-to-Noise;",
        "GLX = Glutamine and glutamate;",
        "NAA = N-acetylaspartate;",
        "CRLB = Cramer-Rao lower bounds")
    ) %>% 
    as_flex_table() %>% 
    bold(part = "header") %>% 
    save_as_docx(
        path = here("outputs", "tables", "mrs_qa_table.docx"),
        pr_section = sect_properties)
    
    
