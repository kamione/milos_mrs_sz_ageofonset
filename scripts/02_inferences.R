# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(ggpubr)
library(glue)
library(rstatix)
library(gtsummary)
library(ggthemes)
library(gt)
library(flextable)
library(officer)

sect_properties <- prop_section(
    page_size = page_size(orient = "portrait",
                          width = 11.7, height = 8.3),
    type = "continuous",
    page_margins = 
        page_mar(
            bottom = 0.5,
            top = 0.5,
            right = 0.5,
            left = 0.5,
        )
)

sect_properties_landscape <- prop_section(
    page_size = page_size(orient = "landscape",
                          width = 8.3, height = 11.7),
    type = "continuous",
    page_margins = 
        page_mar(
            bottom = 0.5,
            top = 0.5,
            right = 0.5,
            left = 0.5,
        )
)

source(here("src", "R", "stat.R"))
filter <- dplyr::filter # resolve the conflict


# Data I/O ---------------------------------------------------------------------
data <- read_rds(here("data", "processed", "data4analysis.rds"))


# 1. Metabolite Group Difference -----------------------------------------------
imaging_list1 <- c("BG Glx" = "bg_glugln",
                   "ACC Glx" = "acc_glugln",
                   "BG NAA" = "bg_naa",
                   "ACC NAA" = "acc_naa")

figs_2groups <- lapply(seq_along(imaging_list1), function(x) {
    
    ggviolin(data,
             x = "group",
             y = imaging_list1[x],
             add = c("jitter", "boxplot"),
             add.params = list(color = "gray30", alpha = 0.8)) +
        scale_y_continuous(labels = scales::label_percent(accuracy = 0.01, scale = 1)) +
        labs(x = "", y = "Concentration Ratio", title = names(imaging_list1[x]))  +
        stat_compare_means(
            method = "t.test",
            label = "p.signif",
            symnum.args <- list(
                cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                symbols = c("***", "**", "*", "ns")
            ),
            comparisons = list(c("FEP", "HC"))
        ) +
        theme_pander() +
        theme(plot.margin = unit(rep(2, 4), "mm"))
})

ggarrange(plotlist = figs_2groups, ncol = 2, nrow = 2, labels = LETTERS[1:4]) %>% 
    ggexport(filename = here("outputs", "figs", "mrs_2groups_fig.pdf"),
             width = 10, height = 7)

y_positions <- c(2.8, 3.9, 2.3, 2.3)

figs_4groups <- lapply(seq_along(imaging_list1), function(x) {
    
    fit <- kruskal.test(as.formula(glue("{imaging_list1[x]} ~ group2")), data = data)
    print(broom::glance(fit))
    
    # get the effect size with 95% CI
    effectsize <- glue("{imaging_list1[x]} ~ group2") %>% 
        as.formula() %>% 
        kruskal_effsize(ci = TRUE, data = data)
    print(effectsize)
    
    step_increase <- ifelse(imaging_list1[x] == "bg_naa", 2.1, 0.3)
    
    # post-hoc analysis
    stat.test <- data %>% 
        t_test(
            formula(glue("{imaging_list1[x]} ~ group2")),
            comparisons = list(c("AOS", "AOS-C"), c("LOS", "LOS-C"), c("AOS", "LOS")),
            p.adjust.method = "fdr",
            var.equal = FALSE,
            detailed = TRUE
        ) %>%
        add_xy_position(step.increase = step_increase)
    print(stat.test %>% select(statistic, df, p.adj, conf.low, conf.high))
    
    data %>% 
        cohens_d(
            formula(glue("{imaging_list1[x]} ~ group2")),
            comparisons = list(c("AOS", "AOS-C"), c("LOS", "LOS-C"), c("AOS", "LOS")),
            var.equal = FALSE
        ) %>% 
        print()
    
    ggviolin(data,
              x = "group2",
              y = imaging_list1[x],
              add = c("jitter", "boxplot"),
              add.params = list(color = "gray30", alpha = 0.8)) +
        scale_y_continuous(labels = scales::label_percent(accuracy = 0.01, scale = 1)) +
        labs(x = "", y = "Concentration Ratio", title = names(imaging_list1[x])) +
        stat_pvalue_manual(stat.test, hide.ns = TRUE) +
        stat_compare_means(label.x = 0.8, label.y = y_positions[x]) +
        theme_pander() +
        theme(plot.margin = unit(c(2, 2, 2, 2), "mm"))
    
})

ggarrange(plotlist = figs_4groups, ncol = 2, nrow = 2, labels = LETTERS[1:4]) %>% 
    ggexport(filename = here("outputs", "figs", "mrs_4groups_fig.pdf"),
             width = 10, height = 7)

# 2. Correlation between Metabolite and Cognitive Performance ------------------
corr_res <- data %>% 
    select(bg_glugln:acc_naa, DS_F:Monotone) %>% 
    select(-c(WC_error, VF_error_count)) %>% 
    rename(
        "BG Glx" = "bg_glugln",
        "BG NAA" = "bg_naa",
        "ACC Glx" = "acc_glugln",
        "ACC NAA" = "acc_naa",
        "Digit Span (Forward)" = "DS_F",
        "Digit Span (Backward)" = "DS_B",
        "Visual Pattern (Correct)" = "Visual_Pattern_correct",
        "Visual Pattern (Longest)" = "Visual_Pattern_Longest_Passed",
        "Logical Memory (Immediate)" = "LM_immediate",
        "Logical Memory (Delayed)" = "LM_delayed",
        "Digit Symbol" = "Digit_symbol",
        "MWCST (Correct)" = "WC_correct",
        "MWCST (Non-Perseverative Error)" = "WC_NP_error",
        "MWCST (Perseverative Error)" = "WC_P_error",
        "MWCST (Category Complete)" = "WC_categories",
        "Verbal Fluency (Count)" = "VF_count",
        "Verbal Fluency (Duplicate)" = "VF_duplicate_count"
    ) %>% 
    psych::corr.test(method = "spearman")

# plot correlation heat map
pdf(here("outputs", "figs", "spearman_corr_fdr.pdf"), height = 10, width = 10)
corrplot::corrplot(
    corr_res$r,
    col = colorRampPalette(c("#0C6291", "#FBFEF9", "#A63446"))(256),
    mar = rep(0, 4),
    method = "square",
    type = "upper",
    p.mat = corr_res$p,
    insig = "blank",
    diag = FALSE,
    tl.col = "#1C1C1C",
    tl.srt = 45,
    rect.col = "lightgrey",
)
dev.off()

# plot scatter plots based on correlation heat map
scatter_1 <- data %>% 
    ggscatter(x = "bg_naa", y = "Visual_Pattern_correct", color = "group2", shape = "group2", size = 3) +
    scale_color_manual(values = c("#268785", "#72636E", "tomato3", "grey30")) +
    labs(x = "BG NAA Concentration Ratio", y = "Visual Pattern (Correct)", color = "", shape = "") +
    #scale_x_continuous(breaks = seq(0.5, 1.8, 0.25), limit = c(0.5, 1.8)) +
    scale_y_continuous(breaks = seq(0, 25, 5), limit = c(-5, 25)) +
    ggthemes::theme_pander() +
    theme(plot.margin = margin(2, 2, 2, 2, "mm")) +
    geom_smooth(method = "lm", color = "grey20", fill = "grey75") +
    annotate(geom = "text",
             x = 0.5,
             y = 24,
             label = glue("italic(r)=={format_decimals(corr_res$r['BG NAA', 'Visual Pattern (Correct)'])}*','~italic(q)=={format_decimals(corr_res$p['BG NAA', 'Visual Pattern (Correct)'])}"),
             color = "gray20",
             hjust = 0,
             parse = TRUE)

scatter_2 <- data %>% 
    ggscatter(x = "bg_naa", y = "WC_NP_error", color = "group2", shape = "group2", size = 3) +
    scale_color_manual(values = c("#268785", "#72636E", "tomato3", "grey30")) +
    labs(x = "BG NAA Concentration Ratio", y = "MWCST Non-Perseverative Errors", color = "", shape = "") +
    #scale_x_continuous(breaks = seq(0.5, 1.8, 0.25), limit = c(0.5, 1.8)) +
    scale_y_continuous(breaks = seq(0, 25, 5), limit = c(-5, 25)) +
    ggthemes::theme_pander() +
    theme(plot.margin = margin(2, 2, 2, 2, "mm")) +
    geom_smooth(method = "lm", color = "grey20", fill = "grey75") +
    annotate(geom = "text",
             x = 1,
             y = 24,
             label = glue("italic(r)=={format_decimals(corr_res$r['BG NAA', 'MWCST (Non-Perseverative Error)'])}*','~italic(q)=={format_decimals(corr_res$p['BG NAA', 'MWCST (Non-Perseverative Error)'])}"),
             color = "gray20",
             hjust = 0,
             parse = TRUE)

scatter_3 <- data %>% 
    ggscatter(x = "acc_naa", y = "WC_NP_error", color = "group2", shape = "group2", size = 3) +
    scale_color_manual(values = c("#268785", "#72636E", "tomato3", "grey30")) +
    labs(x = "ACC NAA Concentration Ratio", y = "MWCST Non-Perseverative Errors", color = "", shape = "") +
    scale_x_continuous(breaks = seq(0.75, 1.8, 0.25)) +
    scale_y_continuous(breaks = seq(0, 25, 5), limit = c(-5, 25)) +
    ggthemes::theme_pander() +
    theme(plot.margin = margin(2, 2, 2, 2, "mm")) +
    geom_smooth(method = "lm", color = "grey20", fill = "grey75") +
    annotate(geom = "text",
             x = 1.32,
             y = 24,
             label = glue("italic(r)=={format_decimals(corr_res$r['ACC NAA', 'MWCST (Non-Perseverative Error)'])}*','~italic(q)=={format_decimals(corr_res$p['ACC NAA', 'MWCST (Non-Perseverative Error)'])}"),
             color = "gray20",
             hjust = 0,
             parse = TRUE)

scatter_combined <- ggarrange(scatter_1, scatter_2, scatter_3, nrow = 1, ncol = 3, common.legend = TRUE)

ggexport(scatter_combined,
         filename = here("outputs", "figs", "scatter_corr.pdf"),
         height = 4,
         width = 12)


comparison_violin <- data %>% 
    filter(group2 %in% c("AOS", "LOS")) %>% 
    select(group2, Visual_Pattern_correct, WC_NP_error) %>% 
    mutate(group2 = droplevels(group2)) %>% 
    rename(
        "Visual Pattern (Correct)" = "Visual_Pattern_correct",
        "MWCST (Non-Perseverative Error)" = "WC_NP_error"
    ) %>% 
    pivot_longer(!group2, names_to = "measures", values_to = "value") %>%
    ggviolin(x = "group2",
             y = "value",
             facet.by = "measures",
             add = c("jitter", "boxplot"),
             add.params = list(color = "gray30", alpha = 0.8)) +
    labs(x = "", y = "") +
    theme(
        strip.background = element_rect(
            color = "white", fill = "white", linetype = "solid"
        )
    ) +
    theme_pander() +
    stat_compare_means(
        method = "t.test",
        label = "p.signif",
        symnum.args <- list(
            cutpoints = c(0, 0.001, 0.01, 0.05, 1),
            symbols = c("***", "**", "*", "ns")
        ),
        comparisons = list(c("LOS", "AOS"))
    )

ggsave(filename = here("outputs", "figs", "aos_cog_violin_comparison.pdf"),
       comparison_violin,
       height = 4.5,
       width = 6)

# 3. Cognition Group Difference ------------------------------------------------
table_cog_fep <- data %>% 
    filter(group2 %in% c("AOS", "LOS")) %>% 
    select(group2, DS_F:Monotone) %>% 
    select(-c(WC_error, VF_error_count)) %>% 
    mutate(group2 = droplevels(group2)) %>% 
    rename(
        "Digit Span (Forward)" = "DS_F",
        "Digit Span (Backward)" = "DS_B",
        "Visual Pattern (Correct)" = "Visual_Pattern_correct",
        "Visual Pattern (Longest)" = "Visual_Pattern_Longest_Passed",
        "Logical Memory (Immediate)" = "LM_immediate",
        "Logical Memory (Delayed)" = "LM_delayed",
        "Digit Symbol" = "Digit_symbol",
        "MWCST (Correct)" = "WC_correct",
        "MWCST (Non-Perseverative Error)" = "WC_NP_error",
        "MWCST (Perseverative Error)" = "WC_P_error",
        "MWCST (Category Complete)" = "WC_categories",
        "Verbal Fluency (Correct)" = "VF_count",
        "Verbal Fluency (Duplicate)" = "VF_duplicate_count"
    ) %>% 
    tbl_summary(
        by = group2,
        missing = "no",
        type = where(is.numeric) ~ "continuous",
        statistic = list(everything() ~ "{mean} ({sd})")
    ) %>% 
    add_p(
        test = everything() ~ "t.test", 
        test.args = all_tests("t.test") ~ list(var.equal = TRUE)
    ) %>% 
    add_q(method = "fdr") %>% 
    add_stat(fns = everything() ~ cohen_d) %>% 
    bold_p(q = TRUE) %>%
    modify_header(add_stat_1 ~ "**Cohen'd (95% CI)**")

table_cog_fep %>% 
    as_gt() %>% 
    gtsave(filename = here("outputs", "tables", "table_cog_fep.html"))

table_cog_fep %>% 
    as_flex_table() %>% 
    bold(part = "header") %>% 
    fontsize(size = 8, part = "all") %>% 
    save_as_docx(
        path = here("outputs", "tables", "table_cog_fep.docx"),
        pr_section = sect_properties
    )

table_cog_all <- data %>% 
    select(group2, DS_F:Monotone) %>% 
    select(-c(WC_error, VF_error_count)) %>% 
    mutate(group2 = droplevels(group2)) %>% 
    rename(
        "Digit Span (Forward)" = "DS_F",
        "Digit Span (Backward)" = "DS_B",
        "Visual Pattern (Correct)" = "Visual_Pattern_correct",
        "Visual Pattern (Longest)" = "Visual_Pattern_Longest_Passed",
        "Logical Memory (Immediate)" = "LM_immediate",
        "Logical Memory (Delayed)" = "LM_delayed",
        "Digit Symbol" = "Digit_symbol",
        "MWCST (Correct)" = "WC_correct",
        "MWCST (Non-Perseverative Error)" = "WC_NP_error",
        "MWCST (Perseverative Error)" = "WC_P_error",
        "MWCST (Categories)" = "WC_categories",
        "Verbal Fluency (Correct)" = "VF_count",
        "Verbal Fluency (Duplicate)" = "VF_duplicate_count"
    ) %>% 
    tbl_summary(
        by = group2,
        missing = "no",
        type = where(is.numeric) ~ "continuous",
        statistic = list(everything() ~ "{mean} ({sd})")
    ) %>% 
    add_p(test = list(everything() ~ "aov")) %>% 
    add_q(method = "fdr") %>% 
    bold_p(q = TRUE) %>%
    add_stat(everything() ~ add_stat_pairwise) %>% 
    modify_column_hide(c("**AOS-C vs. LOS**", "**AOS vs. LOS-C**"))
    
table_cog_all %>% 
    as_gt() %>% 
    gtsave(filename = here("outputs", "tables", "table_cog_all.html"))

table_cog_all %>% 
    as_flex_table() %>% 
    save_as_docx(path = here("outputs", "tables", "table_cog_all.docx"),
                 pr_section = sect_properties_landscape)

# 4. Regression Models ---------------------------------------------------------                                                                                                                                                                                                                                     
data_fep <- data %>% 
    filter(group == "FEP") %>% 
    mutate(group2 = droplevels(group2))

lm_res1 <- data_fep %>%
    filter(!is.na(CPZ_before_scan)) %>% 
    rename("naa" = "bg_naa") %>% 
    lm(Visual_Pattern_correct ~ Age + Gender + Years_of_education + DUP_months + aoo + naa * group2, data = .)
lm_res1_cpz <- data_fep %>% 
    rename("naa" = "bg_naa") %>% 
    lm(Visual_Pattern_correct ~ Age + Gender + Years_of_education + DUP_months + CPZ_before_scan + naa * group2, data = .)
broom::tidy(lm_res1) %>% mutate_if(is.numeric, round, 10)
broom::glance(lm_res1)
report::report(lm_res1)
table_lm_1_org <- lm_res1 %>% 
    tbl_regression(
        label = list(
            aoo ~ "Age of Onset",
            naa ~ "NAA",
            Gender ~ "Sex",
            Years_of_education ~ "Years of Education",
            DUP_months ~ "DUP (months)",
            group2 ~ "Onset Group"
        )
    ) %>% 
    bold_p()
table_lm_1_std <- lm_res1 %>% 
    tbl_regression(
        tidy_fun = tidy_standardize,
        label = list(
            aoo ~ "Age of Onset",
            naa ~ "NAA",
            Gender ~ "Sex",
            Years_of_education ~ "Years of Education",
            DUP_months ~ "DUP (months)",
            group2 ~ "Onset Group"
        )
    ) %>% 
    add_glance_table(include = c(r.squared, adj.r.squared))

table_lm_1 <- tbl_merge(list(table_lm_1_std, table_lm_1_org)) %>% 
    modify_header(estimate_1 ~ glue("**{rlang::expr('\U03B2')}**")) %>% 
    modify_column_hide(columns = c(estimate_2, ci_2)) %>% 
    modify_spanning_header(
        c(estimate_1, ci_1, p.value_2) ~ 
            "**Visual Pattern (Correct) ~ BG**")

lm_res2 <- data_fep %>%
    rename("naa" = "bg_naa") %>% 
    lm(WC_NP_error ~ Age + Gender + Years_of_education + DUP_months + aoo + naa * group2, data = .)
broom::tidy(lm_res2) %>% mutate_if(is.numeric, round, 10)
broom::glance(lm_res2)
report::report(lm_res2)
table_lm_2_org <- lm_res2 %>% 
    tbl_regression(
        label = list(
            aoo ~ "Age of Onset",
            naa ~ "NAA",
            Gender ~ "Sex",
            Years_of_education ~ "Years of Education",
            DUP_months ~ "DUP (months)",
            group2 ~ "Onset Group"
        )
    ) %>% 
    bold_p()

table_lm_2_std <- lm_res2 %>% 
    tbl_regression(
        tidy_fun = tidy_standardize,
        label = list(
            aoo ~ "Age of Onset",
            naa ~ "NAA",
            Gender ~ "Sex",
            Years_of_education ~ "Years of Education",
            DUP_months ~ "DUP (months)",
            group2 ~ "Onset Group"
        )
    ) %>% 
    add_glance_table(include = c(r.squared, adj.r.squared))

table_lm_2 <- tbl_merge(list(table_lm_2_std, table_lm_2_org)) %>% 
    modify_header(estimate_1 ~ glue("**{rlang::expr('\U03B2')}**")) %>% 
    modify_column_hide(columns = c(estimate_2, ci_2)) %>% 
    modify_spanning_header(
        c(estimate_1, ci_1, p.value_2) ~ 
            "**MWCST (Non-Perseverative Error) ~ BG**")

lm_res3 <- data_fep %>%
    rename("naa" = "acc_naa") %>% 
    lm(WC_NP_error ~ Age + Gender + Years_of_education + DUP_months + aoo + naa * group2, data = .)
broom::tidy(lm_res3) %>% mutate_if(is.numeric, round, 10)
broom::glance(lm_res3)
report::report(lm_res3)
table_lm_3_org <- lm_res3 %>% 
    tbl_regression(
        label = list(
            aoo ~ "Age of Onset",
            naa ~ "NAA",
            Gender ~ "Sex",
            Years_of_education ~ "Years of Education",
            DUP_months ~ "DUP (months)",
            group2 ~ "Onset Group"
        )
    ) %>% 
    bold_p()

table_lm_3_std <- lm_res3 %>% 
    tbl_regression(
        tidy_fun = tidy_standardize,
        label = list(
            aoo ~ "Age of Onset",
            naa ~ "NAA",
            Gender ~ "Sex",
            Years_of_education ~ "Years of Education",
            DUP_months ~ "DUP (months)",
            group2 ~ "Onset Group"
        )
    ) %>% 
    add_glance_table(include = c(r.squared, adj.r.squared))


table_lm_3 <- tbl_merge(list(table_lm_3_std, table_lm_3_org)) %>% 
    modify_header(estimate_1 ~ glue("**{rlang::expr('\U03B2')}**")) %>% 
    modify_column_hide(columns = c(estimate_2, ci_2)) %>% 
    modify_spanning_header(
        c(estimate_1, ci_1, p.value_2) ~ 
            "**MWCST (Non-Perseverative Error) ~ ACC**")

# merge table 1, 2 and 3
merged_table <- tbl_merge(
    list(table_lm_1, table_lm_2, table_lm_3),
    tab_spanner = c("**A. Visual Pattern (Correct) ~ BG**", 
                    "**B. MWCST (Non-Perseverative Error) ~ BG**",
                    "**C. MWCST (Non-Perseverative Error) ~ ACC**")
    ) %>% 
    modify_footnote(
        c("label") ~ 
        "DUP = Duration of untreated psychosis;
        AOS = Average-onset schizophrenia;
        LOS = Late-onset schizophrenia") %>% 
    modify_footnote(
        c("estimate_1_1", "estimate_1_2", "estimate_1_3") ~ "Standardized Betas"
    )

merged_table %>% 
    as_gt() %>%
    gtsave(filename = here("outputs", "tables", "regression_models_comparison.html"))

merged_table %>% 
    as_flex_table() %>% 
    autofit() %>% 
    bold(part = "header") %>% 
    save_as_docx(
        path = here("outputs", "tables", "regression_models_comparison.docx"),
        pr_section = sect_properties_landscape)


# only model 1 shows an interaction effect
# plot interaction effect 
interaction_effect <- lm_res1 %>% 
    sjPlot::plot_model(type = "pred", terms = c("naa", "group2")) +
        labs(x = "Basal Ganglia NAA Concentration Ratio",
             y = "Visual Pattern (Correct)",
             color = "",
             fill = "",
             title = "Interaction Effects") +
        scale_x_continuous(labels = scales::label_percent(accuracy = 0.01, scale = 1)) +
        scale_color_manual(values = c("#268785", "tomato3")) +
        scale_fill_manual(values = c("#268785", "tomato3")) +
        theme_pander() +
        theme(plot.margin = margin(2, 2, 2, 2, "mm"),
              legend.position = c(.875, .225),
              legend.background = element_rect(fill = "transparent"))

ggsave(filename = here("outputs", "figs", "interaction_effect.pdf"),
       plot = interaction_effect,
       width = 6, 
       height = 4.5)
