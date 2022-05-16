cohen_d <- function(data, variable, by, ...) {
    d <- effsize::cohen.d(data[[variable]] ~ as.factor(data[[by]]), 
                     conf.level = .95, pooled = TRUE, paired = FALSE, 
                     hedges.correction = TRUE)
    
    # Formatting statistic with CI
    est <- style_sigfig(d$estimate)
    ci <- style_sigfig(d$conf.int) %>% paste(collapse = ", ")
    
    # returning estimate with CI together
    str_glue("{est} ({ci})")
}

add_stat_pairwise <- function(data, variable, by, ...) {
    # calculate pairwise p-values
    pw <- pairwise.t.test(data[[variable]], data[[by]], p.adj = "none")
    
    # convert p-values to list
    index <- 0L
    p.value.list <- list()
    for (i in seq_len(nrow(pw$p.value))) {
        for (j in seq_len(nrow(pw$p.value))) {
            index <- index + 1L
            
            p.value.list[[index]] <- 
                c(pw$p.value[i, j]) %>%
                setNames(glue::glue("**{colnames(pw$p.value)[j]} vs. {rownames(pw$p.value)[i]}**"))
        }
    }
    
    # convert list to data frame
    p.value.list %>% 
        unlist() %>%
        purrr::discard(is.na) %>%
        t() %>%
        as.data.frame() %>%
        # formatting/roundign p-values
        dplyr::mutate(dplyr::across(everything(), style_pvalue))
}


format_decimals <- function(x, decimal_places = 3) {
    format(round(x, decimal_places), nsmall = decimal_places)
}
