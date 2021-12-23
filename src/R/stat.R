cohen_d <- function(data, variable, by, ...) {
    effsize::cohen.d(data[[variable]] ~ as.factor(data[[by]]), 
                     conf.level = .95, pooled = TRUE, paired = FALSE, 
                     hedges.correction = TRUE)$estimate
}