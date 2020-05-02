### Figure 1



## Load needed libraries.
library(dplyr)
library(magrittr)


## Set global parameters.
stddev <- 2
alpha_touse <- 0.05
diff1 <- 12
diff2 <- 2
simpower <- 0.9


## Generate sample sizes.
maxdiff <- 13
diffinc <- 0.01
tbl_ss_by_diff <- dplyr::tibble(
  diff = seq(1, maxdiff, diffinc) %x% rep(1, 3),
  Power = rep(1, length(seq(1, maxdiff, diffinc))) %x% c(0.7, 0.8, 0.9),
  sampsizeraw = apply(
    X = as.matrix(cbind(diff, Power)),
    MARGIN = 1,
    FUN = function(x)
      pwr::pwr.t.test(
        d = x[1] / stddev,
        sig.level = alpha_touse / 2,
        power = x[2],
        type = "one.sample",
        alternative = "two.sided"
      )$n
  ),
  sampsize = ceiling(sampsizeraw)
)
sampsize1 <- (tbl_ss_by_diff %>% dplyr::filter(diff == diff1, Power == simpower))$sampsize
sampsize2 <- (tbl_ss_by_diff %>% dplyr::filter(diff == diff2, Power == simpower))$sampsize


## Generate significance curves.
sampsizes <- c(sampsize1, ceiling(mean(sampsize1, sampsize2)), sampsize2)
diffinc2 <- 0.1
tbl_statsig_by_diff <- dplyr::tibble(
  diff = seq(1, maxdiff, diffinc2) %x% rep(1, length(sampsizes)),
  sampsize = rep(1, length(seq(1, maxdiff, diffinc2))) %x% sampsizes,
  pvalue_alt = 2 * pt(
    diff / (stddev / sqrt(sampsize)),
    df = sampsize - 1,
    lower.tail = FALSE
  )
)
example_statsig <- dplyr::tibble(
  highdiff = diff1,
  highdiff_pvalue = unique(tbl_statsig_by_diff$pvalue_alt[tbl_statsig_by_diff$sampsize == min(sampsizes) & tbl_statsig_by_diff$diff == highdiff]),
  lowdiff = diff2,
  lowdiff_pvalue = unique(tbl_statsig_by_diff$pvalue_alt[tbl_statsig_by_diff$sampsize == max(sampsizes) & tbl_statsig_by_diff$diff == lowdiff])
)
ggp_statsig <- tbl_statsig_by_diff %>%
  dplyr::mutate(`Sample Size` = as.factor(sampsize)) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = diff,
      y = log(log(1/pvalue_alt)),
      group = `Sample Size`,
      color = `Sample Size`
      # linetype = `Sample Size`
    )
  ) +
  ggplot2::theme_classic() +
  ggplot2::scale_x_continuous(breaks = seq(1, maxdiff, 1)) +
  ggplot2::scale_y_continuous(
    limits = c(
      log(log(1/0.1)),
      log(log(1/0.000000001))
    ),
    breaks = c(
      log(log(1/alpha_touse)),
      log(log(1/0.01)),
      log(log(1/0.001)),
      log(log(1/0.0001))
    ),
    labels = c(
      as.character(alpha_touse),
      "0.01",
      "0.001",
      "0.0001"
    )
  ) +
  ggplot2::scale_color_manual(values = c("dark gray", "black")) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(yintercept = log(log(1/alpha_touse)), linetype = "dashed") +
  ggplot2::geom_hline(yintercept = log(log(1/0.01)), linetype = "dotted") +
  ggplot2::geom_hline(yintercept = log(log(1/0.001)), linetype = "dotted") +
  ggplot2::geom_hline(yintercept = log(log(1/0.0001)), linetype = "dotted") +
  ggplot2::geom_point(
    x = example_statsig$highdiff,
    y = log(log(1/example_statsig$highdiff_pvalue)),
    color = "dark gray"
  ) +
  ggplot2::geom_point(
    x = example_statsig$lowdiff,
    y = log(log(1/example_statsig$lowdiff_pvalue)),
    color = "black"
  ) +
  ggplot2::xlab("Estimated Reduction (Days)") +
  ggplot2::ylab("Statistical Significance (P-value)") +
  ggplot2::ggtitle(paste0("Figure 1. Statistical Significance by Estimated Reduction in Days to Recovery"))



### Figures A1 and A2. Scenarios 1 and 2


## Simulate RCTs.
sampnum <- 10000
sampnumplot <- sampnum

# scenario 1
seed1 <- 2004301758
tbl_pvalues1 <- dplyr::tibble(
  sampidx = seq(1, sampnum, 1),
  rsample = lapply(
    seq(1, sampnum, 1),
    function(x) {
      set.seed(seed1 + x)
      rnorm(n = sampsize1, mean = diff1, sd = stddev)
    }
  ),
  estimate = as.numeric(
    lapply(
      rsample,
      function(x) {
        (t.test(x, alternative = "two.sided"))$estimate
      }
    )
  ),
  pvalue = as.numeric(
    lapply(
      rsample,
      function(x) {
        (t.test(x, alternative = "two.sided"))$p.value
      }
    )
  ),
  statsig = (pvalue <= alpha_touse / 2)
) %>%
  dplyr::mutate(
    statsigprop = ifelse(
      statsig == TRUE,
      paste(statsig, "(", round(mean(statsig) * 100), "%)"),
      paste(statsig, "(", round(mean(1-statsig) * 100), "%)")
    )
  )
example1 <- tbl_pvalues1 %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::mutate(sampsize = sampsize1)
ggp_pvalues1 <- tbl_pvalues1 %>%
  dplyr::filter(sampidx <= sampnumplot) %>%
  dplyr::rename(`Statistically Significant` = statsigprop) %>%
  ggplot2::ggplot(ggplot2::aes(x = sampidx, y = log(pvalue), color = `Statistically Significant`)) +
  ggplot2::theme_classic() +
  ggplot2::scale_y_continuous(
    breaks = c(
      log(0.0001),
      log(0.001),
      log(0.01),
      log(alpha_touse / 2),
      log(alpha_touse)
    ),
    labels = c(
      "0.0001",
      "0.001",
      "0.01",
      as.character(alpha_touse / 2),
      as.character(alpha_touse)
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::geom_hline(yintercept = log(alpha_touse / 2)) +
  ggplot2::xlab("Sample Number") +
  ggplot2::ylab("P-value") +
  ggplot2::ggtitle(
    label = paste0("Scenario 1: P-values for ", sampnumplot, " RCTs (n = ", sampsize1, ")"),
    subtitle = paste0("True average reduction in days to recovery = ", diff1, " (std. dev. = ", stddev, "), power = ", round(simpower * 100), "%, alpha = ", alpha_touse)
  )

# scenario 2
seed2 <- round(2004301807 / 100)
tbl_pvalues2 <- dplyr::tibble(
  sampidx = seq(1, sampnum, 1),
  rsample = lapply(
    seq(1, sampnum, 1),
    function(x) {
      set.seed(seed2 + x)
      rnorm(n = sampsize2, mean = diff2, sd = stddev)
    }
  ),
  estimate = as.numeric(
    lapply(
      rsample,
      function(x) {
        (t.test(x, alternative = "two.sided"))$estimate
      }
    )
  ),
  pvalue = as.numeric(
    lapply(
      rsample,
      function(x) {
        (t.test(x, alternative = "two.sided"))$p.value
      }
    )
  ),
  statsig = (pvalue <= alpha_touse / 2)
) %>%
  dplyr::mutate(
    statsigprop = ifelse(
      statsig == TRUE,
      paste(statsig, "(", round(mean(statsig) * 100), "%)"),
      paste(statsig, "(", round(mean(1-statsig) * 100), "%)")
    )
  )
example2 <- tbl_pvalues2 %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::mutate(sampsize = sampsize2)
ggp_pvalues2 <- tbl_pvalues2 %>%
  dplyr::filter(sampidx <= sampnumplot) %>%
  dplyr::rename(`Statistically Significant` = statsigprop) %>%
  ggplot2::ggplot(ggplot2::aes(x = sampidx, y = log(pvalue), color = `Statistically Significant`)) +
  ggplot2::theme_classic() +
  ggplot2::scale_y_continuous(
    breaks = c(
      log(0.0001),
      log(0.001),
      log(0.01),
      log(alpha_touse / 2),
      log(alpha_touse)
    ),
    labels = c(
      "0.0001",
      "0.001",
      "0.01",
      as.character(alpha_touse / 2),
      as.character(alpha_touse)
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::geom_hline(yintercept = log(alpha_touse / 2)) +
  ggplot2::xlab("Sample Number") +
  ggplot2::ylab("P-value") +
  ggplot2::ggtitle(
    label = paste0("Scenario 2: P-values for ", sampnumplot, " RCTs (n = ", sampsize2, ")"),
    subtitle = paste0("True average reduction in days to recovery = ", diff2, " (std. dev. = ", stddev, "), power = ", round(simpower * 100), "%, alpha = ", alpha_touse)
  )