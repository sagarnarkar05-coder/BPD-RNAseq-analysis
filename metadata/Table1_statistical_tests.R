# Table 1 demographic analysis

library(readxl)

# File paths

stats_dir <- "/home/sagar007/RNAseq/Stats"

metadata_file <- file.path(
  stats_dir,
  "Sample_metadata.xlsx"
)

# Read sample metadata

sample_data <- read_excel(
  metadata_file,
  sheet = "Sample metadata"
)

infant_data <- read_excel(
  metadata_file,
  sheet = "Table1 infant data"
)

# Group order

infant_data$Group <- factor(
  infant_data$Group,
  levels = c("BPD", "nonBPD")
)

# Convert binary variables

to_binary <- function(x) {

  x <- trimws(tolower(as.character(x)))

  out <- rep(NA, length(x))

  out[x %in% c("true", "1", "yes", "y")] <- TRUE
  out[x %in% c("false", "0", "no", "n")] <- FALSE

  out
}

infant_data$Fetal_IUGR <- to_binary(
  infant_data$Fetal_IUGR
)

infant_data$Maternal_antenatal_steroid <- to_binary(
  infant_data$Maternal_antenatal_steroid
)

infant_data$Cesarean <- to_binary(
  infant_data$Cesarean
)

infant_data$Chorioamnionitis <- to_binary(
  infant_data$Chorioamnionitis
)

# Welch t-tests

GA_test <- t.test(
  GA_wk ~ Group,
  data = infant_data,
  var.equal = FALSE
)

BW_test <- t.test(
  BW_g ~ Group,
  data = infant_data,
  var.equal = FALSE
)

# Fisher exact test

run_fisher <- function(variable, group) {

  keep <- !is.na(variable) & !is.na(group)

  variable <- variable[keep]
  group <- group[keep]

  if (length(unique(variable)) < 2) {
    return(NA_real_)
  }

  fisher.test(
    table(group, variable)
  )$p.value
}

IUGR_p <- run_fisher(
  infant_data$Fetal_IUGR,
  infant_data$Group
)

Steroid_p <- run_fisher(
  infant_data$Maternal_antenatal_steroid,
  infant_data$Group
)

Cesarean_p <- run_fisher(
  infant_data$Cesarean,
  infant_data$Group
)

Chorio_p <- run_fisher(
  infant_data$Chorioamnionitis,
  infant_data$Group
)

# Summary functions

mean_sd <- function(x, digits = 1) {

  paste0(
    round(mean(x, na.rm = TRUE), digits),
    " ± ",
    round(sd(x, na.rm = TRUE), digits)
  )
}

binary_summary <- function(x) {

  x <- x[!is.na(x)]

  yes <- sum(x)
  total <- length(x)

  paste0(
    yes,
    "/",
    total,
    " (",
    round(100 * yes / total, 1),
    "%)"
  )
}

format_p <- function(p) {

  if (is.na(p)) {
    return("—")
  }

  if (p < 0.001) {
    return("<0.001")
  }

  sprintf("%.3f", p)
}

# Demographic table

bpd <- infant_data$Group == "BPD"
nonbpd <- infant_data$Group == "nonBPD"

table1 <- data.frame(

  Characteristic = c(
    "Infants, n",
    "RNA-seq samples, n",
    "Gestational age, weeks",
    "Birth weight, g",
    "Fetal IUGR, n/N (%)",
    "Maternal antenatal steroid, n/N (%)",
    "Cesarean delivery, n/N (%)",
    "Chorioamnionitis, n/N (%)"
  ),

  BPD = c(
    sum(bpd),
    sum(sample_data$Group == "BPD"),
    mean_sd(infant_data$GA_wk[bpd]),
    mean_sd(infant_data$BW_g[bpd], 0),
    binary_summary(infant_data$Fetal_IUGR[bpd]),
    binary_summary(infant_data$Maternal_antenatal_steroid[bpd]),
    binary_summary(infant_data$Cesarean[bpd]),
    binary_summary(infant_data$Chorioamnionitis[bpd])
  ),

  nonBPD = c(
    sum(nonbpd),
    sum(sample_data$Group == "nonBPD"),
    mean_sd(infant_data$GA_wk[nonbpd]),
    mean_sd(infant_data$BW_g[nonbpd], 0),
    binary_summary(infant_data$Fetal_IUGR[nonbpd]),
    binary_summary(infant_data$Maternal_antenatal_steroid[nonbpd]),
    binary_summary(infant_data$Cesarean[nonbpd]),
    binary_summary(infant_data$Chorioamnionitis[nonbpd])
  ),

  Test = c(
    "—",
    "—",
    "Welch t-test",
    "Welch t-test",
    "Fisher exact test",
    "Fisher exact test",
    "Not tested",
    "Not tested"
  ),

  P_value = c(
    "—",
    "—",
    format_p(GA_test$p.value),
    format_p(BW_test$p.value),
    format_p(IUGR_p),
    format_p(Steroid_p),
    format_p(Cesarean_p),
    format_p(Chorio_p)
  ),

  stringsAsFactors = FALSE
)

# Print table

print(
  table1,
  row.names = FALSE
)

# Save table

write.csv(
  table1,
  file.path(
    stats_dir,
    "Table1_demographic_statistics.csv"
  ),
  row.names = FALSE
)

# Save exact test results

test_results <- data.frame(
  Variable = c(
    "Gestational age",
    "Birth weight",
    "Fetal IUGR",
    "Maternal antenatal steroid",
    "Cesarean delivery",
    "Chorioamnionitis"
  ),
  Test = c(
    "Welch t-test",
    "Welch t-test",
    "Fisher exact test",
    "Fisher exact test",
    "Not tested",
    "Not tested"
  ),
  P_value = c(
    GA_test$p.value,
    BW_test$p.value,
    IUGR_p,
    Steroid_p,
    Cesarean_p,
    Chorio_p
  )
)

write.csv(
  test_results,
  file.path(
    stats_dir,
    "Table1_test_results.csv"
  ),
  row.names = FALSE
)
