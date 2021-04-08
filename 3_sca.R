library(tidyverse)
library(haven)
library(broom)
library(glue)
library(tictoc)
library(magrittr)
library(furrr)
library(estimatr)
library(rsample)
library(here)
library(ggridges)

rm(list = ls())

# 1. Prepare Data ----
load("Data/df_full.Rdata")
load("Data/depvar_sd.Rdata")

scale_var <- function(var, df){
  mean <- depvar_sd$mean[depvar_sd$var == var]
  sd <- depvar_sd$sd[depvar_sd$var == var]
  x <- (df[[var]] - mean)/sd
}

df_sca <- df %>%
  filter(!is.na(Survey_Weight_W8))
df_sca[depvar_sd$var] <- map(depvar_sd$var, scale_var, df_sca)
save(df_sca, file = "Data/df_sca.Rdata")


# 2. Define Specifications ----
blocks <- list(
  MH = c("GHQ_W2", "GHQ_W4"),
  PH = c("GenHealth_W2", "GenHealth_W4", "Disabled"),
  SES = c("NSSEC3", "ParentEduc5", "SES"),
  LOC = c("LOC", "Int_LOC", "Ext_LOC"), 
  IMD = "IMD_W2",
  Dem = c("Female", "Ethnicity"),
  Educ = "Education_W8",
  Behav = c("Bullied", "SchoolAtt", "Risk")
)

get_vars <- function(stub){
  str_subset(names(df), paste0("^", stub))
}

v_combn <- function(block){
  df <- do.call("expand_grid", block)
  
  v_block <- df %>%
    rowwise() %>%
    group_split() %>%
    map(as.character)
  
  all_combn <- map(seq_along(block), ~ map(v_block, combn, .x, FUN = list)) %>%
    rlang::squash() %>%
    tibble() %>%
    unique() %>%
    pull(1)
}

combos <- blocks %>%
  map(~ structure(.x, names = .x)) %>%
  map(~ map(.x, get_vars)) %>%
  map(v_combn)

map_dbl(combos, length)
map_dbl(combos, length) %>%
  prod()

combos$MH <- Filter(function(x) 
  sum(str_count(x, "W2")) == 1 &
    sum(str_count(x, "Factor")) %in% c(0, length(x)),
  combos$MH)
combos$PH <- Filter(function(x)
  sum(str_count(x, "W2") + str_count(x, "Disabled")) == 2,
  combos$PH)
combos$SES <- Filter(function(x) 
  (length(x) == 2 & sum(str_count(x, "^SES")) == 0) |
    identical(x,"SES_W1"),
  combos$SES)
combos$LOC <- Filter(function(x) 
  (length(x) == 1 & sum(str_count(x, "^LOC"))==1) |
    (length(x) == 2 & sum(str_count(x, "^LOC"))==0 &
       sum(str_count(x, "Factor")) %in% c(0,2)),
  combos$LOC)
combos$Dem <- Filter(function(x) length(x)==2,
                     combos$Dem)
combos$Behav <- Filter(function(x) 
  length(x) == 3 &
    sum(str_count(x, "Factor")) %in% c(0,2),
  combos$Behav)

map_dbl(combos, length)
map_dbl(combos, length) %>%
  prod()

sca_specs <- map(combos, ~ map(.x,paste, collapse = " + ")) %>%
  do.call("expand_grid", .) %>%
  unnest(cols = everything()) %>%
  mutate(base_spec = glue("{MH} + {PH} + {Dem} + {SES} + {IMD} + {LOC}")) %>%
  uncount(3, .id = "spec_id") %>%
  mutate(spec = case_when(
    spec_id == 1 ~ glue("{base_spec} + {Educ} + {Behav}"),
    spec_id == 2 ~ glue("{base_spec} + {Educ}"),
    spec_id == 3 ~ glue("{base_spec} + {Behav}")
  )) %>%
  select(spec_id, spec) %>%
  expand_grid(depvar = map(c("GHQ_W8", "Patience", "Height"), get_vars) %>% unlist(),
              unem = get_vars("Unem")) %>%
  filter(!(str_detect(unem, "(FTE|2006)") & str_detect(spec, "_W4")),
         !(str_detect(unem, "(FTE|2006)") & str_detect(spec, "SES"))) %>%
  mutate(form = glue("{depvar} ~ {unem} + {spec}")) %>%
  select(- spec) %>%
  distinct()
save(sca_specs, file = "Data/sca_specs.Rdata")

sca_specs %>% 
  filter(str_detect(depvar, "GHQ")) %>%
  nrow()

set.seed(1)
sca_specs_n <- sca_specs %>%
  group_by(depvar) %>%
  sample_n(20000) %>%
  ungroup()
save(sca_specs_n, file = "Data/sca_specs_n.Rdata")


# 3. Run SCA ----
load("Data/df_sca.Rdata")
load("Data/sca_specs_n.Rdata")

get_lm <- function(form, unem){
  form <- as.formula(form)
  mod <- lm_robust(form, df_sca, 
                   weights = Survey_Weight_W8,
                   se_type = "stata")
  n_unem <- model.frame(form, df_sca)[[unem]] %>%
    table() %>% pluck(2)
  
  coef <- mod %>%
    broom::tidy(conf.int = TRUE) %>%
    select(term, beta = estimate,
           lci = conf.low, uci =  conf.high,
           p = p.value, t = statistic) %>%
    filter(str_detect(term, unem)) %>%
    mutate(n_unem = n_unem) %>%
    as_tibble()
}
save(get_lm, file = "Data/get_lm.Rdata")

set.seed(1)
tic()
plan(multisession)
sca_results <- sca_specs_n %>%
  mutate(mod = future_map2(form, unem, get_lm, .progress = TRUE)) %>%
  unnest(mod)
future:::ClusterRegistry("stop")
toc()
save(sca_results, file = "Data/sca_results.Rdata")


# 4. SCA Inferential Statistics ----
load("Data/df_sca.Rdata")
load("Data/sca_results.Rdata")

boots <- 500
boot_splits  <- bootstraps(df_sca, boots)

sca_inference <- function(boot_specs, positive = TRUE){
  
  sca_bootstrap <- function(boot){
    lm_bootstrap <- function(row, df_boot){
      spec <- boot_specs[row,]
      
      df_boot[[spec$depvar]] <- ifelse(df_boot[[spec$unem]] == levels(df_boot[[spec$unem]])[2], 
                                       df_boot[[spec$depvar]] - spec$beta, df_boot[[spec$depvar]])
      
      form <- as.formula(spec$form)
      mod <- lm_robust(form, df_boot, 
                       weights = Survey_Weight_W8,
                       se_type = "stata") %>%
        broom::tidy() %>%
        filter(str_detect(term, "^Unem_")) %>%
        select(beta = estimate, p = p.value, t = statistic)
    }
    
    df_boot <- analysis(boot_splits$splits[[boot]])
    bootstrap_results <- map_dfr(1:nrow(boot_specs), lm_bootstrap, df_boot)
    if (positive == TRUE){
      bootstrap_results <- bootstrap_results %>%
        mutate(signif = ifelse(p < 0.05 & beta > 0, TRUE, FALSE))
    } else{
      bootstrap_results <- bootstrap_results %>%
        mutate(signif = ifelse(p < 0.05 & beta < 0, TRUE, FALSE))
    }
    bootstrap_results <- bootstrap_results %>%
      summarise(median = median(beta),
                signif = mean(signif),
                t = mean(t))
  }
  
  bootstraps_stats <- function(bootstraps_results){
    if (positive == TRUE){
      sca_stats <- boot_specs %>%
        mutate(signif = ifelse(p < 0.05 & beta > 0, TRUE, FALSE))
    } else{
      sca_stats <- boot_specs %>%
        mutate(signif = ifelse(p < 0.05 & beta < 0, TRUE, FALSE))
    }
    sca_stats <- sca_stats %>%
      summarise(median = median(beta),
                signif = mean(signif),
                t = mean(t)) %>%
      unlist()
    
    if (positive == TRUE){
      bootstrap_stats <- bootstraps_results %>%
        mutate(median = ifelse(median > sca_stats["median"], TRUE, FALSE),
               t = ifelse(t > sca_stats["t"], TRUE, FALSE))
    } else{
      bootstrap_stats <- bootstraps_results %>%
        mutate(median = ifelse(median < sca_stats["median"], TRUE, FALSE),
               t = ifelse(t < sca_stats["t"], TRUE, FALSE))
    }
    bootstrap_stats <- bootstrap_stats %>%
      mutate(signif = ifelse(signif > sca_stats["signif"], TRUE, FALSE)) %>%
      summarise_all(mean)
  }
  
  
  sca_p <- map_dfr(1:nrow(boot_splits), sca_bootstrap) %>%
    bootstraps_stats()
  
  sca_p
}

sample_sca <- function(depvar, seed){
  set.seed(seed)
  sca_results %>%
    filter(str_detect(depvar, !!depvar),
           spec_id == 1) %>%
    sample_n(1000)
}

sca_p <- list()
tic()
sca_p$ghq <- sample_sca("GHQ", 2) %>%
  sca_inference()
toc()

save(sca_p, file = "Data/sca_inference.Rdata")

load("Data/sca_inference.Rdata")
tic()
sca_p$height <- sample_sca("Height", 3) %>%
  sca_inference(positive = FALSE)
sca_p$patience <- sample_sca("Patience", 4) %>%
  sca_inference(positive = FALSE)
toc()

sca_p <- bind_rows(sca_p, .id = "depvar_stub")

save(sca_p, file = "Data/sca_inference.Rdata")


# 5. Plots & Writing ----
load("Data/sca_inference.Rdata")
load("Data/sca_results.Rdata")

get_dense <- function(df){
  df %>%
    group_by(depvar_stub) %>%
    mutate(dense = dense_rank(beta)) %>%
    group_by(depvar_stub, spec_id) %>%
    mutate(spec_dense = percent_rank(beta)) %>%
    ungroup()
}

sca_results <- sca_results %>%
  mutate(unem = str_replace(unem, "_6Months", "_6_Months_Cont"),
         unem = str_replace(unem, "^Unem_", ""),
         unem = str_replace(unem, "m_", "_"),
         signif = ifelse(p <= 0.05, "p <= 0.05", "p > 0.05")) %>%
  separate(depvar, into = c("depvar_stub", "depvar_wave", "depvar_type"),
           fill = "right", extra = "merge", remove = FALSE) %>%
  separate(depvar_type, into = c("depvar_score", "depvar_combine"),
           fill = "right", extra = "merge", remove = FALSE) %>%
  mutate(depvar_combine = ifelse(is.na(depvar_combine) & depvar_stub == "GHQ",
                                 "Score", depvar_combine)) %>%
  separate(unem, into = c("unem_months", "unem_range",  "unem_type", "unem_dur",
                          "unem_start", "unem_end"), fill = "right", convert = TRUE) %>%
  get_dense()

# Plot
sca_plot <- sca_results %>%
  ggplot() +
  aes(x = dense, y = beta, color = signif) +
  facet_wrap(~ depvar_stub, scales = "free_x") +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.3) +
  scale_x_continuous(labels = scales::comma) +
  scale_color_manual(values = c("#1b9e77", "#d95f02")) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  theme_minimal() +
  theme(legend.position = c(0.1, 0.85),
        axis.text.x = element_text(angle = 45),
        panel.spacing.x = unit(2, "lines")) +
  labs(x = "Rank", y = "Standardized Coefficient", color = NULL)
ggsave(filename = "Images/figure_1.png", plot = sca_plot,
       width = 21, height = 9.9, units = "cm", dpi = 800)  

specs <- c("+ Full Controls", "excl. Behaviour", "excl. Education")

sca_spec <- sca_results %>%
  filter(depvar_stub == "GHQ") %>%
  mutate(spec = factor(specs[spec_id], specs)) %>%
  ggplot() +
  aes(x = spec_dense, y = beta, color = spec) +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.3) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_brewer(palette = "Paired") +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  theme_minimal() +
  theme(legend.position = c(0.1, 0.85),
        axis.text.x = element_text(angle = 45),
        panel.spacing.x = unit(2, "lines")) +
  labs(x = "% Rank", y = "Standardized Coefficient", color = NULL)
ggsave(filename = "Images/figure_2.png", plot = sca_spec,
       width = 21, height = 9.9, units = "cm", dpi = 800)  

sca_results %>%
  filter(depvar_stub == "GHQ") %>%
  group_by(spec_id) %>%
  summarise(median = median(beta),
            mean = mean(beta),
            .groups = "drop")

# Paragraph
get_res <- function(var){
  var %>%
    glue(" ~ Unem_6Months") %>%
    as.formula() %>%
    lm_robust(df_sca, Survey_Weight_W8) %>%
    tidy() %>%
    filter(str_detect(term, "Unem")) %>%
    select(estimate, conf.low, conf.high) %>%
    as_tibble() %>%
    mutate(across(everything(), round, 2),
           string = glue("b = {estimate}, 95% CI = {conf.low}, {conf.high}")) %>%
    pull(string)
}
res_bivar <- tibble(outcome = c("GHQ_W8_Likert", "Patience_W8", "Height_W8")) %>%
  mutate(res = map_chr(outcome, get_res)) %>%
  deframe()


res_med <- sca_results %>%
  group_by(depvar_stub) %>%
  summarise(beta = median(beta),
            .groups = "drop") %>%
  deframe() %>%
  round(2)
res_med

res_signif <- sca_results %>%
  group_by(depvar_stub) %>%
  summarise(prop = 100*sum(p <= 0.05)/n(),
            .groups = "drop") %>%
  deframe() %>%
  round(2)
res_signif

res_neg <- sca_results %>%
  group_by(depvar_stub) %>%
  summarise(prop = 100*sum(beta < 0)/n(),
            .groups = "drop") %>%
  deframe() %>%
  round(2)
res_neg

sca_p

paragraph <- "
Using our preferred operationalisations there was a bivariate association
 between youth unemployment and GHQ scores ({res_bivar['GHQ_W8_Likert']})
 and patience ({res_bivar['Patience_W8']}) but not height ({res_bivar['Height_W8']})
 suggesting that the latter is not a good placebo outcome.
 Figure 1 presents the range of (standardized) estimates across model specifications.
 The median effect size for GHQ scores was {res_med['GHQ']} SD.
 {res_signif['GHQ']}% of specifications were statistically significant
 and only {res_neg['GHQ']}% of specifications predicted better mental health among the youth unemployed.
 None of the 500 bootstrap samples produced larger median effect sizes,
 a higher proportion of significant results,
 or larger average z-statistics.
 {res_signif['Height']}% and {res_signif['Patience']}% of specifications reached statistical significance for height and patience, respectively,
 and the median effect sizes were {res_med['Height']} and {res_med['Patience']} SD.
 A small number of specifications reached substantial effect sizes, however.
 More than half of the bootstrap samples produced more extreme results than the original sample.
" %>%
  glue() %>%
  str_replace_all("\n", "")
paragraph


# Unemployment Length Plots
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- sca_results %>%
  filter(str_detect(depvar, "GHQ_W8")) %>%
  mutate(unem_months = glue("{unem_months}+") %>% 
           fct_reorder(unem_months)) %>%
  ggplot() +
  aes(x = beta, y = unem_months) +
  geom_density_ridges(color = "grey40", fill = cbbPalette[3]) +
  geom_vline(xintercept = 0) +
  theme_minimal() +
  labs(x = "Standardized Coefficient",
       y = "Unemployment Duration (Months)")
ggsave(filename = "Images/sca_unem.png", plot = p,
       width = 21, height = 9.9,
       units = "cm", dpi = 300)  
