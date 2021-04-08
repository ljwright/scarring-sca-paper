library(tidyverse)
library(haven)
library(glue)
library(lavaan)
library(Hmisc)
library(officer)
library(flextable)
library(broom)
library(here)
library(estimatr)
# library(psych)

rm(list = ls())

# Load Data
df <- read_dta("Data/Dataset.dta") %>%
  as_factor() %>%
  zap_formats()

mod_lookup <- c(
  single = "One Factor", twofactor = "Two Factor", 
  graetz = "Graetz (1991)", hankins = "Hankins (2008)", 
  hankins2 = "Hankins (2008)\n(External-wording covariance)",
  molina = "Molina et al. (2014)", rodrigo = "Rodrigo et al. (2019)"
)

# 1. GHQ ----
df_ghq <- df %>%
  select(NSID, contains("Weight"), matches("^W\\d_GHQ"))

run_cfa <- function(w, t){
  df_analysis <- df_ghq %>%
    select(NSID, contains(glue("w{w}_ghq_{t}"))) %>%
    drop_na()
  
  neg_1 <- paste0("W", w,"_GHQ_", t,"_", c(2, 5, 6, 9))
  neg_2 <- paste0("W", w,"_GHQ_", t,"_", c(10, 11))
  pos <- paste0("W", w,"_GHQ_", t,"_", c(1, 3, 4, 7, 8, 12))
  neg <- c(neg_1, neg_2)
  items <- c(neg, pos)
  
  neg_cov <- combn(neg, 2, simplify = FALSE) %>%
    map_chr(paste, collapse = " ~~ ") %>%
    paste(collapse = "\n")
  
  mods <- list(
    single = glue("MH =~ {paste(items, collapse = ' + ')}"),
    hankins = glue("MH =~ {paste(items, collapse = ' + ')}; {neg_cov}"),
    molina =  glue("MH =~ {paste(items, collapse = ' + ')}
NEG =~ {paste(neg, collapse = ' + ')}; MH ~~ 0*NEG"),
    rodrigo = glue("MH =~ {paste(items, collapse = ' + ')}
NEG =~ {paste(neg, collapse = ' + ')}
POS =~ {paste(pos, collapse = ' + ')}
MH ~~ 0*NEG; MH ~~ 0*POS; NEG ~~ 0*POS"),
    graetz = glue("NEG1 =~ {paste(neg_1, collapse = ' + ')}
NEG2 =~ {paste(neg_2, collapse = ' + ')}
POS =~ {paste(pos, collapse = ' + ')}"),
    twofactor = glue("NEG =~ {paste(neg, collapse = ' + ')}
POS =~ {paste(pos, collapse = ' + ')}")
    ) %>%
    map(cfa, data = df_analysis, estimator = "DWLS")
  
  nice <- function(x) round(x, 3)

  fit_stats <- function(mod){
    fit <- summary(mod, fit.measures = TRUE) %>%
      pluck("FIT")

    if (!is.null(fit)){
      res <- tibble(measure = c("RMSEA", "CFI"),
                    beta = c(fit[['rmsea']], fit[['cfi']]),
                    lci = c(fit[['rmsea.ci.lower']], NA),
                    uci = c(fit[['rmsea.ci.upper']], NA)) %>%
        mutate(string = ifelse(measure == "RMSEA",
                               glue("{nice(beta)} ({nice(lci)}, {nice(uci)})"),
                               glue("{nice(beta)}")))
    } else(
      res <- tibble(measure = c("RMSEA", "CFI"),
                    beta = NA, lci = NA, uci = NA)
    )
    res %>%
      mutate(w = w, t = t)
  }

  stats <- map_dfr(mods, fit_stats, .id = "mod")

  get_latent <- function(mod){
   lavPredict(mod) %>%
      as.data.frame() %>% as_tibble() %>%
      bind_cols(select(df_analysis, "NSID")) %>%
      mutate(w = w, t = t)
  }

  pred <- map_dfr(mods, get_latent, .id = "mod") %>%
    select(NSID, mod, w, t, MH) %>%
    filter(!is.na(MH))

  list(stats = stats, pred = pred)
  
}

ghq_cfa <- list()
ghq_cfa$mods <- map2(rep(c(2,4,8), each = 3),
             rep(c("Caseness", "Corrected", "Likert"), 3),
             run_cfa)

ghq_cfa$stats <- map_dfr(ghq_cfa$mods, ~ .x$stats) %>%
  mutate(mod_clean = factor(mod_lookup[mod], levels = mod_lookup))

ghq_cfa$tbl <- ghq_cfa$stats %>%
  mutate(w = case_when(w == 2 ~ "Age 14/15", w == 4 ~ "Age 16/17", w == 8 ~ "Age 25")) %>%
  select(w, mod_clean, t, measure, string) %>%
  pivot_wider(names_from = c("t", "measure"), values_from = "string") %>%
  arrange(w, mod_clean) %>%
  flextable() %>%
  merge_v(j = ~ w) %>%
  set_header_labels(w = "Age", mod_clean = "Model",
                    Caseness_RMSEA = "RMSEA (95% CI)", 
                    Caseness_CFI = "CFI",
                    Corrected_RMSEA = "RMSEA (95% CI)",
                    Corrected_CFI = "CFI",
                    Likert_RMSEA = "RMSEA (95% CI)", 
                    Likert_CFI = "CFI") %>%
  border_remove() %>%
  add_header(w = "", mod_clean = "", 
             Caseness_RMSEA = "Caseness Scoring", 
             Caseness_CFI = "Caseness Scoring",
             Corrected_RMSEA = "Corrected Scoring",
             Corrected_CFI = "Corrected Scoring",
             Likert_RMSEA = "Likert Scoring", 
             Likert_CFI = "Likert Scoring") %>%
  fontsize(size = 11, part = "all") %>%
  merge_h(part = "header") %>%
  border_inner_h(border = fp_border(color="gray30", style = "dashed")) %>%
  hline_top(border = fp_border(color="black", width = 2), part = "all") %>%
  hline_bottom(border = fp_border(color="black", width = 2), part = "all") %>%
  fix_border_issues(part = "all") %>%
  align(j=1:2, align="right", part = "all") %>%
  align(j=3:8, align="center", part = "all") %>%
  valign(j = 1, valign = "top") %>%
  autofit()
ghq_cfa$tbl
save_as_docx(ghq_cfa$tbl, path = "Tables/ghq_cfa.docx")


ghq_cfa$weights <- df_ghq %>%
  select(NSID, contains("Weight")) %>%
  pivot_longer(cols = -NSID, names_to = "w", values_to = "weight") %>%
  mutate(w = str_sub(w, -1, -1) %>% as.numeric()) %>%
  filter(!is.na(w))

df_ghq_latent <- map_dfr(ghq_cfa$mods, ~ .x$pred) %>%
  filter(mod == "rodrigo") %>%
  left_join(ghq_cfa$weights, by = c("NSID", "w"))  %>%
  mutate(name = glue("GHQ_W{w}_{str_to_sentence(t)}_Factor")) %>%
  group_by(name) %>%
  mutate(mean = wtd.mean(MH, weights = weight, normwt = TRUE),
         var = wtd.var(MH, weights = weight, normwt = TRUE),
         sd = sqrt(var),
         MH = (MH - mean)/sd) %>%
  ungroup() %>%
  select(NSID, name, MH) %>%
  pivot_wider(names_from = "name", values_from = "MH") %>%
  full_join(select(df, NSID), by = "NSID") %>%
  select(-matches("W2_Likert"))
save(df_ghq_latent, file = "Data/df_ghq_latent.Rdata")


get_alpha <- function(w, t){
  reg <- glue("^W{w}(.*){t}")
  a <- df_ghq %>%
    select(matches(reg)) %>%
    psych::alpha() %>% pluck("total") %>%
    pluck("raw_alpha")
  tibble(w = w, t = t, alpha = a)
}
ghq_alpha <- map2_dfr(rep(c(2,4,8), each = 3),
                      rep(c("Caseness", "Corrected", "Likert"), 3),
                      get_alpha) %>%
  mutate(w = case_when(w == 2 ~ "Age 14/15", w == 4 ~ "Age 16/17", w == 8 ~ "Age 25")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  pivot_wider(names_from = t, values_from = alpha) %>%
  flextable() %>%
  set_header_labels(w = "Age") %>%
  border_remove() %>%
  fontsize(size = 11, part = "all") %>%
  merge_h(part = "header") %>%
  border_inner_h(border = fp_border(color="gray30", style = "dashed")) %>%
  hline_top(border = fp_border(color="black", width = 2), part = "all") %>%
  hline_bottom(border = fp_border(color="black", width = 2), part = "all") %>%
  fix_border_issues(part = "all") %>%
  align(j=1, align="right", part = "all") %>%
  align(j=2:3, align="center", part = "all") %>%
  valign(j = 1, valign = "top") %>%
  autofit()
ghq_alpha
save_as_docx(ghq_alpha, path = "Tables/ghq_alpha.docx")

# 2. Locus of Control ----
# Confirmatory Factor Analysis
df_loc <- df %>%
  select(NSID, Survey_Weight_W2, matches("LOC(.*)Item")) %>%
  drop_na()

int <- paste0("Int_LOC_W2_Item", 1:3)
ext <- paste0("Ext_LOC_W2_Item", 1:3)
items <- c(ext, int)

ext_cov <- combn(ext, 2, simplify = FALSE) %>%
  map_chr(paste, collapse = " ~~ ") %>%
  paste(collapse = "\n")
int_cov <- combn(int, 2, simplify = FALSE) %>%
  map_chr(paste, collapse = " ~~ ") %>%
  paste(collapse = "\n")

loc_cfa <- list()
loc_cfa$mods <- list(
    single = glue("LOC =~ {paste(items, collapse = ' + ')}"),
    hankins = glue("LOC =~ {paste(items, collapse = ' + ')}; {int_cov}"),
    hankins2 = glue("LOC =~ {paste(items, collapse = ' + ')}; {ext_cov}"),
    twofactor = glue("INT =~ {paste(int, collapse = ' + ')}; EXT =~ {paste(ext, collapse = ' + ')}")
    ) %>%
    map(cfa, data = df_loc, estimator = "DWLS")

fit_stats <- function(mod){
  nice <- function(x) round(x, 3)
  
  fit <- summary(mod, fit.measures = TRUE) %>%
    pluck("FIT")
  
  if (!is.null(fit)){
    res <- tibble(measure = c("RMSEA", "CFI"),
                  beta = c(fit[['rmsea']], fit[['cfi']]),
                  lci = c(fit[['rmsea.ci.lower']], NA),
                  uci = c(fit[['rmsea.ci.upper']], NA)) %>%
      mutate(string = ifelse(measure == "RMSEA",
                             glue("{nice(beta)} ({nice(lci)}, {nice(uci)})"),
                             glue("{nice(beta)}")))
  } else(
    res <- tibble(measure = c("RMSEA", "CFI"),
                  beta = NA, lci = NA, uci = NA)
  )
}

loc_cfa$tbl <- map_dfr(loc_cfa$mods, fit_stats, .id = "mod") %>%
  mutate(mod_clean = factor(mod_lookup[mod], levels = mod_lookup)) %>%
  select(mod_clean, measure, string) %>%
  pivot_wider(names_from = "measure", values_from = "string") %>%
  arrange(mod_clean) %>%
  flextable() %>%
  set_header_labels(mod_clean = "Model",
                    RMSEA = "RMSEA (95% CI)") %>%
  border_remove() %>%
  fontsize(size = 11, part = "all") %>%
  hline_top(border = fp_border(color="black", width = 2), part = "all") %>%
  hline_bottom(border = fp_border(color="black", width = 2), part = "all") %>%
  fix_border_issues(part = "all") %>%
  align(j=1, align="right", part = "all") %>%
  align(j=2:3, align="center", part = "all") %>%
  valign(j = 1, valign = "top") %>%
  autofit()
loc_cfa$tbl
save_as_docx(loc_cfa$tbl, path = "Tables/loc_cfa.docx")

get_latent <- function(mod){
  lavPredict(mod) %>%
    as.data.frame() %>% as_tibble() %>%
    bind_cols(select(df_loc, "NSID"))
}

df_loc_latent <- map_dfr(loc_cfa$mods, get_latent, .id = "mod") %>%
  pivot_longer(cols = -c(mod, NSID), names_to = "factor") %>%
  filter(!(mod %in% c("single", "hankins2")), !is.na(value))  %>%
  left_join(select(df_loc, NSID, Survey_Weight_W2), by = c("NSID"))  %>% 
  mutate(name = case_when(factor == "LOC" ~ "LOC_W2_Factor",
                          factor == "INT" ~ "Int_LOC_W2_Factor",
                          factor == "EXT" ~ "Ext_LOC_W2_Factor")) %>%
  group_by(name) %>%
  mutate(mean = wtd.mean(value, weights = Survey_Weight_W2, normwt = TRUE),
         var = wtd.var(value, weights = Survey_Weight_W2, normwt = TRUE),
         sd = sqrt(var),
         value = (value - mean)/sd) %>%
  ungroup() %>%
  select(NSID, Survey_Weight_W2, name, value) %>%
  pivot_wider(names_from = "name", values_from = "value")

wtd_quantiles <- function(x, weight){
  y <- wtd.quantile(x, weights = weight, 
                    normwt = TRUE, probs = c(0.25, 0.5, 0.75))
  tibble(p25 = y[1], p50 = y[2], p75 = y[3])
}
qtile <- df_loc_latent %>%
  select(Survey_Weight_W2, Ext_LOC_W2_Factor, Int_LOC_W2_Factor) %>%
  pivot_longer(cols = c(Ext_LOC_W2_Factor, Int_LOC_W2_Factor)) %>%
  group_by(name) %>%
  summarise(qtile = list(wtd_quantiles(value, Survey_Weight_W2))) %>%
  unnest(qtile) %>% 
  arrange(name)

df_loc_latent <- df_loc_latent %>%
  mutate(
    LOC_W2_Factor_50 = case_when(
      Ext_LOC_W2_Factor > qtile$p50[1] & Int_LOC_W2_Factor > qtile$p50[2] ~ "Internal",
      Ext_LOC_W2_Factor < qtile$p50[1] & Int_LOC_W2_Factor < qtile$p50[2] ~ "External",
      TRUE ~ "Neither") %>% factor(c("External", "Neither", "Internal")),
    LOC_W2_Factor_75 = case_when(
      Ext_LOC_W2_Factor > qtile$p75[1] & Int_LOC_W2_Factor > qtile$p75[2] ~ "Internal",
      Ext_LOC_W2_Factor < qtile$p25[1] & Int_LOC_W2_Factor < qtile$p25[2] ~ "External",
      TRUE ~ "Neither") %>% factor(c("External", "Neither", "Internal"))
  ) %>%
  select(-Survey_Weight_W2) %>%
  full_join(select(df, NSID), by = "NSID")
save(df_loc_latent, file = "Data/df_loc_latent.Rdata")

tidy_cfa <- function(mod){
 mod %>%
    tidy() %>% 
    filter(str_detect(op, "=~")) %>%
    separate(term, into = c("latent", "op", "item"), sep = " ") %>%
    select(latent, item, std.all)
}

loc_lab <- select(df_loc, matches("LOC")) %>%
  map_chr(attr, "label")

loc_cfa$tbl_2 <- map_dfr(loc_cfa$mods, tidy_cfa, .id = "mod") %>%
  filter(mod != "hankins2") %>%
  mutate(std.all = round(std.all, 3),
         label = loc_lab[item], #%>% fct_rev(),
         mod_clean = mod_lookup[mod]) %>%
  select(label, std.all, mod_clean) %>%
  pivot_wider(names_from = "mod_clean", values_from = "std.all") %>%
  arrange(label) %>%
  flextable() %>%
  set_header_labels(label = "Item") %>%
  border_remove() %>%
  fontsize(size = 11, part = "all") %>%
  hline_top(border = fp_border(color="black", width = 2), part = "all") %>%
  hline_bottom(border = fp_border(color="black", width = 2), part = "all") %>%
  fix_border_issues(part = "all") %>%
  align(j=1, align="right", part = "all") %>%
  align(j=2:4, align="center", part = "all") %>%
  valign(j = 1, valign = "top") %>%
  autofit()
loc_cfa$tbl_2
save_as_docx(loc_cfa$tbl_2, path = "Tables/loc_loadings.docx")

# Cronbach's Alpha
get_alpha <- function(vars){
  v <- paste(vars, collapse = " ")
  a <- df_loc %>%
    select(all_of(vars)) %>%
    psych::alpha() %>% pluck("total") %>%
    pluck("raw_alpha")
  tibble(alpha = a, len = length(vars), items = v)
}
loc_alpha <- list()
loc_alpha$df <- map(3:6, ~ combn(items, .x, simplify = FALSE)) %>%
  flatten() %>%
  map_dfr(get_alpha) %>%
  arrange(desc(alpha))
loc_alpha$tbl <- loc_alpha$df %>%
  mutate(n = row_number(),
         observed = "x",
         items = str_split(items, " "),
         alpha = round(alpha, 2)) %>%
  unnest(items) %>%
  mutate(q = loc_lab[items],
         q = str_sub(q, 1, 2)) %>%
  pivot_wider(id_cols = c(n, alpha), names_from = q, 
              values_from = observed, values_fill = "") %>%
  select(alpha, Q1, Q2, Q3, Q4, Q5, Q6) %>%
  flextable() %>%
  set_header_labels(alpha = "Alpha") %>%
  border_remove() %>%
  fontsize(size = 11, part = "all") %>%
  hline_top(border = fp_border(color="black", width = 2), part = "all") %>%
  hline_bottom(border = fp_border(color="black", width = 2), part = "all") %>%
  fix_border_issues(part = "all") %>%
  align(j=1:7, align="center", part = "all") %>%
  valign(j = 1, valign = "top") %>%
  autofit()
loc_alpha$tbl
save_as_docx(loc_alpha$tbl, path = "Tables/loc_alpha.docx")


# Exploratory Factor Analysis
loc_efa <- read_dta("Data/loc_factors.dta") %>%
  mutate(type = ifelse(str_detect(item, "Int"), "Internal LOC", "External LOC"),
         q = loc_lab[item],
         n = str_sub(q, 1, 2))
ggplot(loc_efa) +
  aes(x = Factor1, y = Factor2, color = type, label = n) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_text() +
  guides(label = FALSE, 
         colour = guide_legend(override.aes = list(size=5,
                                                   label = "\U2022"))) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Factor 1", y = "Factor 2", color = "Wording") +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.85))
ggsave(filename = "Images/loc_efa.png",
       width = 21, height = 9.9, units = "cm", dpi = 300)


# 3. School Attitude and Risk Behaviours ----
get_fa <- function(vars){
  f <- df %>%
    select(matches(vars)) %>%
    drop_na() %>%
    as.matrix() %>%
    psych::fa(rotate = 'none', fm = 'pa', max.iter = 1) %>%
    pluck(loadings)
  class(f) <- "matrix"
  
  as_tibble(f, rownames = "var") %>%
    mutate(PA1 = round(PA1, 2)) %>%
    separate(var, c("var", "wave"), sep = "_") %>%
    pivot_wider(names_from = wave, values_from = PA1)
}
behav_efa <- bind_rows(get_fa("SchoolAtt_W"), get_fa("Risk_W")) %>%
  mutate(var = ifelse(var == "Risk", "Risk Behaviours", "Attitude to School")) %>%
  flextable() %>%
  set_header_labels(var = "Variable", W1 = "13/14", 
                    W2 = "14/15", W3 = "15/16") %>%
  border_remove() %>%
  add_header(var = "", W1 = "Age", W2 = "Age", W3 = "Age") %>%
  merge_h(part = "header")  %>%
  fontsize(size = 11, part = "all") %>%
  hline_top(border = fp_border(color="black", width = 2), part = "all") %>%
  hline_bottom(border = fp_border(color="black", width = 2), part = "all") %>%
  fix_border_issues(part = "all") %>%
  align(j=1, align="right", part = "all") %>%
  align(j=2:4, align="center", part = "all") %>%
  valign(j = 1, valign = "top") %>%
  autofit()
behav_efa
save_as_docx(behav_efa, path = "Tables/behaviour_efa.docx")

# 4. Final Dataset ----
df <- df %>%
  select(-matches("LOC(.*)Item"), -matches("^W\\d_GHQ")) %>%
  full_join(df_loc_latent, by = "NSID") %>%
  full_join(df_ghq_latent, by = "NSID")
save(df, file = "Data/df_full.Rdata")


depvar_sd <- df %>%
  select(starts_with("GHQ_W8"), Patience_W8, Height_W8, Survey_Weight_W8) %>%
  pivot_longer(cols = - Survey_Weight_W8, names_to = "var") %>%
  group_by(var) %>%
  summarise(mean = wtd.mean(value, Survey_Weight_W8),
            sd = wtd.var(value, Survey_Weight_W8, normwt = TRUE) %>% sqrt())
save(depvar_sd, file = "Data/depvar_sd.Rdata")


# 5. Negative Control Quality ----
lm_robust(Height_W8 ~ -1 + Female:NSSEC3_W1, df, Survey_Weight_W8) %>%
  tidy() %>%
  select(term, beta = 2, lci = 6, uci = 7) %>%
  separate(term, c("sex", "sep"), sep = ":") %>%
  mutate(sex = str_replace(sex, "Female", ""),
         sep = str_replace(sep, "NSSEC3_W1", "") %>%
           factor(levels(df$NSSEC3_W1))) %>%
  ggplot() +
  aes(x = sep, y = beta, ymin = lci, ymax = uci,
      color = sex, group = sex) +
  geom_line() +
  geom_pointrange() +
  scale_color_brewer(palette = "Set1") +
  labs(x = NULL, y = "Height", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("Images/height_quality.png", height = 9.9, width = 21, units = "cm")


lm_robust(Patience_W8 ~ -1 + Education_W8, df, Survey_Weight_W8) %>%
  tidy() %>%
  select(term, beta = 2, lci = 6, uci = 7) %>%
  mutate(term = str_replace(term, "Education_W8", "") %>%
           factor(levels(df$Education_W8))) %>%
  ggplot() +
  aes(x = term, y = beta, ymin = lci, ymax = uci, color = term) +
  geom_pointrange() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = NULL, y = "Patience", color = NULL) +
  guides(color = FALSE) +
  theme_minimal()
ggsave("Images/patience_quality.png", height = 9.9, width = 21, units = "cm")
