#-------------------------------------------------------------------------------
#
#     When do correlations between spatial sampling bias and a species'
#     environmental niche affect species distribution models
#
#     by David J Baker
#
#-------------------------------------------------------------------------------


#~# Clear environment
rm(list = ls())

#-------------------------------------------------------------------------------
#
# SECTION 1. Setup
#
#-------------------------------------------------------------------------------


#---- Increase memory size
memory.size(1000000)

#---- Set seed
set.seed(17)

#---- Load packages
library(tidyverse)
library(snowfall)
library(raster)
library(Rfast)
library(gtools)
library(rgeos)
library(data.table)
library(mgcv)
library(brms)
library(speedglm)
library(spatstat)
library(Metrics)
library(ecospat)
library(ggthemes)
library(ggforce)
library(gghalves)
library(see)
library(ggh4x)
library(cowplot)
library(scico)
library(gtable)
library(UCSCXenaShiny) # https://github.com/openbiox/UCSCXenaShiny

#---- Path to source functions
source("Scripts/BiasStructureAnalysis.R")

#-------------------------------------------------------------------------------
#
# Simulation models
#
#-------------------------------------------------------------------------------


#----
#---- Simulate environmental variable (e.g. forest cover, elevation etc.)
#----

#- For an entirely reproducible analysis use simulated environmental data, here
#- generated on a 500 x 500 m grid for 5 variables
enviro <- generate.enviro(n = 100, phi = c(1, 0.1, 0.5, 0.1, 1))
#save(enviro, file = "Outputs/enviro_dat_500x500_phi_1_01_05_01_1.Rdata")

#- Or, load pre-simulated data (for speed and reproducability)
load("Outputs/enviro_dat_500x500_phi_1_01_05_01_1.Rdata")

#- Create plots of variables
plot.variables(enviro, "V1", "Outputs/Figures/enviro_v1")
plot.variables(enviro, "V2", "Outputs/Figures/enviro_v2")
plot.variables(enviro, "V3", "Outputs/Figures/enviro_v3")
plot.variables(enviro, "V4", "Outputs/Figures/enviro_v4")
plot.variables(enviro, "V5", "Outputs/Figures/enviro_v5")

#----
#---- Generate biased surface towards key fine-scale variable
#----

#-- Towards key fine-scale variable
# bias_fs_pos <-
#   samp.bias.Landcov(enviro[, c("id", "X", "Y", "V1")], B = 0.6, A = -0.06)
# save(bias_fs_pos, file= "Outputs/Bias_data_files/bias_fs_pos.Rdata")
load("Outputs/Bias_data_files/bias_fs_pos.Rdata")
plot.variables(bias_fs_pos, "wgt", "Outputs/Figures/bias_fs_pos",
               colPal = "scico")

#-- Away from key fine-scale variable
# bias_fs_neg <-
#   samp.bias.Landcov(enviro[, c("id", "X", "Y", "V1")], B = 0.4, A = 0.06)
# save(bias_fs_neg, file= "Outputs/Bias_data_files/bias_fs_neg.Rdata")
load("Outputs/Bias_data_files/bias_fs_neg.Rdata")
plot.variables(bias_fs_neg, "wgt", "Outputs/Figures/bias_fs_neg",
               colPal = "scico")

#-- Towards a key coarse-scale variable
# bias_cs_pos <-
#   samp.bias.Landcov(enviro[,c("id","X","Y","V2")], B = 0.6, A = -0.06)
# save(bias_cs_pos, file= "Outputs/Bias_data_files/bias_cs_pos.Rdata")
load("Outputs/Bias_data_files/bias_cs_pos.Rdata")
plot.variables(bias_fs_neg, "wgt", "Outputs/Figures/bias_cs_pos",
               colPal = "scico")

#-- Away from a key coarse-scale variable
# bias_cs_neg <-
#   samp.bias.Landcov(enviro[,c("id","X","Y","V2")], B = 0.4, A = 0.06)
# save(bias_cs_neg, file= "Outputs/Bias_data_files/bias_cs_neg.Rdata")
load("Outputs/Bias_data_files/bias_cs_neg.Rdata")
plot.variables(bias_cs_neg, "wgt", "Outputs/Figures/bias_cs_neg",
               colPal = "scico")

#-- Towards a neutral coarse-scale variable
# bias_cs_neu <-
#   samp.bias.Landcov(enviro[,c("id","X","Y","V4")], B = 0.6, A = -0.06)
# save(bias_cs_neu, file= "Outputs/Bias_data_files/bias_cs_neu.Rdata")
load("Outputs/Bias_data_files/bias_cs_neu.Rdata")
plot.variables(bias_cs_neu, "wgt", "Outputs/Figures/bias_cs_neu",
               colPal = "scico")

#-- Towards a neutral fine-scale variable
# bias_fs_neu <-
#   samp.bias.Landcov(enviro[,c("id","X","Y","V5")], B = 0.6, A = -0.06)
# save(bias_fs_neu, file= "Outputs/Bias_data_files/bias_fs_neu.Rdata")
load("Outputs/Bias_data_files/bias_fs_neu.Rdata")
plot.variables(bias_fs_neu, "wgt", "Outputs/Figures/bias_fs_neu",
               colPal = "scico")

#----
#----  Species type 1: Low habitat specificity
#----

#--- Simulate generalist species
sp_com_gen <- cbind.data.frame(
  enviro,
  species.occ(
    K = 0.8,
    B = 0.75,
    A = -0.1,
    Var1 = enviro$V1,
    k = 0.8,
    b = 0.65,
    a = 0.2,
    Var2 = enviro$V2,
    Var3 = enviro$V3
  )
)
#save(sp_com_gen,  file = "Outputs/Species_data/sp_com_gen.Rdata")

#-- Or, load pre-simulated species (for reproducability)
load("Outputs/Species_data/sp_com_gen.Rdata")

#- POC map figure
plot.variables(sp_com_gen,
               "poc",
               "Outputs/Figures/sp_com_generalist_POC",
               colPal = "scico")

#- PA map figure
plot.variables(sp_com_gen, "pa", "Outputs/Figures/sp_com_generalist_PA")

#- Number of presences and prev
sum(sp_com_gen$pa) / nrow(enviro)

#--- Simulate specialist species
sp_com_spec <- cbind.data.frame(
  enviro,
  species.occ(
    K = 0.8,
    B = 0.7,
    A = -0.04,
    Var1 = enviro$V1,
    k = 0.8,
    b = 0.65,
    a = 0.2,
    Var2 = enviro$V2,
    Var3 = enviro$V3
  )
)
#save(sp_com_spec,  file = "Outputs/Species_data/sp_com_spec.Rdata")

#-- Or, load pre-simulated species (for reproducability)
load("Outputs/Species_data/sp_com_spec.Rdata")

#- POC map figure
plot.variables(sp_com_spec,
               "poc",
               "Outputs/Figures/sp_com_specialist_POC",
               colPal = "scico")

#- PA map figure
plot.variables(sp_com_spec, "pa", "Outputs/Figures/sp_com_specialist_PA")

#- Number of presences and prev
sum(sp_com_spec$pa) / nrow(enviro)

#-------------------------------------------------------------------------------
#
# Run simulation experiements
#
#-------------------------------------------------------------------------------

#----
#---- Combinations of models that need to be run
#----
xp <- expand.grid(
  spName = c("sp_com_gen", "sp_com_spec"),
  biasType = c(
    "NULL",
    "bias_fs_pos",
    "bias_fs_neg",
    "bias_cs_pos",
    "bias_cs_neg",
    "bias_cs_neu",
    "bias_fs_neu"
  ),
  nPres = 200,
  detProb = c(0.4, 0.7, 1.0)
)
xp$falPosRat <- 0
xp$listInclusion <- 1
xp$effort <- 100

#----
#---- Detection / non-detection / no bias correction
#----
run_exp(xp,
        modType = "detNondet",
        biasCor = FALSE,
        biasMeth = "None")

#----
#---- Detection / non-detection / bias correction / thin majority
#----
run_exp(xp[3:nrow(xp),],
        modType = "detNondet",
        biasCor = TRUE,
        biasMeth = "thin_maj")

#----
#---- Detection / non-detection / bias correction / thin all
#----
run_exp(xp,
        modType = "detNondet",
        biasCor = TRUE,
        biasMeth = "thin_all")

#----
#---- Detection-only / no bias correction
#----
run_exp(xp,
        modType = "detOnly",
        biasCor = FALSE,
        biasMeth = "None")

#----
#---- Detection-only / bias correction (cluster)
#----
run_exp(xp,
        modType = "detOnly",
        biasCor = TRUE,
        biasMeth = "cluster")

#----
#---- Detection-only / bias correction (real)
#----
run_exp(xp,
        modType = "detOnly",
        biasCor = TRUE,
        biasMeth = "real")


#-------------------------------------------------------------------------------
#
#   Analyse Model Results
#
#-------------------------------------------------------------------------------

#----
#----  What are the effects of biases and uncertainty on model performance?
#----  Using rank correlation coefficients
#----

#--- Facet labels for panel plots
det.labs <- c("DetProb = 0.4", "DetProb = 0.7", "DetProb = 1.0")
names(det.labs) <- c("0.4", "0.7", "1")

#--- Precision for x
scaleFUN <- function(x) {
  sprintf("%.2f", x)
}

#--- Load results data
f <-
  list.files("Outputs/Model_results/",
             full.names = TRUE,
             pattern = "_200")
modRes <- lapply(f, function(x)
  get(load(x)))
modRes <- do.call(rbind, modRes) %>%
  dplyr::select(modType,
                biasMeth,
                spName,
                detProb,
                fpr,
                biasType,
                n_sites,
                n_det,
                corP,
                corS) %>%
  filter(fpr == 0) %>%
  droplevels()

#-- Create longer df for correlation statistics
modResCor <- pivot_longer(modRes, cols = c(corP, corS)) %>%
  mutate(biasMeth = factor(
    biasMeth,
    levels = c("None",
               "thin_all",
               "thin_maj",
               "real",
               "cluster")
  ))
#- Correlation coefficients can be negative but need to set there slightly about
#- zero in order to use a beta model.
modResCor$value[modResCor$value <= 0] <- 0.001

#-- Recode factors
#' Don't do this before running the brms models because the pretty names
#' are difficult to handle in the model.
modRes <- modRes %>%
  mutate(
    biasType = fct_recode(
      biasType,
      "No bias" = "NULL",
      "fs(+)" = "bias_fs_pos",
      "fs(-)" =  "bias_fs_neg",
      "cs(+)" = "bias_cs_pos",
      "cs(-)" = "bias_cs_neg",
      "cs(0)" = "bias_cs_neu",
      "fs(0)" = "bias_fs_neu"
    )
  ) %>%
  mutate(spName = fct_recode(spName,
                             "Generalist" = "sp_com_gen",
                             "Specialist" = "sp_com_spec")) %>%
  mutate(
    modType = fct_recode(
      modType,
      "Detection / NonDetection" = "detNondet",
      "Detection / Background" = "detOnly"
    )
  ) %>%
  mutate(
    biasMeth = fct_recode(
      biasMeth,
      "Naive" = "None",
      "Thin all" = "thin_all",
      "Thin maj." = "thin_maj",
      "True bias" = "real",
      "Cluster" = "cluster"
    )
  ) %>%
  mutate(biasMeth = factor(
    biasMeth,
    levels = c("Naive",
               "Thin all",
               "Thin maj.",
               "True bias",
               "Cluster")
  ))

modRes %>% filter(modType == "Detection / Background", biasMeth ==  "Naive") %>%
  group_by(spName, biasType , detProb) %>%
  summarise(mean(corS))

#---
#--- Figure 2 panel
#---

#-- Test significance of bias type relative to the naive models for each panel
ex <- expand.grid(
  spName = c("sp_com_gen", "sp_com_spec"),
  detProb  = c(0.4, 0.7, 1.0),
  fpr = 0,
  modType = c("detNondet", "detOnly"),
  name = c("corS")
)
for (i in 1:nrow(ex)) {
  mdat <- modResCor %>%
    dplyr::filter(detProb == ex[i, "detProb"],
                  fpr == ex[i, "fpr"],
                  spName == ex[i, "spName"],
                  modType == ex[i, "modType"],
                  biasMeth == "None",
                  name == ex[i, "name"])
  m1 <- brm(
    bf(value ~ biasType, phi ~ biasType),
    family = Beta(),
    cores = 2,
    chains = 2,
    data = mdat
  )
  save(
    m1,
    file = paste0(
      "Outputs/brms_models/pl_",
      paste(ex[i, 1],
            ex[i, 2],
            ex[i, 3],
            ex[i, 4],
            "None",
            ex[i, 5], sep = "_"),
      ".Rdata"
    ),
    compress = "xz"
  )
}

#-- Process BRMs results
brms_res <- lapply(1:nrow(ex), function(i) {
  print(i)
  tmp <- get(load(paste0(
    "Outputs/brms_models/pl_",
    paste(ex[i, 1],
          ex[i, 2],
          ex[i, 3],
          ex[i, 4],
          "None",
          ex[i, 5], sep = "_"),
    ".Rdata"
  )))
  tmpSum <- summary(tmp)
  tmpSum <- as.data.frame(tmpSum$fixed[3:8, 3:4])
  names(tmpSum) <- c("l95", "u95")
  tmpSum$Star <- 0
  tmpSum$Star[sign(tmpSum[, 1]) == sign(tmpSum[, 2])] <- 1
  tmpSum$biasType <- row.names(tmpSum)
  tmpSum$biasType <- sub("biasType", "", tmpSum$biasType)
  tmpSum$spName <- ex[i, 1]
  tmpSum$detProb <- ex[i, 2]
  tmpSum$fpr <- ex[i, 3]
  tmpSum$modType <- ex[i, 4]
  tmpSum$biasMeth <- "None"
  tmpSum$name <- ex[i, 5]
  rownames(tmpSum) <- 1:nrow(tmpSum)
  tmpSum
  
})
brms_res <- do.call(rbind, brms_res)
brms_res_df <- brms_res %>%
  mutate(
    biasType = fct_recode(
      biasType,
      "fs(+)" = "bias_fs_pos",
      "fs(-)" =  "bias_fs_neg",
      "cs(+)" = "bias_cs_pos",
      "cs(-)" = "bias_cs_neg",
      "cs(0)" = "bias_cs_neu",
      "fs(0)" = "bias_fs_neu"
    )
  ) %>%
  mutate(spName = fct_recode(spName,
                             "Generalist" = "sp_com_gen",
                             "Specialist" = "sp_com_spec")) %>%
  mutate(
    modType = fct_recode(
      modType,
      "Detection / NonDetection" = "detNondet",
      "Detection / Background" = "detOnly"
    )
  ) %>%
  mutate(biasMeth = fct_recode(biasMeth, "Naive" = "None"))

#-- Fig. 2a detection / non-detection

#- Group means for no bias scenario
fig_2a_means <- modRes %>%
  filter(biasMeth == "Naive",
         biasType == "No bias",
         modType == "Detection / NonDetection",
         fpr == 0) %>%
  group_by(modType, spName, fpr, detProb) %>%
  summarise(none_m = mean(corS))

#- Significant differences from brms models
brms_fig_2a_n <- filter(
  brms_res_df,
  biasMeth == "Naive",
  modType == "Detection / NonDetection",
  fpr == 0,
  Star == 1,
  sign(l95) == -1,
  name == "corS"
)

#- Plot (a)
fig_2a <- modRes %>%
  filter(biasMeth == "Naive",
         modType == "Detection / NonDetection") %>%
  ggplot(aes(
    x = biasType,
    y = corS,
    fill = biasType,
    colour = biasType
  )) +
  geom_half_violin(trim = TRUE,
                   scale = "width",
                   side = "r") +
  geom_hline(aes(yintercept = none_m),
             data = fig_2a_means,
             colour = "grey50") +
  geom_jitter2(
    aes(colour =  biasType),
    position = position_nudge(x = -0.2),
    alpha = 0.2,
    size = 0.75
  ) +
  scale_fill_colorblind() +
  scale_colour_colorblind() +
  geom_text(
    aes(x = biasType,
        y = -Inf,
        label = "*(-) "),
    colour = "red",
    fontface = "bold",
    size = 3,
    hjust = "inward",
    vjust = -0.4,
    data = brms_fig_2a_n
  ) +
  facet_nested(modType + spName ~ detProb,
               scales = "fixed",
               labeller = labeller(detProb = det.labs)) +
  scale_y_continuous(labels = scaleFUN,
                     limits = c(NA, 1)) +
  xlab("Bias type") +
  ylab("Spearman's rank correlation coefficient") +
  theme_bw(14) %+replace%
  theme(
    legend.position =  "none",
    axis.text.x = element_text(angle = 45),
    strip.background = element_rect(color = "Black",
                                    fill = "White")
  ) +
  coord_flip()

#-- Fig. 2b detection / non-detection & FPR = 0.05
#- Group means for no bias scenario
fig_2b_means <- modRes %>%
  filter(biasMeth == "Naive",
         biasType == "None",
         modType == "Detection / Background") %>%
  group_by(modType, spName, detProb) %>%
  summarise(none_m = mean(corS))
#- Significant differences from brms models
brms_fig_2b_n <- filter(
  brms_res_df,
  biasMeth == "Naive",
  modType == "Detection / Background",
  Star == 1,
  sign(l95) == -1,
  name == "corS"
)
brms_fig_2b_p <- filter(
  brms_res_df,
  biasMeth == "Naive",
  modType == "Detection / Background",
  Star == 1,
  sign(l95) == 1,
  name == "corS"
)
#- Plot (b)
fig_2b <- modRes %>%
  filter(biasMeth == "Naive",
         modType == "Detection / Background") %>%
  ggplot(aes(
    x = biasType,
    y = corS,
    fill = biasType,
    colour = biasType
  )) +
  geom_half_violin(trim = TRUE,
                   scale = "width",
                   side = "r") +
  geom_hline(aes(yintercept = none_m),
             data = fig_2b_means,
             colour = "grey50") +
  geom_jitter2(
    aes(colour =  biasType),
    position = position_nudge(x = -0.2),
    alpha = 0.2,
    size = 0.75
  ) +
  scale_fill_colorblind() +
  scale_colour_colorblind() +
  geom_text(
    aes(x = biasType,
        y = -Inf,
        label = "*(-) "),
    colour = "red",
    fontface = "bold",
    size = 3,
    hjust = "inward",
    vjust = -0.4,
    data = brms_fig_2b_n
  ) +
  facet_nested(modType + spName ~ detProb,
               scales = "fixed",
               labeller = labeller(detProb = det.labs)) +
  scale_y_continuous(labels = scaleFUN, limits = c(NA, 1)) +
  xlab("Bias type") +
  ylab("Spearman's rank correlation coefficient") +
  theme_bw(14) %+replace%
  theme(
    legend.position =  "none",
    axis.text.x = element_text(angle = 45),
    strip.background = element_rect(color = "Black",
                                    fill = "White")
  ) +
  coord_flip()


#-- Fig. 2 panel
fig_2 <- plot_grid(
  fig_2a,
  fig_2b,
  nrow = 2,
  labels = c("a)", "b)"),
  align = "hv"
)
#- Save Fig. 2
ggsave(
  filename = "Outputs/Figures/Figure_2ab_corS_Naive.png",
  plot = fig_2,
  width = 8,
  height = 10
)

#----
#---- Summarise the major results of the effects of bias on rank correlation
#----

#---
#--- Detection / Non-detection
#---

#- Mean and range for no bias
modRes %>%
  filter(biasMeth == "Naive",
         biasType == "No bias",
         modType == "Detection / NonDetection") %>%
  group_by(spName) %>%
  summarise(
    corS_m = mean(corS),
    corS_l95 = quantile(corS, prob = 0.025),
    corS_u95 = quantile(corS, prob = 0.975)
  )

#- Test of effect of detection prob for generalist
gdat <- modRes %>%
  filter(
    biasMeth == "Naive",
    biasType == "No bias",
    modType == "Detection / NonDetection",
    spName == "Generalist"
  )
g1 <- brm(
  bf(corS ~ detProb, phi ~ detProb),
  family = Beta(),
  cores = 2,
  chains = 2,
  data = gdat
)
pp_check(g1)
summary(g1)
# Population-Level Effects:
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept         3.40      0.08     3.25     3.56 1.00     1631     1554
# phi_Intercept     4.62      0.26     4.11     5.12 1.00     1780     1204
# detProb           0.55      0.10     0.36     0.74 1.00     1860     1564
# phi_detProb       1.30      0.34     0.65     1.99 1.00     1901     1317

#- Test of effect of detection prob for specialists
sdat <- modRes %>%
  filter(
    biasMeth == "Naive",
    biasType == "No bias",
    modType == "Detection / NonDetection",
    spName == "Specialist"
  )
s1 <- brm(
  bf(corS ~ detProb, phi ~ detProb),
  family = Beta(),
  cores = 2,
  chains = 2,
  data = sdat
)
pp_check(s1)
summary(s1)
# Population-Level Effects:
#                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept         3.42      0.03     3.37     3.48 1.00     1601     1182
# phi_Intercept     5.85      0.25     5.33     6.32 1.00     1723     1100
# detProb           0.17      0.03     0.11     0.23 1.00     1677     1247
# phi_detProb       2.67      0.34     1.99     3.34 1.00     1561     1202


#--- Summarise experiment results

#- No bias vs fs
modRes %>%
  filter(biasMeth == "Naive",
         modType == "Detection / NonDetection") %>%
  group_by(biasType, modType, detProb, spName) %>%
  summarise(corS_m = mean(corS)) %>%
  filter(biasType %in% c("No bias", "fs(-)", "fs(+)"))

#- No bias vs cs
modRes %>%
  filter(biasMeth == "Naive",
         modType == "Detection / NonDetection") %>%
  group_by(biasType, modType, detProb, spName) %>%
  summarise(corS_m = mean(corS)) %>%
  filter(biasType %in% c("No bias", "cs(-)", "cs(+)"))


#---
#--- Detection / background
#---

#- Mean and range for no bias
modRes %>%
  filter(biasMeth == "Naive",
         biasType == "No bias",
         modType == "Detection / Background") %>%
  group_by(spName) %>%
  summarise(
    corS_m = mean(corS),
    corS_l95 = quantile(corS, prob = 0.025),
    corS_u95 = quantile(corS, prob = 0.975)
  )

#- Test of effect of detection prob for generalist
gdat <- modRes %>%
  filter(
    biasMeth == "Naive",
    biasType == "No bias",
    modType == "Detection / Background",
    spName == "Generalist"
  )
g1 <- brm(
  bf(corS ~ detProb, phi ~ detProb),
  family = Beta(),
  cores = 2,
  chains = 2,
  data = gdat
)
pp_check(g1)
summary(g1)
# Population-Level Effects:
#                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept         3.22      0.08     3.06     3.38 1.00     1310     1416
# phi_Intercept     3.82      0.25     3.34     4.32 1.00     1791     1343
# detProb           0.71      0.10     0.52     0.91 1.00     1491     1337
# phi_detProb       2.48      0.34     1.82     3.13 1.00     1867     1358

#- Test of effect of detection prob for specialists
sdat <- modRes %>%
  filter(
    biasMeth == "Naive",
    biasType == "No bias",
    modType == "Detection / Background",
    spName == "Specialist"
  )
s1 <- brm(
  bf(corS ~ detProb, phi ~ detProb),
  family = Beta(),
  cores = 2,
  chains = 2,
  data = sdat
)
pp_check(s1)
summary(s1)
# Population-Level Effects:
#                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept         3.42      0.03     3.37     3.48 1.00     1601     1182
# phi_Intercept     5.85      0.25     5.33     6.32 1.00     1723     1100
# detProb           0.17      0.03     0.11     0.23 1.00     1677     1247
# phi_detProb       2.67      0.34     1.99     3.34 1.00     1561     1202


#--- Summarise experiment results

#- No bias vs fs
modRes %>%
  filter(biasMeth == "Naive",
         modType == "Detection / Background") %>%
  group_by(biasType, modType, detProb, spName) %>%
  summarise(corS_m = mean(corS)) %>%
  filter(biasType %in% c("No bias", "fs(-)", "fs(+)"))

#- No bias vs cs
modRes %>%
  filter(biasMeth == "Naive",
         modType == "Detection / Background") %>%
  group_by(biasType, modType, detProb, spName) %>%
  summarise(corS_m = mean(corS)) %>%
  filter(biasType %in% c("No bias", "cs(-)", "cs(+)"))


#----
#----  Effect of bias correction on rank correlation
#----

#  What are the effects of bias correction model performance under different
#  conditions? Using spearman's rank correlation

ex2 <- expand.grid(
  spName = c("sp_com_gen", "sp_com_spec"),
  detProb  = c(0.4, 0.7, 1.0),
  fpr = 0,
  modType = c("detNondet", "detOnly"),
  name = c("corS")
)

#-- Test significance of bias type relative to the Naive models for each panel
#-- Run model for each panel i = 7
for (i in 1:nrow(ex2)) {
  mdat <- dplyr::filter(modResCor,
                        detProb == ex2[i, "detProb"],
                        fpr == ex2[i, "fpr"],
                        spName == ex2[i, "spName"],
                        modType == ex2[i, "modType"],
                        name == ex2[i, "name"]) %>%
    droplevels()
  m1 <-
    brm(
      bf(value ~ biasType * biasMeth, phi ~ biasType + biasMeth),
      family = Beta(),
      cores = 2,
      chains = 2,
      data = mdat
    )
  save(
    m1,
    file = paste0(
      "Outputs/brms_models/bc_",
      paste(ex2[i, 1],
            ex2[i, 2],
            ex2[i, 3],
            ex2[i, 4],
            ex2[i, 5], sep = "_"),
      ".Rdata"
    ),
    compress = "xz"
  )
  
}

#- Extract BRMs results
brms_tab <- lapply(1:nrow(ex2), function(i) {
  print(i)
  tmp <- get(load(paste0(
    "Outputs/brms_models/bc_",
    paste(ex2[i, 1],
          ex2[i, 2],
          ex2[i, 3],
          ex2[i, 4],
          ex2[i, 5], sep = "_"),
    ".Rdata"
  )))
  tmpSum <-
    as.data.frame(posterior_summary(tmp, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  
  # interactions
  r_id <-
    grep(row.names(tmpSum), pattern = 'b_biasMeth|b_biasTypebias.*:')
  tmpSum <- tmpSum[r_id, c(1, 3:7)]
  names(tmpSum) <- c("Est", "l95", "Q25", "Q50", "Q75", "u95")
  tmpSum$biasType <-
    sub("b_biasMeth|b_biasType", "",  gsub("(:).*", "", row.names(tmpSum)))
  tmpSum$biasType[tmpSum$biasType %in%
                    c("thin_all", "thin_maj", "cluster", "real")] <-
    "None"
  tmpSum$biasMeth <-
    sub("b_biasMeth|biasMeth", "",  gsub(".*(:)", "", row.names(tmpSum)))
  
  # Scenarios
  tmpSum$Star <- 0
  tmpSum$Star[sign(tmpSum[, "l95"]) == sign(tmpSum[, "u95"])] <- 1
  tmpSum$spName <- ex2[i, 1]
  tmpSum$detProb <- ex2[i, 2]
  tmpSum$fpr <- ex2[i, 3]
  tmpSum$modType <- ex2[i, 4]
  tmpSum$metric <- ex2[i, 5]
  row.names(tmpSum) <- 1:nrow(tmpSum)
  tmpSum
  
})
brms_tab <- do.call(rbind, brms_tab)
brms_tab_df <- brms_tab %>%
  mutate(
    biasType = fct_recode(
      biasType,
      "No bias" = "None",
      # "Thin all" = "thin_all",
      # "Thin maj." = "thin_maj",
      "fs(+)" = "bias_fs_pos",
      "fs(-)" = "bias_fs_neg",
      "cs(+)" = "bias_cs_pos",
      "cs(-)" = "bias_cs_neg",
      "cs(0)" = "bias_cs_neu",
      "fs(0)" = "bias_fs_neu"
    )
  ) %>%
  mutate(spName = fct_recode(spName,
                             "Generalist" = "sp_com_gen",
                             "Specialist" = "sp_com_spec")) %>%
  mutate(
    modType = fct_recode(
      modType,
      "Detection / NonDetection" = "detNondet",
      "Detection / Background" = "detOnly"
    )
  ) %>%
  mutate(
    biasMeth = fct_recode(
      biasMeth,
      "Thin all" = "thin_all",
      "Thin maj." = "thin_maj",
      "True bias" = "real",
      "Cluster" = "cluster"
    )
  ) %>%
  mutate(signfx = sign(u95)) %>%
  dplyr::select(biasType, spName, modType, biasMeth, detProb, Star, signfx)

#- Conditional effects
brms_cond_effects <- lapply(1:nrow(ex2), function(i) {
  print(i)
  tmp <- get(load(paste0(
    "Outputs/brms_models/bc_",
    paste(ex2[i, 1],
          ex2[i, 2],
          ex2[i, 3],
          ex2[i, 4],
          ex2[i, 5], sep = "_"),
    ".Rdata"
  )))
  
  # Interactions
  coeff1 <-
    conditional_effects(tmp, robust = FALSE)$'biasType:biasMeth'
  coeff1 <-
    coeff1[, c("biasType", "biasMeth", "estimate__", "lower__", "upper__")]
  names(coeff1) <-
    c("biasType", "biasMeth", "Estimate", "l95", "u95")
  
  # Scenarios
  coeff1$spName <- ex2[i, 1]
  coeff1$detProb <- ex2[i, 2]
  coeff1$fpr <- ex2[i, 3]
  coeff1$modType <- ex2[i, 4]
  coeff1$metric <- ex2[i, 5]
  coeff1
  
})
brms_cond_effects <- do.call(rbind, brms_cond_effects)
brms_cond_effects_df <- brms_cond_effects %>%
  mutate(
    biasType = fct_recode(
      biasType,
      "No bias" = "None",
      "Thin all" = "thin_all",
      "Thin maj." = "thin_maj",
      "Cluster" = "cluster",
      "No bias" = "NULL",
      "fs(+)" = "bias_fs_pos",
      "fs(-)" = "bias_fs_neg",
      "cs(+)" = "bias_cs_pos",
      "cs(-)" = "bias_cs_neg",
      "cs(0)" = "bias_cs_neu",
      "fs(0)" = "bias_fs_neu"
    )
  ) %>%
  mutate(spName = fct_recode(spName,
                             "Generalist" = "sp_com_gen",
                             "Specialist" = "sp_com_spec")) %>%
  mutate(
    modType = fct_recode(
      modType,
      "Detection / NonDetection" = "detNondet",
      "Detection / Background" = "detOnly"
    )
  ) %>%
  mutate(
    biasMeth = fct_recode(
      biasMeth,
      "Uncorrected" = "None",
      "Thin all" = "thin_all",
      "Thin maj." = "thin_maj",
      "True bias" = "real",
      "Cluster" = "cluster"
    )
  )
#- Merge
condFxDf <- left_join(
  brms_cond_effects_df,
  brms_tab_df,
  by = c("biasType", "spName", "modType", "biasMeth", "detProb")
)
condFxDf$Star[condFxDf$Star == 0] <- NA
condFxDf$signfx[is.na(condFxDf$Star)] <- NA
condFxDf$sigLab <- NA
condFxDf$sigLab[condFxDf$signfx == -1] <- "atop('*'^'-')"
condFxDf$sigLab[condFxDf$signfx == 1] <- "atop('*'^'+')"

#- Figure 3
fig3a <-
  filter(condFxDf,
         modType == "Detection / NonDetection") %>%
  ggplot(aes(x = biasType,
             y = Estimate,
             colour = biasMeth)) +
  geom_point2(position = position_dodge(width = 0.75)) +
  geom_errorbar(
    aes(x = biasType,
        ymin = l95,
        ymax = u95),
    position = position_dodge(width = 0.75),
    width = 0.75
  ) +
  scale_colour_colorblind() +
  ylab("Conditional effect") +
  xlab("Bias type") +
  geom_hline(
    yintercept = 1,
    linetype = 1,
    size = 0.25,
    alpha = 1
  ) +
  geom_vline(xintercept = seq(1.5, 6.5, 1), linetype = 3) +
  geom_text(
    aes(
      x = biasType,
      y = 1.01,
      label = sigLab,
      colour = biasMeth
    ),
    fontface = "bold",
    position = position_dodge(width = 0.75),
    size = 5,
    vjust = 0.9,
    hjust = 0.7,
    parse = TRUE,
    show.legend = FALSE
  ) +
  facet_nested(modType + spName ~ detProb,
               scales = "fixed",
               labeller = labeller(detProb = det.labs)) +
  theme_bw(14) %+replace%
  theme(
    legend.position =  c(0.095, 0.11),
    legend.title = element_blank(),
    legend.key.size = unit(0.3, 'cm'),
    legend.background = element_rect(fill = "white",
                                     colour = "transparent"),
    axis.text.x = element_text(angle = 45),
    axis.title.y = element_blank(),
    strip.background = element_rect(color = "Black",
                                    fill = "White"),
    panel.grid = element_blank()
  ) +
  coord_flip()

fig3b <-
  filter(condFxDf,
         modType == "Detection / Background") %>%
  ggplot(aes(x = biasType,
             y = Estimate,
             colour = biasMeth)) +
  geom_point2(position = position_dodge(width = 0.75)) +
  geom_errorbar(
    aes(x = biasType,
        ymin = l95,
        ymax = u95),
    position = position_dodge(width = 0.75),
    width = 0.75
  ) +
  scale_colour_colorblind() +
  scale_y_continuous(
    labels = c("0.2", "0.4", "0.6", "0.8", "1.0"),
    breaks = seq(0.2, 1.0, 0.2)
  ) +
  ylab("Conditional effect") +
  xlab("Bias type") +
  geom_hline(
    yintercept = 1,
    linetype = 1,
    size = 0.25,
    alpha = 1
  ) +
  geom_vline(xintercept = seq(1.5, 6.5, 1), linetype = 3) +
  geom_text(
    aes(
      x = biasType,
      y = 1.075,
      label = sigLab,
      colour = biasMeth
    ),
    fontface = "bold",
    position = position_dodge(width = 0.75),
    size = 5,
    vjust = 0.9,
    hjust = 0.7,
    parse = TRUE,
    show.legend = FALSE
  ) +
  facet_nested(modType + spName ~ detProb,
               scales = "fixed",
               labeller = labeller(detProb = det.labs)) +
  theme_bw(14) %+replace%
  theme(
    legend.position =  c(0.095, 0.11),
    legend.title = element_blank(),
    legend.key.size = unit(0.3, 'cm'),
    legend.background = element_rect(fill = "white",
                                     colour = "transparent"),
    axis.text.x = element_text(angle = 45),
    axis.title.y = element_blank(),
    strip.background = element_rect(color = "Black",
                                    fill = "White"),
    panel.grid = element_blank()
  ) +
  coord_flip()

#- Figure 3 panel
fig3 <-
  plot_grid(
    fig3a,
    fig3b,
    nrow = 2,
    labels = c("a)", "b)"),
    align = "hv"
  )

#- Save figure 3
ggsave(
  filename = "Outputs/Figures/Figure_3_BiasMod_Regression_Panel.png",
  plot = fig3,
  width = 8,
  height = 10
)

#---
#--- Summary statistics for experiments with bias correction
#---

#- "Detection / NonDetection"
modRes %>%
  filter(modType == "Detection / NonDetection") %>%
  group_by(spName, biasMeth, detProb, biasType) %>%
  summarise(
    corS_m = mean(corS),
    corS_l95 = quantile(corS, prob = 0.025),
    corS_u95 = quantile(corS, prob = 0.975)
  ) %>%
  filter(
    biasType %in% c("No bias", "fs(-)", "fs(+)", "cs(-)", "cs(+)"),
    spName == "Generalist",
    detProb == 0.4
  )

#- "Detection / background"
modRes %>%
  filter(modType == "Detection / Background") %>%
  group_by(spName, biasMeth, biasType) %>%
  summarise(
    corS_m = mean(corS),
    corS_l95 = quantile(corS, prob = 0.025),
    corS_u95 = quantile(corS, prob = 0.975)
  ) %>%
  filter(biasType %in% c("No bias", "fs(-)", "fs(+)", "cs(-)", "cs(+)"),
         spName == "Generalist")


#-------------------------------------------------------------------------------
#
#  Technical insights into how structure affects models
#
#-------------------------------------------------------------------------------

biasTypeALL = c("NULL",
                "bias_fs_pos",
                "bias_fs_neg",
                "bias_cs_pos",
                "bias_cs_neg")

biasExGen <- lapply(biasTypeALL, function(biasType) {
  # Load species data
  sp_in <-
    get(load(
      paste0("Outputs/Species_data_files/", "sp_com_gen", ".Rdata")
    ))
  
  # Load bias data
  if (biasType == "NULL") {
    landuse.bias <- rep(1, nrow(sp_in))
  } else {
    landuse.bias <-
      get(load(paste0(
        "Outputs/Bias_data_files/", biasType, ".Rdata"
      )))
    landuse.bias <- landuse.bias$wgt
  }
  
  # Sample distribution
  sp_samp <- species.sampling(
    sp_dat = sp_in,
    bias = landuse.bias,
    nPres = 200,
    detProb = 1,
    fpr = 0,
    listProb = 1
  )
  sp_samp$biasType <- biasType
  sp_samp
  
})
biasExGen <- do.call(rbind, biasExGen)
Fig4a <- biasExGen %>%
  mutate(
    biasType = fct_recode(
      biasType,
      "No bias" = "NULL",
      "fs(+)" = "bias_fs_pos",
      "fs(-)" = "bias_fs_neg",
      "cs(+)" = "bias_cs_pos",
      "cs(-)" = "bias_cs_neg"
    )
  ) %>%
  mutate(biasType = factor(biasType, levels = c(
    "No bias",
    "fs(+)",
    "fs(-)",
    "cs(+)",
    "cs(-)"
  ))) %>%
  filter(det == 1) %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_point2(aes(colour = biasType), alpha = 0.5) +
  geom_mark_ellipse(aes(color = biasType, linetype = biasType), size = 1.5) +
  scale_color_colorblind() +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw(14) %+replace% theme(legend.position = "none") +
  xlab('Fine-scale "habitat" variable') +
  ylab('Coarse-scale "elevation" variable') +
  ggtitle("Generalist")

biasExSpe <- lapply(biasTypeALL, function(biasType) {
  # Load species data
  sp_in <-
    get(load(
      paste0("Outputs/Species_data_files/", "sp_com_spec", ".Rdata")
    ))
  
  # Load bias data
  if (biasType == "NULL") {
    landuse.bias <- rep(1, nrow(sp_in))
  } else {
    landuse.bias <-
      get(load(paste0(
        "Outputs/Bias_data_files/", biasType, ".Rdata"
      )))
    landuse.bias <- landuse.bias$wgt
  }
  
  # Sample distribution
  sp_samp <- species.sampling(
    sp_dat = sp_in,
    bias = landuse.bias,
    nPres = 200,
    detProb = 1,
    fpr = 0,
    listProb = 1
  )
  sp_samp$biasType <- biasType
  sp_samp
  
})
biasExSpe <- do.call(rbind, biasExSpe)
Fig4b <- biasExSpe %>%
  mutate(
    biasType = fct_recode(
      biasType,
      "No bias" = "NULL",
      "fs(+)" = "bias_fs_pos",
      "fs(-)" = "bias_fs_neg",
      "cs(+)" = "bias_cs_pos",
      "cs(-)" = "bias_cs_neg"
    )
  ) %>%
  mutate(biasType = factor(biasType, levels = c(
    "No bias",
    "fs(+)",
    "fs(-)",
    "cs(+)",
    "cs(-)"
  ))) %>%
  filter(det == 1) %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_point2(aes(colour = biasType), alpha = 0.5) +
  geom_mark_ellipse(aes(color = biasType, linetype = biasType), size = 1.5) +
  scale_color_colorblind(name = "Bias type") +
  scale_linetype(name = "Bias type") +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw(14) %+replace% theme(
    legend.position = c(0.2, 0.3),
    legend.key.size = unit(0.3, 'cm'),
    legend.title = element_text(size = 12)
  ) +
  xlab('Fine-scale "habitat" variable') +
  ylab('Coarse-scale "elevation" variable') +
  ggtitle("Specialist")

#- Figure 4 panel
fig4 <-
  plot_grid(
    Fig4a,
    Fig4b,
    ncol = 2,
    labels = c("a)", "b)"),
    align = "hv"
  )

#- Save figure 3
ggsave(
  filename = "Outputs/Figures/Figure_4_SampleSpace_Panel.png",
  plot = fig4,
  width = 10,
  height = 5
)

#-------------------------------------------------------------------------------#
#-------------------------------- END ------------------------------------------#
#-------------------------------------------------------------------------------#