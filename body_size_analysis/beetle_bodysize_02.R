
#### Libraries etc ####

library(tidyverse)
library(ggridges)
library(brms)
library(tidybayes)
library(modelr)
library(patchwork)
library(RColorBrewer)
library(broom.mixed)

set.seed(29)


custom_colors <- c(brewer.pal(9, "Set1")[1],   # Red
                   brewer.pal(9, "Set1")[5], # Orange
                   brewer.pal(9, "Set1")[2])   # Blue


#### Model checking functions ####

check_brms <- function(model,             # brms model
                       integer = FALSE,   # integer response? (TRUE/FALSE)
                       plot = TRUE,       # make plot?
                       ...                # further arguments for DHARMa::plotResiduals 
) {
  
  mdata <- brms::standata(model)
  if (!"Y" %in% names(mdata))
    stop("Cannot extract the required information from this brms model")
  
  dharma.obj <- DHARMa::createDHARMa(
    simulatedResponse = t(brms::posterior_predict(model, ndraws = 1000)),
    observedResponse = mdata$Y, 
    fittedPredictedResponse = apply(
      t(brms::posterior_epred(model, ndraws = 1000, re.form = NA)),
      1,
      mean),
    integerResponse = integer)
  
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  
  invisible(dharma.obj)
  
}


#### Wrangling ####

dat <- read_csv("g3personalitymean.csv") 
breeding <- read_csv("breeding_data.csv")

breeding <- breeding %>% 
  mutate(caretreatment = substr(care,1,1),
         fnum = paste(caretreatment,pair,sep=""))

df_size <- dat %>% 
  rename(fnum = newfamily,
         caretreatment = tr) %>% 
  left_join(breeding %>% select(fnum,broodsize,broodmass, wt))

df_size

df_size %>% 
  select(caretreatment,fnum) %>% 
  distinct() %>% 
  group_by(caretreatment) %>% 
  tally()

df_size <- df_size |> 
  mutate(caretreatment_fct = factor(caretreatment,
                                    levels = c("n","m","f")))

df_size <- df_size |> 
  mutate(caretreatment_fct = factor(case_when(
    caretreatment == "f" ~ "FC",
    caretreatment == "m" ~ "MC",
    caretreatment == "n" ~ "NC"
  ),
  levels = c("FC","MC","NC")))


df_size |> 
  ggplot(aes(y = caretreatment_fct, x = bodysize,
             fill = caretreatment_fct)) +
  geom_density_ridges(alpha = 0.6) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  theme_classic(base_size = 14) +
  coord_cartesian(xlim = c(2.6,6.5)) +
  coord_flip()+
  labs(y = "Care treatment",
       x = "Body size")


df_size |> 
  ggplot(aes(y = caretreatment_fct, x = bodysize,
             fill = caretreatment_fct)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.6, alpha = 0.7) +
  stat_dotsinterval(side = "bottom", scale = 0.6, slab_linewidth = NA, alpha = 0.7) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  theme_classic(base_size = 14)

df_size |> 
  ggplot(aes(y = caretreatment_fct, x = bodysize,
             fill = caretreatment_fct)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.6, alpha = 0.6) +
  stat_dotsinterval(side = "bottom", scale = 0.6, slab_linewidth = NA, alpha = 0.7) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  theme_classic(base_size = 14) +
  coord_flip()

#### Modelling ####

##
# These models are from the 'newSun_data' script
# - Do we need to consider sex interactions with other fixed effects?
#  -> I don't remember if we did that already
##

#### .Model runs ####

## body size with sex
m_bs <- brm(brmsformula(bodysize ~ sex + 
                          caretreatment_fct * scale(broodsize, scale=FALSE) +
                          (0 + caretreatment_fct||fnum),
                        sigma ~ 0 + caretreatment_fct),
            data = as.data.frame(df_size),
            iter=12000,
            warmup=2000,
            chains=6,
            cores=6,
            control=list(adapt_delta=0.99,max_treedepth=15),
            backend="cmdstanr",
            file = "Models/newSun.RDS")


## body size without sex covariate
m_bs2 <- brm(brmsformula(bodysize ~ caretreatment_fct * scale(broodsize, scale=FALSE) +
                           (0 + caretreatment_fct||fnum),
                         sigma ~ 0 + caretreatment_fct),
             data = as.data.frame(df_size),
             iter=12000,
             warmup=2000,
             chains=6,
             cores=6,
             control=list(adapt_delta=0.99,max_treedepth=15),
             backend="cmdstanr",
             file = "Models/newSun_noSex.RDS")


## body size with all interactions with sex
m_bs_ints <- brm(brmsformula(bodysize ~ sex * 
                               caretreatment_fct * scale(broodsize, scale=FALSE) +
                               (0 + caretreatment_fct||fnum),
                             sigma ~ 0 + caretreatment_fct),
                 data = as.data.frame(df_size),
                 iter=12000,
                 warmup=2000,
                 chains=6,
                 cores=6,
                 control=list(adapt_delta=0.99,max_treedepth=15),
                 backend="cmdstanr",
                 file = "Models/newSun_allints.RDS")


#### Modelling - include carcass weight ####

#### .Model runs ####

## body size with sex
m_bs_wt <- brm(brmsformula(bodysize ~ sex + 
                             caretreatment_fct * scale(broodsize, scale=FALSE) +
                             scale(wt) + #caretreatment_fct:scale(wt) +
                             (0 + caretreatment_fct||fnum),
                           sigma ~ 0 + caretreatment_fct),
               data = as.data.frame(df_size),
               iter=12000,
               warmup=2000,
               chains=6,
               cores=6,
               control=list(adapt_delta=0.99,max_treedepth=15),
               backend="cmdstanr",
               file = "Models/newSun_wt.RDS")


## body size without sex covariate
m_bs2_wt <- brm(brmsformula(bodysize ~ caretreatment_fct * scale(broodsize, scale=FALSE) +
                              scale(wt) + #caretreatment_fct:scale(wt) +
                              (0 + caretreatment_fct||fnum),
                            sigma ~ 0 + caretreatment_fct),
                data = as.data.frame(df_size),
                iter=12000,
                warmup=2000,
                chains=6,
                cores=6,
                control=list(adapt_delta=0.99,max_treedepth=15),
                backend="cmdstanr",
                file = "Models/newSun_noSex_wt.RDS")



## This might be stupid
## body size with sex
# m_bs_wt_broodres <- brm(brmsformula(bodysize ~ sex + 
#                                       caretreatment_fct * scale(broodsize, scale=FALSE) +
#                                       scale(wt) + caretreatment_fct:scale(wt) +
#                                       (0 + caretreatment_fct||fnum),
#                                     sigma ~ 0 + caretreatment_fct*scale(broodsize, scale = FALSE)),
#                         data = as.data.frame(df_size),
#                         iter=12000,
#                         warmup=2000,
#                         chains=6,
#                         cores=6,
#                         control=list(adapt_delta=0.99,max_treedepth=15),
#                         backend="cmdstanr",
#                         file = "Models/newSun_wt_broodres.RDS")





#### .Check whether variance terms have strong effects ####

m_bs_nosig <- brm(brmsformula(bodysize ~ sex + 
                                caretreatment_fct * scale(broodsize, scale=FALSE) +
                                (0 + caretreatment_fct||fnum)),
                  #sigma ~ 0 + caretreatment_fct),
                  data = as.data.frame(df_size),
                  iter=12000,
                  warmup=2000,
                  chains=6,
                  cores=6,
                  control=list(adapt_delta=0.99,max_treedepth=15),
                  backend="cmdstanr",
                  file = "Models/newSun_nosig.RDS")

m_bs_notrtvar <- brm(brmsformula(bodysize ~ sex + 
                                   caretreatment_fct * scale(broodsize, scale=FALSE) +
                                   (1|fnum)),
                     #sigma ~ 0 + caretreatment_fct),
                     data = as.data.frame(df_size),
                     iter=12000,
                     warmup=2000,
                     chains=6,
                     cores=6,
                     control=list(adapt_delta=0.99,max_treedepth=15),
                     backend="cmdstanr",
                     file = "Models/newSun_notrtvar.RDS")



m_bs_wt_nosig <- brm(brmsformula(bodysize ~ sex + 
                                caretreatment_fct * scale(broodsize, scale=FALSE) +
                                  scale(wt) + 
                                (0 + caretreatment_fct||fnum)),
                  #sigma ~ 0 + caretreatment_fct),
                  data = as.data.frame(df_size),
                  iter=12000,
                  warmup=2000,
                  chains=6,
                  cores=6,
                  control=list(adapt_delta=0.99,max_treedepth=15),
                  backend="cmdstanr",
                  file = "Models/newSun_wt_nosig.RDS")

m_bs_wt_notrtRE <- brm(brmsformula(bodysize ~ sex + 
                                   caretreatment_fct * scale(broodsize, scale=FALSE) +
                                   scale(wt) + 
                                   (1|fnum),
                     sigma ~ 0 + caretreatment_fct),
                     data = as.data.frame(df_size),
                     iter=12000,
                     warmup=2000,
                     chains=6,
                     cores=6,
                     control=list(adapt_delta=0.99,max_treedepth=15),
                     backend="cmdstanr",
                     file = "Models/newSun_wt_notrtRE.RDS")

m_bs_wt_notrtvar <- brm(brmsformula(bodysize ~ sex + 
                                   caretreatment_fct * scale(broodsize, scale=FALSE) +
                                     scale(wt) + 
                                   (1|fnum)),
                     #sigma ~ 0 + caretreatment_fct),
                     data = as.data.frame(df_size),
                     iter=12000,
                     warmup=2000,
                     chains=6,
                     cores=6,
                     control=list(adapt_delta=0.99,max_treedepth=15),
                     backend="cmdstanr",
                     file = "Models/newSun_wt_notrtvar.RDS")


#### .Test the treatment effects ####

m_bs_testint <- brm(brmsformula(bodysize ~ sex + 
                                  caretreatment_fct + 
                                  scale(broodsize, scale=FALSE) +
                                  (0 + caretreatment_fct||fnum),
                                sigma ~ 0 + caretreatment_fct),
                    data = as.data.frame(df_size),
                    iter=12000,
                    warmup=2000,
                    chains=6,
                    cores=6,
                    control=list(adapt_delta=0.99,max_treedepth=15),
                    backend="cmdstanr",
                    file = "Models/newSun_testint.RDS")

m_bs_testtrt <- brm(brmsformula(bodysize ~ sex + 
                                  #caretreatment_fct + 
                                  scale(broodsize, scale=FALSE) +
                                  (0 + caretreatment_fct||fnum),
                                sigma ~ 0 + caretreatment_fct),
                    data = as.data.frame(df_size),
                    iter=12000,
                    warmup=2000,
                    chains=6,
                    cores=6,
                    control=list(adapt_delta=0.99,max_treedepth=15),
                    backend="cmdstanr",
                    file = "Models/newSun_testtrt.RDS")


m_bs_wt_testint <- brm(brmsformula(bodysize ~ sex + 
                                  caretreatment_fct + 
                                  scale(broodsize, scale=FALSE) +
                                    scale(wt) + 
                                  (0 + caretreatment_fct||fnum),
                                sigma ~ 0 + caretreatment_fct),
                    data = as.data.frame(df_size),
                    iter=12000,
                    warmup=2000,
                    chains=6,
                    cores=6,
                    control=list(adapt_delta=0.99,max_treedepth=15),
                    backend="cmdstanr",
                    file = "Models/newSun_wt_testint.RDS")

m_bs_wt_testtrt <- brm(brmsformula(bodysize ~ sex + 
                                  #caretreatment_fct + 
                                  scale(broodsize, scale=FALSE) +
                                    scale(wt) + 
                                  (0 + caretreatment_fct||fnum),
                                sigma ~ 0 + caretreatment_fct),
                    data = as.data.frame(df_size),
                    iter=12000,
                    warmup=2000,
                    chains=6,
                    cores=6,
                    control=list(adapt_delta=0.99,max_treedepth=15),
                    backend="cmdstanr",
                    file = "Models/newSun_wt_testtrt.RDS")



#### .Model inspection ####

summary(m_bs)
pp_check(m_bs, ndraws = 100)
check_brms(m_bs, integer = FALSE)
conditional_effects(m_bs)



summary(m_bs2)
pp_check(m_bs2, ndraws = 100)
check_brms(m_bs2, integer = FALSE)
conditional_effects(m_bs2)

summary(m_bs_ints)
pp_check(m_bs_ints, ndraws = 100)
check_brms(m_bs_ints, integer = FALSE)
conditional_effects(m_bs_ints)

## Including the interactions doesn't fit that hump

summary(m_bs_wt)
pp_check(m_bs_wt, ndraws = 100)
check_brms(m_bs_wt, integer = FALSE)


#### .Model comparisons ####


## compare models
m_bs <- add_criterion(m_bs, "loo")
m_bs2 <- add_criterion(m_bs2, "loo")
m_bs_ints <- add_criterion(m_bs_ints, "loo")


loo_compare(m_bs, m_bs2, m_bs_ints)
# elpd_diff se_diff
# brms_bodysize_m1  0.0       0.0   
# brms_bodysize_m2  0.0       1.3   
# brms_bodysize_m3 -6.7       4.3 

# Suggests models are not very different from one another
# - variances are basically unchanged across all

## compare models
m_bs_wt <- add_criterion(m_bs_wt, "loo")
m_bs2_wt <- add_criterion(m_bs2_wt, "loo")

loo_compare(m_bs, m_bs2, m_bs_ints,
            m_bs_wt, m_bs2_wt)



m_bs_nosig <- add_criterion(m_bs_nosig, "loo")
m_bs_notrtvar <- add_criterion(m_bs_notrtvar, "loo")


loo_compare(m_bs, m_bs2, m_bs_ints,
            m_bs_nosig, m_bs_notrtvar)


m_bs_testint <- add_criterion(m_bs_testint, "loo")
m_bs_testtrt <- add_criterion(m_bs_testtrt, "loo")

loo_compare(m_bs, m_bs2, 
            m_bs_testint, m_bs_testtrt)


m_bs_wt_testint <- add_criterion(m_bs_wt_testint, "loo")
m_bs_wt_testtrt <- add_criterion(m_bs_wt_testtrt, "loo")

loo_compare(m_bs_wt,
            m_bs_wt_testint, 
            m_bs_wt_testtrt)



m_bs_wt_nosig <- add_criterion(m_bs_wt_nosig, "loo")
m_bs_wt_notrtRE <- add_criterion(m_bs_wt_notrtRE, "loo")
m_bs_wt_notrtvar <- add_criterion(m_bs_wt_notrtvar, "loo")

loo_compare(m_bs_wt,
            m_bs_wt_nosig, 
            m_bs_wt_notrtRE,
            m_bs_wt_notrtvar)

##
# Get p-values if you want?
##

loo_comp <- loo_compare(list(
  no_int = loo(m_bs_wt_testint), 
  no_trt = loo(m_bs_wt_testtrt),
  base = loo(m_bs_wt)))

loo_comp

1 - pnorm(-loo_comp[3,1], loo_comp[2,2])


1 - pnorm(6.4, 4.4)

1 - pnorm(6.4, 4.4)


#### Plots ####

#### .Variances ####



draws_brms_bodysize <- m_bs_wt |> 
  gather_draws(sd_fnum__caretreatment_fctFC,
               sd_fnum__caretreatment_fctMC,
               sd_fnum__caretreatment_fctNC,
               b_sigma_caretreatment_fctFC,
               b_sigma_caretreatment_fctMC,
               b_sigma_caretreatment_fctNC) |>  
  separate(.variable,
           into = c("scale_type",
                    "group_level",
                    "ct",
                    "treatment")) |> 
  mutate(value_var = ifelse(scale_type == "sd",
                            .value^2,
                            exp(.value)^2))

draws_brms_bodysize_nice <- draws_brms_bodysize |> 
  mutate(treatment_nice = str_sub(treatment, -2, -1),
         group_level_nice = case_when(
           group_level == "fnum" ~ "Among-family",
           group_level == "sigma" ~ "Residual"
         )) |> 
  mutate(treatment_nice = fct(treatment_nice, levels = c("FC", "MC", "NC")))



df_bs_wt_drawsummaries <- draws_brms_bodysize_nice |> 
  select(.chain:.draw,
         treatment_nice,
         group_level_nice,
         value_var) |> 
  pivot_wider(names_from = group_level_nice,
              values_from = value_var) |> 
  group_by(treatment_nice) |> 
  median_qi(`Among-family`,
            `Residual`,
            vp = (`Among-family` + `Residual`),
            H2 = (`Among-family`/vp),
            .width = c(.95,.66))

df_bs_wt_drawsummaries

gg_var_H2 <- df_bs_wt_drawsummaries |> 
  select(treatment_nice, contains("H2"),
         .width, .point, .interval) |> 
  ggplot(aes(y = treatment_nice, x = `H2`,
             xmin = `H2.lower`,
             xmax = `H2.upper`,
             colour = treatment_nice)) +
  geom_pointinterval() +
  scale_colour_manual(values = custom_colors, guide = "none") +
  #coord_cartesian(xlim = c(0,1)) +
  xlim(0,1) +
  labs(x = expression(paste("Broad-sense heritability (", H^2, ")")),
       y = "Care treatment") +
  theme_classic(base_size = 18) +
  coord_flip()

gg_var_H2


gg_var_VE <- df_bs_wt_drawsummaries |> 
  select(treatment_nice, contains("Residual"),
         .width, .point, .interval) |> 
  ggplot(aes(y = treatment_nice, x = `Residual`,
             xmin = `Residual.lower`,
             xmax = `Residual.upper`,
             colour = treatment_nice)) +
  geom_pointinterval() +
  scale_colour_manual(values = custom_colors, guide = "none") +
  #coord_cartesian(xlim = c(0,1)) +
  xlim(0,0.32) +
  labs(x = expression(paste("Residual variance")),
       y = "Care treatment") +
  theme_classic(base_size = 18) +
  coord_flip()

gg_var_VF <- df_bs_wt_drawsummaries |> 
  select(treatment_nice, contains("Among-family"),
         .width, .point, .interval) |> 
  ggplot(aes(y = treatment_nice, x = `Among-family`,
             xmin = `Among-family.lower`,
             xmax = `Among-family.upper`,
             colour = treatment_nice)) +
  geom_pointinterval() +
  scale_colour_manual(values = custom_colors, guide = "none") +
  #coord_cartesian(xlim = c(0,1)) +
  xlim(0,0.32) +
  labs(x = expression(paste("Among-family variance")),
       y = "Care treatment") +
  theme_classic(base_size = 18) +
  coord_flip()

gg_var_VP <- df_bs_wt_drawsummaries |> 
  select(treatment_nice, contains("vp"),
         .width, .point, .interval) |> 
  ggplot(aes(y = treatment_nice, x = `vp`,
             xmin = `vp.lower`,
             xmax = `vp.upper`,
             colour = treatment_nice)) +
  geom_pointinterval() +
  scale_colour_manual(values = custom_colors, guide = "none") +
  #coord_cartesian(xlim = c(0,1)) +
  xlim(0,0.32) +
  labs(x = expression(paste("Phenotypic variance")),
       y = "Care treatment") +
  theme_classic(base_size = 18) +
  coord_flip()

gg_var_VF + gg_var_VE + gg_var_VP + 
  plot_layout(nrow = 1,
              # axes = "collect",
              guides = "collect") +
  plot_annotation(tag_levels = 'A')

ggsave("Figures/variances_triplet.pdf",
       width = 18,
       height = 6)

gg_var_VP + gg_var_VF + gg_var_VE + 
  plot_layout(nrow = 1,
              # axes = "collect",
              guides = "collect") +
  plot_annotation(tag_levels = 'A')

ggsave("Figures/variances_triplet-alt.pdf",
       width = 18,
       height = 6)


gg_var_VF + gg_var_VE + gg_var_VP + 
  plot_layout(nrow = 3,
              axes = "collect",
              guides = "collect") +
  plot_annotation(tag_levels = 'A')



#### .Predictions ####

# df_pred_draws <- df_size |> 
#   group_by(caretreatment_fct) |> 
#   modelr::data_grid(.model = m_bs_wt,
#                     broodsize = seq_range(broodsize, by = 1)) |> 
#   add_predicted_draws(m_bs_wt, re_formula = NA)
# 
# df_size |> 
#   group_by(caretreatment_fct) |> 
#   expand(broodsize,
#          sex = "f",
#          wt = mean(wt, na.rm = TRUE))

df_pred_draws <- df_size |> 
  group_by(caretreatment_fct) |> 
  expand(broodsize = seq_range(broodsize, by = 1),
         sex = "f") |> 
  ungroup() |> 
  add_column(wt = mean(df_size$wt, na.rm = TRUE)) |> 
  add_predicted_draws(m_bs_wt, re_formula = NA)



gg_bodysize_broodint <- df_pred_draws |> 
  ggplot(aes(x = broodsize, y = .prediction, 
             color = caretreatment_fct, fill = caretreatment_fct)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .66), alpha = 0.1) +
  stat_lineribbon(aes(y = .prediction), .width = c(0)) +
  geom_point(data = df_size, aes(y = bodysize), alpha = 0.7) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  scale_colour_manual(values = custom_colors, guide = "none") +
  labs(x = "Brood size",
       y = "Body size") +
  theme_classic(base_size = 18)

gg_bodysize_broodint

ggsave("Figures/bodysize_brmspreds.pdf", width = 7, height = 5)

gg_bodysize_broodint_meanhdci <- df_pred_draws |> 
  ggplot(aes(x = broodsize, y = .prediction, 
             color = caretreatment_fct, fill = caretreatment_fct)) +
  stat_lineribbon(aes(y = .prediction), point_interval = "mean_hdci",
                  .width = c(.95, .66), alpha = 0.1) +
  stat_lineribbon(aes(y = .prediction), point_interval = "mean_hdci",
                  .width = c(0)) +
  geom_point(data = df_size, aes(y = bodysize), alpha = 0.7) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  scale_colour_manual(values = custom_colors, guide = "none") +
  labs(x = "Brood size",
       y = "Body size") +
  theme_classic(base_size = 18)

gg_bodysize_broodint_meanhdci

# df_pred_draws_broodmn <- df_size |> 
#   group_by(caretreatment_fct) |> 
#   expand(broodsize = mean(broodsize, na.rm = TRUE),
#          sex = "f") |> 
#   ungroup() |> 
#   add_column(wt = mean(df_size$wt, na.rm = TRUE)) |> 
#   add_predicted_draws(m_bs_wt, re_formula = NA)

df_pred_draws_broodmn <- df_size |> 
  expand(caretreatment_fct,
         broodsize = mean(broodsize, na.rm = TRUE),
         sex = "f",
         wt = mean(wt, na.rm = TRUE)) |> 
  add_predicted_draws(m_bs_wt, re_formula = NA)


# df_pred_draws_broodmn |> 
#   ggplot(aes(x = caretreatment_fct, y = .prediction, 
#              color = caretreatment_fct, fill = caretreatment_fct)) +
#   stat_halfeye(aes(y = .prediction), colour = "black", alpha = 0.7) +
#   scale_fill_manual(values = custom_colors, guide = "none") +
#   scale_colour_manual(values = custom_colors, guide = "none") +
#   labs(x = "Care treatment",
#        y = "Body size") +
#   theme_classic(base_size = 18)

gg_bodysize_broodmean <- df_pred_draws_broodmn |> 
  ggplot(aes(x = caretreatment_fct, y = .prediction, 
             color = caretreatment_fct,
             fill = caretreatment_fct)) +
  #stat_interval(aes(y = .prediction), alpha = 0.2) +
  # geom_point(data = df_size, aes(y = bodysize), alpha = 0.2,
  #            position = position_jitter(width = 0.2)) +
  geom_dots(data = df_size, aes(y = bodysize), 
            alpha = 0.1, colour = NA, dotsize = 0.8,
            layout = "weave", side = "both") +
  geom_pointrange(data = df_size |> 
                    group_by(caretreatment_fct, fnum) |> 
                    summarise(bs_mean = mean(bodysize, na.rm = TRUE),
                              bs_min = min(bodysize, na.rm = TRUE),
                              bs_max = max(bodysize, na.rm = TRUE)),
                  aes(y = bs_mean,
                      ymin = bs_min,
                      ymax = bs_max),
                  position = position_jitter(width = 0.3),
                  alpha = 0.3) +
  stat_pointinterval(aes(y = .prediction),
                     point_size = 4) +
  scale_colour_manual(values = custom_colors, guide = "none") +
  scale_fill_manual(values = custom_colors, guide = "none") +
  labs(x = "Care treatment",
       y = "Body size") +
  theme_classic(base_size = 18)

gg_bodysize_broodmean

ggsave("Figures/bodysize_brmspreds_data.pdf", width = 9, height = 6)


gg_bodysize_broodmean_black <- df_pred_draws_broodmn |> 
  ggplot(aes(x = caretreatment_fct, y = .prediction, 
             color = caretreatment_fct,
             fill = caretreatment_fct)) +
  #stat_interval(aes(y = .prediction), alpha = 0.2) +
  # geom_point(data = df_size, aes(y = bodysize), alpha = 0.2,
  #            position = position_jitter(width = 0.2)) +
  geom_dots(data = df_size, aes(y = bodysize), 
            alpha = 0.1, colour = NA, dotsize = 0.8,
            layout = "weave", side = "both") +
  geom_pointrange(data = df_size |> 
                    group_by(caretreatment_fct, fnum) |> 
                    summarise(bs_mean = mean(bodysize, na.rm = TRUE),
                              bs_min = min(bodysize, na.rm = TRUE),
                              bs_max = max(bodysize, na.rm = TRUE)),
                  aes(y = bs_mean,
                      ymin = bs_min,
                      ymax = bs_max),
                  position = position_jitter(width = 0.3),
                  alpha = 0.5) +
  stat_pointinterval(aes(y = .prediction),
                     point_size = 4, colour = "black") +
  scale_colour_manual(values = custom_colors, guide = "none") +
  scale_fill_manual(values = custom_colors, guide = "none") +
  labs(x = "Care treatment",
       y = "Body size") +
  theme_classic(base_size = 18)

gg_bodysize_broodmean_black

ggsave("Figures/bodysize_brmspreds_data-bl.pdf", width = 9, height = 6)



df_size_summ <- df_size |> 
  group_by(caretreatment_fct, fnum) |> 
  summarise(bs_mean = mean(bodysize, na.rm = TRUE),
            bs_min = min(bodysize, na.rm = TRUE),
            bs_max = max(bodysize, na.rm = TRUE))

df_pred_draws_broodmn |> 
  ggplot(aes(x = caretreatment_fct, y = .prediction, 
             color = caretreatment_fct,
             fill = caretreatment_fct)) +
  stat_interval(aes(y = .prediction), alpha = 0.2) +
  stat_pointinterval(aes(y = .prediction)) +
  # geom_point(data = df_size, aes(y = bodysize), alpha = 0.2,
  #            position = position_jitter(width = 0.2)) +
  geom_dots(data = df_size, aes(y = bodysize), 
            alpha = 0.1, colour = NA, dotsize = 0.8,
            layout = "weave", side = "both") +
  geom_linerange(data = df_size_summ,
                  aes(y = bs_mean,
                      ymin = bs_min,
                      ymax = bs_max),
                  position = position_jitter(width = 0.3),
                  alpha = 0.3) +
  geom_point(data = df_size_summ,
                  aes(y = bs_mean,),
                  position = position_jitter(width = 0.3),
                  alpha = 0.6) +
  scale_colour_manual(values = custom_colors, guide = "none") +
  scale_fill_manual(values = custom_colors, guide = "none") +
  labs(x = "Care treatment",
       y = "Body size") +
  theme_classic(base_size = 18)


#### .Plots together ####

gg_bodysize_broodmean_black +
  gg_bodysize_broodint +
  gg_var_H2 +
  plot_layout(nrow = 1,
              axes = "collect",
              guides = "collect") +
  plot_annotation(tag_levels = 'A')

ggsave("Figures/bodysize_triplet.pdf",
       width = 18,
       height = 6)


gg_bodysize_broodmean_black +
  gg_var_H2 +
  plot_layout(nrow = 1,
              axes = "collect",
              guides = "collect") +
  plot_annotation(tag_levels = 'A')

ggsave("Figures/bodysize_double.pdf",
       width = 12,
       height = 6)

gg_bodysize_broodint

ggsave("Figures/bodysize_int.pdf",
       width = 6,
       height = 6)


#### .Export tables as csv ####

tidy(m_bs_wt, effects = "fixed", conf.method="HPDinterval") |> 
  filter(str_sub(term,1,5) != "sigma") |> 
  write_csv("Tables/Model_summary.csv")

tidy(m_bs_wt, effects = "ran_pars", conf.method="HPDinterval") |> 
  write_csv("Tables/Model_varcomp_sd.csv")
