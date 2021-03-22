## ---- out.width = "50%", echo = FALSE, fig.align = "center", fig.cap = "Dream pet dragon"----------------------------------------------------------------------------------------------
knitr::include_graphics("images/fun_dragon.jpg")


## ----loadlibs_mm, message=FALSE, results='hide', warning=FALSE-------------------------------------------------------------------------------------------------------------------------

library(lmerTest)
library(tidyverse)
library(asreml)
library(MCMCglmm)
library(nadiv)


## ---- out.width = "50%", echo = FALSE, fig.align = "center", fig.cap = "Blue dragon male"----------------------------------------------------------------------------------------------
knitr::include_graphics("images/blue_dragon.jpg")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_dragons <- read.csv("data/dragons.csv")
str(df_dragons)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_dragons <- df_dragons %>%
  mutate(
    body_size_sc = scale(body_size),
    assay_rep_sc = scale(assay_rep, scale = FALSE)
  )


## ---- fig.cap = "Checking assumptions of model lmer_f"---------------------------------------------------------------------------------------------------------------------------------
lmer_f <- lmer(max_speed ~ assay_rep_sc + body_size_sc + (1 | ID),
  data = df_dragons
)
par(mfrow = c(1, 3))
plot(resid(lmer_f, type = "pearson") ~ fitted(lmer_f))
qqnorm(residuals(lmer_f))
qqline(residuals(lmer_f))
hist(residuals(lmer_f))
summary(lmer_f)


## ---- results = "hide", cache = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
rep_flying <- as.data.frame(VarCorr(lmer_f)) %>%
  select(grp, vcov) %>%
  spread(grp, vcov) %>%
  mutate(repeatability = ID / (ID + Residual))
rep_flying


## ---- echo = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::kable(rep_flying,
  digits = 3,
  caption = "Variance components and repeatbility for the maximum
   flying speed of blue dragons"
)


## ---- fig.cap = "Checking assumptions of model lmer_e"---------------------------------------------------------------------------------------------------------------------------------
lmer_e <- lmer(exploration ~ assay_rep_sc + body_size_sc + (1 | ID),
  data = df_dragons
)
par(mfrow = c(1, 3))
plot(resid(lmer_e, type = "pearson") ~ fitted(lmer_e))
qqnorm(residuals(lmer_e))
qqline(residuals(lmer_e))
hist(residuals(lmer_e))
summary(lmer_e)


## ---- results = "hide", cache = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
rep_expl <- as.data.frame(VarCorr(lmer_e)) %>%
  select(grp, vcov) %>%
  spread(grp, vcov) %>%
  mutate(repeatability = ID / (ID + Residual))
rep_expl


## ---- echo = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::kable(rep_flying,
  digits = 3,
  caption = "Variance components and repeatability for exploration
   behaviour of blue dragons"
)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_blups_fe <- merge(
  as.data.frame(ranef(lmer_f)),
  as.data.frame(ranef(lmer_e)),
  by = "grp"
) %>%
  mutate(
    speed = condval.x,
    exploration = condval.y
  )


## ---- fig.cap = "Relation between exploration and flying speed using BLUPs from univariate models"-------------------------------------------------------------------------------------
(cor_blups <- with(df_blups_fe, cor.test(condval.x, condval.y)))
ggplot(df_blups_fe, aes(x = exploration, y = speed)) +
  geom_point() +
  labs(xlab = "Exploration (BLUP)", ylab = "Flying speed (BLUP)") +
  theme_classic()


## ---- fig.cap = "Relation between exploration and flying speed using BLUPs from univariate models including +/- SE as error bars"------------------------------------------------------
ggplot(df_blups_fe, aes(x = exploration, y = speed)) +
  geom_point() +
  geom_linerange(aes(
    xmin = exploration - condsd.x,
    xmax = exploration + condsd.x
  )) +
  geom_linerange(aes(
    ymin = speed - condsd.y,
    ymax = speed + condsd.y
  )) +
  labs(
    xlab = "Exploration (BLUP +/- SE)",
    ylab = "Flying speed (BLUP +/- SE)"
  ) +
  theme_classic()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_dragons <- df_dragons %>%
  mutate(
    ID = as.factor(ID),
    speed_sc = scale(max_speed),
    exploration_sc = scale(exploration)
  )

asr_us <- asreml(cbind(speed_sc, exploration_sc) ~ trait +
  trait:assay_rep_sc + trait:body_size_sc,
random = ~ ID:us(trait),
residual = ~ units:us(trait),
data = df_dragons,
maxiter = 100
)


## ---- fig.cap = "Checking assumptions of model asr_us"---------------------------------------------------------------------------------------------------------------------------------
par(mfrow = c(1, 3))
plot(residuals(asr_us) ~ fitted(asr_us))
qqnorm(residuals(asr_us))
qqline(residuals(asr_us))
hist(residuals(asr_us))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(asr_us, coef = TRUE)$coef.fixed
wald(asr_us, ssType = "conditional", denDF = "numeric")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(asr_us)$varcomp


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
vpredict(asr_us, rep_speed ~ V1 / (V1 + V5))
vpredict(asr_us, rep_expl ~ V3 / (V3 + V7))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
(cor_fe <- vpredict(asr_us, rep_expl ~ V2 / (sqrt(V1 * V3))))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
asr_us <- asreml(cbind(speed_sc, exploration_sc) ~ trait +
  trait:assay_rep_sc + trait:body_size_sc,
random = ~ ID:idh(trait),
residual = ~ units:us(trait),
data = df_dragons,
maxiter = 100
)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
(p_biv <- pchisq(2 * (asr_us$loglik - asr_us$loglik),
  df = 1,
  lower.tail = FALSE
))


## ---- fig.cap = "Correlation estimates (with CI) using 2 different methods"------------------------------------------------------------------------------------------------------------
df_cor <- data.frame(
  Method = c("ASReml", "BLUPs"),
  Correlation = c(as.numeric(cor_fe[1]), cor_blups$estimate),
  low = c(as.numeric(cor_fe[1] - 1.96 * cor_fe[2]), cor_blups$conf.int[1]),
  high = c(as.numeric(cor_fe[1] + 1.96 * cor_fe[2]), cor_blups$conf.int[2])
)
ggplot(df_cor, aes(x = Method, y = Correlation)) +
  geom_point() +
  geom_linerange(aes(ymin = low, ymax = high)) +
  ylim(-1, 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
prior_1ex <- list(
  R = list(V = diag(2), nu = 0.002),
  G = list(G1 = list(
    V = diag(2) * 0.02, nu = 3,
    alpha.mu = rep(0, 2),
    alpha.V = diag(1000, 2, 2)
  ))
)


## ---- message = FALSE, warning = FALSE, cache = TRUE-----------------------------------------------------------------------------------------------------------------------------------
mcmc_us <- MCMCglmm(cbind(speed_sc, exploration_sc) ~ trait - 1 +
    trait:assay_rep_sc +
    trait:body_size_sc,
  random = ~ us(trait):ID,
  rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian"),
  prior = prior_1ex,
  nitt = 420000,
  burnin = 20000,
  thin = 100,
  verbose = FALSE,
  data = df_dragons
)


## ---- fig.cap = "MCMC trace and Posterior distribution of the (co)variance estimates of model mcmc_us"---------------------------------------------------------------------------------
plot(mcmc_us$VCV[, c(1, 2, 4)])
plot(mcmc_us$VCV[, c(5, 6, 8)])


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(mcmc_us)


## ---- fig.cap = "Posterior trace and distribution of the repeatability in flying speed"------------------------------------------------------------------------------------------------
mcmc_prop_f <- mcmc_us$VCV[, 1] /
  (mcmc_us$VCV[, 1] + mcmc_us$VCV[, 5])
plot(mcmc_prop_f)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
posterior.mode(mcmc_prop_f)
HPDinterval(mcmc_prop_f)


## ---- fig.cap = "Posterior trace and distribution of the repeatbility of exploration"--------------------------------------------------------------------------------------------------
mcmc_prop_e <- mcmc_us$VCV[, 4] /
  (mcmc_us$VCV[, 4] + mcmc_us$VCV[, 8])
plot(mcmc_prop_e)
posterior.mode(mcmc_prop_e)
HPDinterval(mcmc_prop_e)


## ---- fig.cap = "Posterior trace and distribution of the correlation between flying speed and exploration"-----------------------------------------------------------------------------
mcmc_cor_fe <- mcmc_us$VCV[, 2] /
  sqrt(mcmc_us$VCV[, 1] * mcmc_us$VCV[, 4])
plot(mcmc_cor_fe)
posterior.mode(mcmc_cor_fe)
HPDinterval(mcmc_cor_fe)


## ---- fig.cap = "Correlation estimates (with CI) using 3 different methods"------------------------------------------------------------------------------------------------------------
df_cor[3, 1] <- "MCMCglmm"
df_cor[3, -1] <- c(posterior.mode(mcmc_cor_fe), HPDinterval(mcmc_cor_fe))

ggplot(df_cor, aes(x = Method, y = Correlation)) +
  geom_point() +
  geom_linerange(aes(ymin = low, ymax = high)) +
  ylim(-1, 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic()


## ---- out.width = "50%", echo = FALSE, fig.align = "center", fig.cap = "A female blue dragon of the West"------------------------------------------------------------------------------
knitr::include_graphics("images/blue_dragon.jpg")

