## ---- eval = TRUE, warning = FALSE, message = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
library(lme4)
library(plyr)
library(tidyverse)
library(patchwork)
library(lattice)
source("data/glmm_funs.R")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dat_tf <- read.csv("data/Banta_TotalFruits.csv")
str(dat_tf)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dat_tf <- mutate(
  dat_tf,
  X = factor(X),
  gen = factor(gen),
  rack = factor(rack),
  amd = factor(amd, levels = c("unclipped", "clipped")),
  nutrient = factor(nutrient, label = c("Low", "High"))
)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
(reptab <- with(dat_tf, table(popu, gen)))


## **Exercise**: this mode of inspection is OK for this data set but might fail for much larger data sets or for more levels of nesting. See if you can think of some other numerical or graphical methods for inspecting the structure of data sets.

## 
## 1. plot(reptab) gives a mosaic plot of the two-way table; examine this, see if you can figure out how to interpret it, and decide whether you think it might be useful

## 2. try the commands colSums(reptab>0) (and the equivalent for rowSums) and figure out what they are telling you.

## 3. Using this recipe, how would you compute the range of number of genotypes per treatment combination?

## 

## ---- fig.cap ="A truly useless plot no one can understand"------------------------------------------------------------------------------------------------------------------------------------------
plot(reptab)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
colSums(reptab > 0)
rowSums(reptab > 0)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
reptab2 <- with(dat_tf, table(paste(amd, nutrient, sep = "_"), gen))
range(reptab2)


## ---- echo = FALSE, fig.cap = "Number of fruits (log + 1) as a function of treatments"---------------------------------------------------------------------------------------------------------------
p1 <- qplot(
    interaction(nutrient, amd),
    log(1 + total.fruits),
    data = dat_tf, geom = "boxplot") +
  facet_wrap(~reg, nrow = 1) +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Boxplot")
p2 <- qplot(
    interaction(nutrient, amd),
    log(1 + total.fruits),
    data = dat_tf) +
  facet_wrap(~reg, nrow = 1) +
  stat_sum() +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Dot plot")
p1 + p2


## **Exercise** generate these plots and figure out how they work before continuing. Try conditioning/faceting on population rather than region: for facet_wrap you might want to take out the nrow=1 specification. If you want try reorder the subplots by overall mean fruit set and/or colour the points according to the region they come from.


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dat_tf <- dat_tf %>%
  mutate(
    gna = interaction(gen, nutrient, amd),
    gna = reorder(gna, total.fruits, mean)
  )


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(dat_tf, aes(x = gna, y = log(1 + total.fruits))) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
grp_means <- with(dat_tf, tapply(total.fruits, list(gna), mean))
summary(grp_means)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
grp_vars <- with(
  dat_tf,
  tapply(
    total.fruits,
    list(gna), var
  )
)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lm1 <- lm(grp_vars ~ grp_means - 1) ## `quasi-Poisson' fit
phi_fit <- coef(lm1)
lm2 <- lm((grp_vars - grp_means) ~ I(grp_means^2) - 1)
k_fit <- 1 / coef(lm2)


## ---- fig.cap = "Graphical evaluation of distribution to use"----------------------------------------------------------------------------------------------------------------------------------------
plot(grp_vars ~ grp_means, xlab = "group means", ylab = "group variances")
abline(c(0, 1), lty = 2)
text(105, 500, "Poisson")
curve(phi_fit * x, col = 2, add = TRUE)
## bquote() is used to substitute numeric values
## in equations with symbols
text(110, 3900,
  bquote(paste("QP: ", sigma^2 == .(round(phi_fit, 1)) * mu)),
  col = 2
)
curve(x * (1 + x / k_fit), col = 4, add = TRUE)
text(104, 7200, paste("NB: k=", round(k_fit, 1), sep = ""), col = 4)
l_fit <- loess(grp_vars ~ grp_means)
mvec <- 0:120
lines(mvec, predict(l_fit, mvec), col = 5)
text(100, 2500, "loess", col = 5)


## ---- fig.cap = "Graphical evaluation of distribution to use with ggplot"----------------------------------------------------------------------------------------------------------------------------
ggplot(
  data.frame(grp_means, grp_vars),
  aes(x = grp_means, y = grp_vars)) +
  geom_point() +
  geom_smooth(
    aes(colour = "Loess"), se = FALSE) +
  geom_smooth(
    method = "lm", formula = y ~ x - 1, se = FALSE,
    aes(colour = "Q_Pois")) +
  stat_function(
    fun = function(x) x * (1 + x / k_fit),
    aes(colour = "Neg_bin")
  ) +
  geom_abline(
    aes(intercept = 0, slope = 1, colour = "Poisson")) +
  scale_colour_manual(
    name = "legend",
    values = c("blue", "purple", "black", "red")) +
  scale_fill_manual(
    name = "legend",
    values = c("blue", "purple", "black", "red")) +
  guides(fill = FALSE)


## ---- fig.cap = "Graphical evaluation of distribution to use including quadratic effect"-------------------------------------------------------------------------------------------------------------
lm3 <- lm(grp_vars ~ I(grp_means)^2 - 1) ## quadratic fit
quad_fit <- coef(lm3)

ggplot(
  data.frame(grp_means, grp_vars),
  aes(x = grp_means, y = grp_vars)) +
  geom_point() +
  geom_smooth(
    method = "lm", formula = y ~ x - 1, se = FALSE,
    aes(colour = "Q_Pois")) +
  stat_function(
    fun = function(x) x * (1 + x / k_fit),
    aes(colour = "Neg_bin")
  ) +
  geom_smooth(
    method = "lm", formula = y ~ I(x^2) - 1, se = FALSE,
    aes(colour = "Quad")) +
  scale_colour_manual(
    name = "legend",
    values = c("blue", "purple", "black")) +
  scale_fill_manual(
    name = "legend",
    values = c("blue", "purple", "black")) +
  guides(fill = FALSE)


## ---- fig.cap = "Fruit production by treatments by population"---------------------------------------------------------------------------------------------------------------------------------------
ggplot(dat_tf, aes(x = amd, y = log(total.fruits + 1), colour = nutrient)) +
  geom_point() +
  ## need to use as.numeric(amd) to get lines
  stat_summary(aes(x = as.numeric(amd)), fun = mean, geom = "line") +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines")) +
  facet_wrap(~popu)


## ---- fig.cap = "Fruit production by genotype by treatments"-----------------------------------------------------------------------------------------------------------------------------------------
ggplot(dat_tf, aes(x = amd, y = log(total.fruits + 1), colour = gen)) +
  geom_point() +
  stat_summary(aes(x = as.numeric(amd)), fun = mean, geom = "line") +
  theme_bw() +
  ## label_both adds variable name ('nutrient') to facet labels
  facet_grid(. ~ nutrient, labeller = label_both)


## ---- fig.cap = "Model coefficients for GLM fits on each genotype"-----------------------------------------------------------------------------------------------------------------------------------
glm_lis <- lmList(
  total.fruits ~ nutrient * amd | gen,
  data = dat_tf,
  family = "poisson")
plot.lmList(glm_lis)


## ---- fig.cap = "Q-Q plots of model coefficients for GLM fits on each genotype"----------------------------------------------------------------------------------------------------------------------
qqmath.lmList(glm_lis)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mp1 <- glmer(total.fruits ~ nutrient * amd +
  rack + status +
  (amd * nutrient | popu) +
  (amd * nutrient | gen),
data = dat_tf, family = "poisson"
)
overdisp_fun(mp1)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mp2 <- update(mp1, . ~ . + (1 | X))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
attr(VarCorr(mp2)$gen, "correlation")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mp3 <- glmer(total.fruits ~ nutrient * amd +
  rack + status +
  (amd + nutrient | popu) +
  (amd + nutrient | gen) + (1 | X),
data = dat_tf, family = "poisson"
)

attr(VarCorr(mp3)$gen, "correlation")
attr(VarCorr(mp3)$popu, "correlation")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mp_obs <- glmer(total.fruits ~ nutrient * amd +
  rack + status +
  (1 | X),
data = dat_tf, family = "poisson"
)


## **Exercise** using update, fit the models with

## 
## 1. clipping variation at both genotype and population levels;

## 2. nutrient variation at both genotype and populations; convince yourself that trying to fit variation in either clipping or nutrients leads to overfitting (perfect correlations).

## 3. Fit the model with only intercept variation at the population and genotype levels, saving it as mp4; show that there is non-zero variance estimated

## 

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mpcli <- update(mp_obs, . ~ . + (amd | gen) + (amd | popu))
VarCorr(mpcli)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mpnut <- update(mp_obs, . ~ . + (nutrient | gen) + (nutrient | popu))
VarCorr(mpnut)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mp4 <- update(mp_obs, . ~ . + (1 | gen) + (1 | popu))
VarCorr(mp4)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
overdisp_fun(mp4)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
VarCorr(mp4)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mp4v1 <- update(mp_obs, . ~ . + (1 | popu)) ## popu only (drop gen)
mp4v2 <- update(mp_obs, . ~ . + (1 | gen)) ## gen only (drop popu)
anova(mp4, mp4v1)
anova(mp4, mp4v2)


## ----simdev_glmm, eval = FALSE, echo = TRUE----------------------------------------------------------------------------------------------------------------------------------------------------------
## simdev <- function() {
##   newdat <- simulate(mp4v1)
##   reduced <- lme4::refit(mp4v1, newdat)
##   full <- lme4::refit(mp4, newdat)
##   2 * (c(logLik(full) - logLik(reduced)))
## }
## 
## set.seed(101)
## ## raply in plyr is a convenient wrapper for repeating the simulation many times
## nulldist <- raply(200, simdev(), .progress = "text")
## ## zero spurious (small) negative values
## nulldist[nulldist < 0 & abs(nulldist) < 1e-5] <- 0
## obsdev <- 2 * c(logLik(mp4) - logLik(mp4v1))





## ---- eval = TRUE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mean(c(nulldist, obsdev) >= obsdev)


## ---- fig.cap = "Distribution of BLUPs for genotypes and populations"--------------------------------------------------------------------------------------------------------------------------------
r1 <- as.data.frame(ranef(mp4, condVar = TRUE, whichel = c("gen", "popu")))
p1 <- ggplot(subset(r1, grpvar == "gen"), aes(y = grp, x = condval)) +
  geom_point() +
  geom_pointrange(
    aes(xmin = condval - condsd * 1.96, xmax = condval + condsd * 1.96)
  ) +
  geom_vline(aes(xintercept = 0, color = "red")) +
  theme_classic() +
  theme(legend.position = "none")
p2 <- ggplot(subset(r1, grpvar == "popu"), aes(y = grp, x = condval)) +
  geom_point() +
  geom_pointrange(
    aes(xmin = condval - condsd * 1.96, xmax = condval + condsd * 1.96)
  ) +
  geom_vline(aes(xintercept = 0, color = "red")) +
  theme_classic() +
  theme(legend.position = "none")
p1 + p2


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
(dd_aic <- dfun(drop1(mp4)))
(dd_lrt <- drop1(mp4, test = "Chisq"))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mp5 <- update(mp4, . ~ . - amd:nutrient)
anova(mp5, mp4)


## ---- fig.cap = "Q-Q plot of BLUPs from model mp5"---------------------------------------------------------------------------------------------------------------------------------------------------
r5 <- as.data.frame(ranef(mp5))
ggplot(data = r5, aes(sample = condval)) +
  geom_qq() + geom_qq_line() +
  facet_wrap(~ grpvar) +
  theme_classic()


## **Exercise**

## 
## - Re-do the analysis with region as a fixed effect.

## - Re-do the analysis with a one-way layout as suggested above

## 

## ---- out.width = "50%", echo = FALSE, fig.align = "center", fig.cap = "A GLMM character"------------------------------------------------------------------------------------------------------------
knitr::include_graphics("images/Thorn.png")

