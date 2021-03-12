## ----loadlibs, message=FALSE, results='hide', warning=FALSE------------------------------------------------------------------------------------------------------------------------------------------

library(lmerTest)
library(tidyverse)
library(rptR)


## ---- out.width="50%", echo = FALSE, fig.align="center", fig.cap = "The superb unicorn of the Scottish Highlands"------------------------------------------------------------------------------------
knitr::include_graphics("images/unicorn.png")


## ----load_data_aggr, message=FALSE, results='hide', warning=FALSE------------------------------------------------------------------------------------------------------------------------------------

unicorns <- read.csv("data/unicorns_aggression.csv")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
str(unicorns)
summary(unicorns)


## ----mod1--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

m_1 <- lmer(aggression ~ 1 +  (1 | ID), data = unicorns)



## ----mod1_summary------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

summary(m_1)



## ----mod1_ranova-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ranova(m_1)



## ----mod1_plot, fig.cap = "Fitted values vs residuals for a simple mixed model of unicorn aggression"------------------------------------------------------------------------------------------------

plot(m_1)





## ----rplotaggr, eval = TRUE, warning = FALSE, fig.cap = "Unicorn aggressivity as a function of opponent size when fighting for sweets", fig.align='center'-------------------------------------------
ggplot(unicorns, aes(x = opp_size, y = aggression)) +
  geom_jitter(
    alpha = 0.5,
    width = 0.05
  ) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  labs(
    x = "Opponent size (SD)",
    y = "Aggression"
  ) +
  theme_classic()


## ----mean_aggr, echo=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

unicorns %>%
  group_by(opp_size) %>%
  summarise(mean_aggr = mean(aggression)) %>%
  knitr::kable(digits = 2)



## ----mod2--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

m_2 <- lmer(aggression ~ opp_size  + block + (1 | ID), data = unicorns)



## ---- fig.cap = "Checking residuals have Gaussian distribution"--------------------------------------------------------------------------------------------------------------------------------------
par(mfrow = c(1, 2)) # multiple graphs in a window
qqnorm(residuals(m_2)) # a q-q plot
qqline(residuals(m_2))
hist(resid(m_2)) # are the residuals roughly Gaussian?


## ---- fig.cap = "Residuals by fitted values for model m_2 to check homoscedasticity"-----------------------------------------------------------------------------------------------------------------
plot(m_2)


## ---- fig.cap = "Checking random effects are gaussian"-----------------------------------------------------------------------------------------------------------------------------------------------
# extracting blups
r1 <- as.data.frame(ranef(m_2, condVar = TRUE))

par(mfrow = c(1, 2))
hist(r1$condval)
qqnorm(r1$condval)
qqline(r1$condval)


## ----mod2_plots, eval=TRUE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
r1 <- r1[order(r1$condval), ] # sorting the BLUPs
ggplot(r1, aes(y = grp, x = condval)) +
  geom_point() +
  geom_pointrange(
    aes(xmin = condval - condsd * 1.96, xmax = condval + condsd * 1.96)
  ) +
  geom_vline(aes(xintercept = 0, color = "red")) +
  theme_classic() +
  theme(legend.position = "none")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(m_2)


## **Try tweaking the fixed part of your model:**

## 
## - What happens if you add more fixed effects? Try it!

## - Could focal body size also matter? If so, should you rescale before adding it to the model?

## - Should you add interactions (e.g. block:opp_size)?

## - Should you drop non-significant fixed effects?

## 

## **Having changed the fixed part of your model, do the variance estimates change at all?**

## 
## * Is among-individual variance always estimated as zero regardless of fixed effects?

## * Is among-individual variance significant with some fixed effets structures but not others?

## 

## ----mod3--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

m_3 <- lmer(
  aggression ~ opp_size + scale(body_size)
              + scale(assay_rep, scale = FALSE) + block
              + (1 | ID),
  data = unicorns)
summary(m_3)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
0.02538 / (0.02538 + 0.58048)


## Do you see where I took the numbers ?


## ----mod3_Rcalc--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
v_id <- VarCorr(m_3)$ID[1, 1]
v_r <- attr(VarCorr(m_3), "sc")^2
r_man <- v_id / (v_id + v_r)
r_man


## ---- echo = TRUE, results =  TRUE, warning = FALSE, message = FALSE---------------------------------------------------------------------------------------------------------------------------------
r_rpt <- rptGaussian(
  aggression ~ opp_size + block + (1 | ID),
  grname = "ID", data = unicorns)
r_rpt


## ----mod_size_wrong----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

lmer_size <- lmer(body_size ~ block + (1 | ID),
                data = unicorns)



## **What might you conclude, and why would this be foolish?**


## **How can we do it properly?**


## ----mod_size_right----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

unicorns2 <- unicorns[unicorns$assay_rep == 1, ]
lmer_size2 <- lmer(body_size ~   block +  (1 | ID),
                data = unicorns2)



## ---- out.width = "20%", echo = FALSE, fig.align = "center", fig.cap = "The superb unicorn"----------------------------------------------------------------------------------------------------------
knitr::include_graphics("images/unicorn.png")

