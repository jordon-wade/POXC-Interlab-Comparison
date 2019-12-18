library(dabestr)
library(ggplot2)
library(lme4)
library(emmeans)
library(tidyr)
library(car)
library(plyr)
library(dplyr)
library(raster)
library(Rmisc)

### LOAD IN DATA ###
## Raw POXC data ##
POXC.dat <- read.csv("Wade Margenot/POXC interlab comparison/Manuscript(s) Things/Published data and code/POXC.comparison_published.csv", header=TRUE)
POXC.dat[1] <- lapply(POXC.dat[1], as.numeric)
POXC.dat[2:7] <- lapply(POXC.dat[2:7], as.factor)
POXC.dat[8:10] <- lapply(POXC.dat[8:10], as.numeric)
POXC.dat$Detect <- ifelse(POXC.dat$POXC == POXC.dat$POXC.raw, "1", "0")
POXC.dat$Detect <- replace_na(POXC.dat$Detect, 0)
POXC.data <- as.data.frame(POXC.dat[,c(2:8,10:11)])
str(POXC.dat) # check structure, "Detect" should be chr

## Raw soil characterization data ##
Soil.data <- read.csv("Wade Margenot/POXC interlab comparison/Manuscript(s) Things/Published data and code/POXC.comparison_soil data_published.csv", header=TRUE)

## Wrangling analytical variability Data ##
POXC.CV <- as.data.frame(POXC.data %>% dplyr::group_by(Soil.Order, Soil.ID, Lab.ID, Mass.class, Sieve.Size) %>% dplyr::summarise(avg=mean(POXC), stdev=sd(POXC), CV.pct=cv(POXC), na.rm=TRUE))
POXC.CV$log.CV.pct <- log(POXC.CV$CV.pct)

## Equation 2: used to calculate all estimated marginal means ##
# Deriving it: should we use Soil.Order or Soil.ID? #
CV.sources1.1 <- lm(log.CV.pct ~ Lab.ID + Soil.Order*Mass.class*Particle.Size, data=POXC.CV)
plot(resid(CV.sources1))
shapiro.test(resid(CV.sources1))
Anova(CV.sources1, type=2)
CV.sources1.1 <- lm(log.CV.pct ~ Lab.ID + Soil.Order + Mass.class + Particle.Size + Soil.Order:Mass.class:Particle.Size, data=POXC.CV)
plot(resid(CV.sources1.1))
shapiro.test(resid(CV.sources1.1))
Anova(CV.sources1.1, type=2)
# final model w/ Soil.Order
CV.sources1.2 <- lm(log.CV.pct ~ Lab.ID + Soil.Order + Mass.class + Particle.Size + Soil.Order:Mass.class:Particle.Size, data=POXC.CV)
plot(resid(CV.sources1.2))
shapiro.test(resid(CV.sources1.2))
Anova(CV.sources1.2, type=2)

CV.sources2 <- lm(log.CV.pct ~ Lab.ID + Soil.ID*Mass.class*Particle.Size, data=POXC.CV)
plot(resid(CV.sources2))
shapiro.test(resid(CV.sources2))
Anova(CV.sources2, type=2)
CV.sources2.1 <- lm(log.CV.pct ~ Lab.ID + Soil.ID + Mass.class + Particle.Size + Soil.ID:Particle.Size, data=POXC.CV)
plot(resid(CV.sources2.1))
shapiro.test(resid(CV.sources2.1))
Anova(CV.sources2.1, type=2)

anova(CV.sources1.2, CV.sources2.1)
library(lmtest)
lrtest(CV.sources1.2, CV.sources2.1)
detach("package:lmtest", unload=TRUE)
AIC(CV.sources1.2, CV.sources2.1)
# Soil.ID is the more parsimonious model

# THIS IS THE ONE WE USE FOR EVERYTHING! #
CV.sources <- lm(log.CV.pct ~ Lab.ID + Soil.ID + Mass.class + Sieve.Size + Soil.ID:Sieve.Size, data=POXC.CV)
plot(resid(CV.sources))
shapiro.test(resid(CV.sources)) # W > 0.95
plot(CV.sources) # check visually

### TABLES ###
## Table 1: Physicochemical properties ##
Table.1 <- Soil.data[,c(1:3)]
Table.1$SOC.pct <- round(Soil.data$SOC.pct, 1)
Table.1$C.N <- round(Soil.data$C.N, 1)
Table.1$pH <- Soil.data$pH
Table.1$Clay.mass <- round(Soil.data$Clay.pct*10, 0)
Table.1$Sand.mass <- round(Soil.data$Sand.pct*10, 0)
Table.1$CEC <- round(Soil.data$CEC, 1)
Table.1

## Table 2: Central tendencies/deviation of absolute values and analytical variability ##
AbsoluteValues.summary <- as.data.frame(POXC.data %>% dplyr::group_by(Soil.ID, Mass.class, Sieve.Size) %>% dplyr::summarize(avg=mean(POXC, na.rm=TRUE), med = median(POXC, na.rm=TRUE), abs.dev=mad(POXC, na.rm=TRUE)))
AbsoluteValues.summary <- AbsoluteValues.summary %>% mutate_at(vars(avg, med, abs.dev), funs(round(., 0))) # values for central tendency measures

CV.summary <- as.data.frame(POXC.CV %>% dplyr::group_by(Soil.ID, Mass.class, Sieve.Size) %>% dplyr::summarize(avg=mean(CV.pct, na.rm=TRUE)))
CV.summary # analytical variability values

## Table 3: Detection ##
Detection <- as.data.frame(POXC.data %>% dplyr::group_by(Sieve.Size, Mass.class, .drop=TRUE) %>% dplyr::count(Detect) %>% spread(key=Detect, value=n), na.rm=TRUE)
Detection$Pct.detectable <- 100*(Detection$`1`/(Detection$`0` + Detection$`1`))
Detection %>% mutate_at(vars(Pct.detectable), funs(round(., 1)))

## Table 4: Sources of variation ##
Anova(CV.sources, type=2)
summary(CV.sources)

## Table 5: conversions between treatments ##
# Wrangle data #
AbsoluteValues.summary.wide <- pivot_wider(AbsoluteValues.summary, names_from=c("Sieve.Size", "Mass.class"), values_from=c("avg", "med", "abs.dev"))
# only inlcude soils SOC < 10% #
AbsoluteValues.summary.wide.reduced <- AbsoluteValues.summary.wide[-c(8,9,34,35,36),]
# Calculate equations row by row #
Conversion.2mm <- lm(AbsoluteValues.summary.wide.reduced$`avg_<2mm_2.5g` ~ AbsoluteValues.summary.wide.reduced$`avg_<2mm_0.75g`, data=AbsoluteValues.summary.wide.reduced)
round(Conversion.2mm$coef, 2) # regression coefficients
sqrt(mean(Conversion.2mm$residuals^2)) # RMSE
summary(Conversion.2mm)$r.squared #R^2

Conversion.0.5mm <- lm(AbsoluteValues.summary.wide.reduced$`avg_<0.5mm_2.5g` ~ AbsoluteValues.summary.wide.reduced$`avg_<0.5mm_0.75g`, data=AbsoluteValues.summary.wide.reduced)
round(Conversion.0.5mm$coef, 2) # regression coefficients
sqrt(mean(Conversion.0.5mm$residuals^2)) # RMSE
summary(Conversion.0.5mm)$r.squared #R^2

Conversion.0.75g <- lm(AbsoluteValues.summary.wide.reduced$`avg_<2mm_0.75g` ~ AbsoluteValues.summary.wide.reduced$`avg_<0.5mm_0.75g`, data=AbsoluteValues.summary.wide.reduced)
round(Conversion.0.75g$coef, 2) # regression coefficients
sqrt(mean(Conversion.0.75g$residuals^2)) # RMSE
summary(Conversion.0.75g)$r.squared #R^2

Conversion.2.5g <- lm(AbsoluteValues.summary.wide.reduced$`avg_<2mm_2.5g` ~ AbsoluteValues.summary.wide.reduced$`med_<0.5mm_2.5g`, data=AbsoluteValues.summary.wide.reduced)
round(Conversion.2.5g$coef, 2) # regression coefficients
sqrt(mean(Conversion.2.5g$residuals^2)) # RMSE
summary(Conversion.2.5g)$r.squared #R^2

# for soils > 10% SOC #
Conversion.0.75g.all <- lm(AbsoluteValues.summary.wide$`avg_<2mm_0.75g` ~ AbsoluteValues.summary.wide$`avg_<0.5mm_0.75g`, data=AbsoluteValues.summary.wide)
round(Conversion.0.75g.all$coef, 2) # regression coefficients
sqrt(mean(Conversion.0.75g.all$residuals^2)) # RMSE
summary(Conversion.0.75g.all)$r.squared #R^2

### FIGURES ###
## Figure 1a: absolute differences due to sieve size ##
SoilSummary.sieve <- as.data.frame(POXC.data %>% dplyr::group_by(Soil.ID, Sieve.Size) %>% dplyr::summarise(avg=mean(POXC, na.rm=TRUE)))
SoilSummary.sieve.wide <- spread(SoilSummary.sieve, key="Sieve.Size", value="avg") %>% cbind(Soil.data)
SoilSummary.sieve.wide <- SoilSummary.sieve.wide[,-4] # delete column with duplicate data
SoilSummary.sieve.wide$diff <- SoilSummary.sieve.wide$`<0.5mm`- SoilSummary.sieve.wide$`<2mm`
summary(SoilSummary.sieve.wide$diff)

Sieve.diff.boxplot <- ggplot(SoilSummary.sieve.wide, aes(x= "", y=diff)) + geom_boxplot(width = .6, fill="#EBEBEB") + coord_flip() + theme1 + geom_jitter(shape=16, position=position_jitter(0.05), size=2.5) + labs(y= expression("Difference (<0.5mm \u2013 <2.0mm)")) + stat_summary(fun.y=mean, geom="point", shape=23, size=4, fill="white") + scale_y_continuous(limits=c(-200, 1000))
Sieve.diff.boxplot #export as PNG at 700x200

## Figure 1b: absolute differences due to soil mass ##
SoilSummary.mass <- as.data.frame(na.omit(POXC.data %>% dplyr::group_by(Soil.ID, Mass.class) %>% dplyr::summarise(avg=mean(POXC, na.rm=TRUE))))
SoilSummary.mass.wide <- spread(SoilSummary.mass, key="Mass.class", value="avg") %>% cbind(Soil.data)
SoilSummary.mass.wide$diff <- SoilSummary.mass.wide$`0.75g`- SoilSummary.mass.wide$`2.5g`
SoilSummary.mass.wide <- SoilSummary.mass.wide[-c(8,9,34,35,36), -4] # get rid of duplicate column and rows without both soil masses
summary(SoilSummary.mass.wide$diff)

Mass.diff.boxplot <- ggplot(SoilSummary.mass.wide, aes(x= "", y=diff)) + geom_boxplot(width = .6, fill="#EBEBEB") + coord_flip() + theme1 + geom_jitter(shape=16, position=position_jitter(0.05), size=2.5) + labs(y= expression("Difference (0.75g \u2013 2.50g)")) + stat_summary(fun.y=mean, geom="point", shape=23, size=4, fill="white") + scale_y_continuous(limits=c(-50, 650))
Mass.diff.boxplot #export as PNG at 700x200

# While we're at it, % differences for both... #
SoilSummary.sieve.wide$Pct.chg <- 100*(SoilSummary.sieve.wide$diff/SoilSummary.sieve.wide$`<0.5mm`)
summary(SoilSummary.sieve.wide$Pct.chg)

SoilSummary.mass.wide$Pct.chg <- 100*(SoilSummary.mass.wide$diff/SoilSummary.mass.wide$`2.5g`)
summary(SoilSummary.mass.wide$Pct.chg)


## Figure 2: Detection rates by SOC ##
Detection.data <- as.data.frame(POXC.data %>% dplyr::group_by(Soil.ID, Sieve.Size, Mass.class, .drop=TRUE) %>% dplyr::count(Detect) %>% spread(key=Detect, value=n), na.rm=TRUE)
Detection.data$Pct.detectable <- 100*(Soil.prop$`1`/(Soil.prop$`0` + Soil.prop$`1`))
Detection.data_tmp <- Detection.data %>% dplyr::select(-c(4,5)) %>% pivot_wider(names_from=c(Sieve.Size, Mass.class), values_from=Pct.detectable)
Detection.data.avg <- cbind(Detection.data_tmp, Soil.data)
Detection.data.avg <- Detection.data.avg[, -6]
Detection.data.avg
det.plot_a <- ggplot(data=Detection.data.avg, aes(x=SOC.pct, y=`<2mm_0.75g`)) + geom_point(size=2) + theme1 + ylim(60, 100) + geom_hline(yintercept = 95.7, size = 1, lty=5, color = "#D55E00")
det.plot_b <- ggplot(data=Detection.data.avg, aes(x=SOC.pct, y=`<0.5mm_0.75g`)) + geom_point(size=2) + theme1 + ylim(60, 100) + geom_hline(yintercept = 96.0, size = 1, lty=5, color = "#D55E00")
det.plot_c <- ggplot(data=Detection.data.avg, aes(x=SOC.pct, y=`<2mm_2.5g`)) + geom_point(size=2) + theme1 + ylim(60, 100) + xlim(0, 10) + geom_hline(yintercept = 99.1, size = 1, lty=5, color = "#D55E00")
det.plot_d <- ggplot(data=Detection.data.avg, aes(x=SOC.pct, y=`<0.5mm_2.5g`)) + geom_point(size=2) + theme1 + ylim(60, 100) + xlim(0, 10) + geom_hline(yintercept = 99.5, size = 1, lty=5, color = "#D55E00")
det.plot_a # export image at 600x430
det.plot_b # export image at 600x430
det.plot_c # export image at 600x430
det.plot_d # export image at 600x430


## Figure 3 (top): probability distribution ##
# But first... data wrangling #
CV.lab <- as.data.frame(emmeans(CV.sources, ~ Lab.ID))
CV.lab$mean.bt <- exp(CV.lab$emmean) # backtransform to natural units

summary(CV.lab$emmean)
summary(CV.lab$mean.bt)

lab.boxplot <- ggplot(CV.lab, aes(x= "", y=mean.dt)) + geom_boxplot(width = .6, fill="#EBEBEB") + coord_flip() + theme1 + geom_jitter(shape=16, position=position_jitter(0.05), size=2.5) + labs(y= expression("Coefficient of Variation"~"(%)")) + stat_summary(fun.y=mean, geom="point", shape=23, size=5, fill="white") + scale_y_continuous(limits=c(0,25)) + geom_hline(yintercept = 6.49, size = 1, lty=2, color = "#D55E00") + geom_hline(yintercept = 7.70, size = 1, lty=2, color = "#56B4E9")
lab.boxplot # export as PNG at 450 x 180

## Figure 3 (bottom): boxplot ##
lab.densityplot <- ggplot(CV.lab, aes(x=mean.dt)) + geom_density(fill="#999999", colour=NA, alpha=.15) + geom_line(stat="density") + xlim(0, 25) + theme1 + geom_vline(xintercept = 6.49, size = 1, lty=2, color = "#D55E00") + geom_vline(xintercept = 7.70, size = 1, lty=2, color = "#56B4E9")
lab.densityplot # export 350 x 250


## Figure 4a: soil mass estimation stats ##
CV.mass.data <- as.data.frame(emmeans(CV.sources, ~ Mass.class | Lab.ID | Soil.ID | Sieve.Size))
CV.mass.data <- CV.mass.data[,-c(6:9)]
CV.mass.data$CV.mean.dt <- exp(CV.mass.data$emmean)
CV.mass <- dabest(CV.mass.data, x=Mass.class, y=CV.mean.dt, func=median, idx=c("2.5g", "0.75g"), seed=12345)
plot(CV.mass, rawplot.markersize = .6, palette= "Dark2") # export at 750x453
print(CV.mass) # mean difference = 6.49 (5.55, 9.34) so 2.5g is 6.49% LESS variable on average than 0.75g

## Figure 4b: sieve sive estimation stats ##
CV.sieve.data <- as.data.frame(emmeans(CV.sources, ~ Sieve.Size | Lab.ID | Soil.ID | Mass.class))
CV.sieve.data <- CV.sieve.data[,-c(6:9)]
CV.sieve.data$CV.mean.dt <- exp(CV.sieve.data$emmean)
CV.sieve <- dabest(CV.sieve.data, x=Sieve.Size, y=CV.mean.dt, func=median, idx=c("<2mm", "<0.5mm"), seed=12345)
plot(CV.sieve, rawplot.markersize = 0.6, palette= "Dark2") # export at 750x453
print(CV.sieve) # mean difference = -1.78 (-2.57, -0.906) so <0.5mm is 1.78% LESS variable on average than <2mm


## Figure 5: particle size x soil interaction ##
# Data wrassling #
interaction.data <- as.data.frame(emmeans(CV.sources,  ~  Soil.ID | Sieve.Size))
interaction.data.long <- as.data.frame(interaction.data[,c(1:3)])
interaction.data.wide <- pivot_wider(data=interaction.data.long, names_from=Sieve.Size, values_from=emmean)
interaction.data.wide$diff <- interaction.data.wide$`<2mm` - interaction.data.wide$`<0.5mm`
interaction.data.wide$diff.backtransform <- exp(interaction.data.wide$`<2mm`) - exp(interaction.data.wide$`<0.5mm`)
interaction.data.wide$SOC.pct <- Soil.data$SOC.pct
interaction.data <- gather(key="Sieve.Size", value="POXC", '<0.5mm', '<2mm', data=interaction.data.wide)
interaction.data$mean.bt <- exp(interaction.data$POXC) # backtransform data

interaction.data$Soil.id <- factor(interaction.data$Soil.ID, levels=levels(rev(interaction.data$Soil.ID)))
interaction.data$SOC <- as.factor(interaction.data$SOC.pct)

interaction.dotplot <- ggplot(interaction.data, aes(mean.bt, Soil.ID)) + geom_line(aes(group = Soil.ID), size=1, alpha=0.65) + geom_point(aes(color = Sieve.Size), size=2.5, alpha=0.65) + scale_y_discrete(name="", limits = rev(levels(interaction.data$Soil.ID))) + theme1 + theme(legend.key=element_blank(), legend.title = element_blank()) + geom_vline(xintercept = 8.86, linetype="dotted")
interaction.dotplot #export as PNG at 650 x 700

summary(abs(interaction.data$diff.backtransform))
summary(interaction.data$diff.backtransform)

## Figure 6a: analytical variability of 0.75g at <2mm by SOC ##
# Wrangle data #
Soil.data$SOC2 <- (Soil.data$SOC.pct)^2
Soil.CV.data <- cbind(Soil.data, CV.data.wide)
Soil.CV.data <- Soil.CV.data[,-11]
colnames(Soil.CV.data)[colnames(Soil.CV.data)=="<0.5mm_0.75g"] <- "CV.trans_.5mm.75g"
colnames(Soil.CV.data)[colnames(Soil.CV.data)=="<0.5mm_2.5g"] <- "CV.trans_.5mm.2.5g"
colnames(Soil.CV.data)[colnames(Soil.CV.data)=="<2mm_0.75g"] <- "CV.trans_2mm.75g"
colnames(Soil.CV.data)[colnames(Soil.CV.data)=="<2mm_2.5g"] <- "CV.trans_2mm.2.5g"
Soil.CV.data$CV_.5mm.75g <- exp(Soil.CV.data$CV.trans_.5mm.75g)
Soil.CV.data$CV_.5mm.2.5g <- exp(Soil.CV.data$CV.trans_.5mm.2.5g)
Soil.CV.data$CV_2mm.75g <- exp(Soil.CV.data$CV.trans_2mm.75g)
Soil.CV.data$CV_2mm.2.5g <- exp(Soil.CV.data$CV.trans_2mm.2.5g)
Soil.CV.data <- as.data.frame(Soil.CV.data)
# what kind of relationship is it? Linear/quadratic/exponential?

# Compare model fits... <2mm first (1), then <0.5mm (2) #

SOC.CV.lm1 <- glm(CV_2mm.75g ~ SOC.pct, family=gaussian, data=Soil.CV.data)
SOC.CV.quad1 <- glm(CV_2mm.75g ~ SOC2+ SOC.pct,family=gaussian, data=Soil.CV.data)
SOC.CV.exp1 <- glm(CV_2mm.75g ~ log(SOC.pct),family=gaussian, data=Soil.CV.data)

summary(SOC.CV.lm1)
summary(SOC.CV.quad1)
summary(SOC.CV.exp1)
AIC(SOC.CV.lm, SOC.CV.quad, SOC.CV.exp1, SOC.CV.exp2)
sqrt(mean(residuals(SOC.CV.lm1)^2))
sqrt(mean(residuals(SOC.CV.quad1)^2))
sqrt(mean(residuals(SOC.CV.exp1)^2))
library(lmtest)
lrtest(SOC.CV.lm, SOC.CV.quad, SOC.CV.exp1)
detach("package:lmtest", unload=TRUE)

library(rcompanion)
compareGLM(SOC.CV.lm1, SOC.CV.quad1, SOC.CV.exp1)
detach("package:rcompanion", unload=TRUE)
# exponential is best fit

SOC.CV.lm2 <- glm(CV_.5mm.75g ~ SOC.pct, family=gaussian, data=Soil.CV.data)
SOC.CV.quad2 <- glm(CV_.5mm.75g ~ SOC2+ SOC.pct,family=gaussian, data=Soil.CV.data)
SOC.CV.exp2 <- glm(CV_.5mm.75g ~ log(SOC.pct),family=gaussian, data=Soil.CV.data)

summary(SOC.CV.lm2)
summary(SOC.CV.quad2)
summary(SOC.CV.exp2)
AIC(SOC.CV.lm2, SOC.CV.quad2, SOC.CV.exp2)
sqrt(mean(residuals(SOC.CV.lm2)^2))
sqrt(mean(residuals(SOC.CV.quad2)^2))
sqrt(mean(residuals(SOC.CV.exp2)^2))
library(lmtest)
lrtest(SOC.CV.lm2, SOC.CV.quad2, SOC.CV.exp2)
detach("package:lmtest", unload=TRUE)

library(rcompanion)
compareGLM(SOC.CV.lm2, SOC.CV.quad2, SOC.CV.exp2)
detach("package:rcompanion", unload=TRUE)
# apparently, exponential is best for both!

# Get model coefficient and R2
summary(SOC.CV.exp1)
summary(SOC.CV.exp2)

# make predicted line for graphs 
xmin <- min(Soil.CV.data$SOC.pct)
xmax <- max(Soil.CV.data$SOC.pct)
predicted <- data.frame(SOC.pct=seq(xmin, xmax, length.out=100))
predicted$predvals1 <- predict(SOC.CV.exp1, predicted)
predicted$predvals2 <- predict(SOC.CV.exp2, predicted)
predicted <- as.data.frame(predicted)

SOC.CV.plot_a <- ggplot(data=Soil.CV.data, aes(x=SOC.pct, y=CV_2mm.75g)) + geom_point(size=3) + geom_line(data=predicted, aes(x=SOC.pct, y=predvals1), size=1) + theme1
SOC.CV.plot_a # export as PNG at 700x431

## Figure 6b: analytical variability of 0.75g at <0.5mm by SOC ##
SOC.CV.plot_b <- ggplot(data=Soil.CV.data, aes(x=SOC.pct, y=CV_.5mm.75g)) + geom_point(size=3) + geom_line(data=predicted, aes(x=SOC.pct, y=predvals2), size=1) + theme1
SOC.CV.plot_b # export as PNG at 700x431

### SUPPLEMENTARY DATA ###
## Figure S1 and Table S2: analytical variability by soil order ##
CV.by.order <- emmeans(CV.sources1.2,  ~  Soil.Order) # calculate using Soil.Order "best" model
CLD(CV.by.order, details = FALSE, sort = TRUE,, alpha = 0.05, reversed=TRUE, Letters=letters)

order.boxplot <- ggplot(CV.order, aes(x=fct_reorder(Soil.Order, emmean, .desc=TRUE), y=emmean, fill=Soil.Order)) + geom_boxplot() + theme1 + scale_fill_brewer(palette="Paired", guide=FALSE) + geom_jitter(size = 1.5, shape = 16, alpha=0.5, width=0.15, stroke=0)
order.boxplot

## Table S3: Sample size, given +/- 20% of the overall mean for a variety of confidence intervals ##
# equation based on (http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Power/BS704_Power_print.html) for one sample with continuous outcome #
CV.func <- function(x){sd(x, na.rm=TRUE)*100/(mean(x, na.rm=TRUE))}
# N.XX = number of samples for a given confidence interval XX
N.99 <- function(x){((6.635776)*((CV.func(x))^2))/(400)} # 99% CI
N.95 <- function(x){((3.8416)*((CV.func(x))^2))/(400)} # 95% CI
N.90 <- function(x){((2.706025)*((CV.func(x))^2))/(400)} # 90% CI
N.80 <- function(x){((1.643524)*((CV.func(x))^2))/(400)} # 80% CI

N.99(CV.data.low$CV_.5mm.75g)
N.95(CV.data.low$CV_.5mm.75g)
N.90(CV.data.low$CV_.5mm.75g)
N.80(CV.data.low$CV_.5mm.75g)
N.99(CV.data.med$CV_.5mm.75g)
N.95(CV.data.med$CV_.5mm.75g)
N.90(CV.data.med$CV_.5mm.75g)
N.80(CV.data.med$CV_.5mm.75g)
N.99(CV.data.high$CV_.5mm.75g)
N.95(CV.data.high$CV_.5mm.75g)
N.90(CV.data.high$CV_.5mm.75g)
N.80(CV.data.high$CV_.5mm.75g)


### IN-TEXT CALCULATIONS AKA "data not showns" (in order of appearance) ###
# Range of CV values #
round(summary(POXC.CV$CV.pct), 2)

# 95% confidence interval for analytical variability #
summary(POXC.CV$log.CV.pct)
CI(na.omit(POXC.CV$log.CV.pct, ci=0.95)) # backtransform these values

# median absolute difference #
summary(AbsoluteValues.summary$abs.dev)
CI(AbsoluteValues.summary$abs.dev, ci=0.95)

# Reader method comparison #
Method.df <- as.data.frame(c("cuvette", "cuvette", "plate reader", "cuvette", "plate reader", "cuvette", "plate reader", "plate reader", "plate reader", "cuvette", "plate reader", "cuvette"))
POXC.method.data <- cbind(CV.lab, Method.df)
names(POXC.method.data)[8] <- "Reader.Method"
Method.CV <- lm(mean.dt ~ Reader.Method, data=POXC.method.data)
Anova(Method.CV,type=3)

# Magnitude of difference between sieve sizes and SOC #
int.lm <- lm(abs(diff.backtransform) ~ SOC.pct, data= interaction.data)
summary(int.lm)

## CODE SNIPPETS ##
# Theme for graphs cuz... aesthetics #
theme1 <- theme(plot.background=element_blank(), panel.grid.major=element_line(color="#999999", size=0.15), panel.grid.minor=element_line(color="#CCCCCC", size=0.15), panel.border=element_rect(color="#000000", size=0.5, fill=NA), panel.background=element_blank())
