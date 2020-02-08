################################################################################
##                                                                            ##
##                       SleepCare Trial Power Analysis                       ##
##                                                                            ##
################################################################################

## set location for SLeepWell trial dataset used for analysis to determine
## parameter estimates for power simulations
if (Sys.info()[["user"]] == "jwile") {
  loc.data <- "g:/My Drive/CBTI+Light project/Recruitment and data collection/"
} else {
  loc.data <- "g:/My Drive/CBTI+Light project/Recruitment and data collection/"
}

## load various relevant packages

## graphics packages
library(scales)
library(RColorBrewer)
library(grid)
library(ggplot2)
library(cowplot)
library(visreg)
library(ggthemes)
library(viridis)

## tools and miscellaneous packages
library(devtools)
library(testthat)
library(mvtnorm)
library(car)
library(extraoperators)
library(JWileymisc)

## data related packages
library(reshape2)
library(xtable)
library(chron)
library(zoo)
library(data.table)

## modelling packages
library(MASS)
library(MplusAutomation)

theme_set(theme_cowplot())
Sys.setenv(TZ = "Australia/Melbourne")

## load SleepWell Trial data
d.all <- readRDS(file.path(loc.data, "Data", "SleepWell_Surveys.RDS"))

SEback <- function(x) 100 - ((100 - x)^2)

## define function to run monte carlo simulations for power analysis
SleepCare_MonteCarlo <- function(swap, modelout = "power_sim.inp",
                                 TITLE = "SleepCare Monte Carlo;") {
MONTECARLO <-  "
NAMES ARE Y1 Y2 Y3 dstrat2 dstrat3 dstrat4
 dstrat5 dstrat6 dstrat7 dstrat8
 blt cbt bxc;
NOBSERVATIONS = %nobs%;
NREPS = %nreps%;
SEED = 58459;
PATMISS = Y1(%y1miss%) Y2(%y2miss%) Y3(%y3miss%);
PATPROBS = 1;
CUTPOINTS = dstrat2 (%strat2cut%) dstrat3 (%strat3cut%) dstrat4 (%strat4cut%)
 dstrat5 (%strat5cut%) dstrat6 (%strat6cut%) dstrat7 (%strat7cut%)
 dstrat8 (%strat8cut%)
 blt (0) cbt(0) bxc (1.098612);
RESULTS = sleepcare_sim_estimate.dat;
"
MODELPOPULATION = "
[Y1@%y1m% Y2@%y2m% Y3@%y3m%];
Y1 - Y3@%resvar%;
dstrat2-bxc@1;
int BY Y1@1 Y2@1 Y3@1;
slope BY Y1@0 Y2@0.5 Y3@1;
int ON
  dstrat2@%istrat2% dstrat3@%istrat3% dstrat4@%istrat4%
  dstrat5@%istrat5% dstrat6@%istrat6% dstrat7@%istrat7%
  dstrat8@%istrat8%
  blt@%iblt% cbt@%icbt% bxc@%ibxc% ;
slope ON blt@%sblt% cbt@%scbt% bxc@%sbxc% ;
[int@%intm% slope@%slopem%];
int@%intv% slope@%slopev%;
int WITH slope@%iscov%;
"
MODEL <- "
[Y1 - Y3@0]; ! fix residual intercepts to 0
Y1-Y3*%resvar% (resvar); ! homogenous residual variances
int BY Y1@1 Y2@1 Y3@1;
slope BY Y1@0 Y2@0.5 Y3@1;
int ON
  dstrat2@%istrat2% dstrat3@%istrat3% dstrat4@%istrat4%
  dstrat5@%istrat5% dstrat6@%istrat6% dstrat7@%istrat7%
  dstrat8@%istrat8%
  blt@0 cbt@0 bxc@0 ;
slope ON blt*%sblt% cbt*%scbt% bxc*%sbxc% (s1-s3);
[int*%intm% slope*%slopem%] (int1-int2);
int*%intv% slope*%slopev%;
int WITH slope*%iscov%;
"

for (v in names(swap)) {
  MONTECARLO <- gsub(v, swap[[v]], MONTECARLO)
  MODELPOPULATION <- gsub(v, swap[[v]], MODELPOPULATION)
  MODEL <- gsub(v, swap[[v]], MODEL)
}

  msimfit <- mplusModeler(
    mplusObject(
  TITLE = TITLE,
  ANALYSIS = "ESTIMATOR = ML; PROCESSORS = 4;",
  MONTECARLO = MONTECARLO,
  MODELPOPULATION = MODELPOPULATION,
  MODEL = MODEL,
  MODELCONSTRAINT = "
NEW(ss1 ss2 ss3 ss4 mblt mcbt);
ss1 = int2;
ss2 = ss1 + s1;
ss3 = ss1 + s2;
ss4 = ss1 + s3;
mblt = s1 + (s3/2);
mcbt = s2 + (s3/2);
"), modelout = modelout, run=TRUE)

est <- matrix(as.numeric(
  scan("sleepcare_sim_estimate.dat", what = character()
       ) %s!in% as.character(seq(from = 1, to = swap[["%nreps%"]]))),
  nrow = 100, byrow=TRUE)

  return(list(Model = msimfit, Estimates = est))
}


################################################################################
##                                                                            ##
##              Power Analysis for Insomnia Severity Index (ISI)              ##
##                                                                            ##
################################################################################

isidat <- reshape(d.all[, .(ID, condition, baseline_isi_high, baseline_stage_high, ISI_Total, Survey)],
        v.names = "ISI_Total",
        timevar = "Survey",
        idvar = c("ID", "condition", "baseline_isi_high", "baseline_stage_high"),
        direction = "wide")[!is.na(ISI_Total.Baseline)]
setnames(isidat, names(isidat), c("id", "condition", "isihigh", "stagehigh",
                                  "ISI0", "ISI1", "ISI2", "ISI3", "ISI4"))
isidat[, dcond := as.integer(condition == "CBT+")]
isidat[, dstrat1 := as.integer(isihigh == "ISI <= 7" & stagehigh == "Stage <= 2")]
isidat[, dstrat2 := as.integer(isihigh == "ISI <= 7" & stagehigh == "Stage >= 3")]
isidat[, dstrat3 := as.integer(isihigh == "ISI >= 8" & stagehigh == "Stage <= 2")]
isidat[, dstrat4 := as.integer(isihigh == "ISI >= 8" & stagehigh == "Stage >= 3")]

m <- mplusObject(
  TITLE = "Estimates for Power Analysis",
  ANALYSIS = "ESTIMATOR = ML;!ESTIMATOR = BAYES; ALGORITHM=GIBBS(RW); BITER = 100000 (20000);",
  MODEL = "
  [ISI1 - ISI3@0]; ! fix residual intercepts to 0
  ISI1* ISI2* ISI3* (resvar); ! homogenous residual variances
  int BY ISI1@1 ISI2@1 ISI3@1;
  slope BY ISI1@0 ISI2@0.5 ISI3@1;
  int ON dstrat2 dstrat3 dstrat4 dcond@0 (i1 - i4);
  slope ON dcond (s4);
  [int* slope*] (rm1 - rm2);
  int* slope* (rv1-rv2);
  int WITH slope*0;
",
OUTPUT = "STDYX;",
usevariables = c("ISI1", "ISI2", "ISI3", "dstrat2", "dstrat3", "dstrat4", "dcond"),
rdata = as.data.frame(isidat))
mfit <- mplusModeler(m, dataout = "power_isi.dat", run = TRUE)
coef(mfit)

powerGRIDISI <- as.data.table(expand.grid(
  NUM = seq(from = 100, to = 300, by = 10),
  resvar = 8.066 * c(1, 1.1),
  intvar = 11.105 * c(1, 1.1),
  blt = -5 - -2,
  cbt = c(-6 - -2, -5 - -2),
  bxc = c(-3, 0, +3)))

table(rowSums(powerGRIDISI[,.(blt, cbt, bxc)]))
table(rowSums(powerGRIDISI[,.(blt, cbt, bxc)])-2)

#### STRATA info
## 1: PM - 1to3 - low ISI
## 2: PM - 1to3 - High ISI
## 3: PM -   4  - low ISI
## 4: PM -   4  - High ISI
## 5: MH - 1to3 - low ISI
## 6: MH - 1to3 - High ISI
## 7: MH -   4  - low ISI
## 8: MH -   4  - High ISI

## assume 60 v 40 split on PM vs MH
## assume a 75 v 25 split on high v low ISI
## assume a 75 v 25 split on 1to3 v 4 (cancer stage)
## use this to calculate probability in each strata

if (FALSE) {
isipower <- lapply(1:nrow(powerGRIDISI), function(i) {
  SleepCare_MonteCarlo(swap = c(
  "%nobs" = powerGRIDISI$NUM[i], "%nreps%" = 500,
  "%y1miss%" = 0, "%y2miss%" = .25, "%y3miss%" = .35,
  "%strat2cut%" = round(qlogis(.6 * .75 * .25, lower.tail=FALSE), 3),
  "%strat3cut%" = round(qlogis(.6 * .25 * .75, lower.tail=FALSE), 3),
  "%strat4cut%" = round(qlogis(.6 * .25 * .25, lower.tail=FALSE), 3),
  "%strat5cut%" = round(qlogis(.4 * .75 * .75, lower.tail=FALSE), 3),
  "%strat6cut%" = round(qlogis(.4 * .75 * .25, lower.tail=FALSE), 3),
  "%strat7cut%" = round(qlogis(.4 * .25 * .75, lower.tail=FALSE), 3),
  "%strat8cut%" = round(qlogis(.4 * .25 * .25, lower.tail=FALSE), 3),
  "%y1m%" = 0, "%y2m%" = 0, "%y3m%" = 0,
  "%resvar%" = powerGRIDISI$resvar[i],
  "%istrat2%" = 5,
  "%istrat3%" = 1,
  "%istrat4%" = 6,
  "%istrat5%" = 0,
  "%istrat6%" = 5,
  "%istrat7%" = 1,
  "%istrat8%" = 6,
  "%iblt%" = 0, "%icbt%" = 0, "%ibxc%" = 0,
  "%sblt%" = powerGRIDISI$blt[i],
  "%scbt%" = powerGRIDISI$cbt[i],
  "%sbxc%" = powerGRIDISI$bxc[i],
  "%intm%" = 8.469, "%slopem%" = -2,
  "%intv%" = powerGRIDISI$intvar[i], "%slopev%" = 6.765, "%iscov%" = 2.314))
  })
saveRDS(isipower, file = "isipower.RDS", compress = "xz")
} else {
  isipower <- readRDS("isipower.RDS")
}

powerGRIDISI[, ResidualVariance := factor(resvar, levels = sort(unique(resvar)),
                                          labels = c("Matching SleepWell", "SleepWell + 10%"))]
powerGRIDISI[, InterceptVariance := factor(intvar, levels = sort(unique(intvar)),
                                           labels = c(
                                             "Intercept Variance\nMatching SleepWell",
                                             "Intercept Variance\nSleepWell + 10%"))]
powerGRIDISI[, CBT := factor(cbt, levels = c(-4, -3),
                             labels = c(
                               "CBT Simple Effect\n(High, -4 ISI)",
                               "CBT Simple Effect\n(MID, -3 ISI)"))]
powerGRIDISI[, InteractionEffect := factor(bxc, levels = c(-3, 0, 3),
                                     labels = c("BLT x CBT (Neg, -3 ISI)",
                                                "BLT x CBT (None, +0 ISI)",
                                                "BLT x CBT (Pos, +3 ISI)"))]

powerGRIDISI$powerSBLT <- unlist(lapply(isipower, function(m) {
  m$Model$results$parameters$unstandardized[17, "pct_sig_coef"]
}))
powerGRIDISI$powerSCBT <- unlist(lapply(isipower, function(m) {
  m$Model$results$parameters$unstandardized[18, "pct_sig_coef"]
}))
powerGRIDISI$powerSBXC <- unlist(lapply(isipower, function(m) {
  m$Model$results$parameters$unstandardized[19, "pct_sig_coef"]
}))
powerGRIDISI$powerMBLT <- unlist(lapply(isipower, function(m) {
  m$Model$results$parameters$unstandardized[35, "pct_sig_coef"]
}))
powerGRIDISI$powerMCBT <- unlist(lapply(isipower, function(m) {
  m$Model$results$parameters$unstandardized[36, "pct_sig_coef"]
}))


table(powerGRIDISI[, .(MBLT = blt + bxc/2)])
table(powerGRIDISI[, .(MBLT = blt + bxc/2 - 2)])

table(rowSums(powerGRIDISI[, c("blt", "cbt", "bxc")])-2)
table(rowSums(powerGRIDISI[, c("blt"), drop=FALSE])-2)


powerISIPlot <- function(y = "", Title = "", data = powerGRIDISI) {
  ggplot(data, aes_string(x = "NUM", y = y, colour = "ResidualVariance",
                          linetype = "InteractionEffect")) +
    stat_smooth(se=FALSE, method = "loess") +
    scale_x_continuous(
      "Total Participants Randomised",
      breaks = seq(from = 100, to = 300, by = 20)) +
    scale_y_continuous("Power (Monte Carlo Simulation)",
                       breaks = seq(0, 1, by = .1),
                       labels = paste0(seq(0, 100, by = 10), "%")) +
    scale_colour_manual(values = c(
                          "Matching SleepWell" = "black",
                          "SleepWell + 10%" = "grey60")) +
    scale_linetype_manual(values = c(
                            "BLT x CBT (Neg, -3 ISI)" = 1,
                            "BLT x CBT (None, +0 ISI)" = 2,
                            "BLT x CBT (Pos, +3 ISI)" = 3)) +
    geom_hline(yintercept = .8) +
    geom_vline(xintercept = c(210)) +
    facet_grid(CBT ~ InterceptVariance, scale = "free_y") +
    ggtitle(label = Title,
            subtitle = paste(
              "MID = Minimally Important Difference (Change of 3 in ISI).",
              "BLT Simple Effect held at MID, -3 ISI.",
              sep = "\n")) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.key.width = unit(1.5, "cm"))
}

powerISIPlot("powerMBLT", Title = "Power for the Main Effect of BLT on ISI")
powerISIPlot("powerMCBT", Title = "Power for the Main Effect of CBT on ISI")
powerISIPlot("powerSBXC", Title = "Power for the BLT x CBT Interaction on ISI")


ggplot(powerGRIDISI, aes(NUM, powerMBLT, linetype = ResidualVariance,
                         colour = InteractionEffect)) +
  stat_smooth(se=FALSE, method = "loess") +
  geom_hline(yintercept = .8) +
  geom_vline(xintercept = c(220)) +
  facet_grid(CBT ~ InterceptVariance)

ggplot(powerGRIDISI, aes(NUM, powerMCBT, linetype = ResidualVariance,
                         colour = Interaction)) +
  stat_smooth(se=FALSE, method = "loess") +
  geom_hline(yintercept = .8) +
  geom_vline(xintercept = c(220)) +
  facet_grid(CBT ~ InterceptVariance)

ggplot(powerGRIDISI, aes(NUM, powerSBXC, linetype = ResidualVariance,
                         colour = Interaction)) +
  stat_smooth(se=FALSE, method = "loess") +
  geom_hline(yintercept = .8) +
  geom_vline(xintercept = c(220)) +
  facet_grid(CBT ~ InterceptVariance)

################################################################################
##                                                                            ##
##                Power Analysis for Fatigue Symptoms (PROMIS)                ##
##                                                                            ##
################################################################################

fadat <- reshape(d.all[, .(ID, condition, baseline_isi_high, baseline_stage_high, FA_TScore, Survey)],
        v.names = "FA_TScore",
        timevar = "Survey",
        idvar = c("ID", "condition", "baseline_isi_high", "baseline_stage_high"),
        direction = "wide")[!is.na(FA_TScore.Baseline)]

setnames(fadat, names(fadat), c("id", "condition", "isihigh", "stagehigh",
                                  "FA0", "FA1", "FA2", "FA3", "FA4"))
fadat[, dcond := as.integer(condition == "CBT+")]
fadat[, dstrat1 := as.integer(isihigh == "ISI <= 7" & stagehigh == "Stage <= 2")]
fadat[, dstrat2 := as.integer(isihigh == "ISI <= 7" & stagehigh == "Stage >= 3")]
fadat[, dstrat3 := as.integer(isihigh == "ISI >= 8" & stagehigh == "Stage <= 2")]
fadat[, dstrat4 := as.integer(isihigh == "ISI >= 8" & stagehigh == "Stage >= 3")]


mfa <- mplusObject(
  TITLE = "Estimates for Power Analysis",
  ANALYSIS = "ESTIMATOR = ML;!ESTIMATOR = BAYES; ALGORITHM=GIBBS(RW); BITER = 100000 (20000);",
  MODEL = "
  [FA1 - FA3@0]; ! fix residual intercepts to 0
  FA1*30 FA2*30 FA3*30 ;!(resvar); ! homogenous residual variances
  int BY FA1@1 FA2@1 FA3@1;
  slope BY FA1@0 FA2*0.5 FA3@1;
  int ON dstrat2 dstrat3 dstrat4 dcond@0 (i1 - i4);
  slope ON dcond (s4);
  [int* slope*] (rm1 - rm2);
  int*30 slope*19 (rv1-rv2);
  int WITH slope*0;
",
OUTPUT = "STDYX;",
usevariables = c("FA1", "FA2", "FA3", "dstrat2", "dstrat3", "dstrat4", "dcond"),
rdata = as.data.frame(fadat))
mfafit <- mplusModeler(mfa, dataout = "power_fa.dat", run = TRUE)
coef(mfafit)


powerGRIDFA <- as.data.table(expand.grid(
  NUM = seq(from = 100, to = 300, by = 10),
  resvar = round(mean(c(9.2, 34.8, 21.9)) * c(1, 1.1), 1),
  intvar = round(28.7 * c(1, 1.1), 1),
  blt = c(-5.2 - -.2, -4.2 - -.2),
  cbt = c(-4.2 - -.2),
  bxc = c(-4, 0, +4)))

table(rowSums(powerGRIDFA[,.(blt, cbt, bxc)]))
table(rowSums(powerGRIDFA[,.(blt, cbt, bxc)])-.2)

#### STRATA info
## 1: PM - 1to3 - low ISI
## 2: PM - 1to3 - High ISI
## 3: PM -   4  - low ISI
## 4: PM -   4  - High ISI
## 5: MH - 1to3 - low ISI
## 6: MH - 1to3 - High ISI
## 7: MH -   4  - low ISI
## 8: MH -   4  - High ISI

## assume 60 v 40 split on PM vs MH
## assume a 75 v 25 split on high v low ISI
## assume a 75 v 25 split on 1to3 v 4 (cancer stage)
## use this to calculate probability in each strata

if (FALSE) {
  fapower <- lapply(1:nrow(powerGRIDFA), function(i) {
    SleepCare_MonteCarlo(swap = c(
                           "%nobs" = powerGRIDFA$NUM[i], "%nreps%" = 500,
                           "%y1miss%" = 0, "%y2miss%" = .25, "%y3miss%" = .35,
                           "%strat2cut%" = round(qlogis(.6 * .75 * .25, lower.tail=FALSE), 3),
                           "%strat3cut%" = round(qlogis(.6 * .25 * .75, lower.tail=FALSE), 3),
                           "%strat4cut%" = round(qlogis(.6 * .25 * .25, lower.tail=FALSE), 3),
                           "%strat5cut%" = round(qlogis(.4 * .75 * .75, lower.tail=FALSE), 3),
                           "%strat6cut%" = round(qlogis(.4 * .75 * .25, lower.tail=FALSE), 3),
                           "%strat7cut%" = round(qlogis(.4 * .25 * .75, lower.tail=FALSE), 3),
                           "%strat8cut%" = round(qlogis(.4 * .25 * .25, lower.tail=FALSE), 3),
                           "%y1m%" = 0, "%y2m%" = 0, "%y3m%" = 0,
                           "%resvar%" = powerGRIDFA$resvar[i],
                           "%istrat2%" = 2,
                           "%istrat3%" = 4,
                           "%istrat4%" = 6,
                           "%istrat5%" = 0,
                           "%istrat6%" = 2,
                           "%istrat7%" = 4,
                           "%istrat8%" = 6,
                           "%iblt%" = 0, "%icbt%" = 0, "%ibxc%" = 0,
                           "%sblt%" = powerGRIDFA$blt[i],
                           "%scbt%" = powerGRIDFA$cbt[i],
                           "%sbxc%" = powerGRIDFA$bxc[i],
                           "%intm%" = 52.2, "%slopem%" = -.2,
                           "%intv%" = powerGRIDFA$intvar[i], "%slopev%" = 24.3, "%iscov%" = -6.5),
                         modelout = "power_fa_sim.inp", TITLE = "SleepCare FA Power;")
  })
  saveRDS(fapower, file = "fapower.RDS", compress = "xz")
} else {
  fapower <- readRDS("fapower.RDS")
}

powerGRIDFA[, ResidualVariance := factor(resvar, levels = sort(unique(resvar)),
                                          labels = c("Matching SleepWell", "SleepWell + 10%"))]
powerGRIDFA[, InterceptVariance := factor(intvar, levels = sort(unique(intvar)),
                                          labels = c("Intercept Variance\nMatching SleepWell",
                                                     "Intercept Variance\nSleepWell + 10%"))]
powerGRIDFA[, BLT := factor(blt, levels = c(-5, -4),
                            labels = c(
                              "BLT Simple Effect\n(High, -5 Fatigue)",
                              "BLT Simple Effect\n(MID, -4 Fatigue)"))]
powerGRIDFA[, InteractionEffect := factor(bxc, levels = c(-4, 0, 4),
                                     labels = c("BLT x CBT (Neg, -4 Fatigue)",
                                                "BLT x CBT (None, +0 Fatigue)",
                                                "BLT x CBT (Pos, +4 Fatigue)"))]

powerGRIDFA$powerSBLT <- unlist(lapply(fapower, function(m) {
  m$Model$results$parameters$unstandardized[17, "pct_sig_coef"]
}))
powerGRIDFA$powerSCBT <- unlist(lapply(fapower, function(m) {
  m$Model$results$parameters$unstandardized[18, "pct_sig_coef"]
}))
powerGRIDFA$powerSBXC <- unlist(lapply(fapower, function(m) {
  m$Model$results$parameters$unstandardized[19, "pct_sig_coef"]
}))
powerGRIDFA$powerMBLT <- unlist(lapply(fapower, function(m) {
  m$Model$results$parameters$unstandardized[35, "pct_sig_coef"]
}))
powerGRIDFA$powerMCBT <- unlist(lapply(fapower, function(m) {
  m$Model$results$parameters$unstandardized[36, "pct_sig_coef"]
}))

powerFAPlot <- function(y = "", Title = "", data = powerGRIDFA) {
  ggplot(data, aes_string(x = "NUM", y = y, colour = "ResidualVariance",
                          linetype = "InteractionEffect")) +
    stat_smooth(se=FALSE, method = "loess") +
    scale_x_continuous(
      "Total Participants Randomised",
      breaks = seq(from = 100, to = 300, by = 20)) +
    scale_y_continuous("Power (Monte Carlo Simulation)",
                       breaks = seq(0, 1, by = .1),
                       labels = paste0(seq(0, 100, by = 10), "%")) +
    scale_colour_manual(values = c(
                          "Matching SleepWell" = "black",
                          "SleepWell + 10%" = "grey60")) +
    scale_linetype_manual(values = c(
                            "BLT x CBT (Neg, -4 Fatigue)" = 1,
                            "BLT x CBT (None, +0 Fatigue)" = 2,
                            "BLT x CBT (Pos, +4 Fatigue)" = 3)) +
    geom_hline(yintercept = .8) +
    geom_vline(xintercept = c(210)) +
    facet_grid(BLT ~ InterceptVariance, scale = "free_y") +
    ggtitle(label = Title,
            subtitle = paste(
              "MID = Minimally Important Difference (Change of 4 in Fatigue).",
              "CBT Simple Effect held at MID, -4 Fatigue.",
              sep = "\n")) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.key.width = unit(1.5, "cm"))
}


################################################################################
##                                                                            ##
##               Save Plots of Power Analysis for ISI and Fatigue             ##
##                                                                            ##
################################################################################

pdf("SleepCare_power_plots.pdf", width = 9, height = 12)
powerISIPlot("powerMBLT", Title = "Power for the Main Effect of BLT on ISI")
powerISIPlot("powerMCBT", Title = "Power for the Main Effect of CBT on ISI")
powerISIPlot("powerSBXC", Title = "Power for the BLT x CBT Interaction on ISI")
powerFAPlot("powerMBLT", Title = "Power for the Main Effect of BLT on Fatigue")
powerFAPlot("powerMCBT", Title = "Power for the Main Effect of CBT on Fatigue")
powerFAPlot("powerSBXC", Title = "Power for the BLT x CBT Interaction on Fatigue")
dev.off()

png("SleepCare_MainEffect_BLT_ISI.png", width = 9, height = 12, units = "in", res = 300)
powerISIPlot("powerMBLT", Title = "Power for the Main Effect of BLT on ISI")
dev.off()

png("SleepCare_MainEffect_CBT_ISI.png", width = 9, height = 12, units = "in", res = 300)
powerISIPlot("powerMCBT", Title = "Power for the Main Effect of CBT on ISI")
dev.off()

png("SleepCare_Interaction_ISI.png", width = 9, height = 12, units = "in", res = 300)
powerISIPlot("powerSBXC", Title = "Power for the BLT x CBT Interaction on ISI")
dev.off()

png("SleepCare_MainEffect_BLT_Fatigue.png", width = 9, height = 12, units = "in", res = 300)
powerFAPlot("powerMBLT", Title = "Power for the Main Effect of BLT on Fatigue")
dev.off()

png("SleepCare_MainEffect_CBT_Fatigue.png", width = 9, height = 12, units = "in", res = 300)
powerFAPlot("powerMCBT", Title = "Power for the Main Effect of CBT on Fatigue")
dev.off()

png("SleepCare_Interaction_Fatigue.png", width = 9, height = 12, units = "in", res = 300)
powerFAPlot("powerSBXC", Title = "Power for the BLT x CBT Interaction on Fatigue")
dev.off()

