##### Load packages #####
library(readxl)
library(tidyverse)
library(zoo)
library(kableExtra)
library(reshape2)
library(ggpubr)
library(caret)
library(rms)

library(glmnet)
library(boot)
library(cowplot)

library(survival)
library(survminer)
library(sva)

library(gtsummary)

#### PCA table function #####
#### Input: data frame for table, column index, PCA summaries for both hemispheres
#### Output: ggplot object
pca.table <- function(table.tmp, ind, a) {
  tmp.a <- summary(a)
  table.tmp[1:10,ind] <- c(tmp.a$importance['Proportion of Variance',c(1:10)])
  return(table.tmp)
}

#### Function to calculate lambda for LASSO regression #####
#### Input: summary from cross-validation procedure, tolerance value
#### Output: optimized lambda

get_lambda <- function(fit, tol = .5) {
  error <- fit$cvm[fit$lambda == fit$lambda.min]
  sd <- fit$cvsd[fit$lambda == fit$lambda.min]
  tolerance <- error + tol*sd
  max(fit$lambda[fit$cvm < tolerance])
}

#### Function to compute LASSO coefficients #####
#### Input: input data, subject subset (determined by bootstrapping function)
#### Output: coefficients

lasso.coef <- function(data, indices) {
  d <- data[indices,]
  fit <- glmnet(x = d[,c(3:ncol(d))], y = d[,c(1,2)], 
                    alpha = 1, family = "cox", lambda = get_lambda(lasso.fit.cv, 0.2), standardize = TRUE)
  return(coef(fit)[,1])
}

#### PCA variance plot function #####
#### Input: PCA results, title
#### Output: ggplot object

pca_plot <- function(plot.data, title){
  
  tmp.plot <- ggplot(plot.data, aes(x = x, y = value, group = variable)) + 
    geom_line(aes(color = variable), size = 1.1) +
    geom_point(aes(color = variable), size = 3) +
    scale_x_continuous(name ="First ten principal components", breaks = c(1:10),
                limits=c(1,10)) +
    scale_y_continuous(name="% of variance explained", breaks = c(0,25,50,75), limits=c(0, 75)) +
    theme_bw() +
    scale_color_brewer(palette="Set1") +
    ggtitle(title) + 
    labs(color = 'Cortical measure') +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 16),
         axis.title.y = element_text(size = 16),
         legend.text = element_text(size = 14),
         legend.title = element_text(size = 16),
         axis.text.x = element_text(size = 14),
         axis.text.y = element_text(size = 14)) +
    theme(legend.position="bottom")
  return(tmp.plot)
}

#### PCA loading plot function #####
#### Input: loading values for domain in each hemisphere, title
#### Output: ggplot object

loading_plot <- function(domain.data.lh, title){

  load_lh <- domain.data.lh$rotation[,1]
  lh_best <- as.data.frame(abs(load_lh[names(sort(abs(load_lh), decreasing = TRUE)[1:10])]))
  lh_best$area <- rownames(lh_best) %>%   strsplit( '^[^_]*?[_][^_]*?(*SKIP)(*F)|_' , perl=TRUE) %>%   sapply( "[", 1 )
  colnames(lh_best) <- c('Loading', 'Area')
  lh_best$Area <- factor(lh_best$Area, levels = lh_best$Area)
  
  left_loading <- ggplot(lh_best, aes(x=Loading, y=Area))+
    geom_bar(stat="identity", width=0.7, fill="steelblue")+
    theme_minimal() +
    scale_y_discrete(name ="Cortical segments", limits=rev, labels = ) +
    scale_x_continuous(name="Factor loading onto PC1", breaks = c(0, 0.1, 0.2, 0.3), limits=c(0,0.3)) +
    ggtitle(title) +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
       axis.title.x = element_text(size = 16),
       axis.title.y = element_text(size = 16),
       legend.text = element_text(size = 14),
       legend.title = element_text(size = 16),
       axis.text.x = element_text(size = 14),
       axis.text.y = element_text(size = 14)) 
  
  return(left_loading)
}

#### Cox model function #####
#### Input: predictors for Cox model, row index for table, tables, additional predictors (optional)
#### Output: summary and overfitting table

cox_model <- function(cox.bl, i, cox.table, overfit.table, ...){

  x = list(...)
  # Get bootstrapped model
  bl.conf <- bootcov(cox.bl, B=1000, pr=TRUE)
  # Extract regression parameters for each predictor
  bl.an <- anova(bl.conf)
  # Store regression summary
  bl.summ <- summary(bl.conf, pre.cond = 3)
  # Internal validation using bootstrapping
  bl.val <- validate(cox.bl, method="boot", B=1000, bw=FALSE)
  # Store overfitting values
  overfit.table[i,] <- c(bl.val['Dxy', 'optimism'], bl.val['R2','optimism'], bl.val['Slope','optimism'])
  
  if (length(x) == 0){
    # Store regression results
    cox.table[i,] <- c(bl.summ[2,4], bl.an[1,3], bl.summ[4,4], bl.an[2,3], bl.summ[6,4], bl.an[3,3], bl.summ[8,4], bl.an[3,3], 0, 0, 0, 0, LR.orig, 0, 1, bl.val['Dxy',   'index.corrected']/2 + 0.5)
    return(list(cox.table, overfit.table))

  }else if (nrow(bl.summ) == 10){
        bl.p <- lrtest(x[[1]], bl.conf)
  
    cox.table[i,] <- c(bl.summ[2,4], bl.an[1,3], bl.summ[4,4], bl.an[2,3], bl.summ[8,4], bl.an[3,3],bl.summ[10,4], bl.an[3,3], bl.summ[6,4], bl.an[4,3],0,0, cox.bl$stats['Model L.R.'], bl.p$stats['P'], 1-(x[[2]]/cox.bl$stats['Model L.R.']), bl.val['Dxy', 'index.corrected']/2 + 0.5)
        return(list(cox.table, overfit.table))
  } else {
        bl.p <- lrtest(x[[1]], bl.conf)
  
    cox.table[i,] <- c(bl.summ[2,4], bl.an[1,3], bl.summ[4,4], bl.an[2,3], bl.summ[10,4], bl.an[3,3],bl.summ[12,4], bl.an[3,3], bl.summ[6,4], bl.an[4,3], bl.summ[8,4], bl.an[5,3], cox.bl$stats['Model L.R.'], bl.p$stats['P'], 1-(x[[2]]/cox.bl$stats['Model L.R.']), bl.val['Dxy', 'index.corrected']/2 + 0.5)
        return(list(cox.table, overfit.table))
  }
  
}

#### Calibration plot function #####
#### Input: predictors for Cox model, timepoint
#### Output: ggplot object

calibration_plot <- function(S, time.to.pace, gaf, caarms.cd, caarms.tc, pre.cond, pc1.thick, pc2.thick, t){
  
  set.seed(123)

  dd <- datadist(time.to.pace, gaf, caarms.cd, caarms.tc, pre.cond, pc1.thick, pc2.thick)
  options(datadist='dd')
  
  bl <- cph(S ~ time.to.pace + gaf + pre.cond , x=TRUE, y=TRUE, surv = TRUE, time.inc = t*12)
  bl.caarms <- cph(S ~ time.to.pace + gaf + pre.cond  + caarms.tc + caarms.cd, x=TRUE, y=TRUE, surv = TRUE, time.inc = t*12)
  bl.thick <- cph(S ~ time.to.pace + gaf + pre.cond  + pc1.thick + pc2.thick,  x=TRUE, y=TRUE, surv=TRUE, time.inc = t*12)
  
  bl.cal <- calibrate(bl.thick, method = "boot", u= t*12, m=10, B=1000, bw=FALSE)
  bl.caarms.cal <- calibrate(bl.caarms, method = "boot", u= t*12, m=10, B=1000, bw=FALSE)
  bl.base.cal <- calibrate(bl, method = "boot", u= t*12, m=10, B=1000, bw=FALSE)
  
  
  pred <- as.data.frame(1 - bl.cal[,'pred'])
  colnames(pred) <- c('x')
  
  cal.correct <- as.data.frame(1 - bl.cal[,'calibrated.corrected'])
  colnames(cal.correct) <- c('Value')
  cal.correct$Type <- 'MRI'
  
  pred2 <- as.data.frame(1 - bl.base.cal[,'pred'])
  colnames(pred2) <- c('x')
  
  cal.base.correct <- as.data.frame(1 - bl.base.cal[,'calibrated.corrected'])
  colnames(cal.base.correct) <- c('Value')
  cal.base.correct$Type <- 'clinical'
  
  pred3 <- as.data.frame(1 - bl.caarms.cal[,'pred'])
  colnames(pred3) <- c('x')
  
  cal.caarms.correct <- as.data.frame(1 - bl.caarms.cal[,'calibrated.corrected'])
  colnames(cal.caarms.correct) <- c('Value')
  cal.caarms.correct$Type <- 'with CAARMS'
  
  cal_plot_data <- rbind(cbind(pred,cal.correct), cbind(pred2, cal.base.correct), cbind(pred3, cal.caarms.correct))

  cal.plot <- ggplot(cal_plot_data, aes(x = x, y = Value, group = Type)) + 
         geom_line(aes(color = Type), size = 1.1) +
         scale_x_continuous(name =paste("Predicted ", as.character(t),"-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.2,0.4,0.6,0.8), limits=c(0, 0.9)) +
         scale_y_continuous(name=paste("Observed ", as.character(t),"-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.2,0.4,0.6,0.8), limits=c(0, 0.9)) +
         theme_bw() +
         labs(color = '') +
         geom_abline(aes(intercept = 0, slope = 1, color='black'))	+
          scale_color_manual(values = c("black" = "black", "clinical" = "#1b9e77",  "with CAARMS" = "#d95f02", "MRI" = "#7570b3"),
                             labels = c("ideal",  "base", "base + CAARMS subscales", "base + MRI")) +
         theme(plot.title = element_text(size = 18, hjust = 0.5),
               axis.title.x = element_text(size = 16),
               axis.title.y = element_text(size = 16),
               legend.text = element_text(size = 14),
               legend.title = element_text(size = 16),
               axis.text.x = element_text(size = 14),
               axis.text.y = element_text(size = 14)) +
         theme(legend.position="bottom")
  
  return(cal.plot)
  
  
}

calibration_plot_all <- function(S, time.to.pace, gaf, caarms.cd, caarms.tc, pre.cond){
  
  set.seed(123)

  dd <- datadist(time.to.pace, gaf, caarms.cd, caarms.tc, pre.cond)
  options(datadist='dd')
  
  bl <- cph(S ~ time.to.pace + gaf + pre.cond + caarms.cd + caarms.tc , x=TRUE, y=TRUE, surv = TRUE, time.inc = 2*12)
  
  bl.base.cal <- calibrate(bl, method = "boot", u= 2*12, m=10, B=1000, bw=FALSE)
  
  
  pred <- as.data.frame(1 - bl.base.cal[,'pred'])
  colnames(pred) <- c('x')
  
  cal.base.correct <- as.data.frame(1 - bl.base.cal[,'calibrated.corrected'])
  colnames(cal.base.correct) <- c('Value')
  cal.base.correct$Type <- paste("year",as.character(2), sep = " ")
  
  cal_plot_data <- cbind(pred, cal.base.correct)
  
  for (t in c(4,8,12)){
  
    bl <- cph(S ~ time.to.pace + gaf + pre.cond + caarms.cd + caarms.tc , x=TRUE, y=TRUE, surv = TRUE, time.inc = t*12)
    
    bl.base.cal <- calibrate(bl, method = "boot", u= t*12, m=10, B=1000, bw=FALSE)
    
    pred <- as.data.frame(1 - bl.base.cal[,'pred'])
    colnames(pred) <- c('x')
    
    cal.base.correct <- as.data.frame(1 - bl.base.cal[,'calibrated.corrected'])
    colnames(cal.base.correct) <- c('Value')
    if (t == 12){
      cal.base.correct$Type <- paste("zear",as.character(t), sep = " ")
    }else{
      cal.base.correct$Type <- paste("year",as.character(t), sep = " ")
    }
    
    cal_plot_data <- rbind(cal_plot_data, cbind(pred, cal.base.correct))
  
  }
  
  cal.plot <- ggplot(cal_plot_data, aes(x = x, y = Value, group = Type)) + 
         geom_line(aes(color = Type), size = 1.1) +
         scale_x_continuous(name =paste("Predicted Probability of Transition to Psychosis",sep = ""), breaks = c(0.2,0.4,0.6,0.8), limits=c(0, 0.9)) +
         scale_y_continuous(name=paste("Observed Probability of Transition to Psychosis",sep = ""), breaks = c(0.2,0.4,0.6,0.8), limits=c(0, 0.9)) +
         theme_bw() +
         labs(color = '') +
         geom_abline(aes(intercept = 0, slope = 1, color="black"))	+
         scale_color_manual(values = c( "black" = "black", "year 2" = "#a6cee3", "year 4" = "#1f78b4", "year 8" = "#b2df8a", "zear 12" = "#33a02c"),
                             labels = c( "ideal", "year 2", "year 4", "year 8", "year 12"),
                                     guide = guide_legend(override.aes = list(pch = c(16, NA, NA, NA, NA), linetype = c(1, 1, 1, 1, 1)))) +
         theme(plot.title = element_text(size = 18, hjust = 0.5),
               axis.title.x = element_text(size = 16),
               axis.title.y = element_text(size = 16),
               legend.text = element_text(size = 14),
               legend.title = element_text(size = 16),
               axis.text.x = element_text(size = 14),
               axis.text.y = element_text(size = 14)) +
         theme(legend.position="bottom")
  
  return(cal.plot)
  
  
}

#### LOAD DATA #####
# Load clinical measures
PACE_bl <- read_excel("path\to\baseline\data")

# Load follow-up data
PACE_fup <- read_excel("path\to\follow-up\data")
# Load MRI Quality Control data
PACE_MRIQC <- read_excel("path\to\MRI_QC\data")
# Load look-up table MRI ID and PACE 400 ID
caseid_MRI <- read_excel("path\to\MRI\data")

# Load MRI thickness data
PACE_MRIthick <- read.delim("path\to\cortical_thickness\data",header = TRUE, sep = ",")
# Load MRI curvature data
PACE_MRIcurve <- read.delim("path\to\cortical_curvature\data",header = TRUE, sep = ",")
# Load MRI surface data
PACE_MRIsurf <- read.delim("path\to\cortical_surface\data",header = TRUE, sep = ",")
#Load MRI volume data
PACE_MRIvol <- read.delim("path\to\cortical_volume\data",header = TRUE, sep = ",")

# Load neurocognition data
PACE_Cog <- read_excel("path\to\cognition\data")

# Load brain age values
PACE_brainage <- read.delim("path\to\brainage\data",header = TRUE, sep = ",")

#### MERGE TABLES ####

# Remove follow-up gaf and qlst from follow-up data (not needed in analysis, avoids two gaf and qlst columns after merger with baseline table)
PACE_fup <- PACE_fup[,!(names(PACE_fup) %in% c('gaf', 'qlst'))] 

# Merge baseline with follow-up according to PACE 400 id
PACE_clinical <- list(PACE_bl, PACE_fup) %>% reduce(inner_join, by="caseid")

PACE_all <- PACE_clinical # store data for cumulative hazard analysis

#### PREPARE CLINICAL DATA ####

# Remove subjects with not gaf or timepace information
PACE_clinical <- subset(PACE_clinical, !is.na(timepace))
PACE_clinical <- subset(PACE_clinical, !is.na(gaf))
PACE_clinical <- subset(PACE_clinical, !is.na(tc))
PACE_clinical <- subset(PACE_clinical, !is.na(cd))


# Calculate days until transition
PACE_clinical$transdays <- as.Date(as.character(PACE_clinical$enddate), format = "%Y-%m-%d") - as.Date(as.character(PACE_clinical$startdate), format = "%Y-%m-%d")

PACE_clinical$transmonths <- (as.yearmon(as.character(PACE_clinical$enddate), format = "%Y-%m-%d") - as.yearmon(as.character(PACE_clinical$startdate), format = "%Y-%m-%d"))*12

# Calculate UHR categories
PACE_clinical$pre.cond <- PACE_clinical$blips + 2*PACE_clinical$atten + 4*PACE_clinical$vulner
# Remove subject with no UHR category
PACE_clinical <- subset(PACE_clinical, !is.na(pre.cond))
# Code UHR category according to dominant condition (BLIPS -> Attenuated -> Family History)
PACE_clinical[PACE_clinical$pre.cond == 3 | PACE_clinical$pre.cond == 5  | PACE_clinical$pre.cond == 7, c('pre.cond')] <- 1
PACE_clinical[PACE_clinical$pre.cond == 6 , c('pre.cond')] <- 2
PACE_clinical[PACE_clinical$pre.cond == 4 , c('pre.cond')] <- 3

#### PREPARE NEUROCOGNITION TABLE ####

# Select neurocognition variables (age adjusted scales for Ward's subtests Arithmetic and Digit Symbol Coding)
PACE_neurocog <- PACE_Cog[,c('caseid','Ax1CodingASS','Ax1ArithASS', 'Ax1RAVLT_tot3')]

# Remove rows from MRI tables with less than 50%
PACE_neurocog <- PACE_neurocog[which(rowMeans(!is.na(PACE_neurocog)) > 0.5), ]
# Set subjects with missing Digit Symbol Coding or Arithmetic scores to NA
PACE_neurocog[PACE_neurocog$Ax1CodingASS < 0 & !is.na(PACE_neurocog$Ax1CodingASS), 'Ax1CodingASS'] <- NA
PACE_neurocog[PACE_neurocog$Ax1ArithASS < 0 & !is.na(PACE_neurocog$Ax1ArithASS), 'Ax1ArithASS'] <- NA
PACE_neurocog[PACE_neurocog$Ax1RAVLT_tot3 < 0 & !is.na(PACE_neurocog$Ax1RAVLT_tot3), 'Ax1RAVLT_tot3'] <- NA
# Remove subjects with NA values
PACE_neurocog <- subset(PACE_neurocog, !is.na(Ax1CodingASS))
PACE_neurocog <- subset(PACE_neurocog, !is.na(Ax1ArithASS))
PACE_neurocog <- subset(PACE_neurocog, !is.na(Ax1RAVLT_tot3))
# Merge neurocogntion data with clinical data according to PACE 400 ID
PACE_neurocog <- list(PACE_clinical, PACE_neurocog) %>% reduce(inner_join, by="caseid")

# Remove duplicated rows
PACE_neurocog <- PACE_neurocog[!duplicated(PACE_neurocog$caseid), ]

#### PREPARE NEUROIMAGING TABLE ####

## Merge caseid and MRI quality data according to MRI ID
mriQC <- merge(PACE_MRIQC,caseid_MRI, by.x = "BaselineMRI", by.y = "BaselineMRI")

# Remove rows from MRI tables with less than 99% data prior to merging
PACE_MRIthick <- PACE_MRIthick[which(rowMeans(!is.na(PACE_MRIthick)) > 0.99), ]
PACE_MRIcurve <- PACE_MRIcurve[which(rowMeans(!is.na(PACE_MRIcurve)) > 0.99), ]
PACE_MRIsurf <- PACE_MRIsurf[which(rowMeans(!is.na(PACE_MRIsurf)) > 0.99), ]
PACE_MRIvol <- PACE_MRIvol[which(rowMeans(!is.na(PACE_MRIvol)) > 0.99), ]

# Join all MRI tables according to PACE 400 id
PACE_MRI <- list(PACE_MRIthick, PACE_MRIcurve, PACE_MRIsurf, PACE_MRIvol) %>% reduce(inner_join, by=c("caseid","ID"))

# Merge clinical and all three MRI tables (Quality Control, MRI values, and brainage)
PACE_MRI <- merge(PACE_MRI, mriQC, by.x = "caseid", by.y = "caseid")
PACE_MRI <- merge(PACE_MRI, PACE_brainage, by.x = "caseid", by.y = "caseid")
PACE_MRI <- merge(PACE_clinical, PACE_MRI, by.x = "caseid", by.y = "caseid")

# Only select subjects with good MRI quality 
PACE_MRI <- subset(PACE_MRI, QC == 1)

# Remove duplicated rows
PACE_MRI <- PACE_MRI[!duplicated(PACE_MRI$caseid), ]

# Encode MRI site variable

PACE_MRI[PACE_MRI$MRIsite == "RCH" & !is.na(PACE_MRI$MRIsite), "MRIsite"] <- 0
PACE_MRI[PACE_MRI$MRIsite == "RMH" & !is.na(PACE_MRI$MRIsite), "MRIsite"] <- 1
PACE_MRI[PACE_MRI$MRIsite == "RCH (Orig 1.5)" & !is.na(PACE_MRI$MRIsite), "MRIsite"] <- 0
PACE_MRI[PACE_MRI$MRIsite == "BRI" & !is.na(PACE_MRI$MRIsite), "MRIsite"] <- 3

PACE_MRI$MRIsite_cat <- as.numeric(PACE_MRI$MRIsite == 3)

#### PREPARE TABLES FOR PCA ANALYSIS ####
PACE_MRIthick <- merge(PACE_MRI[,c('caseid','MRIsite')], PACE_MRIthick, by.x = 'caseid')
PACE_MRIcurve <- merge(PACE_MRI[,c('caseid','MRIsite')], PACE_MRIcurve, by.x = 'caseid')
PACE_MRIsurf <- merge(PACE_MRI[,c('caseid','MRIsite')], PACE_MRIsurf, by.x = 'caseid')
PACE_MRIvol <- merge(PACE_MRI[,c('caseid','MRIsite')], PACE_MRIvol, by.x = 'caseid')


##### Create descriptive table  #####

PACE_MRI_desc <- PACE_MRI %>% select(caseid, ttmt, sex, age_bl, timepace, vulner, blips, atten, bprst, sanst, gaf, qlst, tc, pa, cd, transtat, fuptime, transdays, sofas, corrected_gap)

PACE_MRI_desc$Ax1ArithASS <- NA
PACE_MRI_desc$Ax1CodingASS <- NA
PACE_MRI_desc$Ax1RAVLT_tot3 <- NA
PACE_MRI_desc$group <- "MRI"

PACE_neurocog_desc <- PACE_neurocog %>% select(caseid, ttmt, sex, age_bl, timepace, vulner, blips, atten, bprst, sanst, gaf, qlst, tc, pa, cd, transtat, fuptime, transdays, sofas, Ax1CodingASS, Ax1ArithASS, Ax1RAVLT_tot3)

PACE_neurocog_desc$corrected_gap <- NA
PACE_neurocog_desc$group <- "Cognition"

PACE_table <- rbind(PACE_MRI_desc, PACE_neurocog_desc)

# Calculate UHR categories
PACE_table$uhr_group <- PACE_table$blips + 2*PACE_table$atten + 4*PACE_table$vulner

# Remove subject with no UHR category
PACE_table <- subset(PACE_table, !is.na(uhr_group))

# Code UHR category according to dominant condition (BLIPS -> Attenuated -> Family History)
PACE_table[PACE_table$uhr_group == 3 | PACE_table$uhr_group == 5  | PACE_table$uhr_group == 7, c('uhr_group')] <- 1
PACE_table[PACE_table$uhr_group == 6 , c('uhr_group')] <- 2
PACE_table[PACE_table$uhr_group == 4 , c('uhr_group')] <- 3

PACE_table$transdays <- PACE_table$transdays * PACE_table$transtat

table_demog <- PACE_table %>% 
  select(group, ttmt, sex, age_bl, timepace, uhr_group, bprst, sanst, gaf, qlst, tc, pa, cd, transtat, fuptime, transdays, sofas, Ax1CodingASS, Ax1ArithASS, Ax1RAVLT_tot3, corrected_gap ) %>% # keep only columns of interest
  mutate(
    ttmt = factor(ttmt, labels = c("Intervention", "Standard")), 
    transtat = factor(transtat, labels = c("No Transition", "Transition")), 
    sex = factor(sex, labels = c("Male", "Female")),
    uhr_group = factor(uhr_group, labels = c("any BLIPS", "APS or APS+Trait", "Trait"))
  ) %>% 
  tbl_summary(     
    by = group,                                               # stratify entire table by outcome
    statistic = list(all_continuous() ~ "{mean}+/-{sd} ({N_nonmiss})",        # stats and format for continuous columns
                     all_categorical() ~ "{p}% ({n}/{N}) "),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = list(c(ttmt, transtat, sex) ~ "dichotomous",
                  c(age_bl, timepace, bprst, sanst, gaf, qlst, tc, pa, cd, fuptime, transdays, sofas, Ax1CodingASS, Ax1ArithASS, Ax1RAVLT_tot3, corrected_gap) ~ "continuous",
                  c(uhr_group) ~ "categorical"),
    value = list(sex ~ "Female",
                 transtat ~ "Transition",
                 ttmt ~ "Intervention"),
    label  = list(                                              # display labels for column names
      age_bl ~ "Age at baseline (years)",
      sex ~ "Gender",
      ttmt ~ "Treatment",
      uhr_group ~ "UHR subgroup",
      timepace ~ "Time between symptom onset and first contact with PACE (days)",
      tc ~ "CAARMS Disorders of Thought Content, severity",
      pa ~ "CAARMS Perceptual Abnormalities, severity",
      cd ~ "CAARMS Conceptual Disorganisation, severity",
      bprst ~ "BPRS Total",
      sanst ~ "SANS Total",
      gaf ~ "GAF",
      qlst ~ "QLS Total",
      Ax1CodingASS ~ "Coding",
      Ax1ArithASS ~ "Arithmetic",
      Ax1RAVLT_tot3 ~ "RAVLT total",
      corrected_gap ~ "Brain Age Gap",
      transtat ~ "Transition status",
      transdays ~ "Time to transition (days)",
      fuptime ~ "Follow-up time (days)",
      sofas ~ "SOFAS at follow-up"
    ),
    missing = "no", # don't list missing data separately;
  ) %>%
  modify_header(label = "Variable") %>% # update the column header
  bold_labels()

table_demog %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "U:/PostDoc/Research/Publications/Own/PACE_400_Cox/demog_table.docx")

##### Plot cumulative hazard functions  #####

PACE_all$transyear <- (as.yearmon(as.character(PACE_all$enddate), format = "%Y-%m-%d") - as.yearmon(as.character(PACE_all$startdate), format = "%Y-%m-%d"))

# Calculate UHR categories
PACE_all$uhr_group <- PACE_all$blips + 2*PACE_all$atten + 4*PACE_all$vulner

# Remove subject with no UHR category
PACE_all <- subset(PACE_all, !is.na(uhr_group))

# Code UHR category according to dominant condition (BLIPS -> Attenuated -> Family History)
PACE_all[PACE_all$uhr_group == 3 | PACE_all$uhr_group == 5  | PACE_all$uhr_group == 7, c('uhr_group')] <- 1
PACE_all[PACE_all$uhr_group == 6 , c('uhr_group')] <- 2
PACE_all[PACE_all$uhr_group == 4 , c('uhr_group')] <- 3

fit <- survfit(Surv(transyear, transtat) ~ uhr_group, data = PACE_all)

cumhaz_plot <- ggsurvplot(fit,
          conf.int = FALSE,
          pval = TRUE,
          pval.method = TRUE,
          log.rank.weights = "n",
          pval.coord = c(12.5,0.1),
          pval.method.coord = c(12.5,0.15),
          pval.method.size = 5,
          xlab = "Time in years",   # customize X axis label.
          break.time.by = 5,     # break X axis in time intervals by 200.
          ylab = "Cumulative hazard function for transitioning to FEP",
          legend.labs = c("BLIPS", "Attenuated psychosis", "Vulnerability"), 
          risk.table.col = "strata", # Change risk table color by groups
          ggtheme = theme_bw(), # Change ggplot2 theme
          palette = c("#1b9e77", "#d95f02", "#7570b3"), # check colorbrewer2.org for color palettes
          fun = "cumhaz")

cumhaz_plot


#####Principal component analysis  #####
# Each structural map contains 34 variables each linked to a region of the brain. To reduce the dimension of the feature space, we applied principal component analysis (PCA). 

# Split each domain in left and right hemisphere
PACE_MRIthick.caseid <- PACE_MRIthick$caseid
PACE_MRIthick.bl <- PACE_MRIthick[,c(2,4:(ncol(PACE_MRIthick)-1))]
PACE_MRIthick.bl <- subset(PACE_MRIthick.bl, select = -c(lh_MeanThickness_thickness))

PACE_MRIcurve.caseid <- PACE_MRIcurve$caseid
PACE_MRIcurve.bl <- PACE_MRIcurve[,c(2,4:ncol(PACE_MRIcurve))]

PACE_MRIsurf.caseid <- PACE_MRIsurf$caseid
PACE_MRIsurf.bl <- PACE_MRIsurf[,c(2,4:(ncol(PACE_MRIsurf)-1))]
PACE_MRIsurf.bl <- subset(PACE_MRIsurf.bl, select = -c(lh_WhiteSurfArea_area))

PACE_MRIvol.caseid <- PACE_MRIvol$caseid
PACE_MRIvol.bl <- PACE_MRIvol[,c(2,4:ncol(PACE_MRIvol))]




##### Harmonize MRI data based on scanner site  #####

# Split each domain in left and right hemisphere
PACE_MRIthick.bl.harm <- ComBat(t(PACE_MRIthick.bl[,c(2:ncol(PACE_MRIthick.bl))]), PACE_MRIthick.bl$MRIsite, mod = PACE_MRI$transtat)
PACE_MRIthick.bl.harm <- as.data.frame(t(PACE_MRIthick.bl.harm))

PACE_MRIcurve.bl.harm <- ComBat(t(PACE_MRIcurve.bl[,c(2:ncol(PACE_MRIcurve.bl))]), PACE_MRIcurve.bl$MRIsite, mod = PACE_MRI$transtat)
PACE_MRIcurve.bl.harm <- as.data.frame(t(PACE_MRIcurve.bl.harm))

PACE_MRIsurf.bl.harm <- ComBat(t(PACE_MRIsurf.bl[,c(2:ncol(PACE_MRIsurf.bl))]), PACE_MRIsurf.bl$MRIsite, mod = PACE_MRI$transtat)
PACE_MRIsurf.bl.harm <- as.data.frame(t(PACE_MRIsurf.bl.harm))

PACE_MRIvol.bl.harm <- ComBat(t(PACE_MRIvol.bl[,c(2:ncol(PACE_MRIvol.bl))]), PACE_MRIvol.bl$MRIsite, mod = PACE_MRI$transtat)
PACE_MRIvol.bl.harm <- as.data.frame(t(PACE_MRIvol.bl.harm))

#### Calculate PCA for each domain and hemisphere ####

# Thickness
PACE_MRIthick.pcabl.harm <- prcomp(PACE_MRIthick.bl.harm, center = TRUE,scale. = TRUE, rank.=10)
# Select the first thickness PCA for each hemisphere
MRI.thick <- as.data.frame(cbind(PACE_MRIthick.caseid, PACE_MRIthick.pcabl.harm$x[,'PC1'], PACE_MRIthick.pcabl.harm$x[,'PC2']))
# Rename selected PCAs
colnames(MRI.thick) <- c('caseid', 'thickness_bl_harm1', 'thickness_bl_harm2')


# Volume
PACE_MRIvol.pcabl.harm  <- prcomp(PACE_MRIvol.bl.harm, center = TRUE,scale. = TRUE, rank.=10)
# Select the first volume PCA for each hemisphere
MRI.vol <- as.data.frame(cbind(PACE_MRIvol.caseid, PACE_MRIvol.pcabl.harm$x[,'PC1'], PACE_MRIvol.pcabl.harm$x[,'PC2']))
# Rename selected PCAs
colnames(MRI.vol) <- c('caseid', 'volume_bl_harm1', 'volume_bl_harm2')


# Surface area
PACE_MRIsurf.pcabl.harm  <- prcomp(PACE_MRIsurf.bl.harm, center = TRUE,scale. = TRUE, rank.=10)
# Select the first surface area PCA for each hemisphere
MRI.surf <- as.data.frame(cbind(PACE_MRIsurf.caseid, PACE_MRIsurf.pcabl.harm$x[,'PC1'], PACE_MRIsurf.pcabl.harm$x[,'PC2']))
# Rename selected PCAs
colnames(MRI.surf) <- c('caseid', 'surface_bl_harm1', 'surface_bl_harm2')


# Curvature
PACE_MRIcurve.pcabl.harm  <- prcomp(PACE_MRIcurve.bl.harm, center = TRUE,scale. = TRUE, rank.=10)
# Select the first curvature PCA for each hemisphere
MRI.curve <- as.data.frame(cbind(PACE_MRIcurve.caseid, PACE_MRIcurve.pcabl.harm$x[,'PC1'], PACE_MRIcurve.pcabl.harm$x[,'PC2']))
# Rename selected PCAs
colnames(MRI.curve) <- c('caseid', 'curvature_bl_harm1', 'curvature_bl_harm2')

# Merge all selected MRI principal components
MRI.all <- merge(MRI.thick, MRI.vol, by.x = "caseid", by.y = "caseid")
MRI.all <- merge(MRI.all, MRI.curve, by.x = "caseid", by.y = "caseid")
MRI.all <- merge(MRI.all, MRI.surf, by.x = "caseid", by.y = "caseid")

# Merge MRI commponents with PACE data
PACE_MRI <- merge(PACE_MRI, MRI.all, by.x = "caseid", by.y = "caseid")


##### Create PCA table  #####

pc <- data.frame(Thickness=double(),
                 Volume=double(), 
                 Surface=double(), 
                 Curve=double()) 

pc <- pca.table(pc, 1, PACE_MRIthick.pcabl.harm)
pc <- pca.table(pc, 4, PACE_MRIcurve.pcabl.harm)
pc <- pca.table(pc, 2, PACE_MRIvol.pcabl.harm)
pc <- pca.table(pc, 3, PACE_MRIsurf.pcabl.harm)
pc <- round(pc, 2)
row.names(pc) <- c(paste('PC', c(1:10), sep=''))


##### Show PCA table #####

pc %>%
  kbl(caption = "Table 1: Proportion of variance [%] included in first 10 principal components of each structural map") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  pack_rows("Bilateral", 1, 10)


##### Boxplots of principal components categorised by MRI site #####


# plot_thickness <- ggplot(PACE_MRI, aes(x=factor(MRIsite), y=thickness_bl, fill=as.factor(transtat)))+
#   geom_boxplot() +
#   labs(fill='Transition status') + 
#   ggtitle("Original values") +
#   scale_y_continuous(breaks=c(seq(-10,10,1)), limits = c(-10, 10))+
#   labs(x="MRI site", y="Thickness left hemisphere PCA score") +
#   scale_x_discrete(labels=c("RCH","RMH", "BRI")) +
#   scale_fill_brewer(labels=c('No Transition', 'Transition'), palette="Dark2") + 
#   theme_minimal() +
#   theme(plot.title = element_text(size=15, hjust = 0.5),
#         axis.text.x = element_text(size = 10),
#         axis.title.x = element_text(size = 15),
#         axis.text.y = element_text(size = 10),
#         axis.title.y = element_text(size =15),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 10)) 
# 
# plot_thickness_harm <- ggplot(PACE_MRI, aes(x=factor(MRIsite), y=thickness_bl_harm1, fill=as.factor(transtat)))+
#   geom_boxplot() +
#   labs(fill='Transition status') + 
#   ggtitle("Harmonized values") +
#   scale_y_continuous(breaks=c(seq(-10,10,1)), limits = c(-10, 10))+
#   labs(x="MRI site", y="Thickness left hemisphere PCA score") +
#   scale_x_discrete(labels=c("RCH","RMH", "BRI")) +
#   scale_fill_brewer(labels=c('No Transition', 'Transition'), palette="Dark2") + 
#   theme_minimal() +
#   theme(plot.title = element_text(size=15, hjust = 0.5),
#         axis.text.x = element_text(size = 10),
#         axis.title.x = element_text(size = 15),
#         axis.text.y = element_text(size = 10),
#         axis.title.y = element_text(size = 15),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 10)) 
# 
# plot_thickness
# plot_thickness_harm
# 
# ggarrange(plot_thickness, plot_thickness_harm, ncol = 2, legend = "bottom", common.legend = TRUE)


PACE_MRIthick.bl.harm$MRIsite <- PACE_MRIthick$MRIsite
PACE_MRIthick.bl.harm$transtat <-   PACE_MRI$transtat

plot_thickness2 <- ggplot(PACE_MRI, aes(x=factor(MRIsite), y=rh_rostralmiddlefrontal_thickness, fill=as.factor(transtat)))+
  geom_boxplot() +
  labs(fill='Transition status') + 
  ggtitle("Original values") +
  scale_y_continuous(breaks=c(seq(-10,10,1)), limits = c(2, 3.5))+
  labs(x="MRI site", y="Thickness right rostralmiddlefrontal") +
  scale_x_discrete(labels=c("RCH","RMH", "BRI")) +
  scale_fill_brewer(labels=c('No Transition', 'Transition'), palette="Dark2") + 
  theme_minimal() +
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size =15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) 

plot_thickness_harm2 <- ggplot(PACE_MRIthick.bl.harm, aes(x=factor(MRIsite), y=rh_rostralmiddlefrontal_thickness, fill=as.factor(transtat)))+
  geom_boxplot() +
  ggtitle("Harmonized values") +
  labs(fill='Transition status') + 
  scale_y_continuous(breaks=c(seq(-10,10,1)), limits = c(2, 3.5))+
  labs(x="MRI site", y="Thickness right rostralmiddlefrontal") +
  scale_x_discrete(labels=c("RCH","RMH","BRI")) +
  scale_fill_brewer(labels=c('No Transition', 'Transition'), palette="Dark2") + 
  theme_minimal() +
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) 


ggarrange(plot_thickness2, plot_thickness_harm2, ncol = 2, legend = "bottom", common.legend = TRUE)


##### Visualize PCA results #####

pc_bl <- pc[1:10,]*100
pc_bl$x <- 1:10
pc_bl_plot <- melt(pc_bl, id.vars = "x")

pca_left_plot <- pca_plot(pc_bl_plot, 'Bilateral')
pca_left_plot

##### Thickness loading #####

thick.loading <- loading_plot(PACE_MRIthick.pcabl.harm, 'Thickness')
vol.loading <- loading_plot(PACE_MRIvol.pcabl.harm,  'Volume')
curve.loading <- loading_plot(PACE_MRIcurve.pcabl.harm, 'Mean curvature')
surf.loading <- loading_plot(PACE_MRIsurf.pcabl.harm,  'Area')

ggarrange(thick.loading, vol.loading, curve.loading, surf.loading, ncol = 2, nrow = 2)

##### Cox regression analysis #####

PACE_MRI[,names(PACE_MRIthick.bl.harm)] <- PACE_MRIthick.bl.harm

PACE_MRI[,names(PACE_MRIcurve.bl.harm)] <- PACE_MRIcurve.bl.harm

PACE_MRI[,names(PACE_MRIsurf.bl.harm)] <- PACE_MRIsurf.bl.harm

PACE_MRI[,names(PACE_MRIvol.bl.harm)] <- PACE_MRIvol.bl.harm


# Calculate time to event and event status
dtime <- PACE_MRI$transmonths
event <- PACE_MRI$transtat


# Preparation for imputation
w <- transcan(~ timepace + gaf + cd + tc + ttmt + pre.cond + surface_bl_harm1 + surface_bl_harm2 + thickness_bl_harm1 + thickness_bl_harm2 + volume_bl_harm1 + volume_bl_harm2 + curvature_bl_harm1 + curvature_bl_harm2 + lh_paracentral_thickness + rh_paracentral_thickness + rh_superiortemporal_thickness + lh_fusiform_thickness  + corrected_gap + age_bl, imputed=TRUE, data = PACE_MRI, pl=FALSE, pr = FALSE)

# create variables out of columns
attach(PACE_MRI)

# Imputation (theoretically not needed)
time.to.pace <- log(impute(w, timepace, data = PACE_MRI))
gaf <- impute(w, gaf, data = PACE_MRI)
caarms.tc <- impute(w, tc, data = PACE_MRI)
caarms.cd <- impute(w, cd, data = PACE_MRI)
ttmt <- impute(w, ttmt, data = PACE_MRI)
pre.cond <- impute(w, pre.cond, data = PACE_MRI)
pre.cond <- relevel(factor(pre.cond),ref='3')
pc1.surf <- impute(w, surface_bl_harm1, data = PACE_MRI)
pc2.surf <- impute(w, surface_bl_harm2, data = PACE_MRI)
pc1.thick <- impute(w, thickness_bl_harm1, data = PACE_MRI)
pc2.thick <- impute(w, thickness_bl_harm2, data = PACE_MRI)
pc1.vol <- impute(w, volume_bl_harm1, data = PACE_MRI)
pc2.vol <- impute(w, volume_bl_harm2, data = PACE_MRI)
pc1.curve <- impute(w, curvature_bl_harm1, data = PACE_MRI)
pc2.curve <- impute(w, curvature_bl_harm2, data = PACE_MRI)
para_left <- impute(w, lh_paracentral_thickness, data = PACE_MRI)
para_right <- impute(w, rh_paracentral_thickness, data = PACE_MRI)
superior_right <- impute(w, lh_fusiform_thickness, data = PACE_MRI)
fusi_left <- impute(w, curvature_bl_harm2, data = PACE_MRI)
brain_age <- impute(w, corrected_gap, data = PACE_MRI)
age <- impute(w, age_bl, data = PACE_MRI)
# Preperation for cox model fit
dd <- datadist(time.to.pace, gaf, caarms.cd, caarms.tc, ttmt, pre.cond, pc1.surf, pc2.surf, pc1.thick, pc2.thick, pc1.vol, pc2.vol, pc1.curve, pc2.curve, para_left, para_right, superior_right, fusi_left, brain_age, age)
options(datadist='dd')
units(dtime) <- 'Month'

# Define Surv model
S <- Surv(dtime, event)


##### Cox model using clinical and neuroimaging variables #####

set.seed(12)
# Create empty data frames for regression results table
cox.table <- data.frame(hr1=double(),
                  p1=double(),
                  hr2=double(),
                 p2=double(),
                 hr3=double(),
                 p3=double(),
                 hr4=double(),
                 p4=double(),
                 hr5=double(),
                 p5=double(),
                 hr6=double(),
                 p6=double(),
                 LR_chisq=double(),
                 p = double(),
                 Fraction=double(),
                 cIndex=double()) 

# Create empty data frames for overfitting results table
overfit.table <- data.frame(delta_Dxy=double(),
                 delta_R2=double(), 
                 delta_slope=double())

cox.bl <- cph(S ~ time.to.pace + gaf + pre.cond , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
# Get bootstrapped model
bl.conf <- bootcov(cox.bl, B=1000, pr=TRUE)
LR.orig <- cox.bl$stats['Model L.R.']
cox.res <- cox_model(cox.bl, 1, cox.table, overfit.table)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond  + caarms.tc, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 2, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond  + caarms.cd, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 3, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + pc1.surf , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 4, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + pc1.curve , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 5, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + pc1.vol , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 6, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + pc1.thick , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 7, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + para_left , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 8, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + para_right , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 9, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + superior_right , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 10, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + fusi_left , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 11, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)


cox.table <- cox.res[[1]]
overfit.table <- cox.res[[2]]

cox.table <- round(cox.table, digits = 3)
overfit.table <- round(overfit.table, digits = 3)

row.names(overfit.table) <- c('Clinical','Clinical + UTC', 'Clinical + DC','Clinical + MRI surface','Clinical + MRI curve', 'Clinical + MRI volume', 'Clinical + MRI thickness', 'Clinical + Left Paracentral', 'Clinical + Right Paracentral',  'Clinical + Right Superiortemporal', 'Clinical + Left fusiform')

row.names(cox.table) <-  c('Clinical','Clinical + UTC', 'Clinical + DC','Clinical + MRI surface','Clinical + MRI curve', 'Clinical + MRI volume', 'Clinical + MRI thickness', 'Clinical + Left Paracentral', 'Clinical + Right Paracentral',  'Clinical + Right Superiortemporal', 'Clinical + Left fusiform')

colnames(cox.table) <- c('HR timepace', 'P timepace', 'HR gaf','P gaf', 'HR UHR','P UHR', 'HR2 UHR', 'P UHR', 'HR 4', 'P 4', 'HR 5', 'P 5', 'LR', 'LR P', 'Fraction of new information', 'C index')

colnames(overfit.table) <- c('Optimism Dxy','Optimism R2','Optimism Slope')

#Table lists the Cox model fit measures (Somer's D~xy, R)

cox.table %>%
  kbl(caption = "Table 2: Comparison of model fit measures between nested models using Caarms, cognition, and structural imaging data") %>%
  kable_classic(full_width = F, html_font = "Cambria")

##### Cox model using clinical and brainage variables #####

set.seed(12)
# Create empty data frames for regression results table
cox.table <- data.frame(hr1=double(),
                  p1=double(),
                  hr2=double(),
                 p2=double(),
                 hr3=double(),
                 p3=double(),
                 hr4=double(),
                 p4=double(),
                 hr5=double(),
                 p5=double(),
                 hr6=double(),
                 p6=double(),
                 LR_chisq=double(),
                 p = double(),
                 Fraction=double(),
                 cIndex=double()) 

# Create empty data frames for overfitting results table
overfit.table <- data.frame(delta_Dxy=double(),
                 delta_R2=double(), 
                 delta_slope=double())

cox.bl <- cph(S ~ time.to.pace + gaf + pre.cond , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
# Get bootstrapped model
bl.conf <- bootcov(cox.bl, B=1000, pr=TRUE)
LR.orig <- cox.bl$stats['Model L.R.']
cox.res <- cox_model(cox.bl, 1, cox.table, overfit.table)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond  + caarms.tc, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 2, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + brain_age + age , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 3, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + age , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 4, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + brain_age , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 5, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)


cox.table <- cox.res[[1]]
overfit.table <- cox.res[[2]]

cox.table <- round(cox.table, digits = 3)
overfit.table <- round(overfit.table, digits = 3)

row.names(overfit.table) <- c('Clinical','Clinical + UTC','Clinical + Brain Age + Age','Clinical + Age', 'Clinical + Brain Age')

row.names(cox.table) <- c('Clinical','Clinical + UTC','Clinical + Brain Age + Age','Clinical + Age', 'Clinical + Brain Age')

colnames(cox.table) <- c('HR timepace', 'P timepace', 'HR gaf','P gaf', 'HR UHR','P UHR', 'HR2 UHR', 'P UHR', 'HR 4', 'P 4', 'HR 5', 'P 5', 'LR', 'LR P', 'Fraction of new information', 'C index')

colnames(overfit.table) <- c('Optimism Dxy','Optimism R2','Optimism Slope')


# Table lists the Cox model fit measures (Somer's D~xy, R)


cox.table %>%
  kbl(caption = "Table 2: Comparison of model fit measures between nested models using Caarms or Brain Age Gap") %>%
  kable_classic(full_width = F, html_font = "Cambria")

##### Cox model using clinical and neuroimaging combined variables #####

set.seed(12)
# Create empty data frames for regression results table
cox.table <- data.frame(hr1=double(),
                  p1=double(),
                  hr2=double(),
                 p2=double(),
                 hr3=double(),
                 p3=double(),
                 hr4=double(),
                 p4=double(),
                 hr5=double(),
                 p5=double(),
                 hr6=double(),
                 p6=double(),
                 LR_chisq=double(),
                 p = double(),
                 Fraction=double(),
                 cIndex=double()) 

# Create empty data frames for overfitting results table
overfit.table <- data.frame(delta_Dxy=double(),
                 delta_R2=double(), 
                 delta_slope=double())

cox.bl <- cph(S ~ time.to.pace + gaf + pre.cond , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
# Get bootstrapped model
bl.conf <- bootcov(cox.bl, B=1000, pr=TRUE)
LR.orig <- cox.bl$stats['Model L.R.']
cox.res <- cox_model(cox.bl, 1, cox.table, overfit.table)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond  + caarms.tc + caarms.cd, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 2, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + pc1.vol + pc1.thick , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 3, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + pc1.thick + pc2.thick , x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 4, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + para_left  +  para_right, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 5, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + para_left  +  fusi_left, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 6, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)


cox.table <- cox.res[[1]]
overfit.table <- cox.res[[2]]

cox.table <- round(cox.table, digits = 3)
overfit.table <- round(overfit.table, digits = 3)

row.names(overfit.table) <- c('Clinical','Clinical + UTC','Clinical + Volume + Thickness','Clinical + Thickness','Clinical + Para left and right', 'Clinical + para + fusiform')

row.names(cox.table) <- c('Clinical','Clinical + UTC','Clinical + Volume + Thickness','Clinical + Thickness','Clinical + Para left and right', 'Clinical + para + fusiform')

colnames(cox.table) <- c('HR timepace', 'P timepace', 'HR gaf','P gaf', 'HR UHR','P UHR', 'HR2 UHR', 'P UHR', 'HR 4', 'P 4', 'HR 5', 'P 5', 'LR', 'LR P', 'Fraction of new information', 'C index')

colnames(overfit.table) <- c('Optimism Dxy','Optimism R2','Optimism Slope')


# Table lists the Cox model fit measures (Somer's D~xy, R)

cox.table %>%
  kbl(caption = "Table 2: Comparison of model fit measures between nested models using Caarms or Brain Age Gap") %>%
  kable_classic(full_width = F, html_font = "Cambria")

##### Calibration neuroimaging model #####

year2_plot <- calibration_plot(S, time.to.pace, gaf, caarms.cd, caarms.tc, pre.cond, para_left, fusi_left, 2)

year2_plot


##### Cox Model using clinical variables and neurocognition variables #####

dtime <- PACE_neurocog$transmonths
event <- PACE_neurocog$transtat

rm(time.to.pace, gaf, ttmt, pre.cond, caarms.tc, caarms.cd)

# Preperation for imputation
w <- transcan(~ timepace + gaf +  ttmt + pre.cond + tc + cd + Ax1CodingASS + Ax1ArithASS + Ax1RAVLT_tot3, imputed=TRUE, data = PACE_neurocog, pl=FALSE, pr = FALSE)

# create variables out of columns
attach(PACE_neurocog)

# Imputation
time.to.pace <- log(impute(w, timepace, data = PACE_neurocog))
gaf <- impute(w, gaf, data = PACE_neurocog)
ttmt <- impute(w, ttmt, data = PACE_neurocog)
pre.cond <- impute(w, pre.cond, data = PACE_neurocog)
cog1 <- impute(w, Ax1CodingASS, data = PACE_neurocog)
cog2 <- impute(w, Ax1ArithASS, data = PACE_neurocog)
cog3 <- impute(w, Ax1RAVLT_tot3, data = PACE_neurocog)
caarms.tc <- impute(w, tc, data = PACE_neurocog)
caarms.cd <- impute(w, cd, data = PACE_neurocog)

pre.cond <- relevel(factor(pre.cond),ref=3)
# Preperation for cox model fit
dd <- datadist(time.to.pace, gaf, ttmt, pre.cond, cog1, cog2, cog3, caarms.cd, caarms.tc)
options(datadist='dd')
units(dtime) <- 'Month'

# Define Surv model
S <- Surv(dtime, event)


# Create empty data frames for regression results table
cox.table <- data.frame(hr1=double(),
                  p1=double(),
                  hr2=double(),
                 p2=double(),
                 hr3=double(),
                 p3=double(),
                 hr4=double(),
                 p4=double(),
                 hr5=double(),
                 p5=double(),
                 hr6=double(),
                 p6=double(),
                 LR_chisq=double(),
                 p = double(),
                 Fraction=double(),
                 cIndex=double()) 

# Create empty data frames for overfitting results table
overfit.table <- data.frame(delta_Dxy=double(),
                 delta_R2=double(), 
                 delta_slope=double())

set.seed(12)

cox.bl <- cph(S ~ time.to.pace + gaf + pre.cond, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
# Get bootstrapped model
bl.conf <- bootcov(cox.bl, B=1000, pr=TRUE)
LR.orig <- cox.bl$stats['Model L.R.']
cox.res <- cox_model(cox.bl, 1, cox.table, overfit.table)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + caarms.tc, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 2, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + caarms.cd, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 3, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + cog1, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 4, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + cog2, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 5, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + cog3, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 6, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

cox.table <- cox.res[[1]]
overfit.table <- cox.res[[2]]

cox.table <- round(cox.table, digits = 3)
overfit.table <- round(overfit.table, digits = 3)

row.names(overfit.table) <- c('Clinical','Clinical + thought content', 'Clinical + disorganised','Clinical + Coding', 'Clinical + Arithmetic', 'Clinical + RAVLT')

row.names(cox.table) <- c('Clinical','Clinical + thought content', 'Clinical + disorganised','Clinical + Coding', 'Clinical + Arithmetic', 'Clinical + RAVLT')

colnames(cox.table) <- c('HR timepace', 'P timepace', 'HR gaf','P gaf', 'HR UHR','P UHR', 'HR2 UHR', 'P UHR', 'HR 4', 'P 4', 'HR 5', 'P 5', 'LR', 'LR P', 'Fraction of new information', 'C index')

colnames(overfit.table) <- c('Optimism Dxy','Optimism R2','Optimism Slope')


# Table lists the Cox model fit measures (Somer's D~xy, R)

cox.table %>%
  kbl(caption = "Table 2: Comparison of model fit measures between nested models using Caarms, and cognition") %>%
  kable_classic(full_width = F, html_font = "Cambria")

##### Cox Model using clinical variables and combined neurocognition variables #####

# Create empty data frames for regression results table
cox.table <- data.frame(hr1=double(),
                  p1=double(),
                  hr2=double(),
                 p2=double(),
                 hr3=double(),
                 p3=double(),
                 hr4=double(),
                 p4=double(),
                 hr5=double(),
                 p5=double(),
                 hr6=double(),
                 p6=double(),
                 LR_chisq=double(),
                 p = double(),
                 Fraction=double(),
                 cIndex=double()) 

# Create empty data frames for overfitting results table
overfit.table <- data.frame(delta_Dxy=double(),
                 delta_R2=double(), 
                 delta_slope=double())

set.seed(12)

cox.bl <- cph(S ~ time.to.pace + gaf + pre.cond, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
# Get bootstrapped model
bl.conf <- bootcov(cox.bl, B=1000, pr=TRUE)
LR.orig <- cox.bl$stats['Model L.R.']
cox.res <- cox_model(cox.bl, 1, cox.table, overfit.table)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + caarms.tc + caarms.cd, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 2, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + cog1 + cog2, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 3, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + cog2 + cog3, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 4, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

bl.multi <- cph(S ~ time.to.pace + gaf + pre.cond + cog3 + cog1, x=TRUE, y=TRUE, surv = TRUE, time.inc = 10*12)
cox.res <- cox_model(bl.multi, 5, cox.res[[1]], cox.res[[2]], bl.conf, LR.orig)

cox.table <- cox.res[[1]]
overfit.table <- cox.res[[2]]

cox.table <- round(cox.table, digits = 3)
overfit.table <- round(overfit.table, digits = 3)

row.names(overfit.table) <- c('Clinical','Clinical + Caarms','Clinical + Coding/Arithmetic', 'Clinical + Arithmetic/RAVLT', 'Clinical + RAVLT/Coding')

row.names(cox.table) <- c('Clinical','Clinical + Caarms','Clinical + Coding/Arithmetic', 'Clinical + Arithmetic/RAVLT', 'Clinical + RAVLT/Coding')

colnames(cox.table) <- c('HR timepace', 'P timepace', 'HR gaf','P gaf', 'HR UHR','P UHR', 'HR2 UHR', 'P UHR', 'HR 4', 'P 4', 'HR 5', 'P 5', 'LR', 'LR P', 'Fraction of new information', 'C index')

colnames(overfit.table) <- c('Optimism Dxy','Optimism R2','Optimism Slope')

# Table lists the Cox model fit measures (Somer's D~xy, R)

cox.table %>%
  kbl(caption = "Table 2: Comparison of model fit measures between nested models using Caarms and cognition") %>%
  kable_classic(full_width = F, html_font = "Cambria")



#### Calibration neurocognition model ####

year2_plot <- calibration_plot(S, time.to.pace, gaf, caarms.cd, caarms.tc, pre.cond, cog1, cog2, 2)

year2_plot


