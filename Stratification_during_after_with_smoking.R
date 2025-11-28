#Author: Christina Chatsatourian 

# During-After_with_smoking_PRS_stratification_model

library(survival)

# Adjust `tstop` to avoid issues where `tstart >= tstop`
# Set `tstop` to the minimum of `age` and `bc_age` (censor at BC diagnosis if applicable)
filtered_data$tstop <- pmin(filtered_data$age, filtered_data$bc_age, na.rm = TRUE)

# Filter out rows where `tstop` is less than or equal to 0 (ensuring valid `tstop` values)
filtered_data <- filtered_data[filtered_data$tstop > 0, ]  # Ensure no negative or zero `tstop`

# Initialize the timed data frame with age as the primary time scale
# Adjust tstop to censor at breast cancer diagnosis (bc_age)
t_data <- tmerge(
  data1 = filtered_data, 
  data2 = filtered_data, 
  id = ID, 
  tstart = 0,                        # Start time (age at birth)
  tstop = pmin(age, bc_age, na.rm = TRUE)  # Stop time (censor at BC diagnosis if applicable)
)

#  Add the event status for BC diagnosis
t_data <- tmerge(
  data1 = t_data, 
  data2 = t_data, 
  id = ID, 
  bc_diagnosis = event(bc_age)       # Track BC event 
)

# View the result
head(t_data)

# Add time-dependent smoking statuses
# Current smoking period
smoke_data <- tmerge(
  data1 = t_data,
  data2 = filtered_data,
  id = ID,
  smoking_cur_start = tdc(age_current_start) # Smoking starts (current smokers)
)

# Former smoking period starts
smoke_data <- tmerge(
  data1 = smoke_data,
  data2 = filtered_data,
  id = ID,
  smoking_for_start = tdc(age_former_start)  # Smoking starts (former smoker)
)

# Former smoking period stops
smoke_data <- tmerge(
  data1 = smoke_data,
  data2 = filtered_data,
  id = ID,
  smoking_for_stop = tdc(age_former_stopp)   # Smoking stops (transition to former smoker)
)

# Compute the smoking status
smoke_data$smoking <- with(smoke_data, {
  # Start with everyone as non-smoking (0)
  status <- rep(0, nrow(smoke_data))  # Explicitly set default to 0
  
  # Current smoker (1) from age_current_start
  status[smoking_cur_start == 1] <- 1
  
  # Former smoker (2) between age_former_start and age_former_stopp
  status[smoking_for_start == 1] <- 1 # Mark as current smoker during former smoking period
  status[smoking_for_stop == 1] <- 2 # Transition to former smoker after stopping
  
  # Return updated smoking status
  status
})

# Define non- smoking (0) as the reference group
smoke_data$SmokingStatus <- factor(
  smoke_data$smoking,  # Keep the original `smoking` variable (with NA values)
  levels = c(0, 1, 2),  # Specify the levels: 0 = "non-smoking", 1 = "current smoker", 2 = "former smoker"
  labels = c("non-smoking", "current smoker", "former smoker")  # Assign descriptive labels
)

# View the resulting time-varying dataset
head(smoke_data)

# Add time-dependent hrt statuses (Ever/Never user)
hrt_data <- tmerge(
  data1 = smoke_data,
  data2 = filtered_data,
  id = ID,
  exposed_to_hrt = tdc(hrt_start)        # Change status to "ever-users" (status = 1) at the age of initiation (hrt_start)
)

head(hrt_data)

# Add a time-dependent hrt.status column
hrt_data$hrt.status <- factor(
  ifelse(hrt_data$exposed_to_hrt == 1, "ever", "never"),
  levels = c("never", "ever")  # Set "never" as the reference
)

# Inspect the first rows to verify the output
head(hrt_data)

# Add time-dependent meno status (Had yes/no)
meno_data <- tmerge(
  data1 = hrt_data,
  data2 = filtered_data,
  id = ID,
  meno_stat = tdc(age_meno)        # Change status to "yes-they've had menopause" (status = 1) at the age of  (age_meno)
)

head(meno_data)

# Add a time-dependent meno_status column
meno_data$meno_status <- factor(
  ifelse(meno_data$meno_stat == 1, "yes", "no"),
  levels = c("no", "yes")  # Set "no" as the reference
)

head(meno_data)

# Add time-dependent hyster status (Had yes/no)
hyster_data <- tmerge(
  data1 = meno_data,
  data2 = filtered_data,
  id = ID,
  hyster_stat = tdc(age_hyster)        # Change status to "yes-they've had hyster" (status = 1) at the age of  (age_hyster)
)

head(hyster_data)

# Add a time-dependent hyster_status column
hyster_data$hyster_status <- factor(
  ifelse(hyster_data$hyster_stat == 1, "yes", "no"),
  levels = c("no", "yes")  # Set "no" as the reference
)

head(hyster_data)


# Initialize the timed data frame with age as the primary time scale
t_data_oc <- tmerge(
  data1 = filtered_data, 
  data2 = filtered_data, 
  id = ID, 
  tstart = 0,                       # Start time (age at birth)
  tstop = pmin(age, bc_age, na.rm = TRUE)  # Stop time (censor at BC diagnosis if applicable)
)

# Define time-varying covariates for OC exposure (initiation and discontinuation)
t_data_oc <- tmerge(
  data1 = t_data_oc,         # Original data frame with time and event information
  data2 = t_data_oc,           # Data frame with exposure data
  id = ID,                # Merge based on 'id'
  exposed_to_oc_start = tdc(start_oc)   # OC initiation time (start)
)

# Add OC discontinuation as a time-varying covariate
t_data_oc <- tmerge(
  data1 = t_data_oc, 
  data2 = t_data_oc, 
  id = ID, 
  exposed_to_oc_stop = tdc(stop_oc + 2)   # OC discontinuation time (age stopp + 2years lag)
)

# Create the exposure group column based on exposure status (during, after, never)
t_data_oc$exposure_group <- with(t_data_oc, 
                                 ifelse(exposed_to_oc_start == 1 & exposed_to_oc_stop == 0, 1,    # During exposure (initiated but not stopped)
                                        ifelse(exposed_to_oc_start == 0 & exposed_to_oc_stop == 0, 0,    # No exposure (never used OC)
                                               ifelse(exposed_to_oc_start == 1 & exposed_to_oc_stop == 1, 2,    # After exposure (initiated and stopped)
                                                      NA))))  # This case should ideally never happen

# Verify that the 'exposure_group' column is properly created
head(t_data_oc)


# Step 4: Combine OC exposure with BC data
t_data_combined <- tmerge(
  data1 = hyster_data, 
  data2 = t_data_oc, 
  id = ID, 
  exposure_status = tdc(tstart, exposure_group)  # Merge OC status
)

head(t_data_combined)

# Stratification part
# Create `genetic_risk_group`
t_data_combined$genetic_risk_group <- ifelse(t_data_combined$PRS_decile == 10, "high", "low")

# Define `fin_status`(= combination of oc exposure status and genetic status) with concise numeric labels
t_data_combined$fin_status <- with(t_data_combined, 
                                   ifelse(genetic_risk_group == "low" & exposure_status == 0, 0,  # 0: Never OC, Low Risk
                                          ifelse(genetic_risk_group == "low" & exposure_status == 1, 1,  # 1: During OC, Low Risk
                                                 ifelse(genetic_risk_group == "high" & exposure_status == 0, 2,  # 2: Never OC, High Risk
                                                        ifelse(genetic_risk_group == "high" & exposure_status == 1, 3,  # 3: During OC, High Risk
                                                               ifelse(genetic_risk_group == "low" & exposure_status == 2, 4,  # 4: After OC, Low Risk
                                                                      ifelse(genetic_risk_group == "high" & exposure_status == 2, 5, NA)))))))

# Convert `fin_status` to a factor
t_data_combined$fin_status <- as.factor(t_data_combined$fin_status)

# Filter data for "Never OC" vs. "During OC"
t_data_during <- subset(t_data_combined, fin_status %in% c(0, 1, 2, 3))

# Filter data for "Never OC" vs. "After OC"
t_data_after <- subset(t_data_combined, fin_status %in% c(0, 4, 2, 5))

# Fit the Cox model for "During OC use" comparison
cox_model_during <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ fin_status + SmokingStatus + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob, data = t_data_during)

# Summarize results
summary(cox_model_during)

# Original number of participants in during
length(which(!(unique(t_data_during$ID)%in%cox_model_during$na.action)))

# ~~Extra test~~~
### To compare if there's a significant difference between fin_status2 (reference group) and fin_status3, you only need to relevel the factor for fin_status, setting fin_status2 as the reference category.
# Relevel fin_status to set fin_status2 as the reference group
t_data_during$fin_status <- relevel(t_data_during$fin_status, ref = "2")

# Fit the Cox model for "During OC use" comparison
cox_model_during <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ fin_status + SmokingStatus + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob, data = t_data_during)

# Summarize results
summary(cox_model_during)
# ~~End of Extra test~~

# Fit the Cox model for "After OC use" comparison
cox_model_after <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ fin_status + SmokingStatus + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob, data = t_data_after)

# Summarize results
summary(cox_model_after)

# Original number of participants in after
length(which(!(unique(t_data_after$ID)%in%cox_model_after$na.action)))

# ~~Extra test~~~
### To compare if there's a significant difference between fin_status2 (reference group) and fin_status5, you only need to relevel the factor for fin_status, setting fin_status2 as the reference category.
# Relevel fin_status to set fin_status2 as the reference group
t_data_after$fin_status <- relevel(t_data_after$fin_status, ref = "2")

# Fit the Cox model for "After OC use" comparison
cox_model_after <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ fin_status + SmokingStatus + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob, data = t_data_after)

# Summarize results
summary(cox_model_after)
# ~~End of Extra test~~
