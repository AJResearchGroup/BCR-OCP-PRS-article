#Author: Christina Chatsatourian 

# PRS_stratification_model_Ever-Never_without_smoking

library(survival)

# Adjust `tstop` to avoid issues where `tstart >= tstop`
# Set `tstop` to the minimum of `age` and `bc_age` (censor at BC diagnosis if applicable)
filtered_data$tstop <- pmin(filtered_data$age, filtered_data$bc_age, na.rm = TRUE)

# Filter out rows where `tstop` is less than or equal to 0 (ensuring valid `tstop` values) (No rows were filtered out)
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

# Add the event status for BC diagnosis
t_data <- tmerge(
  data1 = t_data, 
  data2 = t_data, 
  id = ID, 
  bc_diagnosis = event(bc_age)       # Track BC event 
)


# Add time-dependent hrt statuses (Ever/Never user)
hrt_data <- tmerge(
  data1 = t_data,
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

# Add OC exposure as a time-varying covariate with 2 possible statuses (never user and ever user)
t_data_oc <- tmerge(
  data1 = t_data_oc, 
  data2 = t_data_oc, 
  id = ID, 
  exposed_to_oc = tdc(start_oc)        # Change status to "ever-users" (status = 1) at the age of initiation (start_oc)
)

head(t_data_oc)

# Step 4: Combine OC exposure with BC data
t_data_combined <- tmerge(
  data1 = hyster_data, 
  data2 = t_data_oc, 
  id = ID, 
  oc_exposure_status = tdc(tstart, exposed_to_oc)  # Merge OC status
)

head(t_data_combined)

# Stratification part

# Create genetic risk groups
t_data_combined$genetic_risk_group <- ifelse(t_data_combined$PRS_decile == 10, "high", "low")

# Define `fin_status`(= combination of oc exposure status and genetic status) with concise numeric labels
t_data_combined$fin_status <- with(t_data_combined, 
                                   ifelse(oc_exposure_status == 0 & genetic_risk_group == "low", 0,  # Never user - low risk group
                                          ifelse(oc_exposure_status == 1 & genetic_risk_group == "low", 1,  # Ever user - low risk group
                                                 ifelse(oc_exposure_status == 0 & genetic_risk_group == "high", 2,  # Never user - high risk group
                                                        3))))  # Ever user - high risk group

# Convert `fin_status` to a factor for Cox regression
t_data_combined$fin_status <- as.factor(t_data_combined$fin_status)

# Fit the Cox regression model for ever vs never
fit <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ fin_status + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob + number_children, data = t_data_combined)

# Summarize results
summary(fit)

# Original number of participants
length(which(!(unique(t_data_combined$ID)%in%fit$na.action)))

# ~~Extra test~~~
# Relevel fin_status to set fin_status2 as the reference group
t_data_combined$fin_status <- relevel(t_data_combined$fin_status, ref = "2")

# Fit the Cox regression model for ever vs never
fit <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ fin_status + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob + number_children, data = t_data_combined)

# Summarize results
summary(fit)
# ~~End of Extra test~~
