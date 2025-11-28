#Author: Christina Chatsatourian 

# During-After_without_smoking
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



# Add a time-dependent hrt.status column
hrt_data$hrt.status <- factor(
  ifelse(hrt_data$exposed_to_hrt == 1, "ever", "never"),
  levels = c("never", "ever")  # Set "never" as the reference
)


# Add time-dependent meno status (Had yes/no)
meno_data <- tmerge(
  data1 = hrt_data,
  data2 = filtered_data,
  id = ID,
  meno_stat = tdc(age_meno)        # Change status to "yes-they've had menopause" (status = 1) at the age of  (age_meno)
)


# Add a time-dependent meno_status column
meno_data$meno_status <- factor(
  ifelse(meno_data$meno_stat == 1, "yes", "no"),
  levels = c("no", "yes")  # Set "no" as the reference
)


# Add time-dependent hyster status (Had yes/no)
hyster_data <- tmerge(
  data1 = meno_data,
  data2 = filtered_data,
  id = ID,
  hyster_stat = tdc(age_hyster)        # Change status to "yes-they've had hyster" (status = 1) at the age of  (age_hyster)
)


# Add a time-dependent hyster_status column
hyster_data$hyster_status <- factor(
  ifelse(hyster_data$hyster_stat == 1, "yes", "no"),
  levels = c("no", "yes")  # Set "no" as the reference
)


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
  exposed_to_oc_start = tdc(start_oc)   # OC initiation time (start_oc)
)

# Add OC discontinuation as a time-varying covariate
t_data_oc <- tmerge(
  data1 = t_data_oc, 
  data2 = t_data_oc, 
  id = ID, 
  exposed_to_oc_stop = tdc(stop_oc + 2)   # OC discontinuation time (stop_oc age + 2years lag period)
)

# Create the exposure group column based on exposure status (during, after, never)
t_data_oc$exposure_group <- with(t_data_oc, 
                                 ifelse(exposed_to_oc_start == 1 & exposed_to_oc_stop == 0, 1,    # During exposure (initiated but not stopped)
                                        ifelse(exposed_to_oc_start == 0 & exposed_to_oc_stop == 0, 0,    # No exposure (never used OC)
                                               ifelse(exposed_to_oc_start == 1 & exposed_to_oc_stop == 1, 2,    # After exposure (initiated and stopped)
                                                      NA))))  # This case should ideally never happen



# Combine OC exposure with BC data
t_data_combined <- tmerge(
  data1 = hyster_data, 
  data2 = t_data_oc, 
  id = ID, 
  exposure_status = tdc(tstart, exposure_group)  # Merge OC status
)


## Define the reference and target groups for Cox analysis
# The reference group is the "never exposed" group (exposure_group == 0)
# The target groups are the "during exposure" (exposure_group == 1) and "after exposure" (exposure_group == 2)

# Create a variable to denote group categories: 0 = never exposed, 1 = during exposure, 2 = after exposure
t_data_combined$oc_group <- with(t_data_combined,
                                 ifelse(exposure_status == 0, "Unexposed",   # Never exposed (reference group)
                                        ifelse(exposure_status == 1, "During Exposure",  # During exposure
                                               ifelse(exposure_status == 2, "After Exposure", NA))))

# Filter data to only include "Unexposed" and "During Exposure" groups
t_data_during <- subset(t_data_combined, oc_group %in% c("Unexposed", "During Exposure"))

# Define the reference and target groups for during
t_data_during$oc_group <- as.factor(t_data_during$oc_group)

t_data_during$oc_group <- relevel(t_data_during$oc_group, ref = "Unexposed")

levels(t_data_during$oc_group)

# Filter data to only include "Unexposed" and "After Exposure" groups
t_data_after <- subset(t_data_combined, oc_group %in% c("Unexposed", "After Exposure"))

# Define the reference and target groups for after
t_data_after$oc_group <- as.factor(t_data_after$oc_group)

t_data_after$oc_group <- relevel(t_data_after$oc_group, ref = "Unexposed")

levels(t_data_after$oc_group)


# Fit the Cox model comparing "During Exposure" vs "Unexposed"
cox_model_during <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ oc_group + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob, data = t_data_during)

summary(cox_model_during)

#~~Check for missing data~~
# Original number of participants
length(which(!(unique(t_data_combined$ID)%in%cox_model_during$na.action)))

# Original number of participants in during
length(which(!(unique(t_data_during$ID)%in%cox_model_during$na.action)))

# Number of rows excluded due to missing data
length(cox_model_during$na.action)  # This gives the number of rows removed

# Number of unique IDs excluded due to missing data
removed_ids <- unique(t_data_during$ID[cox_model_during$na.action])
length(removed_ids)  # This gives the number of unique individuals removed

# See from which columns specifically are the missing data from
# Define the columns used in order of priority
columns_used <- c("bmi", "tdi", "AgeMenarch", "yob")  # Ordered by priority

# Get the unique IDs that were removed from the Cox model
removed_ids <- unique(t_data_during$ID[cox_model_during$na.action])

# Extract only the rows that correspond to these removed IDs
removed_rows <- t_data_during[t_data_during$ID %in% removed_ids, ]

# Initialize an empty dataframe to track which variable caused exclusion
id_assignment <- data.frame(ID = removed_ids, AssignedVar = NA, stringsAsFactors = FALSE)

# Loop through each variable in priority order and assign the first missing variable
for (var in columns_used) {
  # Find IDs with missing values in the current variable and that haven't been assigned yet
  ids_missing_var <- unique(removed_rows$ID[is.na(removed_rows[[var]])])
  
  # Assign only if the ID has not already been assigned a missing variable
  id_assignment$AssignedVar[id_assignment$ID %in% ids_missing_var & is.na(id_assignment$AssignedVar)] <- var
}

# Count the number of unique IDs assigned to each missing variable
missing_by_variable <- table(id_assignment$AssignedVar, useNA = "no")  # Exclude NAs (not assigned)

# Print results
cat("Total number of unique IDs excluded:", length(removed_ids), "\n")
cat("Number of unique IDs excluded due to specific missing variables:\n")
print(missing_by_variable)
#~~end of check~~

# Fit the Cox model comparing "After Exposure" vs "Unexposed"
cox_model_after <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ oc_group + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob, data = t_data_after)

summary(cox_model_after)

#~~Check for missing data~~
# Original number of participants
length(which(!(unique(t_data_combined$ID)%in%cox_model_after$na.action)))

# Original number of participants in during
length(which(!(unique(t_data_after$ID)%in%cox_model_after$na.action)))

# Number of rows excluded due to missing data
length(cox_model_after$na.action)  # This gives the number of rows removed

# Number of unique IDs excluded due to missing data
removed_ids <- unique(t_data_after$ID[cox_model_after$na.action])
length(removed_ids)  # This gives the number of unique individuals removed

# See from which columns specifically are the missing data from
# Define the columns used in order of priority
columns_used <- c("bmi", "tdi", "AgeMenarch", "yob")  # Ordered by priority

# Get the unique IDs that were removed from the Cox model
removed_ids <- unique(t_data_after$ID[cox_model_after$na.action])

# Extract only the rows that correspond to these removed IDs
removed_rows <- t_data_after[t_data_after$ID %in% removed_ids, ]

# Initialize an empty dataframe to track which variable caused exclusion
id_assignment <- data.frame(ID = removed_ids, AssignedVar = NA, stringsAsFactors = FALSE)

# Loop through each variable in priority order and assign the first missing variable
for (var in columns_used) {
  # Find IDs with missing values in the current variable and that haven't been assigned yet
  ids_missing_var <- unique(removed_rows$ID[is.na(removed_rows[[var]])])
  
  # Assign only if the ID has not already been assigned a missing variable
  id_assignment$AssignedVar[id_assignment$ID %in% ids_missing_var & is.na(id_assignment$AssignedVar)] <- var
}

# Count the number of unique IDs assigned to each missing variable
missing_by_variable <- table(id_assignment$AssignedVar, useNA = "no")  # Exclude NAs (not assigned)

# Print results
cat("Total number of unique IDs excluded:", length(removed_ids), "\n")
cat("Number of unique IDs excluded due to specific missing variables:\n")
print(missing_by_variable)
#~~end of check~~

# 1) PRS filtering 
# Load necessary libraries
library(dplyr)

# Ensure no missing values in PRS_BC 
filtered_data <- filtered_data %>%
  filter(!is.na(prs))


# Fit the models below

cox_model_during_2 <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ oc_group + prs + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob, data = t_data_during)

summary(cox_model_during_2)

#~~Check for missing data~~
# Original number of participants
length(which(!(unique(t_data_combined$ID)%in%cox_model_during_2$na.action)))

# Original number of participants in during
length(which(!(unique(t_data_during$ID)%in%cox_model_during_2$na.action)))

# Number of rows excluded due to missing data
length(cox_model_during_2$na.action)  # This gives the number of rows removed

# Number of unique IDs excluded due to missing data
removed_ids <- unique(t_data_during$ID[cox_model_during_2$na.action])
length(removed_ids)  # This gives the number of unique individuals removed

# See from which columns specifically are the missing data from
# Define the columns used in order of priority
columns_used <- c("bmi", "tdi", "AgeMenarch", "yob")  # Ordered by priority

# Get the unique IDs that were removed from the Cox model
removed_ids <- unique(t_data_during$ID[cox_model_during_2$na.action])

# Extract only the rows that correspond to these removed IDs
removed_rows <- t_data_during[t_data_during$ID %in% removed_ids, ]

# Initialize an empty dataframe to track which variable caused exclusion
id_assignment <- data.frame(ID = removed_ids, AssignedVar = NA, stringsAsFactors = FALSE)

# Loop through each variable in priority order and assign the first missing variable
for (var in columns_used) {
  # Find IDs with missing values in the current variable and that haven't been assigned yet
  ids_missing_var <- unique(removed_rows$ID[is.na(removed_rows[[var]])])
  
  # Assign only if the ID has not already been assigned a missing variable
  id_assignment$AssignedVar[id_assignment$ID %in% ids_missing_var & is.na(id_assignment$AssignedVar)] <- var
}

# Count the number of unique IDs assigned to each missing variable
missing_by_variable <- table(id_assignment$AssignedVar, useNA = "no")  # Exclude NAs (not assigned)

# Print results
cat("Total number of unique IDs excluded:", length(removed_ids), "\n")
cat("Number of unique IDs excluded due to specific missing variables:\n")
print(missing_by_variable)
#~~end of check~~

cox_model_after_2 <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ oc_group + prs + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob, data = t_data_after)

summary(cox_model_after_2)

#~~Check for missing data~~
# Original number of participants
length(which(!(unique(t_data_combined$ID)%in%cox_model_after_2$na.action)))

# Original number of participants in during
length(which(!(unique(t_data_after$ID)%in%cox_model_after_2$na.action)))

# Number of rows excluded due to missing data
length(cox_model_after_2$na.action)  # This gives the number of rows removed

# Number of unique IDs excluded due to missing data
removed_ids <- unique(t_data_after$ID[cox_model_after_2$na.action])
length(removed_ids)  # This gives the number of unique individuals removed

# See from which columns specifically are the missing data from
# Define the columns used in order of priority
columns_used <- c("bmi", "tdi", "AgeMenarch", "yob")  # Ordered by priority

# Get the unique IDs that were removed from the Cox model
removed_ids <- unique(t_data_after$ID[cox_model_after_2$na.action])

# Extract only the rows that correspond to these removed IDs
removed_rows <- t_data_after[t_data_after$ID %in% removed_ids, ]

# Initialize an empty dataframe to track which variable caused exclusion
id_assignment <- data.frame(ID = removed_ids, AssignedVar = NA, stringsAsFactors = FALSE)

# Loop through each variable in priority order and assign the first missing variable
for (var in columns_used) {
  # Find IDs with missing values in the current variable and that haven't been assigned yet
  ids_missing_var <- unique(removed_rows$ID[is.na(removed_rows[[var]])])
  
  # Assign only if the ID has not already been assigned a missing variable
  id_assignment$AssignedVar[id_assignment$ID %in% ids_missing_var & is.na(id_assignment$AssignedVar)] <- var
}

# Count the number of unique IDs assigned to each missing variable
missing_by_variable <- table(id_assignment$AssignedVar, useNA = "no")  # Exclude NAs (not assigned)

# Print results
cat("Total number of unique IDs excluded:", length(removed_ids), "\n")
cat("Number of unique IDs excluded due to specific missing variables:\n")
print(missing_by_variable)
#~~end of check~~

cox_model_during_3 <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ oc_group * prs + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob, data = t_data_during)

summary(cox_model_during_3)

#~~Check for missing data~~
# Original number of participants
length(which(!(unique(t_data_combined$ID)%in%cox_model_during_3$na.action)))

# Original number of participants in during
length(which(!(unique(t_data_during$ID)%in%cox_model_during_3$na.action)))

# Number of rows excluded due to missing data
length(cox_model_during_3$na.action)  # This gives the number of rows removed

# Number of unique IDs excluded due to missing data
removed_ids <- unique(t_data_during$ID[cox_model_during_3$na.action])
length(removed_ids)  # This gives the number of unique individuals removed

# See from which columns specifically are the missing data from
# Define the columns used in order of priority
columns_used <- c("bmi", "tdi", "AgeMenarch", "yob")  # Ordered by priority

# Get the unique IDs that were removed from the Cox model
removed_ids <- unique(t_data_during$ID[cox_model_during_3$na.action])

# Extract only the rows that correspond to these removed IDs
removed_rows <- t_data_during[t_data_during$ID %in% removed_ids, ]

# Initialize an empty dataframe to track which variable caused exclusion
id_assignment <- data.frame(ID = removed_ids, AssignedVar = NA, stringsAsFactors = FALSE)

# Loop through each variable in priority order and assign the first missing variable
for (var in columns_used) {
  # Find IDs with missing values in the current variable and that haven't been assigned yet
  ids_missing_var <- unique(removed_rows$ID[is.na(removed_rows[[var]])])
  
  # Assign only if the ID has not already been assigned a missing variable
  id_assignment$AssignedVar[id_assignment$ID %in% ids_missing_var & is.na(id_assignment$AssignedVar)] <- var
}

# Count the number of unique IDs assigned to each missing variable
missing_by_variable <- table(id_assignment$AssignedVar, useNA = "no")  # Exclude NAs (not assigned)

# Print results
cat("Total number of unique IDs excluded:", length(removed_ids), "\n")
cat("Number of unique IDs excluded due to specific missing variables:\n")
print(missing_by_variable)
#~~end of check~~

cox_model_after_3 <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ oc_group * prs + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob, data = t_data_after)

summary(cox_model_after_3)

#~~Check for missing data~~
# Original number of participants
length(which(!(unique(t_data_combined$ID)%in%cox_model_after_3$na.action)))

# Original number of participants in during
length(which(!(unique(t_data_after$ID)%in%cox_model_after_3$na.action)))

# Number of rows excluded due to missing data
length(cox_model_after_3$na.action)  # This gives the number of rows removed

# Number of unique IDs excluded due to missing data
removed_ids <- unique(t_data_after$ID[cox_model_after_3$na.action])
length(removed_ids)  # This gives the number of unique individuals removed

# See from which columns specifically are the missing data from
# Define the columns used in order of priority
columns_used <- c("bmi", "tdi", "AgeMenarch", "yob")  # Ordered by priority

# Get the unique IDs that were removed from the Cox model
removed_ids <- unique(t_data_after$ID[cox_model_after_3$na.action])

# Extract only the rows that correspond to these removed IDs
removed_rows <- t_data_after[t_data_after$ID %in% removed_ids, ]

# Initialize an empty dataframe to track which variable caused exclusion
id_assignment <- data.frame(ID = removed_ids, AssignedVar = NA, stringsAsFactors = FALSE)

# Loop through each variable in priority order and assign the first missing variable
for (var in columns_used) {
  # Find IDs with missing values in the current variable and that haven't been assigned yet
  ids_missing_var <- unique(removed_rows$ID[is.na(removed_rows[[var]])])
  
  # Assign only if the ID has not already been assigned a missing variable
  id_assignment$AssignedVar[id_assignment$ID %in% ids_missing_var & is.na(id_assignment$AssignedVar)] <- var
}

# Count the number of unique IDs assigned to each missing variable
missing_by_variable <- table(id_assignment$AssignedVar, useNA = "no")  # Exclude NAs (not assigned)

# Print results
cat("Total number of unique IDs excluded:", length(removed_ids), "\n")
cat("Number of unique IDs excluded due to specific missing variables:\n")
print(missing_by_variable)
#~~end of check~~
