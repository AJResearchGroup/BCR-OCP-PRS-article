#Author: Christina Chatsatourian 

## Ever-Never_without_smoking
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

# Add OC exposure as a time-varying covariate with 2 possible statuses (never user and ever user)
t_data_oc <- tmerge(
  data1 = t_data_oc, 
  data2 = t_data_oc, 
  id = ID, 
  exposed_to_oc = tdc(start_oc)        # Change status to "ever-users" (status = 1) at the age of initiation (start_oc)
)


# Step 4: Combine OC exposure with BC data
t_data_combined <- tmerge(
  data1 = hyster_data, 
  data2 = t_data_oc, 
  id = ID, 
  oc_exposure_status = tdc(tstart, exposed_to_oc)  # Merge OC status
)


# Add oc_use column in the dataset
t_data_combined$oc_use <- factor(
  ifelse(t_data_combined$oc_exposure_status == 1, "ever", "never"),
  levels = c("never", "ever")  # Set "never" as the reference
)


# Fit Cox regression model
fit <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ oc_use + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob + number_children, data = t_data_combined)

summary(fit)

#~~Check for missing data~~
# Original number of participants
length(which(!(unique(t_data_combined$ID)%in%fit$na.action)))

# Number of rows excluded due to missing data
length(fit$na.action)  # This gives the number of rows removed

# Number of unique IDs excluded due to missing data
removed_ids <- unique(t_data_combined$ID[fit$na.action])
length(removed_ids)  # This gives the number of unique individuals removed

# See from which columns specifically are the missing data from
# Define the columns used in order of priority
columns_used <- c("bmi", "tdi", "AgeMenarch", "yob", "number_children")  # Ordered by priority

# Get the unique IDs that were removed from the Cox model
removed_ids <- unique(t_data_combined$ID[fit$na.action])

# Extract only the rows that correspond to these removed IDs
removed_rows <- t_data_combined[t_data_combined$ID %in% removed_ids, ]

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


# PRS filtering 
library(dplyr)

# Ensure no missing values in PRS_BC 
filtered_data <- filtered_data %>%
  filter(!is.na(prs))

# Fit the models below

fit_2 <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ oc_use + prs + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob + number_children, data = t_data_combined)

summary(fit_2)

#~~Check for missing data~~
# Original number of participants
length(which(!(unique(t_data_combined$ID)%in%fit_2$na.action)))

# Number of rows excluded due to missing data
length(fit_2$na.action)  # This gives the number of rows removed

# Number of unique IDs excluded due to missing data
removed_ids <- unique(t_data_combined$ID[fit_2$na.action])
length(removed_ids)  # This gives the number of unique individuals removed

# See from which columns specifically are the missing data from
# Define the columns used in order of priority
columns_used <- c("bmi", "tdi", "AgeMenarch", "yob", "number_children")  # Ordered by priority

# Get the unique IDs that were removed from the Cox model
removed_ids <- unique(t_data_combined$ID[fit_2$na.action])

# Extract only the rows that correspond to these removed IDs
removed_rows <- t_data_combined[t_data_combined$ID %in% removed_ids, ]

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


fit_3 <- coxph(Surv(tstart, tstop, bc_diagnosis) ~ oc_use * prs + hrt.status + meno_status + hyster_status + bmi + tdi + AgeMenarch + yob + number_children, data = t_data_combined)

summary(fit_3)

#~~Check for missing data~~
# Original number of participants
length(which(!(unique(t_data_combined$ID)%in%fit_3$na.action)))

# Number of rows excluded due to missing data
length(fit_3$na.action)  # This gives the number of rows removed

# Number of unique IDs excluded due to missing data
removed_ids <- unique(t_data_combined$ID[fit_3$na.action])
length(removed_ids)  # This gives the number of unique individuals removed

# See from which columns specifically are the missing data from
# Define the columns used in order of priority
columns_used <- c("bmi", "tdi", "AgeMenarch", "yob", "number_children")  # Ordered by priority

# Get the unique IDs that were removed from the Cox model
removed_ids <- unique(t_data_combined$ID[fit_3$na.action])

# Extract only the rows that correspond to these removed IDs
removed_rows <- t_data_combined[t_data_combined$ID %in% removed_ids, ]

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
