library(foreign)

# Read-in the data
filename <- 'IUIO2004-2007_eng.sav'
path <- getwd()
data <- read.spss(paste0(path, '\\', filename), to.data.frame=TRUE, use.value.labels=F)

# Save values encodings
values <- c('The sentence is incorrect, not true', 'The sentence is rather irrelevant', 'You have no opinion on this matter, or you are indifferent', 'The sentence is rather accurate', 'The sentence is very accurate, true')

# Import libs
library(dplyr)
library(psych)
library(ggplot2)
library(rela)

# Examine the data
summary(data)

# Transform year
get_year <- function(year) {
  if (nchar(as.character(year)) == 2) {
    return(as.numeric(paste0('19', as.character(year))))
  } else {
    return(year)
  }
}

data$p70a <- sapply(data$p70a, get_year)

# Filter the data
data <- data %>% filter_at(vars(-nr, -p69, -p70a, -DATA), all_vars(.<8))
data <- data %>% filter_at(vars(-nr, -p69, -p70a, -DATA), all_vars(.>0))

# Assume 'hard to say' == 'no opinion'
data[data == 7] = 3

# Drop zero variance column and order labels
data <- data %>% select(-IUIO34)
attributes(data)$variable.labels <- attributes(data)$variable.labels[-35]

# Export data
write.csv(data, paste0(path, '\\', 'data_clean_HMWRK2.csv'))

# Re-examine the data
summary(data)

# Drop ID and demographic columns
data <- data %>% select(-nr, -p69, -p70a)
attributes(data)$variable.labels <- attributes(data)$variable.labels[-c(1, 36, 35)]

# Iteration 2: Drop variables with low communalities
to_drop <- c('IUIO2', 'IUIO3', 'IUIO4', 'IUIO7', 'IUIO8', 'IUIO10', 'IUIO11', 'IUIO12', 'IUIO16', 
             'IUIO18', 'IUIO22', 'IUIO23', 'IUIO24', 'IUIO26', 'IUIO28', 'IUIO33', 'IUIO14',
             'IUIO20', 'IUIO32', 'IUIO27', 'IUIO30', 'IUIO19', 'IUIO31', 'IUIO29')
data <- data %>% select(-to_drop)
attributes(data)$variable.labels <- attributes(data)$variable.labels[-c(2, 3, 4, 7, 8, 10,
                                                                        11, 12, 16, 18, 22, 23, 
                                                                        24, 26, 28, 33, 14, 20, 32, 
                                                                        27, 30, 19, 31, 29)]

# Split studies
study_1 <- data %>% filter(DATA == 1)
study_2 <- data %>% filter(DATA == 2)
study_3 <- data %>% filter(DATA == 3)

# Drop study_id column
study_1 <- study_1 %>% select(-DATA)
attributes(study_1)$variable.labels <- attributes(study_1)$variable.labels[-c(13)]

study_2 <- study_2 %>% select(-DATA)
attributes(study_2)$variable.labels <- attributes(study_2)$variable.labels[-c(13)]

study_3 <- study_3 %>% select(-DATA)
attributes(study_3)$variable.labels <- attributes(study_3)$variable.labels[-c(13)]

# Produce correlation matrices
cor_1 <- data.frame(cor(study_1))
cor_2 <- data.frame(cor(study_2))
cor_3 <- data.frame(cor(study_3))

# Perform KMO
paf(as.matrix(study_1))$KMO
paf(as.matrix(study_2))$KMO
paf(as.matrix(study_3))$KMO

# Perform Bartlett
cortest.bartlett(study_1)
cortest.bartlett(study_2)
cortest.bartlett(study_3)

# Check MSA values
paf(as.matrix(study_1))$MSA
paf(as.matrix(study_2))$MSA
paf(as.matrix(study_3))$MSA

# Get scree plots
scree(study_1)
scree(study_2)
scree(study_3)

# Get parallel plots
fa.parallel(study_1, fm = 'pa', fa = 'fa', n.iter = 50, SMC = TRUE, main='Study 1')
fa.parallel(study_2, fm = 'pa', fa = 'fa', n.iter = 50, SMC = TRUE, main='Study 2')
fa.parallel(study_3, fm = 'pa', fa = 'fa', n.iter = 50, SMC = TRUE, main='Study 3')

# Get solutions
solution_1 <- fa(r = cor_1, nfactors = 3, rotate = "oblimin", fm = "pa", SMC=F) 
solution_2 <- fa(r = cor_2, nfactors = 3, rotate = "oblimin", fm = "pa", SMC=F)
solution_3 <- fa(r = cor_3, nfactors = 3, rotate = "oblimin", fm = "pa", SMC=F)

# Get structure matrices
str_1 = unclass(solution_1$Structure)
str_2 = unclass(solution_2$Structure)
str_3 = unclass(solution_3$Structure)

round(str_1, 3)
round(str_2, 3)
round(str_3, 3)

# Check communalities
comm_1 = solution_1$communality
comm_2 = solution_2$communality
comm_3 = solution_3$communality

comm_1 <- data.frame(comm_1)
comm_2 <- data.frame(comm_2)
comm_3 <- data.frame(comm_3)

comms <- data.frame(comm_1, comm_2, comm_3)

comms

# Define congruence coefficient function
congr_coef <- function(x1, x2) {
  nmrtr = sum(x1*x2)
  dnmntr = sqrt(sum(x1**2) * sum(x2**2))
  nmrtr / dnmntr
  }

# Compute congruence factors
congr_coef(str_1[,1], str_2[,1])
congr_coef(str_1[,2], str_2[,2])
congr_coef(str_1[,3], str_2[,3])
congr_coef(str_1[,1], str_2[,2])
congr_coef(str_1[,1], str_2[,3])
congr_coef(str_1[,2], str_2[,1])
congr_coef(str_1[,2], str_2[,3])
congr_coef(str_1[,3], str_2[,1])
congr_coef(str_1[,3], str_2[,2])

congr_coef(str_1[,1], str_3[,1])
congr_coef(str_1[,2], str_3[,2])
congr_coef(str_1[,3], str_3[,3])
congr_coef(str_1[,1], str_3[,2])
congr_coef(str_1[,1], str_3[,3])
congr_coef(str_1[,2], str_3[,1])
congr_coef(str_1[,2], str_3[,3])
congr_coef(str_1[,3], str_3[,1])
congr_coef(str_1[,3], str_3[,2])

congr_coef(str_2[,1], str_3[,1])
congr_coef(str_2[,2], str_3[,2])
congr_coef(str_2[,3], str_3[,3])
congr_coef(str_2[,1], str_3[,2])
congr_coef(str_2[,1], str_3[,3])
congr_coef(str_2[,2], str_3[,1])
congr_coef(str_2[,2], str_3[,3])
congr_coef(str_2[,3], str_3[,1])
congr_coef(str_2[,3], str_3[,2])
