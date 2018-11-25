# Load data 
data <- read.csv2('C:\\Users\\aleksander.molak\\Documents\\Personal\\ADEM\\KOB-MES_pooled data v eng.csv')

# Import libs
library(dplyr)
library(psych)
library(ggplot2)
library(rela)

# Filter
data <- data %>% filter(A_STUDY == 6)

# Select vars
data <- data %>% select(A_P5_C, starts_with('WA'))
data <- data %>% select(-starts_with('WAB'))

# Get labels
new_names <- c('A_P5_C', 'ambitious', 'warm', 'delicate', 'humorous', 'loving', 'competent', 'pretty', 'mild', 'responsible', 'brave', 'caring', 'industrious', 'handsome', 'strong', 'firm', 'subtle', 'emotional', 'cute', 'sensitive', 'tall', 'athletic')
names(data) <- new_names

# Pick woman types
data_grp_2 <- data %>% filter(A_P5_C == 2)
data_grp_3 <- data %>% filter(A_P5_C == 3)

# Drop ID col
data_grp_2 <- data_grp_2 %>% select(-A_P5_C)
data_grp_3 <- data_grp_3 %>% select(-A_P5_C)

# Produce correlation matrices
cor_2 <- data.frame(cor(data_grp_2))
cor_3 <- data.frame(cor(data_grp_3))

# Perform KMO
paf(as.matrix(data_grp_2))$KMO
paf(as.matrix(data_grp_3))$KMO

# Get scree plots
scree(data_grp_2)
scree(data_grp_3)

# Get solutions
solution_2 <- fa(r = cor_2, nfactors = 3, rotate = "varimax", fm = "pa") 
solution_3 <- fa(r = cor_3, nfactors = 3, rotate = "varimax", fm = "pa")

# Get Pattern Matrices
pttrn_2 <- data.frame(unclass(solution_2$loadings))
pttrn_3 <- data.frame(unclass(solution_3$loadings))

# Get max loadings for both solutions
max_load_2 <- apply(pttrn_2, 1, FUN=max)
max_load_3 <- apply(pttrn_3, 1, FUN=max)
max_v2_vec <- as.matrix(max_load_2)
max_v3_vec <- as.matrix(max_load_3)

# Get congruence coef
fa.congruence(max_v2_vec, max_v3_vec)

# Define congruence coefficient function
congr_coef <- function(x1, x2) {
  nmrtr = sum(x1*x2)
  dnmntr = sqrt(sum(x1**2) * sum(x2**2))
  nmrtr / dnmntr
}

# Get congruence coef
congr_coef(max_v2_vec, max_v3_vec)
