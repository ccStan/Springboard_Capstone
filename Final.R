# Read in and Examine Data Set #
protein <- read.csv("~/R_projects/Data/Capstone Project/Protein Tertiary Structure.csv")
head(protein, 20)
str(protein)
length(protein$RMSD)

# Replace Column Names for ease of typing and conform with naming standards #
protein_names_new <- c('RMSD', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9')
protein_names_new
names(protein) <- protein_names_new
head(protein)


# Create an internal reference table (for humans) for column headings and their defintions #
# Note: two typos from original were carried throught to maintain consistency #

# Create two vectors, variable and definition then cbind into a descrition table #
variable <- c('f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9')
definition <- c('Total surface area', 'Non polar exposed area',
                             'Fractional area of exposed non polar residue',
                             'Fractional area of exposed non polar part of residue',
                             'Molecular mass wieghted exposed area', 
                             'Average deviation from standard exposed area of residue',
                             'Euclidean distance', 'Secondary structure penalty',
                             'Spacial Distribution constraints(N,K, Value)')
descriptions_table <- cbind(variable, definition)
descriptions_table

# Check for missing values and visually inspect length #
protein1 <- na.omit(protein)
length(protein$RMSD)
length(protein1$RMSD)
# Both are of length 45,730 #


# Collect sample2 as per sample 1, but with fewer points, for quick, exploratory plots #
# Use sample size rather than alpha = for dense data #
sample1 <- protein[sample(1:nrow(protein), "1000", replace = FALSE), ]
sample2 <- protein[sample(1:nrow(protein), "100", replace = FALSE), ]

# Get descriptive statistics on Columns #
# Select statistics most useful for this analysis #
df.describe <- describe(protein)
df.describe
df.describe_select <- select(df.describe, mean, sd, skew, kurtosis, median, range)
df.describe_select
# Note the very high kurtosis for f7, Euclidian distance #

# Examine f7 distribution with histograms #
histogram(x = sample1$f7, breaks = 50)
histogram(x = sample2$f7, breaks = 50)
# Note outliers at 8,000 out to 24,000 A.  Small % so ignore #
# Look at f7 over the entire data set #
histogram(x = protein$f7, breaks = 100)
# Very "sharp", centered distribution with a few outliers, so ignore for now #

# START ANALYSIS #

# Plot all variable two-dimensional conbinations and display in a grid #
# Use random samples to simplify and clarify #
plot(sample1) # too dense
plot(sample2)
# Note: high degree of correlation of all "surface" variables, Euclidian distance correlations
# always in a tight band #
# Note: RMSD gives a bifurcated plot with two distinct areas; two factor interaction?! #

# Plot RMSD for further examination #
qplot(x = f1, y = RMSD, data = sample1, geom = 'point')
qplot(x = f1, y = RMSD, data = sample2, geom = 'point')

# Check for anomlies in records from Sample1 where RMSD  = 0 #
# Two ways for 0, RMSD actually is zero, missing value which leads to a zero value #
sample1_zero <- subset(sample1, RMSD == 0)
sample1_zero


# Now Check descriptive statistics #
sample1_zero_describe <- describe(sample1_zero)
sample1_zero_describe_select <- select(sample1_zero_describe, mean, sd, skew, kurtosis, median, range)
sample1_zero_describe_select
# Kurtosis is now "normal" which makes sense since RMSD = 0.  Notice NAN #

# Do same for entire data base, protein
protein_zero <- subset(protein, RMSD == 0)
head(protein_zero)
protein_zero_describe <- describe(protein_zero)
protein_zero_describe_select <- select(sample1_zero_describe, mean, sd, skew, kurtosis, median, range)
protein_zero_describe_select
# Note: same results as sample 1 #

# Get correlation matrix in preparation for regression #
cor(protein[c('RMSD', 'f1', 'f2', 'f3', 'f4', 'f5', 
              'f6', 'f7', 'f8', 'f9')])

# General Linear Regression #
# Start with factors which have higher correlation with dependent variable RMSD #
# Adjusted R-squared copied to Excel spreadsheet #
rmsd1 <- lm(RMSD ~ f3, data = protein)
summary(rmsd1)

rmsd2 <- lm(RMSD ~ f3 + f4, data = protein)
summary(rmsd2)

rmsd3 <- lm(RMSD ~ f3 + f4 + f2, data = protein)
summary(rmsd3)

rmsd4 <- lm(RMSD ~ f3 + f4 + f2 + f9, data = protein)
summary(rmsd4)


rmsd5 <- lm(RMSD ~ f4 + f2, data = protein)
summary(rmsd5)

# General Linear Regression with non-linear Terms #
rmsd6 <- lm(RMSD ~ f2^2, data = protein)
summary(rmsd6)

rmsd7 <- lm(RMSD ~ log(f3), data = protein)
summary(rmsd7)

rmsd8 <- lm(RMSD ~ f4 + f2 + f4*f2, data = protein)
summary(rmsd8)

rmsd9 <- lm(RMSD ~ f3 + f3*f2 + f3*f4, data = protein)
summary(rmsd9)

rmsd10 <- lm(RMSD ~ f3 + f7 + f3*f7, data = protein)
summary(rmsd10)
