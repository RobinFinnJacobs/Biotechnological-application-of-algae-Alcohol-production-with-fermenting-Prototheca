# Statistics 
The codes for the statistical analysis are provided here

## 01 ANOVA fermentation experiment 
```r
#Import CSV file with semicolon as separator
data <- read.csv("~/Desktop/Tabelle_Anova_fermentation.csv", sep = ";", stringsAsFactors = FALSE)

#Replace commas with full stops and convert relevant columns into numerical variables
data$pH_first_day <- as.numeric(sub(",", ".", data$pH_first_day))
data$pH_last_day <- as.numeric(sub(",", ".", data$pH_last_day))
data$pH_diff <- as.numeric(sub(",", ".", data$pH_diff))
data$Vol_proz_last_day <- as.numeric(sub(",", ".", data$Vol_proz_last_day))
data$Vol_proz_diff <- as.numeric(sub(",", ".", data$Vol_proz_diff))

#Convert the column ‘Approach’as a factor
data$Approach <- as.factor(data$Approach)

#ANOVA for the difference in pH value (pH_diff)
anova_pH <- aov(pH_diff ~ Approach, data = data)
summary(anova_pH)

#ANOVA for the difference in volume (Vol_proz_diff)
anova_Vol <- aov(Vol_proz_diff ~ Approach, data = data)
summary(anova_Vol)

# Tukey post hoc test for pH value differences
tukey_pH <- TukeyHSD(anova_pH)
print(tukey_pH)

# Tukey post hoc test for volume percentage differences
tukey_Vol <- TukeyHSD(anova_Vol)
print(tukey_Vol)
```
## 02 ANOVA alcohol tolerance experiment 
```r
#Import CSV file 
data <- read.csv("~/Desktop/Alkoholtoleranz.csv", sep = ";", dec = ",")

#Change data to factors 
data$Organismus <- as.factor(data$Organismus)
data$Alkkonz <- as.factor(data$Alkkonz)
data$Zeit..h. <- as.factor(data$Zeit..h.)

#ANOVA with three-factorial analysis
anova_alktoleranzdreifaktor <- aov(Coverage.. ~ Organismus * Alkkonz * Zeit..h., data = data)
summary(anova_alktoleranzdreifaktor)

#calculate residuals 
residuals <- residuals(anova_alktoleranzdreifaktor)

#Shapiro-Wilk test for normal distribution of residuals  
shapiro_test <- shapiro.test(residuals)
print(shapiro_test)

#Levene test for the homogeneity of variances
install.packages("car")  
library(car)
  
levene_test <- leveneTest(Coverage.. ~ Organismus * Alkkonz * Zeit..h., data = data)
print (levene_test)

# Tukey-HSD-Test 
tukey_test <- TukeyHSD(anova_alktoleranzdreifaktor)
print(tukey_test)
```
