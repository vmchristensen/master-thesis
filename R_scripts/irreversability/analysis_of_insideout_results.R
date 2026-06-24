library(tidyverse)
library(ggplot2)

# Load data
irr <- read.csv('INSIDEOUT_features_final_tau1_m40_leftright.csv')
demo <- read.csv('demographics_cleaned_filtered.csv')

# Clean subject ID format to match
irr$SubjectID <- gsub('sub-', '', irr$SubjectID)

# Merge — keeping your existing group columns
df <- merge(irr, demo[, c('ID', 'LI_GEN', 'HEMISPHERIC','HEMISPHERIC_2')], 
            by.x = 'SubjectID', by.y = 'ID')

# Compute IrrLI
df$IrrLI <- (df$IrrLeft - df$IrrRight) / (df$IrrLeft + df$IrrRight)

# ---- ANALYSIS 1: Correlation IrrLI vs LI_GEN ----
cor_result <- cor.test(df$IrrLI, df$LI_GEN, method = 'pearson')
print(cor_result)

# ---- ANALYSIS 2: Scatterplot ----
ggplot(df, aes(x = IrrLI, y = LI_GEN)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', se = TRUE, color = 'blue') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray') +
  labs(title = 'IrrLI vs LI_GEN',
       x = 'Irreversibility Lateralization Index',
       y = 'Language Lateralization Index (LI_GEN)') +
  theme_minimal()

# ---- ANALYSIS 3: IrrLeft and IrrRight separately vs LI_GEN ----
cor_left  <- cor.test(df$IrrLeft,  df$LI_GEN, method = 'pearson')
cor_right <- cor.test(df$IrrRight, df$LI_GEN, method = 'pearson')
cat('IrrLeft  vs LI_GEN: r =', round(cor_left$estimate, 3),
    'p =', round(cor_left$p.value, 4), '\n')
cat('IrrRight vs LI_GEN: r =', round(cor_right$estimate, 3),
    'p =', round(cor_right$p.value, 4), '\n')

# Remove subjects with empty/NA hemispheric values
df_clean <- df %>% 
  filter(HEMISPHERIC %in% c('Atypical', 'Typical')) %>%
  filter(HEMISPHERIC_2 %in% c('Left', 'Right', 'Bilateral'))

# Check counts
table(df_clean$HEMISPHERIC)
table(df_clean$HEMISPHERIC_2)

# ---- ANALYSIS 4: IrrLI typical vs atypical ----
ggplot(df_clean, aes(x = HEMISPHERIC, y = IrrLI, fill = HEMISPHERIC)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(title = 'IrrLI by language lateralization (typical vs atypical)',
       x = '', y = 'Irreversibility LI') +
  theme_minimal() +
  theme(legend.position = 'none')

wilcox.test(IrrLI ~ HEMISPHERIC, data = df_clean)

# ---- ANALYSIS 5: IrrLI by left right bilateral ----
ggplot(df_clean, aes(x = HEMISPHERIC_2, y = IrrLI, fill = HEMISPHERIC_2)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(title = 'IrrLI by language lateralization group',
       x = '', y = 'Irreversibility LI') +
  theme_minimal() +
  theme(legend.position = 'none')

kruskal.test(IrrLI ~ HEMISPHERIC_2, data = df_clean)
pairwise.wilcox.test(df_clean$IrrLI, df_clean$HEMISPHERIC_2, 
                     p.adjust.method = 'fdr')

# ---- ANALYSIS 6: IrrLeft and IrrRight by group ----
df_long <- df_clean %>%
  pivot_longer(cols = c(IrrLeft, IrrRight),
               names_to = 'Hemisphere',
               values_to = 'Irreversibility')

ggplot(df_long, aes(x = HEMISPHERIC_2, y = Irreversibility, fill = Hemisphere)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = 'IrrLeft and IrrRight by language lateralization group',
       x = '', y = 'Irreversibility') +
  theme_minimal()

# ---- ANALYSIS 7: Contingency table ----
df_clean$Irr_group <- ifelse(df_clean$IrrLI >  0.05, 'IrrLeft',
                             ifelse(df_clean$IrrLI < -0.05, 'IrrRight', 'IrrBilateral'))

cont_table <- table(df_clean$HEMISPHERIC_2, df_clean$Irr_group)
print(cont_table)
chisq.test(cont_table)

# ---- ANALYSIS 8: Mean IrrLI per group ----
df_clean %>%
  group_by(HEMISPHERIC_2) %>%
  summarise(
    n        = n(),
    mean_IrrLI  = mean(IrrLI, na.rm = TRUE),
    sd_IrrLI    = sd(IrrLI,   na.rm = TRUE),
    mean_IrrLeft  = mean(IrrLeft,  na.rm = TRUE),
    mean_IrrRight = mean(IrrRight, na.rm = TRUE)
  )