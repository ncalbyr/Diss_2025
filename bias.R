##### (50 simulation) DENSITY estimation BIAS of each method
library(grid)
library(gridExtra)
# 1. LT2D
bias_lt2d <- (mean_density_lt2d-35)/35*100
bias_lt2d

# 2. Double-Observer/NO mismatch/NO movement
bias_do_nm_none <- (mean_density_nm_none-35)/35*100
bias_do_nm_none

# 3. Double-Observer/NO mismatch/Avoidant movement
bias_do_nm_avoid <- (mean_density_nm_avoid-35)/35*100
bias_do_nm_avoid

# 4. Double-Observer/NO mismatch/Random movement
bias_do_nm_random <- (mean_density_nm_random-35)/35*100
bias_do_nm_random

# 5. Double-Observer/YES mismatch/No movement
bias_do_ym_none <- (mean_density_ym_none-35)/35*100
bias_do_ym_none

# 6. Double-Observer/YES mismatch/Avoidant movement
bias_do_ym_avoid <- (mean_density_ym_avoid-35)/35*100
bias_do_ym_avoid

# 7. Double-Observer/YES mismatch/No movement
bias_do_ym_random <- (mean_density_ym_random-35)/35*100
bias_do_ym_random


# 7 labels for 7 methods
bias_df <- data.frame(
  Method = c(
    "LT2D",
    "Double-Observer / No mismatch / No movement",
    "Double-Observer / No mismatch / Avoidant movement",
    "Double-Observer / No mismatch / Random movement",
    "Double-Observer / Yes mismatch / No movement",
    "Double-Observer / Yes mismatch / Avoidant movement",
    "Double-Observer / Yes mismatch / Random movement"
    ),
  Percent_Bias = c(
    bias_lt2d,
    bias_do_nm_none,
    bias_do_nm_avoid,
    bias_do_nm_random,
    bias_do_ym_none,
    bias_do_ym_avoid,
    bias_do_ym_random
  ))

# Generate better visualization
png("bias.png", width = 800, height = 300)
grid.table(bias_df)
dev.off()
