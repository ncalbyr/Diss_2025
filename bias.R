##### (50 simulation) DENSITY estimation BIAS of each method
library(grid)
library(gridExtra)
# 1. LT2D
bias_lt2d <- (mean_density_lt2d-35)/35*100
bias_lt2d

# 2. Double-Observer/NO mismatch/NO movement
bias_do_nm_none <- (mean_density_none-35)/35*100
bias_do_nm_none

# 3. Double-Observer/NO mismatch/Avoidant movement
bias_do_

# 3. Double-Observer/YES movement/NO mismatch
# (random movement)
bias_do_ym_nm <- (mean_density_ym_nm-35)/35*100
bias_do_ym_nm

# 4. Double Observer/YES movement/YES mismatch
# (random movement)
bias_do_ym_ym <- (mean_density_ym_ym-35)/35*100
bias_do_ym_nm

bias_df <- data.frame(
  Method = c(
    "LT2D",
    "Double-Observer / NO movement / NO mismatch",
    "Double-Observer / YES movement / NO mismatch",
    "Double-Observer / YES movement / YES mismatch"
  ),
  Percent_Bias = c(
    bias_lt2d,
    bias_do_nm_nm,
    bias_do_ym_nm,
    bias_do_ym_ym
  ))

# Generate better visualization
png("bias.png", width = 800, height = 300)
grid.table(bias_df)
dev.off()
