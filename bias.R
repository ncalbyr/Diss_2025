##### (50 simulation) DENSITY estimation BIAS of each method
library(grid)
library(gridExtra)
# 1. LT2D
bias_lt2d <- (mean_density_lt2d-35)/35*100
bias_lt2d

# 2. Double-Observer/NO movement/NO mismatch
bias_do_nm_nm <- (mean_density_none-35)/35*100
bias_do_nm_nm

chapman_none[]
# 3. Double-Observer/YES movement/NO mismatch
# (random movement)
bias_do_ym_nm <- (mean(chapman_df_ym_nm)-35)/35*100
bias_do_ym_nm

# 4. Double Observer/YES movement/YES mismatch
# (random movement)
bias_do_ym_ym <- (mean(chapman_density_ym_ym)-35)/35*100
bias_do_ym_nm

bias <- list(LT2D = bias_lt2d,
      no_movement_no_mismatch = bias_do_nm_nm,
      yes_movement_no_mismatch = bias_do_ym_nm,
      yes_movement_yes_mismatch = bias_do_ym_ym)

# Generate better visualization
png("bias.png", width = 600, height = 250)
grid.table(head(bias, 6), rows = NULL)
dev.off()
