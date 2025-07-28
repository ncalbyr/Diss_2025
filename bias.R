##### (50 simulation) Estimation BIAS of each method
# 1. LT2D
bias_lt2d <- (mean_density_lt2d-35)/35*100
bias_lt2d

# 2. Double-Observer/NO movement/NO mismatch
bias_do_nm_nm <- (mean(lt2d_density)-35)/35*100
bias_do_nm_nm

# 3. Double-Observer/YES movement/NO mismatch
bias_do_ym_nm <- (mean(chapman_density_ym_nm)-35)/35*100
bias_do_ym_nm

# 4. Double Observer/YES movement/YES mismatch
bias_do_ym_ym <- (mean(chapman_density_ym_ym)-35)/35*100
bias_do_ym_nm
