##### (50 simulation) Estimation BIAS of each method
# 1. LT2D
bias_lt2d <- (mean(lt2d_density)-35)/35*100
# 2. Double-Observer/NO movement/NO mismatch
bias_do_nm_nm <- (mean(lt2d_density)-35)/35*100
# 3. Double-Observer/YES movement/NO mismatch
bias_do_ym_nm <- (mean(lt2d_density)-35)/35*100
# 4. Double Observer/YES movement/YES mismatch
bias_do_ym_ym <- (mean(lt2d_density)-35)/35*100

