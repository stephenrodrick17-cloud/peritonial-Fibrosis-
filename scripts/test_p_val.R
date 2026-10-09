f_val <- 5.4321
df_denom <- 35.6713
p_f <- 1 - pf(f_val, df1 = 1, df2 = df_denom)
t_val <- sqrt(f_val)
p_t <- 2 * pt(-abs(t_val), df = df_denom)
cat(sprintf("P from F: %.15f\nP from t: %.15f\nDiff:     %.2e\n", p_f, p_t, abs(p_f - p_t)))
