source("util_IRT_sim.R")

# table 1
summary_results <- summarize_all_designs("../../results/figure_4/r_high_10_ell_1", seeds = 1:2000)
generate_table(summary_results, output_file = "../../tables/table_1.tex", digits = 4)

# supp table 1
summary_results <- summarize_all_designs("../../results/figure_4/r_high_10_ell_0", seeds = 1:2000)
generate_table(summary_results, output_file = "../../tables/supp_table_1.tex", digits = 4)

# supp table 2
summary_results <- summarize_all_designs("../../results/figure_4/r_high_10_ell_2", seeds = 1:2000)
generate_table(summary_results, output_file = "../../tables/supp_table_2.tex", digits = 4)