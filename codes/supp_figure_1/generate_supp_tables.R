source("util_IRT_sim.R")

# table 1
summary_results <- summarize_all_designs("../../results/supp_figure_1/r_high_10_ell_0", seeds = 1:2000)
generate_table(summary_results, output_file = "../../tables/supp_table_3.tex", digits = 4)

# supp table 1
summary_results <- summarize_all_designs("../../results/supp_figure_1/r_high_10_ell_1", seeds = 1:2000)
generate_table(summary_results, output_file = "../../tables/supp_table_4.tex", digits = 4)

# supp table 2
summary_results <- summarize_all_designs("../../results/supp_figure_1/r_high_10_ell_2", seeds = 1:2000)
generate_table(summary_results, output_file = "../../tables/supp_table_5.tex", digits = 4)

# supp table 6 - ell = 0
summary_results <- summarize_all_designs("../../results/supp_table_6/unequal_ell_0", seeds = 1:2000)
generate_table(summary_results, output_file = "../../tables/supp_table_6_ell_0.tex", digits = 4)

# supp table 6 - ell = 1
summary_results <- summarize_all_designs("../../results/supp_table_6/unequal_ell_1", seeds = 1:2000)
generate_table(summary_results, output_file = "../../tables/supp_table_6_ell_1.tex", digits = 4)