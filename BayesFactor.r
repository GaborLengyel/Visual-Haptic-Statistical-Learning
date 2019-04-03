# Loading the data
data = read.csv("data_for_R.csv", header=FALSE)
colnames(data) <- c("exp1_vis", "exp2_vis","exp1_hap","exp2_hap","exp1_con","exp2_con","raw1","raw2","contrl1","contrl2")
fprintf <- function(...) cat(sprintf(...))

# Importing bayesfactor package
library(BayesFactor)
rscale = 0.707 # setting the prior

# bayes factors for the main results of experiment 1
t_v = t.test(na.omit(data$exp1_vis),mu=0.5)
t_h = t.test(na.omit(data$exp1_hap))
bf_v = ttestBF(na.omit(data$exp1_vis),mu=0.5,r=rscale)
bf_h = ttestBF(na.omit(data$exp1_hap),r=rscale)
print("Main results for experiment 1:")
fprintf('Haptic pulling performance: Rho= %2.2f (%2.2f-%2.2f)  with t(%i)=%2.2f  p=%e, Bayes factor = %2.2f\n', mean(na.omit(data$exp1_hap)), t_h$conf.int[1],t_h$conf.int[2], t_h$parameter, t_h$statistic, t_h$p.value, as.vector(bf_h))
fprintf('Visual familiarity performance: fraction correct= %2.2f (%2.2f-%2.2f)  with t(%i)=%2.2f p=%e, Bayes factor = %2.2f\n', mean(na.omit(data$exp1_vis)), t_v$conf.int[1],t_v$conf.int[2], t_v$parameter, t_v$statistic, t_v$p.value, as.vector(bf_v))

# bayes factors for the main results of experiment 2
t_v = t.test(na.omit(data$exp2_vis),mu=0.5)
t_h = t.test(na.omit(data$exp2_hap))
bf_v = ttestBF(na.omit(data$exp2_vis),mu=0.5,r=rscale)
bf_h = ttestBF(na.omit(data$exp2_hap),r=rscale)
print("Main results for experiment 2:")
fprintf('Haptic pulling performance: Rho= %2.2f (%2.2f-%2.2f)  with t(%i)=%2.2f  p=%e, Bayes factor = %2.2f\n', mean(na.omit(data$exp2_hap)), t_h$conf.int[1],t_h$conf.int[2], t_h$parameter, t_h$statistic, t_h$p.value, as.vector(bf_h))
fprintf('Visual familiarity performance: fraction correct= %2.2f (%2.2f-%2.2f)  with t(%i)=%2.2f p=%e, Bayes factor = %2.2f\n', mean(na.omit(data$exp2_vis)), t_v$conf.int[1],t_v$conf.int[2], t_v$parameter, t_v$statistic, t_v$p.value, as.vector(bf_v))

# bayes factors for the comparison of the two experiments
t_v = t.test(na.omit(data$exp1_vis),na.omit(data$exp2_vis))
t_h = t.test(na.omit(data$exp1_hap),na.omit(data$exp2_hap))
bf_v = ttestBF(na.omit(data$exp1_vis),na.omit(data$exp2_vis),r=rscale)
bf_h = ttestBF(na.omit(data$exp1_hap),na.omit(data$exp2_hap),r=rscale)
print("Comparing performance across experiments:")
fprintf('Comparing the Haptic pulling performance: t(%2.2f)=%2.2f  p=%e, Bayes factor = %2.2f\n', t_h$parameter, t_h$statistic, t_h$p.value, as.vector(1/bf_h))
fprintf('Comparing the Visual familiarity performance: t(%2.2f)=%2.2f p=%e, Bayes factor = %2.2f\n', t_v$parameter, t_v$statistic, t_v$p.value, as.vector(bf_v))

# bayes factors for the within-subject consistency analysis
t_ = t.test(c(na.omit(data$exp1_con),na.omit(data$exp2_con)))
bf_ = ttestBF(c(na.omit(data$exp1_con),na.omit(data$exp2_con)),r=rscale)
print("Within-subject consistency:")
fprintf('Correlation between pulling force and familiarity performance in individual scenes: Rho= %2.2f (%2.2f-%2.2f)  with t(%i)=%2.2f  p=%e, Bayes factor = %2.2f\n', mean(c(na.omit(data$exp1_con),na.omit(data$exp2_con))), t_$conf.int[1],t_$conf.int[2], t_$parameter, t_$statistic, t_$p.value, as.vector(bf_))

# bayes factors for the explicit knowledge analysis
c_r = cor.test(na.omit(data$raw1),na.omit(data$raw2), method = "pearson", conf.level = 0.95)
c_c = cor.test(na.omit(data$contrl1),na.omit(data$contrl2), method = "pearson", conf.level = 0.95)
bf_cr = correlationBF(na.omit(data$raw1),na.omit(data$raw2))
bf_cc = correlationBF(na.omit(data$contrl1),na.omit(data$contrl2))
print("Controlling for explicit knowledge:")
fprintf('Raw correlation (joint, with separate intercepts): Rho= %2.2f (%2.2f-%2.2f) p=%e, Bayes factor = %2.2f\n', c_r$estimate, c_r$conf.int[1],c_r$conf.int[2], c_r$p.value, as.vector(bf_cr))
fprintf('Correlation after controlling for explicitness in both training and test (joint, with separate intercepts): Rho= %2.2f (%2.2f-%2.2f) p=%e, Bayes factor = %2.2f\n', c_c$estimate, c_c$conf.int[1],c_c$conf.int[2], c_c$p.value, as.vector(bf_cc))
