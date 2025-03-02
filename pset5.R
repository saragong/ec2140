library(tidyverse)

#### Prep dataframe

setwd('/Users/saragong/Documents/ec2140')

prep_for_hypothesis_tests <- function(df) {
  
  df <- df %>%
    mutate(
      whole_sample = 'All',
      gender = ifelse(x14_male_dummy == 1, 'Men', 'Women'),
      married = ifelse(x4_married == 1, 'Married', 'Unmarried'),
      age = case_when(
        x8_age2225==1 ~ '22-25',
        x9_age2629==1 ~ '26-29',
        x10_age3035==1 ~ '30-35',
        x11_age3644==1 ~ '36-44',
        x12_age4554==1 ~ '45-54'
      )
    ) %>%
    filter(!is.na(age)) %>%
    mutate(
      gender_by_married_by_age = paste(age, married, gender)
    ) %>%
    select(-starts_with('x')) %>%
    pivot_longer(
      cols = c('whole_sample', 'gender', 'married', 'age', 'gender_by_married_by_age'),
      names_to = 'category',
      values_to = 'subgroup'
    ) 
  
  return(df)
}

perform_hypothesis_tests <- function(df) {
  
  hypothesis_tests_df <- map_dfr(
    df %>% pull(category) %>% unique(),
    function(category) {
      
      map_dfr(
        df %>% filter(category == !!category) %>% pull(subgroup) %>% unique(),
        function(subgroup) {
          
          # Filter to a subgroup
          subgroup_df <- df %>%
            filter(category == !!category, subgroup == !!subgroup) %>%
            mutate(
              `ZY` = y_income*z_assigned_to_training,
              `(1-Z)Y` = y_income*(1-z_assigned_to_training),
              `Z`=z_assigned_to_training
            ) 
          
          n <- nrow(subgroup_df)
          
          c1_hat <- mean(subgroup_df$ZY)
          c2_hat <- mean(subgroup_df$`(1-Z)Y`)
          c3_hat <- mean(subgroup_df$Z)
          
          # Calculate diff in means
          Delta_hat <- (c1_hat/c3_hat) - (c2_hat/(1-c3_hat))
          
          # Calculate Jacobian
          F_hat <- matrix(
            c(
              1/c3_hat,
              1/(1-c3_hat),
              -(c1_hat/(c3_hat^2))+(c2_hat/(1-c3_hat)^2)
            ),
            ncol=1)
          
          Omega_hat <- subgroup_df %>% select(ZY, `(1-Z)Y`, Z) %>% cov()
          
          se <- sqrt(t(F_hat) %*% Omega_hat %*% F_hat)
          
          # Under the null, Delta_hat / se ~ N(0,1)
          t_stat <- {sqrt(n) * Delta_hat/se} %>% as.numeric() 
          
          # Two-sided test
          p_value <- 2*(1-pnorm(abs(t_stat)))
          
          # Return a row
          tribble(
            ~category, ~subgroup, ~t_stat, ~p_value,
            category, subgroup, t_stat, p_value
          )
          
        }
      )
    }
  )
  
  return(hypothesis_tests_df)
  
}


#### Question 3

alpha <- 0.05

raw_df <- read.csv('problem_set_5_data/Simplified JTPA Data.csv') %>%
  janitor::clean_names()

# Manual t tests
hypothesis_tests_df <- raw_df %>%
  prep_for_hypothesis_tests() %>%
  perform_hypothesis_tests()

# Off-the-shelf t test functions just for comparison

hypothesis_tests_off_the_shelf_df <- raw_df %>%
  prep_for_hypothesis_tests() %>%
  group_by(category, subgroup) %>%
  # Reverse group order so that it matches with the above
  mutate(
    z_assigned_to_training = factor(z_assigned_to_training, levels=c(1,0))
  ) %>%
  # Calculate t tests
  summarise(
    t_stat_off_the_shelf = t.test(
      y_income ~ z_assigned_to_training, conf.level = 1-alpha, alternative = 'two.sided'
    )$statistic,
    p_value_off_the_shelf = t.test(
      y_income ~ z_assigned_to_training, conf.level = 1-alpha, alternative = 'two.sided'
    )$p.value
  ) 

hypothesis_tests_df <- inner_join(hypothesis_tests_df, hypothesis_tests_off_the_shelf_df)

cat('Avg. difference between t-stats from manual vs. off the shelf calculation', mean(abs(hypothesis_tests_df$t_stat - hypothesis_tests_df$t_stat_off_the_shelf)))

# Looks reasonable! It's 0.001

# Make plot
lapply(
  c('whole_sample', 'gender', 'married', 'age', 'gender_by_married_by_age'),
  function(category) {
    
    hypothesis_tests_df %>%
      filter(category == !!category) %>%
      ggplot(., 
             aes(x=t_stat, y=p_value, label = subgroup)
      ) +
      facet_wrap(~category) +
      geom_point(alpha=0.5) + 
      ggrepel::geom_text_repel(size=1.75, segment.alpha=0.5) +
      geom_vline(xintercept=0, linetype = 'dashed') +
      geom_vline(xintercept=qnorm(alpha/2), linetype = 'dotted') +
      geom_vline(xintercept=qnorm(1-(alpha/2)), linetype = 'dotted') +
      xlim(min(hypothesis_tests_df$t_stat, qnorm(alpha/2)), max(hypothesis_tests_df$t_stat, qnorm(1-(alpha/2)))) +
      ylim(min(hypothesis_tests_df$p_value, qnorm(alpha/2)), max(hypothesis_tests_df$p_value, qnorm(1-(alpha/2)))) +
      theme_minimal() +
      theme(
        legend.title = element_blank(),
        axis.title = element_text(size=8)
      ) 
    
  }
) %>%
  gridExtra::grid.arrange(grobs=., ncol=1)

# For which subgroups do you get t-statistics significantly
# different from zero at the 5\% level (based on a two-sided t-test,
# with critical value 1.96), and what are the t-stats for these groups?

hypothesis_tests_df %>%
  filter(p_value < alpha) %>%
  print()

#### Question 4
# If you instead Bonferroni correct for the number of hypotheses
# tested, how many significant results do you get in part (c)?

n_hypotheses <- nrow(hypothesis_tests_df)

alpha_bonferroni <- alpha/n_hypotheses

hypothesis_tests_df %>%
  filter(p_value < alpha_bonferroni) %>%
  print()

# None.

#### Question 5

n_sims <- 10^3

bootstrapped_hypothesis_tests_df <- map_dfr(
  seq(1:n_sims), function(i) {
    
    cat('simulation', i, '\n')
    
    raw_df %>%
      # Resample
      slice_sample(n=nrow(.), replace = T) %>%
      # Permute instrument to create null distribution
      mutate(z_assigned_to_training = sample(z_assigned_to_training)) %>%
      # Perform test
      prep_for_hypothesis_tests() %>%
      perform_hypothesis_tests() %>%
      mutate(sim=i, .before=category)
    
  }
)

#### Question 6
# For the simulation design described in part (e), what
# fraction of null hypotheses are rejected by the 5\% t-test on average
# across simulations draws $s$? In what fraction of draws $s$ is at
# least one null hypothesis rejected?

bootstrapped_hypothesis_tests_df %>%
  mutate(is_false_positive = ifelse(p_value < alpha, 1, 0)) %>%
  summarise(
    any_false_positives = max(is_false_positive),
    frac_false_positives = mean(is_false_positive),
    .by=sim
  ) %>%
  summarise(
    family_wise_error_rate = mean(any_false_positives),
    false_discovery_rate = mean(frac_false_positives)
  )

#### Question 7
# If we instead Bonferroni correct, in what fraction of
# draws $s$ is at least one null hypothesis rejected? Why is it possible
# for this number to exceed 5\%?
  
bootstrapped_hypothesis_tests_df %>%
  mutate(is_false_positive = ifelse(p_value < alpha_bonferroni, 1, 0)) %>%
  summarise(
    any_false_positives = max(is_false_positive),
    frac_false_positives = mean(is_false_positive),
    .by=sim
  ) %>%
  summarise(
    family_wise_error_rate = mean(any_false_positives),
    false_discovery_rate = mean(frac_false_positives)
  )

#### Question 8
# For the same approach as in part (e), compute the 95th
# percentile of the max of the absolute t-statistics (i.e. for each
# draw $s$, take the max of absolute t-statistics across subgroups,
# and take the 95th percentile of these). How does this compare to the
# Bonferroni critical value?

# Critical value from Bonferroni t tests
cat('Critical value from t stat null distribution using Bonferroni correction', qnorm(1-(alpha_bonferroni/2)))

# Critical value from max t
bootstrapped_hypothesis_tests_df %>%
  summarise(
    max_t = max(abs(t_stat)),
    .by = sim
  ) %>%
  pull(max_t) %>%
  quantile(1-(alpha/2)) %>%
  cat('Critical value from max-t stat null distriburtion', .)

#### Question 9
# Suppose we instead limited ourselves to look at t-statistics
# for the overall sample, men, women, married individuals, unmarried
# individuals, and the five age groups, but did not include interactions.
# For which groups do we obtain significant results after Bonferroni
# correction (in the original data)? Does this suggest anything to you?

filtered_hypothesis_tests_df <- hypothesis_tests_df %>% # original dataframe, not bootstrapped
  filter(category != 'gender_by_married_by_age')

n_hypotheses <- nrow(filtered_hypothesis_tests_df)

alpha_bonferroni <- alpha/n_hypotheses

hypothesis_tests_df %>%
  filter(p_value < alpha_bonferroni) %>%
  print()


