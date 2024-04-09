plot_missingnes_by_marker <- function(geno_obj, bin_size = 0.01){
  gi <- geno_obj$geno
  snp_rates <- adegenet::propTyped(gi, by = 'loc')
  
  # Binning missing rate by groups
  missing_df <- data.frame(MissingRate = snp_rates) %>%
    mutate(
      Bin = cut(MissingRate,
                breaks = seq(min(MissingRate),
                             max(MissingRate) + bin_size,
                             by = bin_size),
                right = F))
  
  
  result <- missing_df %>%
    group_by(Bin) %>%
    summarize(missing_groups_count = n()) %>%
    mutate(
      cummulative_snps = cumsum(missing_groups_count)
    )
  
  result$snp_missing_group <- sapply(result$Bin, get_midpoint)
  
  
  out <- ggplot(data = result, aes(x = 1-snp_missing_group, y = cummulative_snps)) +
    geom_point(size = 1, color = "#0F9D58") + theme_bw(base_size = 9) +
    labs(x = "Missing rate by marker", y = "Number of SNPs") +
    scale_x_continuous(labels = scales::percent)+
    scale_y_continuous(labels = addUnits)
  return(out)
}

plot_MAF <- function(geno_obj,bin_size = 0.01){
  gi <- geno_obj$geno
  snp_rates <- minorAllele(gi)
  
  # Binning missing rate by groups of 1%
  
  maf_df <- data.frame(MinorFreq = snp_rates) %>%
    mutate(
      Bin = cut(MinorFreq,
                breaks = seq(min(MinorFreq),
                             max(MinorFreq) + bin_size,
                             by = bin_size),
                right = F)
    )
  
  result <- maf_df %>%
    group_by(Bin) %>%
    summarize(maf_groups_count = n()) %>%
    mutate(
      cummulative_snps = sum(maf_groups_count) - cumsum(maf_groups_count)
    )
  
  result$snp_maf_group <- sapply(result$Bin, get_midpoint)
  
  out <- ggplot(data = result, aes(x = snp_maf_group, y = cummulative_snps)) +
    geom_point(size = 1, color = "#0F9D58") + theme_bw(base_size = 9) +
    labs(x = "MAF by marker", y = "Number of SNPs") +
    scale_x_continuous(labels = scales::percent)+
    scale_y_continuous(labels = addUnits)
  
  return(out)
}

plot_overall_missingness <- function(geno_obj){
  
  gi <- geno_obj$geno
  
  snp_rates <- adegenet::propTyped(gi, by = 'loc')
  
  n_samples <- adegenet::nInd(gi)
  
  
  n_samp_genotyped <- snp_rates * n_samples
  
  missing_df <- data.frame(NS = as.integer(n_samp_genotyped))
  
  counts <- missing_df %>%
    group_by(NS) %>%
    summarize(fNS = n()) %>%
    mutate(c_calls = NS * fNS,
           cumulative_snps = sum(fNS) - cumsum(fNS),
           cummulative_calls = sum(c_calls) - cumsum(c_calls),
           p_miss = cummulative_calls/(sum(fNS)*n_samples)) %>%
    arrange(p_miss)
  
  lab <- findInterval(seq(.1,1,.1), counts$p_miss)
  lab <-  lab[! (duplicated(lab) | duplicated(lab, fromLast = T))]
  
  indicators <- data.frame(NS = counts[lab,'NS'],
                           p_miss = counts[lab,'p_miss'],
                           cumulative_snps = counts[lab,'cumulative_snps']) %>%
    mutate(
      labC = paste(NS,scales::percent(p_miss), sep = ' < NS\n'),
      labP = paste(NS,format(cumulative_snps, big.mark = ','), sep = ' < NS\n')
    )
  

  a <- ggplot(data = counts, aes(x = NS, y = p_miss)) +
    geom_point(size = 1, color = "#0F9D58") + theme_bw(base_size = 9) +
    geom_point(data = indicators, color = '#F4B400', size = 1.5) +
    labs(x = "Number of accessions genotyped (NS)", y = "Percentage of missing data") +
    geom_text_repel(data = indicators, aes(x = NS, y = p_miss, label = labP), size = 2) +
    scale_y_continuous(labels = scales::percent)
  
  b <- ggplot(data = counts, aes(x = NS, y = cumulative_snps)) +
    geom_point(size = 1, color = "#0F9D58") + theme_bw(base_size = 9) +
    geom_point(data = indicators, color = '#F4B400', size = 1.5) +
    geom_text_repel(data = indicators, aes(x = NS, y = cumulative_snps, label = labC), size = 2) +
    labs(x = "Number of accessions genotyped (NS)", y = "Number of variants") +
    scale_y_continuous(labels = addUnits)
  
  
  c <- ggarrange(a, b, ncol = 2, nrow = 1)
  return(c)
}
