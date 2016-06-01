

print_BX_scatter <- function(tumor_file, normal_file, output_dir) {
  
  if (!file.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  SV_list <- unique(read.table(tumor_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V7)
  
  tumor_tag <- read.table(tumor_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
    group_by(V6, V7, V8) %>% summarize(count = n()) %>%
    spread(V8, count) %>% select(BX = V6, SV = V7, BP1 = `1`, BP2 = `2`) %>%
    filter(!is.na(BP1)) %>% filter(!is.na(BP2))
  
  normal_tag <- read.table(normal_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
    group_by(V6, V7, V8) %>% summarize(count = n()) %>%
    spread(V8, count) %>% select(BX = V6, SV = V7, BP1 = `1`, BP2 = `2`) %>%
    filter(!is.na(BP1)) %>% filter(!is.na(BP2))
  
  
  for (i in 1:length(SV_list)) {
    
    T1 <- tumor_tag %>% filter(SV == SV_list[i])
    T1$type <- rep("tumor", nrow(T1))
    T1$BP1[T1$BP1 >= 50] <- 50
    T1$BP2[T1$BP2 >= 50] <- 50
    
    N1 <- normal_tag %>% filter(SV == SV_list[i])
    N1$type <- rep("control", nrow(N1))
    N1$BP1[N1$BP1 >= 50] <- 50
    N1$BP2[N1$BP2 >= 50] <- 50
    
    ggplot(rbind(T1, N1), aes(x = BP1, y = BP2, colour = type)) +
      geom_jitter(width = 1, height = 1, alpha = 0.7) +
      theme_bw() +
      xlim(c(0, 55)) + ylim(c(0, 55)) +
      ggtitle(SV_list[i])
    
    ggsave(paste(output_dir, "/", SV_list[i], ".png", sep = ""), width = 8, height = 8)
    
  }
  
  
}



print_BX_scatter("Maekawa-1T.BX.tag.txt", "Maekawa-1N.BX.tag.txt", "Maekawa-1")
print_BX_scatter("Maekawa-3T.BX.tag.txt", "Maekawa-3N.BX.tag.txt", "Maekawa-3")


