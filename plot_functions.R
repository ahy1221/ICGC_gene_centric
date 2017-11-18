library(ggplot2)
library(ggrepel)
library(cowplot)
ggplot(dna.VS.rna , aes(x = dna, y = rna)) +
  geom_point( aes( shape = status),  size = 1.2, alpha = .6, color = brewer.pal(3, "Set1")[2])  +
  geom_point( data = dna.VS.rna %>% filter(drivers == "PCAWG_drivers"), color = brewer.pal(3, "Set1")[1]) +
  geom_text_repel(data = dna.VS.rna %>% filter(status != "Others"), aes(label = symbol, color = drivers))+
  #scale_color_manual(values = brewer.pal(3, "Set1")[c(2,1)]) +
  xlab("DNA alteration (Variants + CN)") + ylab("RNA alteration (ASE+altProm+exprOutlier+isSplice+fusion)") +
  theme_minimal()





p1 <- ggplot(for.plot, aes(x = variants.score, y = ase.score)) +
  geom_point(aes( color = driver, size =ase.v.total.hits ), alpha = .5) +
  geom_text_repel(data = for.plot %>%  filter(ase.score > 3e-04 & variants.score == 0 ),  aes(label = symbol), col = cols[1]) +
  geom_text_repel(data = for.plot %>%  filter(ase.score < .5e-04 & variants.score > .8e-04 ),  aes(label = symbol), col = cols[2]) +
  geom_text_repel(data = for.plot %>%  filter(ase.score > .8e-04 & variants.score > .5e-04 ),  aes(label = symbol), col = cols[3]) +
  scale_color_brewer(palette = "Set1") + 
  ylab("ASE score [1/(p.value) * hits]") +  xlab("Variants score[1/(p.value) * hits] ") +
  theme_minimal(base_size = 16)

p2 <- ggplot(for.plot, aes(x = variants.score, y = fusion.score)) +
  geom_point(aes( color = driver, size =fusion.v.total.hits ), alpha = .5) +
  geom_text_repel(data = for.plot %>%  filter(fusion.score > 4e-05 & variants.score <= .1e-04 ),  aes(label = symbol), col = cols[1]) +
  geom_text_repel(data = for.plot %>%  filter(fusion.score <= 1e-04 & variants.score > .8e-04 ),  aes(label = symbol), col = cols[2]) +
  geom_text_repel(data = for.plot %>%  filter(fusion.score > 4e-05 & variants.score > .8e-04 ),  aes(label = symbol), col = cols[3]) +
  scale_color_brewer(palette = "Set1") + 
  ylab("Fusion score [1/(p.value) * hits]") +  xlab("Variants score[1/(p.value) * hits] ") +
  theme_minimal(base_size = 16)
p3 <- ggplot(for.plot, aes(x = variants.score, y = expr_outlier.score)) +
  geom_point(aes( color = driver, size =expr_outlier.v.total.hits ), alpha = .5)  +
  geom_text_repel(data = for.plot %>%  filter(expr_outlier.score > 1.2e-04 & variants.score <= .1e-04 ),  aes(label = symbol), col = cols[1]) +
  geom_text_repel(data = for.plot %>%  filter(expr_outlier.score <= 1e-05 & variants.score > .8e-04 ),  aes(label = symbol), col = cols[2]) +
  geom_text_repel(data = for.plot %>%  filter(expr_outlier.score > 5e-05 & variants.score > .8e-04 ),  aes(label = symbol), col = cols[3]) +
  scale_color_brewer(palette = "Set1") + 
  ylab("Fusion score [1/(p.value) * hits]") +  xlab("Variants score[1/(p.value) * hits] ") +
  theme_minimal(base_size = 16)

p4 <- ggplot(for.plot, aes(x = variants.score, y = isSplice.score)) +
  geom_point(aes( color = driver, size =isSplice.v.total.hits ), alpha = .5)  +
  geom_text_repel(data = for.plot %>%  filter(isSplice.score > 5e-04 & variants.score <= .1e-04 ),  aes(label = symbol), col = cols[1]) +
  geom_text_repel(data = for.plot %>%  filter(isSplice.score <= .5e-05 & variants.score > 1e-04 ),  aes(label = symbol), col = cols[2]) +
  geom_text_repel(data = for.plot %>%  filter(isSplice.score > 1e-04 & variants.score > .5e-04 ),  aes(label = symbol), col = cols[3]) +
  scale_color_brewer(palette = "Set1") + 
  ylab("isSplice score [1/(p.value) * hits]") +  xlab("Variants score[1/(p.value) * hits] ") +
  theme_minimal(base_size = 16)


ggplot(count.liver.gather, aes(x = log10(v.hits +1) , y =log10(hits +1 ))) +
  geom_point( size = 2) +
  geom_point(data = filter(count.liver.gather, drivers == "driver"), color = "red") +
  geom_text_repel(data = filter(count.liver.gather, drivers == "driver"), aes(label = symbol)) +
  facet_grid( .~alteration) 