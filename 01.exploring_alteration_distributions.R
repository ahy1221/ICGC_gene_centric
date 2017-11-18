rm(list = ls())
#source("../pcawg_gene_fusions/aux/")
source("./R/functions.R") # Analysis functions
#source("./gene_sets.R") # Analysis used gene sets


#======= Library and settring
library(tidyr)
library(dplyr)
library(RColorBrewer)
col.set1 = c(brewer.pal(9,"Set1"), brewer.pal(8, "Dark2"))

load("./data/gene.centric.rds")

alterations <- get_all_alterations(gt, alterations = names(gt)[6:14], alteration.names = names(gt)[6:14])

#==== Summary donor hits
count.aggregate.donors<- lapply(alterations, function(x) {setNames(s<-Rfast::rowsums(x), rownames(x))}) %>% 
    do.call("cbind", .) %>% 
    as.data.frame()

count.aggregate.donors[1:4, 1:4]
########################
# Look distribution
#########################
par(mfrow = c(2, 2))

plot(log10(count.aggregate.donors$variants + 1), log10(count.aggregate.donors$ase_all + 1))
plot(log10(count.aggregate.donors$variants + 1), log10(count.aggregate.donors$fusion + 1))
plot(log10(count.aggregate.donors$variants + 1), log10(count.aggregate.donors$isSplice + 1))
plot(log10(count.aggregate.donors$variants + 1), log10(count.aggregate.donors$expr_outlier + 1))

par(mfrow = c(3, 2))
hist(log10(count.aggregate.donors$fusion + 1))
hist(count.aggregate.donors$ase)
hist(count.aggregate.donors$isSplice)
hist(log10(count.aggregate.donors$expr_outlier + 1))
hist(count.aggregate.donors$variants)
hist(log10(count.aggregate.donors$cn + 1))



####################
#Liver
###################
donor.liver <- donor2cancer.use[donor2cancer.use$histology_abbreviation == "Liver-HCC", "icgc_donor_id"]

expr.tpm.liver = expr.tpm.use[,donor.liver]
cn.liver = cn[, donor.liver]
variants.liver = variants[,donor.liver]
ase.liver = ase[,donor.liver]
fusion.liver = fusion[, donor.liver]
alt_prom.liver = alt_prom[,donor.liver]
expr_outlier.liver = expr_outlier[, donor.liver]
isSplice.liver = isSplice[, donor.liver]








for.plot <- data.frame(p.val.ase = p.val.ase, p.val.variants = p.val.variants, 
                       p.val.fusion = p.val.fusion, p.val.expr_outlier = p.val.expr_outlier,
                       p.val.alt_prom = p.val.alt_prom, p.val.isSplice = p.val.isSplice,
                       ase.hits = raw.hits.ase[names(p.val.ase)], 
                       variants.hits = raw.hits.variants[names(p.val.variants)],
                       fusion.hits = raw.hits.fusion[names(p.val.fusion)],
                       expr_outlier.hits = raw.hits.expr_outlier[names(p.val.expr_outlier)],
                       alt_prom.hits = raw.hits.alt_prom[names(p.val.alt_prom)],
                       isSplice.hits = raw.hits.isSplice[names(p.val.isSplice)],
                       symbol = names(p.val.ase))
#save.image("resample.done.Rda")
for.plot$driver <- ifelse(for.plot$symbol %in% HCC.driver, "PCAWG_Driver", "Passenger" )


for.plot <- mutate(for.plot, 
       ase.score = (1/(p.val.ase + 1e5) * ase.hits),
       fusion.score = (1/(p.val.fusion + 1e5) * fusion.hits),
       expr_outlier.score = (1/(p.val.expr_outlier + 1e5) * expr_outlier.hits),
       alt_prom.score = (1/(p.val.alt_prom + 1e5) * alt_prom.hits),
       isSplice.score = (1/(p.val.alt_prom + 1e5) * isSplice.hits),
       variants.score = (1/(p.val.variants +1e5) * variants.hits))
for.plot <- mutate(for.plot,
  ase.v.total.hits = ase.hits + variants.hits,
  fusion.v.total.hits = fusion.hits + variants.hits,
  expr_outlier.v.total.hits = expr_outlier.hits + variants.hits,
  isSplice.v.total.hits = isSplice.hits + variants.hits,
  alt_prom.v.total.hits = alt_prom.hits + variants.hits)


library(RColorBrewer)

cols <- brewer.pal(9, "Spectral")



plot_grid(plotlist = list(p1, p2))


count.liver.gather <- gather(count.liver, key = "alteration", value = "hits", cn.hits, ase.hits, expr_outlier.hits, alt_prom.hits,fusion.hits )
count.liver.gather <- count.liver.gather %>% 
  mutate( drivers = ifelse(symbol %in% liver.driver.TCGA, "driver", "passenger")) 







count.aggregate.donors %>% 
  gather(key = "alteration", value = "hits", variants, cn, ase, expr_outlier, isSplice, fusion, alt_prom ) %>% 
  group_by(alteration, hits) %>% 
  summarise( hit.freq = n())

p.liver <- count.liver %>% 
  gather(key = "alteration", value = "hits", variants, cn, ase, expr_outlier, isSplice, fusion, alt_prom ) %>% 
  group_by(alteration, hits) %>% 
  summarise( hit.freq = n())  %>% 
  ggplot( aes( x= hits, y = hit.freq)) +
  geom_point(size = pts)+
  geom_smooth() +
  facet_grid( alteration ~ ., scale = "free_y" ) +
  theme_minimal(base_size = 16)

p.total <- count.aggregate.donors %>% 
  gather(key = "alteration", value = "hits", variants, cn, ase, expr_outlier, isSplice, fusion, alt_prom ) %>% 
  group_by(alteration, hits) %>% 
  summarise( hit.freq = n())  %>% 
  ggplot( aes( x= hits, y = hit.freq)) +
  geom_point(size = pts)+
  geom_smooth() +
  facet_grid( alteration ~ ., scale = "free_y" ) +
  theme_minimal(base_size = 16)



