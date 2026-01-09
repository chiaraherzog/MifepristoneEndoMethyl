# Colours
pals <- khroma::color('batlow')(8)
# Set treatment colours
cols_trt <- pals[c(1, 6)]


load("out/data.Rdata")
data <- data |> 
  dplyr::mutate(num = gsub("BRCAmut_", "", id),
                trt = gsub(" ", "", trt),
                trt = factor(trt, levels = c("VitaminB", "Mifepristone")))
rownames(data) <- data$id

# Get data in
counts <- data.table::fread("~/Dropbox/data/mifepristone/endometrium/brca/rnaseq/Partek_BRCA_MIMI_SS2_Split_by_attribute_BRCA.txt", sep = '\t')
counts <- as.data.frame(counts) |> tibble::column_to_rownames('Gene Symbol') |> dplyr::select(-c(1:5))
colnames(counts) <- paste0("BRCAmut_", colnames(counts))
counts <- counts[,data$id]
# identical(colnames(counts), data$id)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = data,
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
logcountsMife <- log2(counts(dds, normalized=TRUE) + 1 )

# Gene sets
library(msigdbr)
librar
## Estrogen
h = msigdbr(species = "Homo sapiens",category = "H")
all <- msigdbr(species = 'Homo sapiens')

es_early <- h |> 
  dplyr::filter(gs_id == 'M5906') |> dplyr::filter(gene_symbol %in% rownames(logcountsMife))

es_late <- h |> 
  dplyr::filter(gs_id == 'M5907') |> dplyr::filter(gene_symbol %in% rownames(logcountsMife))

pi3k <- h |> 
  dplyr::filter(gs_id == 'M5923') |> dplyr::filter(gene_symbol %in% rownames(logcountsMife))

ec <- all |> dplyr::filter(gs_id == 'M19877') |> dplyr::filter(gene_symbol %in% rownames(logcountsMife))

identical(data$id, colnames(logcountsMife))

rownames(logcountsMife)
data$es_early <- colSums(logcountsMife[es_early$gene_symbol,])
data$es_late <- colSums(logcountsMife[es_late$gene_symbol,])
data$pi3k <- colSums(logcountsMife[unique(pi3k$gene_symbol),])
data$ec <- colSums(logcountsMife[unique(ec$gene_symbol),])


# Adjuted indices
adjust_age_fib <- function (data, indices){ 
  
  for (i in indices) {
    data[[paste0(i, "_adj")]] <- numeric(nrow(data))
    
    tmp <- data[data$type == "Control", ]
    fit <- lm(tmp[[i]] ~ age + hepidish_Fib, data = tmp)
    data[, paste0(i, "_adj")] <- data[[i]] - as.numeric(predict(fit, newdata = data))
  }
  
  return(data)
}

library(ggtext)

data2 <- data |>
  dplyr::mutate(type = ifelse(time == 'pre', "Control", 'Post'))|> 
  adjust_age_fib(indices = c('es_early', 'es_late', 'pi3k', 'ec'))

a <- data2 |> ggplot(aes(x = time,
                   y = es_early)) +
  geom_boxplot(outlier.shape = NA,
               aes(fill = trt),
               alpha = 0.6) +
  geom_line(aes(group = participant_id)) +
  facet_wrap(~trt) +
  ggpubr::stat_compare_means(paired = T,
                             label = 'p.format',
                             label.y.npc = 0.95) +
  theme_bw() +
  theme(axis.title = element_markdown(),
        legend.position = 'none') +
  labs(x = '',
       y = 'Σ log2(expr)<br><i>unadjusted</i> for cell composition',
       title = 'early response to estrogen (Hallmark M5906)') +
  scale_fill_manual(values = cols_trt,
                    aesthetics = c('colour', 'fill'))

b <- data2 |> ggplot(aes(x = time,
                         y = es_early_adj)) +
  geom_boxplot(outlier.shape = NA,
               aes(fill = trt),
               alpha = 0.6) +
  geom_line(aes(group = participant_id)) +
  facet_wrap(~trt) +
  ggpubr::stat_compare_means(paired = T,
                             label = 'p.format',
                             label.y.npc = 0.95) +
  theme_bw() +
  theme(axis.title = element_markdown(),
        legend.position = 'none') +
  labs(x = '',
       y = 'Σ log2(expr)<br><i>adjusted</i> for cell composition') +
  scale_fill_manual(values = cols_trt,
                    aesthetics = c('colour', 'fill'))

c <- data2 |> ggplot(aes(x = time,
                         y = es_late)) +
  geom_boxplot(outlier.shape = NA,
               aes(fill = trt),
               alpha = 0.6) +
  geom_line(aes(group = participant_id)) +
  facet_wrap(~trt) +
  ggpubr::stat_compare_means(paired = T,
                             label = 'p.format',
                             label.y.npc = 0.95) +
  theme_bw() +
  theme(axis.title = element_markdown(),
        legend.position = 'none') +
  labs(x = '',
       y = 'Σ log2(expr)<br><i>unadjusted</i> for cell composition',
       title = 'late response to estrogen (Hallmark M5907)') +
  scale_fill_manual(values = cols_trt,
                    aesthetics = c('colour', 'fill'))

d <- data2 |> ggplot(aes(x = time,
                         y = es_late_adj)) +
  geom_boxplot(outlier.shape = NA,
               aes(fill = trt),
               alpha = 0.6) +
  geom_line(aes(group = participant_id)) +
  facet_wrap(~trt) +
  ggpubr::stat_compare_means(paired = T,
                             label = 'p.format',
                             label.y.npc = 0.95) +
  theme_bw() +
  theme(axis.title = element_markdown(),
        legend.position = 'none') +
  labs(x = '',
       y = 'Σ log2(expr)<br><i>adjusted</i> for cell composition') +
  scale_fill_manual(values = cols_trt,
                    aesthetics = c('colour', 'fill'))

e <- data2 |> ggplot(aes(x = time,
                         y = pi3k)) +
  geom_boxplot(outlier.shape = NA,
               aes(fill = trt),
               alpha = 0.6) +
  geom_line(aes(group = participant_id)) +
  facet_wrap(~trt) +
  ggpubr::stat_compare_means(paired = T,
                             label = 'p.format',
                             label.y.npc = 0.95) +
  theme_bw() +
  theme(axis.title = element_markdown(),
        legend.position = 'none') +
  labs(x = '',
       y = 'Σ log2(expr)<br><i>unadjusted</i> for cell composition',
       title = 'PI3K/AKT/mTOR signalling (Hallmark M5923)') +
  scale_fill_manual(values = cols_trt,
                    aesthetics = c('colour', 'fill'))

f <- data2 |> ggplot(aes(x = time,
                    y = pi3k_adj)) +
  geom_boxplot(outlier.shape = NA,
               aes(fill = trt),
               alpha = 0.6) +
  geom_line(aes(group = participant_id)) +
  facet_wrap(~trt) +
  ggpubr::stat_compare_means(paired = T,
                             label = 'p.format',
                             label.y.npc = 0.95) +
  theme_bw() +
  theme(axis.title = element_markdown(),
        legend.position = 'none') +
  labs(x = '',
       y = 'Σ log2(expr)<br><i>adjusted</i> for cell composition') +
  scale_fill_manual(values = cols_trt,
                    aesthetics = c('colour', 'fill'))

g <- data2 |> ggplot(aes(x = time,
                    y = ec)) +
  geom_boxplot(outlier.shape = NA,
               aes(fill = trt),
               alpha = 0.6) +
  geom_line(aes(group = participant_id)) +
  facet_wrap(~trt) +
  ggpubr::stat_compare_means(paired = T,
                             label = 'p.format',
                             label.y.npc = 0.95) +
  theme_bw() +
  theme(axis.title = element_markdown(),
        legend.position = 'none') +
  labs(x = '',
       y = 'Σ log2(expr)<br><i>unadjusted</i> for cell composition',
       title = 'KEGG endometrial cancer (M19877)') +
  scale_fill_manual(values = cols_trt,
                    aesthetics = c('colour', 'fill'))

h <- data2 |> ggplot(aes(x = time,
                         y = ec_adj)) +
  geom_boxplot(outlier.shape = NA,
               aes(fill = trt),
               alpha = 0.6) +
  geom_line(aes(group = participant_id)) +
  facet_wrap(~trt) +
  ggpubr::stat_compare_means(paired = T,
                             label = 'p.format',
                             label.y.npc = 0.95) +
  theme_bw() +
  theme(axis.title = element_markdown(),
        legend.position = 'none') +
  labs(x = '',
       y = 'Σ log2(expr)<br><i>adjusted</i> for cell composition') +
  scale_fill_manual(values = cols_trt,
                    aesthetics = c('colour', 'fill'))

design <- '
AABB
CCDD
EEFF
GGHH'

plot <- (a+b+c+d+e+f+g+h) + plot_layout(design = design) + plot_annotation(tag_levels = 'a',title = 'Gene set analysis')

cairo_pdf(filename = 'expression_sets.pdf', width = 9, height = 16)
print(plot)
dev.off()


# Source data figure S7:

out <- a$data|> 
  dplyr::select(-type) |> 
  tidyr::pivot_longer(pi3k:ec_adj,
                      names_to = 'set',
                      values_to = 'value') |> 
  dplyr::group_by(participant_id) |> 
  dplyr::mutate(participant = dplyr::cur_group_id()) |>
  dplyr::ungroup() |> 
  dplyr::select(participant, set, value, trt, time) |> 
  dplyr::arrange(set, trt, time)
write.table(out, file = 'SourceData/s7.txt', quote = F, sep = '\t', row.names = F)
