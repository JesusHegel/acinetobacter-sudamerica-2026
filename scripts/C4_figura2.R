library(brms); library(ggplot2); library(dplyr); library(tidyr)

VARS <- c("sin_carbapenemasa","ST2","blaOXA23","blaOXA72")
ETIQ <- c(sin_carbapenemasa = "Sin carbapenemasa adquirida",
          ST2 = "Linaje ST2", blaOXA23 = "blaOXA-23", blaOXA72 = "blaOXA-72")

d <- read.delim("~/abaumannii/resultados/modelo_entrada.tsv")
todo <- data.frame()

for (v in VARS) {
  m <- readRDS(paste0("~/abaumannii/resultados/modelo/m_", v, ".rds"))
  sub <- droplevels(subset(d, variable == v))
  nd <- data.frame(pais = sub$pais, proyecto = sub$proyecto, n = 1)
  nd$`pais:proyecto` <- paste(sub$pais, sub$proyecto, sep = "_")
  fit <- fitted(m, newdata = nd, scale = "response", summary = TRUE)
  todo <- rbind(todo, data.frame(
    variable = v, pais = sub$pais, proyecto = sub$proyecto,
    n = sub$n, obs = sub$pct/100,
    est = fit[,"Estimate"], lo = fit[,"Q2.5"], hi = fit[,"Q97.5"]))
}

todo$variable <- factor(ETIQ[todo$variable], levels = ETIQ)
ord <- todo %>% group_by(pais) %>% summarise(m = mean(est)) %>% arrange(desc(m))
todo$pais <- factor(todo$pais, levels = ord$pais)
todo <- todo %>% group_by(variable, pais) %>% arrange(est) %>%
  mutate(id = paste0(pais, "_", row_number())) %>% ungroup()

p <- ggplot(todo, aes(x = est, y = reorder(id, est))) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = "grey55", linewidth = .4) +
  geom_point(aes(size = n, colour = pais)) +
  facet_grid(pais ~ variable, scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  scale_size_continuous(range = c(1.3, 4.2), name = "Genomas") +
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  labs(x = "Proporción estimada (media posterior e intervalo creíble del 95 %)", y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
        strip.background = element_rect(fill = "grey93", colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave("~/abaumannii/resultados/figuras/Fig2_estratificacion.pdf", p, width = 9, height = 7)
ggsave("~/abaumannii/resultados/figuras/Fig2_estratificacion.png", p, width = 9, height = 7, dpi = 300)
cat("figuras escritas en resultados/figuras/\n\n")

cat("=== amplitud del rango estimado por pais y variable ===\n")
todo %>% group_by(variable, pais) %>%
  summarise(np = n(), min = min(est), max = max(est), rango = max(est)-min(est), .groups="drop") %>%
  filter(np >= 2) %>% arrange(variable, desc(rango)) %>%
  mutate(across(c(min,max,rango), ~sprintf("%.3f", .))) %>% as.data.frame() %>% print(row.names=FALSE)
