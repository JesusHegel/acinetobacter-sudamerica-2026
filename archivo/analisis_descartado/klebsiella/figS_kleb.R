library(brms); library(ggplot2); library(dplyr)
VARS <- c("sin_carbapenemasa","CC258","blaKPC","blaNDM")
ETIQ <- c(sin_carbapenemasa="Sin carbapenemasa adquirida", CC258="Complejo clonal 258",
          blaKPC="blaKPC", blaNDM="blaNDM")
d <- read.delim("~/abaumannii/tmp_kleb/modelo_entrada_kleb.tsv")
PAIS <- c(Brazil="Brasil", Peru="Perú", Chile="Chile", Colombia="Colombia",
          Argentina="Argentina", Ecuador="Ecuador", Uruguay="Uruguay", Paraguay="Paraguay")
todo <- data.frame()
for (v in VARS) {
  m <- readRDS(paste0("~/abaumannii/tmp_kleb/modelo/m_",v,".rds"))
  sub <- droplevels(subset(d, variable==v))
  nd <- data.frame(pais=sub$pais, proyecto=sub$proyecto, n=1)
  nd$`pais:proyecto` <- paste(sub$pais, sub$proyecto, sep="_")
  fit <- fitted(m, newdata=nd, scale="response", summary=TRUE)
  todo <- rbind(todo, data.frame(variable=v, pais=as.character(sub$pais), n=sub$n,
                 est=fit[,"Estimate"], lo=fit[,"Q2.5"], hi=fit[,"Q97.5"]))
}
todo$variable <- factor(ETIQ[todo$variable], levels=ETIQ)
todo$pais <- PAIS[todo$pais]
ord <- todo %>% count(pais, name="np") %>% arrange(desc(np))
todo$pais <- factor(todo$pais, levels=ord$pais)
todo <- todo %>% group_by(variable,pais) %>% arrange(est) %>%
  mutate(id=paste0(pais,"_",row_number())) %>% ungroup()
p <- ggplot(todo, aes(x=est, y=reorder(id, est))) +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, colour="grey55", linewidth=.35) +
  geom_point(aes(size=n, colour=pais)) +
  facet_grid(pais ~ variable, scales="free_y", space="free_y", switch="y") +
  scale_x_continuous(labels=scales::percent_format(accuracy=1), limits=c(0,1)) +
  scale_size_continuous(range=c(1,4), name="Genomas") +
  scale_colour_brewer(palette="Dark2", guide="none") +
  labs(x="Proporción estimada (media posterior e intervalo creíble del 95 %)", y=NULL) +
  theme_bw(base_size=9) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.text.y.left=element_text(angle=0, hjust=1, size=8),
        strip.background=element_rect(fill="grey93", colour=NA),
        panel.grid.minor=element_blank(), legend.position="bottom")
dir.create("~/abaumannii/resultados/figuras", showWarnings=FALSE, recursive=TRUE)
ggsave("~/abaumannii/resultados/figuras/FigS1_klebsiella.pdf", p, width=9, height=8)
ggsave("~/abaumannii/resultados/figuras/FigS1_klebsiella.png", p, width=9, height=8, dpi=300)
cat("figura escrita\n")
cat("proyectos por pais:\n"); print(table(todo$pais)/4)
