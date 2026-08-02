library(brms); library(cmdstanr)
options(brms.backend="cmdstanr", mc.cores=4)
d <- read.delim("~/abaumannii/tmp_kleb/modelo_entrada_kleb.tsv")
d$pais <- factor(d$pais); d$proyecto <- factor(d$proyecto)
VARS <- c("sin_carbapenemasa","CC258","blaKPC","blaNDM")
dir.create("~/abaumannii/tmp_kleb/modelo", showWarnings=FALSE)
tab <- data.frame()
for (v in VARS) {
  cat("\n=====", v, "=====\n")
  sub <- droplevels(subset(d, variable==v))
  m <- brm(positivos | trials(n) ~ 1 + (1|pais) + (1|pais:proyecto),
           family=binomial(), data=sub,
           prior=c(prior(normal(0,1.5), class=Intercept),
                   prior(normal(0,1.5), class=sd)),
           chains=4, iter=6000, warmup=3000, cores=4, seed=42,
           refresh=0, silent=2, control=list(adapt_delta=0.99, max_treedepth=12))
  saveRDS(m, paste0("~/abaumannii/tmp_kleb/modelo/m_",v,".rds"))
  vc <- VarCorr(m, summary=FALSE)
  sp <- vc$pais$sd[,1]; sy <- vc$`pais:proyecto`$sd[,1]
  vpc <- sy^2/(sy^2+sp^2)
  div <- sum(subset(nuts_params(m), Parameter=="divergent__")$Value)
  tab <- rbind(tab, data.frame(variable=v,
    sd_pais=sprintf("%.2f [%.2f-%.2f]", median(sp), quantile(sp,.025), quantile(sp,.975)),
    sd_proy=sprintf("%.2f [%.2f-%.2f]", median(sy), quantile(sy,.025), quantile(sy,.975)),
    VPC=sprintf("%.3f [%.3f-%.3f]", median(vpc), quantile(vpc,.025), quantile(vpc,.975)),
    P_mayor_05=sprintf("%.3f", mean(vpc>0.5)), div=div,
    rhat=sprintf("%.3f", max(rhat(m), na.rm=TRUE))))
  cat(sprintf("  sd pais %.2f | sd proy %.2f | VPC %.3f | div %d\n",
              median(sp), median(sy), median(vpc), div))
}
cat("\n\n===== RESUMEN K. pneumoniae =====\n")
print(tab, row.names=FALSE)
write.table(tab, "~/abaumannii/tmp_kleb/modelo/varianzas_kleb.tsv", sep="\t", quote=FALSE, row.names=FALSE)
