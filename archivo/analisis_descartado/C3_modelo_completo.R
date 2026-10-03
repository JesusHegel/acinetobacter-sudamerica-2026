library(brms); library(cmdstanr)
options(brms.backend = "cmdstanr", mc.cores = 4)

d <- read.delim("~/abaumannii/resultados/modelo_entrada.tsv")
d$pais <- factor(d$pais); d$proyecto <- factor(d$proyecto)
VARS <- c("sin_carbapenemasa", "ST2", "blaOXA23", "blaOXA72")
res <- list(); tab <- data.frame()

for (v in VARS) {
  cat("\n========== ", v, " ==========\n")
  sub <- droplevels(subset(d, variable == v))
  m <- brm(positivos | trials(n) ~ 1 + (1 | pais) + (1 | pais:proyecto),
           family = binomial(), data = sub,
           prior = c(prior(normal(0, 1.5), class = Intercept),
                     prior(normal(0, 1.5), class = sd)),
           chains = 4, iter = 6000, warmup = 3000, cores = 4,
           seed = 42, refresh = 0, silent = 2,
           control = list(adapt_delta = 0.99, max_treedepth = 12))
  res[[v]] <- m
  saveRDS(m, paste0("~/abaumannii/resultados/modelo/m_", v, ".rds"))

  vc <- VarCorr(m, summary = FALSE)
  sp <- vc$pais$sd[, 1]; sy <- vc$`pais:proyecto`$sd[, 1]
  vpc <- sy^2 / (sy^2 + sp^2)
  np <- nuts_params(m)
  div <- sum(subset(np, Parameter == "divergent__")$Value)
  rh <- max(rhat(m), na.rm = TRUE)

  tab <- rbind(tab, data.frame(
    variable = v,
    sd_pais = sprintf("%.2f [%.2f-%.2f]", median(sp), quantile(sp,.025), quantile(sp,.975)),
    sd_proy = sprintf("%.2f [%.2f-%.2f]", median(sy), quantile(sy,.025), quantile(sy,.975)),
    VPC     = sprintf("%.3f [%.3f-%.3f]", median(vpc), quantile(vpc,.025), quantile(vpc,.975)),
    P_VPC_mayor_0.5 = sprintf("%.3f", mean(vpc > 0.5)),
    divergentes = div, rhat_max = sprintf("%.3f", rh)))
  cat(sprintf("  sd pais %.2f | sd proy %.2f | VPC %.3f | div %d | rhat %.3f\n",
              median(sp), median(sy), median(vpc), div, rh))
}

cat("\n\n===== TABLA RESUMEN =====\n")
print(tab, row.names = FALSE)
write.table(tab, "~/abaumannii/resultados/modelo/varianzas.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n\n===== EFECTOS ALEATORIOS POR PROYECTO (para la Figura 2) =====\n")
ef <- data.frame()
for (v in VARS) {
  r <- ranef(res[[v]])$`pais:proyecto`[, , "Intercept"]
  sub <- droplevels(subset(d, variable == v))
  key <- paste(sub$pais, sub$proyecto, sep = "_")
  ef <- rbind(ef, data.frame(variable = v, grupo = rownames(r),
                             est = r[, "Estimate"], lo = r[, "Q2.5"], hi = r[, "Q97.5"]))
}
write.table(ef, "~/abaumannii/resultados/modelo/efectos_proyecto.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  escrito: resultados/modelo/efectos_proyecto.tsv (", nrow(ef), "filas )\n")
