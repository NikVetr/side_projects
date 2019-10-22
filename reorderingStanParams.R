library(rstan)

t <- "
beta[1] = sigbR0j
beta[2] = sigbR0i
beta[3] = bR0b
beta[4] = sigbR1j
beta[5] = sigbR1i
beta[6] = bR1b
beta[7] = A0b
beta[8] = sigA0j
beta[9] = sigA0i
beta[10] = A1b
beta[11] = sigA1j
beta[12] = sigA1i
beta[13] = tunerExp2
beta[14] = tunerExp3
beta[15:1784] = A0j
beta[1785:1790] = A0i
beta[1791:3560] = A1j
beta[3561:3566] = A1i
beta[3567:5336] = bR0j
beta[5337:5342] = bR0i
beta[5343:7112] = bR1j
beta[7113:7118] = bR1i"

cat(paste0(sapply(1:22, function(x) strsplit(strsplit(t, "\n")[[1]][2:23], " = ")[[x]][2]), " = ", sapply(1:22, function(x) strsplit(strsplit(t, "\n")[[1]][2:23], " = ")[[x]][1]), ";\n"))
