dat <- read.table("EUR.QC.het", header=TRUE)

m <- mean(dat$F)
s <- sd(dat$F)

valid <- subset(dat, F <= m + 3*s & F >= m - 3*s)

write.table(valid[, c(1, 2)], "EUR.valid.sample",
            quote=FALSE, row.names=FALSE)
