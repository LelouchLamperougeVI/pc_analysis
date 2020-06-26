f <- "/mnt/storage/rrr_magnum/M2/EE001/2018_12_18/2018_12_18_1.abf"
abf <- readABF(f)

d <- as.data.frame(abf)

plot(d[1:100, ])