source('readsToWindows.R')

filename = '/home/evansj/me/tools/fourSig/Beta-major-prom-FL.short.tab'
data <- readTab(filename, onlyMappable=FALSE)
DATA.CHR = 1
DATA.READ.COUNT <- 3
chr = 1
window.size = 2
chr.data <- data[data[DATA.CHR] == chr,]
chr.reads <- chr.data[, DATA.READ.COUNT]
windows <- calc.windows(data=chr.reads, window.size=window.size)
print(windows)