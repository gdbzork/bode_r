source('correlate.R')
fns = commandArgs(TRUE)
if (length(fns) != 4) {
  print("Wrong args...")
  quit()
}
afn = fns[1]
bfn = fns[2]
outRaw = fns[3]
outNor = fns[4]

a = loadSumms(afn)
b = loadSumms(bfn)
pdf(outRaw)
tcor(a,b,afn,bfn)
dev.off()
pdf(outNor)
tcorn(a,b,afn,bfn)
dev.off()
