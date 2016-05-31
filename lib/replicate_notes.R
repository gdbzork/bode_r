hsa_do388_map = match(hsa_names,hsa_do388$name)
hsa_do388_map
colnames(hsa_do388) = colnames()
colnames(hsa_do388) = colnames
hsa_do388_map = match(hsa_names,hsa_do388$name)
hsa_do388_map
hsa_do388_filt = hsa_do388[hsa_do388_map,]
dim(hsa_do388_filt)
hsa_do579 = read.table('trna_hsa_counts_do579.txt')
colnames(hsa_do579) = colnames
hsa_do579_map = match(hsa_names,hsa_do579$name)
hsa_do579_map == hsa_do388_map
hsa_do579_filt = hsa_do579[hsa_do579_map,]
hsa_replicates = cbind(hsa_do388_filt$left+hsa_do388_filt$overlap+hsa_do388_filt$right,hsa_do388_filt$left+hsa_do579_filt$overlap+hsa_do579_filt$right,as.vector(hsa_do579_filt$name))
hsa_replicates[1:10,]
colnames(hsa_replicates) = c('do388','do579','name')
hsa_replicates[1:10,]
as.data.frame(hsa_replicates)
write.table(hsa_replicates,file="hsa_replicates.txt",quote=FALSE,row.names=FALSE,sep='\t')
history()
ls()
write.table(cfa_do411_filt,file='trna_cfa_counts_do411_final.txt',quote=FALSE,sep='\t',row.names=FALSE)
write.table(cfa_do652_filt,file='trna_cfa_counts_do652_final.txt',quote=FALSE,sep='\t',row.names=FALSE)
write.table(hsa_do388_filt,file='trna_hsa_counts_do388_final.txt',quote=FALSE,sep='\t',row.names=FALSE)
write.table(hsa_do579_filt,file='trna_hsa_counts_do579_final.txt',quote=FALSE,sep='\t',row.names=FALSE)


countup = function(x) { x$left+x$overlap+x$right }
