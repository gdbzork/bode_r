#1 Start R

#2 load a file into a data frame
a_d = read.table('limma_A-D_all.txt',header=TRUE,sep='\t',quote='')
a_n = read.table('limma_A-N_all.txt',header=TRUE,sep='\t',quote='')

#3 Keep the ones with a p-value < 0.05:
a_d_sig = a_d[a_d$adj.P.Val < 0.05,]
a_n_sig = a_n[a_d$adj.P.Val < 0.05,]
# that gives us lists of the cases where A is significant with respect to
# D and N.

# Or greater than 0.05:
y = x[x$adj.P.Val >= 0.05,]

#4 Create a couple of functions:

subtractList = function(list1,list2) {
  names1 = list1$ProbeId
  names2 = list2$ProbeId
  only1 = !(names1 %in% names2)
  return(list1[only1,])
}

intersectList = function(list1,list2) {
  names1 = list1$ProbeId
  names2 = list2$ProbeId
  inboth = names1 %in% names2
  return(list1[inboth,])
}

#5 Use them to find the list which is significant with respect to both:

a_d_n = intersectList(a_d_sig,a_n_sig)
# this gives you the ones significant to both.  Then repeat with other data
# tables similarly.

