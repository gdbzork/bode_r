filterGo <- function(goset,gocat) {
  flag <- sapply(goset,function(x) x$Ontology == gocat)
  filt <- goset[flag]
  goid <- sapply(filt,function(x) x$GOID)
  onto <- sapply(filt,function(x) x$Ontology)
  return(unique(data.frame(goid,onto,stringsAsFactors=F)))
}

uniprot2go <- function(tbl,map,gocat) {
  xmap <- revmap(map)
  eg <- lapply(tbl$uniprot,function(x) org.Hs.egGO[[xmap[[x]]]])
  eg2 <- lapply(eg,filterGo,gocat)
  return(eg2)
}

loadUniprot <- function(fn) {
  df <- read.table(fn,stringsAsFactors=F,header=F)
  colnames(df) = c('uniprot','symbol','score')
  return(df)
}

init <- function(lib) {
  library(lib)
  library(GO.db)
  library(AnnotationDbi)
}

loadGoa <- function(fn) {
  goa = read.table(fn,header=F,stringsAsFactors=F,sep="\t",quote="")
  colnames(goa) = c('db','uniprot','symbol','qualifier','goid','db_ref',
                    'evidence','with','aspect','name','synonym','object_type',
                    'taxon','date','source','empty','splice_form')
  goa = goa[goa$qualifier == "",]
  return(goa)
}

uniprot2go <- function(uniprot,goa,gocat) {
  x = lapply(uniprot$uniprot,
             function(a) unique(goa$goid[goa$uniprot==a&goa$aspect==gocat]))
  names(x) = uniprot$uniprot
  return(x)
}
