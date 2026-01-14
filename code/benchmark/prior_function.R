get_prior = function(sc, sample = NULL, cell_type = NULL, meta_sc = NULL, filter_pd = T) {
  
  if(is.null(meta_sc)) meta_sc = data.frame(sample = sample, cell_type = cell_type)
  meta_sc$sample = as.character(meta_sc$sample)
  sample = unique(meta_sc[, c('sample')])
  
  cell_type = sort(unique(meta_sc$cell_type))
  K = length(cell_type)
  
  cts = array(NA, dim = c(nrow(sc), length(sample), K))
  rownames(cts) = rownames(sc)
  colnames(cts) = sample
  dimnames(cts)[[3]] = cell_type
  for(j in sample) {
    for(k in dimnames(cts)[[3]]) {
      id = which(meta_sc$sample == j & meta_sc$cell_type == k)
      if(length(id) > 0) cts[,j,k] = rowMeans(sc[, id, drop = F])
    }
    id = which(colMeans(is.na(cts[,j,])) == 0)
    cts[,j,id] = log2(cpm(cts[,j,id]) + 1) # make it log2 CPM + 1
  }
  
  cov = array(NA, dim = c(nrow(cts), K, K))
  rownames(cov) = rownames(sc)
  colnames(cov) = dimnames(cov)[[3]] = cell_type
  for(i in 1:nrow(cov)) {
    cov[i,,] = cov(cts[i,,], use = 'pairwise')
  }
  
  if(filter_pd) {
    gene_pd = apply(cov, 1, is.positive.definite)
    print(paste(sum(1-gene_pd), 'genes are filtered out because cell-type covariance matrix is not positive-definite (PD);'))
    print('The filtering can be disabled by setting filter_pd = FALSE. Note the prior cell-type covariance matrix for each gene is required to be PD.')
  } else gene_pd = 1:nrow(cov)
  cov <- cov[gene_pd,,]
  
  profile = matrix(NA, nrow(sc), K)
  rownames(profile) = rownames(sc)
  colnames(profile) = cell_type
  for(i in cell_type) {
    profile[,i] = log2(cpm(rowMeans(sc[, meta_sc$cell_type == i])) + 1)
  }
  
  return(list(profile = profile[gene_pd,], covariance = cov)) # ctsExp = cts, 
}
