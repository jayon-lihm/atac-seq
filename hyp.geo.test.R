hyp_geo_test_selected <- function(df, x, y, colname_list=FC_id_list, xname="", yname="", lowB, upB){
  ## Test for x vs y & y vs x
  ## Example:
  ## hyp_geo_test_selected(df=f10_all.res, x="AMYGDALA", y="CORTEX", colname_list=FC_id_list, lowB=min(num_FC, num_CTRL), upB=max(num_FC, num_CTRL))
  if (xname==""){ xname <- x }
  if (yname==""){ yname <- y }
  hyp_res.df <- df[,c("name2", "chrom")]
  x_colnames <- colname_list[grep(xname, colname_list)]
  y_colnames <- colname_list[grep(yname, colname_list)]
  num_x <- length(grep(x, colname_list))
  num_y <- length(grep(y, colname_list))
  hyp_res.df[, paste(xname, "_1", sep="")] <- rowSums(df[,x_colnames]) ## Number of 1s for x condition
  hyp_res.df[, paste(xname, "_0", sep="")] <- num_x - hyp_res.df[, paste(xname, "_1", sep="")] ## Number of 0s for x condition
  hyp_res.df[, paste(yname, "_1", sep="")] <- rowSums(df[,y_colnames]) ## Number of 1s for y condition
  hyp_res.df[, paste(yname, "_0", sep="")] <- num_y - hyp_res.df[, paste(yname, "_1", sep="")] ## Number of 0s for y condition
  
  colnames_in_xy <- c(paste(xname, "_1", sep=""), paste(xname, "_0", sep=""), paste(yname, "_1", sep=""), paste(yname, "_0", sep="")) ## Enrichment in X
  colnames_in_yx <- c(paste(yname, "_1", sep=""), paste(yname, "_0", sep=""), paste(xname, "_1", sep=""), paste(xname, "_0", sep="")) ## Enrichment in Y
  
  num_hit <- hyp_res.df[, paste(xname, "_1", sep="")] + hyp_res.df[, paste(yname, "_1", sep="")] ## Total number of 1's
  hyp_res.df2 <- hyp_res.df[which(num_hit >= lowB & num_hit<=upB),] ## Select genes within the range of frequency

  ## Test enrichment in X
  hyp_res.df2$xy_p.value <- apply(as.matrix(hyp_res.df2[,colnames_in_xy]), 1, function(y){ x=y[1]; # number of white balls drawn (=number of amygdala with 1)
                                                                                           m=y[1]+y[2]; # number of white balls=number of amygdala
                                                                                           n=y[3]+y[4]; # number of black balls=number of cortex
                                                                                           k=y[1]+y[3]; # number of draws=number of 1
                                                                                           return( phyper(x, m, n, k, lower.tail=FALSE) + dhyper(x, m, n, k) ) ## p(X>x) + p(X==x)
                                                                                         })
  
  hyp_res.df2$xy_p.value_adj <- p.adjust(hyp_res.df2$xy_p.value, method="BH")

  ## Test enrichment in Y
  hyp_res.df2$yx_p.value <- apply(as.matrix(hyp_res.df2[,colnames_in_yx]), 1, function(y){ x=y[1];
                                                                                           m=y[1]+y[2];
                                                                                           n=y[3]+y[4];
                                                                                           k=y[1]+y[3];
                                                                                           return( phyper(x, m, n, k, lower.tail=FALSE) + dhyper(x, m, n, k) ) ## p(X>x) + p(X==x)
                                                                                         })
  
  hyp_res.df2$yx_p.value_adj <- p.adjust(hyp_res.df2$yx_p.value, method="BH")
  return(hyp_res.df2)
}

