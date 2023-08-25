# FileName="fig2a2"
# i=1
# expr=fread(paste0("OutputData/", FileName, "_count_", i, ".tsv"))
# cell_feature=as.data.frame(fread(paste0("OutputData/", FileName, "_meta_",i,".tsv")))
#

multicell=function(expr, cell_feature, NoSpot=500) {
  cell_loc=cell_feature[,c("x.loc", "y.loc")]
  xrange=range(cell_loc[,1])
  yrange=range(cell_loc[,2])
  rmax=max(xrange, yrange)
  m=round(sqrt(NoSpot))

  r <- raster(ncols=m, nrows=m,xmn=xrange[1], xmx=xrange[2], ymn=yrange[1], ymx=yrange[2])
  spot.idx=cellFromXY(r, cbind(cell_loc[,1], cell_loc[,2]))
  #
  expr2=sapply(1:nrow(expr), function(f) tapply(as.numeric(expr[f,]), spot.idx, sum))
  expr3=t(expr2)
  rownames(expr3)=rownames(expr)
  colnames(expr3)=paste0("Spot", 1:ncol(expr3))

  # Spot Coordinates
  spot_col=colFromCell(r, 1:nrow(expr2))
  spot_row=rowFromCell(r, 1:nrow(expr2))
  spot_coordinates=xyFromCell(r, 1:nrow(expr2))+
    matrix(runif(2*nrow(expr2), -rmax/NoSpot/100, rmax/NoSpot/100), ncol=2)
  spot_loc=data.frame(Spot=colnames(expr3), col=spot_col, row=spot_row, spot_coordinates)


  # Spot Region
  if (is.null(cell_feature$region)==F) {
    dat1=data.frame(spot.idx, region=cell_feature$region) %>%
      group_by(spot.idx, region) %>%
      mutate(count=1) %>%
      summarise(abundance = sum(count)) %>%
      pivot_wider(names_from = region, values_from = abundance, values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-spot.idx) %>%
      rowwise() %>%
      mutate(region = names(.)[which.max(c_across(everything()))])

    spot_loc=data.frame(spot_loc, region=dat1$region)
  }

  # Spot's Cell Type Count
  if (is.null(cell_feature$annotation)==F) {
    dat1=data.frame(spot.idx, cell_feature) %>%
      dplyr::select(spot.idx, annotation) %>%
      group_by(spot.idx, annotation) %>%
      mutate(count=1) %>%
      summarise(abundance = sum(count)) %>%
      pivot_wider(names_from = annotation, values_from = abundance, values_fill = 0)
    spot_loc=data.frame(spot_loc, dat1[,-1])
  }

  return(list(count=expr3, spot_feature=spot_loc))
}
