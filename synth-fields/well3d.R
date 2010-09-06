
createCubeFrame <- function(dimX, dimY, dimZ) {
  zz <- rep(0:(dimZ-1), dimX * dimY)
  yy <- rep(rep(0:(dimY-1), each=dimZ), dimX)
  xx <- rep(0:(dimX-1), each=dimZ*dimY)
  data.frame(x = xx, y = yy, z = zz)
}

prepareFieldFrame <- function(df) {
  df$isElectrode <- 0
  df$V <- 0
  df
}

writeField <- function(df, file) {
  write.table(df[c("x", "y", "z", "isElectrode", "V")], file,
              col.names=F, row.names=F)
}

createSinWell <- function(dimX, dimY, dimZ) {
  d <- prepareFieldFrame(createCubeFrame(dimX, dimY, dimZ))
  d$r <- with(d, sqrt((x - dimX/2)^2 + (y - dimY/2)^2 + (z - dimZ/2)^2))
  d$rArg <- d$r * (2 * pi / (min(dimX, dimY, dimZ) + 1))
  d$V <- (1 - cos(d$rArg))
  d
}

# writeField(createSinWell(100, 100, 100), "../work/sin-well-100x100x100.txt")
