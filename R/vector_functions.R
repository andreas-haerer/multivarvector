#theta function
f.theta <- function(table.of.diff.in.means, unique.id, select.col){
  angle.matrix.radians <- matrix(nrow = length(unique.id), ncol = length(unique.id))
  for(i in 1:length(unique.id)){
    for(j in 1:length(unique.id)){
      angle.matrix.radians[i,j] <- round(acos(cor(x = t(table.of.diff.in.means[i,select.col]), t(table.of.diff.in.means[j,select.col]), use = "na.or.complete")), 3) #
    }
  }

  rownames(angle.matrix.radians) <- unique.id
  colnames(angle.matrix.radians) <- unique.id
  angle.matrix.degrees <- round(angle.matrix.radians*(180/pi), 3)
  angle.output <- list(angle.matrix.radians, angle.matrix.degrees)
  names(angle.output) <- c("theta.radians", "theta.degrees")
  return(angle.output)
}

#meanL function
f.meanL <- function(table.of.diff.in.means, unique.id, select.col){
  get.vectlength <- function(vect){
    return(sqrt(sum(vect^2, na.rm = TRUE)))
  }
  table.of.diff.in.means.t <- table.of.diff.in.means[,select.col]
  length.of.vector.by.group <- apply(table.of.diff.in.means.t, 1 ,get.vectlength) #1 signifies by row. apply the function to every row in the input dataframe
  length.diff.matrix <- matrix(nrow = length(unique.id), ncol = length(unique.id))
  for(i in 1:length(unique.id)){
    for(j in 1:length(unique.id)){
      length.diff.matrix[i,j] <-mean(round(c(length.of.vector.by.group[i], length.of.vector.by.group[j]), 3))
    }
  }
  rownames(length.diff.matrix) <- unique.id
  colnames(length.diff.matrix) <- unique.id

  length.of.vector.by.group <- as.data.frame(length.of.vector.by.group)
  length.of.vector.by.group$wshd <- unique.id
  length.output <- list(length.diff.matrix, length.of.vector.by.group)
  return(length.output)
}

diff <- function(x) {x-lag(x)}

#Vector analysis function
multivariate.vector <- function(infile){

  infile2 <- infile[,-c(1:3)]

  infile_mean <- aggregate(infile2[,2:length(colnames(infile2))], list (group = infile2$group), mean)

  infile_mean2 <- infile_mean %>% add_column(col1 = NA, .after = "group")
  infile_mean2 <- infile_mean2 %>% add_column(col2 = NA, .after = "group")

  colnames(infile_mean2)[2] <- colnames(infile)[1]
  colnames(infile_mean2)[3] <- colnames(infile)[2]

  infile_mean2[,2] <- sapply(strsplit(infile_mean2$group,"_"), `[`, 1)
  infile_mean2[,3] <- sapply(strsplit(infile_mean2$group,"_"), `[`, 2)

  infile_mean2[,1] <- as.factor(infile_mean2[,1])
  infile_mean2[,2] <- as.factor(infile_mean2[,2])
  infile_mean2[,3] <- as.factor(infile_mean2[,3])

  infile_mean3 <- infile_mean2[,-c(1,3)]

  pcoas <- colnames(infile_mean3[2:length(colnames(infile_mean3))])

  differences <- infile_mean3 %>% group_by(infile_mean3[1]) %>% mutate_at(pcoas, function(x) (x- shift(x)))

  differences <- as.data.frame(differences)
  differences2 <- differences[,-1]
  i <- (rowSums(differences2,na.rm=T) !=0)
  differences3 <- differences2[i,]

  unique.id <- c(unique(infile_mean2[,2]))

  #estimate theta
  table.of.diff.in.means <- differences3

  select.col <- c(1:length(colnames(differences3)))
  x.theta <- f.theta(table.of.diff.in.means, unique.id, select.col)

  cat("Angles between vectors\n")
  print(x.theta$theta.degrees)

  cat("\n")

  #estimate mean length
  x.mean_length <- f.meanL(table.of.diff.in.means, unique.id, select.col)
  cat("Mean vector lengths\n")
  print(x.mean_length[[1]])
}
