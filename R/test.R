
alpha <- 0.1
parents <- matrix(0.0, 6,6)
rownames(parents) <- colnames(parents) <- c("A", "B", "C", "AB", "BC", "R")
parents["A","AB"] <- 1.0
parents["C","BC"] <- 1.0
parents["B","AB"] <- alpha
parents["B","BC"] <- 1.0 - alpha
parents["AB","R"] <- 1.0
parents["BC","R"] <- 1.0

parents