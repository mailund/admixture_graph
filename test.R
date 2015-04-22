
parents <- matrix(FALSE, 8, 8)
rownames(parents) <- colnames(parents) <- c("A", "B", "C", "AB", "BC", "ABC", "R", "O")
parents["A","AB"] <- TRUE
parents["C","BC"] <- TRUE
parents["B","AB"] <- TRUE
parents["B","BC"] <- TRUE
parents["AB","ABC"] <- TRUE
parents["BC","ABC"] <- TRUE
parents["ABC","R"] <- TRUE
parents["O","R"] <- TRUE

weights <- matrix("", 8, 8)
rownames(weights) <- colnames(weights) <- c("A", "B", "C", "AB", "BC", "ABC", "R", "O")
weights["B","AB"] <- weights["AB","B"] <- 'a'
weights["B","BC"] <- weights["BC","B"] <- '(1-a)'


graph <- agraph(parents, weights)
f4(graph, "O", "A", "B", "C")
f3(graph, "O", "A", "C")
f3(graph, "O", "A", "B")