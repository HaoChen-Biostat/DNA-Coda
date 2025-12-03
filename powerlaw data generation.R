#power-law degree distribution 

generate_networks <- function(p = 20, gamma0 = 2, s = p*0.2, rho=0.2, r=0.05) {
  
  while (TRUE) {
    # Create a scale-free network following a power-law distribution
    g <- sample_pa(p, gamma0, directed = FALSE)
    adj_matrix <- as_adjacency_matrix(g)
    regular_matrix <- as.matrix(adj_matrix)
    
    # Generate random numbers from a uniform distribution
    a1 <- regular_matrix * sample(c(-1, 1), p * p, replace = TRUE) * runif(p * p, 0.5, 0.7)
    a1 <- a1 / 2##each row was divided by 2 when p = 20, 3 when p = 40, 4 when p = 60/80, 5 when p = 100.
    diag(a1) <- 0
    a1 <- (a1 + t(a1)) / 2
    
    # Calculate the degree of nodes and find the two nodes with the highest degree
    deg <- degree(g)
    top_hubs <- which(deg %in% head(sort(deg, decreasing = TRUE), 2))
    change_edge1 <- order(abs(a1[top_hubs[1], ]), decreasing = TRUE)[1:s]
    change_edge2 <- order(abs(a1[top_hubs[2], ]), decreasing = TRUE)[1:s]
    
    a2 <- a1
    a2[top_hubs[1], change_edge1] <- -1 * a2[top_hubs[1], change_edge1]
    a2[top_hubs[2], change_edge2] <- -1 * a2[top_hubs[2], change_edge2]
    a2[change_edge1, top_hubs[1]] <- -1 * a2[change_edge1, top_hubs[1]]
    a2[change_edge2, top_hubs[2]] <- -1 * a2[change_edge2, top_hubs[2]]
    
    diag(a1) <- 1
    diag(a2) <- 1
    
    diff_network <- a2 - a1
    if (sum(diff_network != 0) == 2 * s) {  # Check the number of non-zero elements
      break
    }
  }
  return(list(a1 = a1, a2 = a2, 
              diff_network = diff_network))
}
