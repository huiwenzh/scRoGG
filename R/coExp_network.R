#' Construct coexpression network and sub-community detection.
#' @author Huiwen Zheng
#' @param edge_dat An edge list dataframe. Recommend to use the RS score generate from robustness function but can be apply to any correlation measurements.
#' @param min_size Minimum size of networks to include in the communty detection.
#' @param n_networks Number of sub-communities to visualise.
#' @usage coExp_network(dat1, min_size=2, n_network=2)
#' @export
coExp_network <- function(edge_dat,min_size=1,n_networks=4){
  #min_max <- function(z){(z-min(z))/(max(z)-min(z))}
  colnames(edge_dat) <- c('from','to','weight')
  #edge_dat$weight <- sign(edge_dat$weight)*min_max(abs(edge_dat$weight))
  edge_dat$weight <- sign(edge_dat$weight)*sqrt(abs(edge_dat$weight))
  G <- igraph::graph_from_data_frame(edge_dat,directed = F)
  G_test <- igraph::induced_subgraph(
    G,   igraph::V(G)[ave(1:igraph::vcount(G), igraph::membership(igraph::components(G)), FUN = length) > min_size])

  m_all <- c()
  res <- c(seq(0.0005,0.5, by=0.005))
  for (i in res){
    c1 = igraph::cluster_leiden(G_test,objective_function ="CPM",resolution_parameter =res[which.max(m_all)] ,n_iterations =100)
    m_all = c(m_all,igraph::modularity(G_test, c1$membership))

  }
  c1 <- igraph::cluster_leiden(G_test,objective_function ="CPM",resolution_parameter =res[which.max(m_all)]  ,n_iterations =100)

  graph_list <- list()
  graph_list[["all"]] <- G_test
  for (n in 1:n_networks){
    selnodes <- igraph::V(G_test)[c1$membership==names(table(c1$membership)[order(table(c1$membership),decreasing = T)])[n]]
    selegoG <- igraph::induced_subgraph(G_test,selnodes)
    igraph::edge_density(selegoG, loops=F)
    selegoG1 <- selegoG
    igraph::E(selegoG1)$weight <- abs(igraph::E(selegoG)$weight)
    l <- igraph::layout.fruchterman.reingold(selegoG1)
    plot(selegoG, layout=l,    # === vertex
         vertex.color =rgb(0.6,0.6,0.6,0.6),          # Node color
         vertex.frame.color = "white",                 # Node border color
         vertex.shape="circle",                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
         vertex.size=13,                               # Size of the node (default is 15)
         vertex.size2=NA,                              # The second size of the node (e.g. for a rectangle)
         # Character vector used to label the nodes
         vertex.label.color="black",

         vertex.label.font=1,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
         vertex.label.cex=0.7,                           # Font size (multiplication factor, device-dependent)
         vertex.label.dist=0,                          # Distance between the label and the vertex
         vertex.label.degree=0 ,                       # The position of the label in relation to the vertex (use pi)

         # === Edge
         edge.color="#cc6600",                           # Edge color
         edge.width=E(selegoG)$weight*10,                                 # Edge width, defaults to 1
         edge.arrow.size=1,                            # Arrow size, defaults to 1
         edge.lty=1,
         edge.curved=0  ,                          # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
    )
    graph_list[[paste0("sub_",n)]] <- selegoG
  }
  graph_list

}
