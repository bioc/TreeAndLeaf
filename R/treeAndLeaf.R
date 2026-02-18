#' Layout a TreeAndLeaf diagram.
#'
#' This function tranforms hclust and phylo objects into tree-and-leaf 
#' igraph objects.
#'
#' @param obj An object of class 'hclust' or 'phylo'.
#' 
#' @return A tree-and-leaf igraph object.
#'
#' @seealso \code{\link{formatTree}}
#' @seealso \code{\link[stats:hclust]{hclust}}
#' @seealso \code{\link[ape:as.phylo]{as.phylo}}
#' @seealso \code{\link[RedeR:addGraph]{addGraph}}
#' @seealso \code{\link[RedeR:relax]{relax}}
#'
#' @examples
#' library(RedeR)
#' hc <- hclust(dist(USArrests), "ave")
#' tal <- treeAndLeaf(hc)
#' 
#' \dontrun{
#' startRedeR()
#' addGraphToRedeR(tal)
#' }
#'
#' @importFrom igraph edge_betweenness remove.edge.attribute is.igraph V E 'V<-' 'E<-' 
#' @importFrom ape as.phylo as.igraph.phylo
#' @importFrom RedeR RedPort
#' @export

treeAndLeaf <- function(obj){
    tal.checks(name="obj", para=obj)
    tal <- .setTaL(obj)
    return(tal)
}
#-------------------------------------------------------------------------------
.setTaL <- function(obj){
    if ("phylo" %in% class(obj)){
        gg <- ape::as.igraph.phylo(obj, directed=FALSE)
        coords <- .get.tree.coords(obj, use.edge.length=FALSE)
        phylo <- obj
    } else if("ggtree" %in% class(obj)){
        lpar <- obj$layers[[1]]$stat_params$layout
        vpar <- c("daylight","ape","fan","equal_angle")
        if(!(lpar %in% vpar) ){
            stop("Please, use one of the 'ggtree' layouts: ",
                paste0(vpar, collapse = ", "),"!", call.=FALSE)
        }
        phylo <- ape::as.phylo(obj$data)
        gg <- ape::as.igraph.phylo(phylo, directed=FALSE)
        coords <- .get.tree.coords(phylo, use.edge.length=FALSE)
        coords$layout[,] <- as.matrix(obj$data[,c("x","y")])
    } else {
        phylo <- ape::as.phylo(obj)
        gg <- ape::as.igraph.phylo(phylo, directed=FALSE)
        coords <- .get.tree.coords(phylo, use.edge.length=FALSE)
    }
    if(!is.null(phylo$edge.length)) E(gg)$edgeLength <- phylo$edge.length
    if(!all(V(gg)$name%in%rownames(coords$layout)))
        stop("An unexpected error occurred while transforming the graph object!")
    #--- set layout
    layout <- coords$layout[V(gg)$name,]
    V(gg)$x <- layout[,"x"]
    V(gg)$y <- layout[,"y"]
    V(gg)$isLeaf <- V(gg)$name%in%coords$phylo$tip.label
    gg$centralVertex <- coords$centralVertex
    #--- set edge weights
    if(!is.null(E(gg)$weight)) gg <- remove.edge.attribute(gg, "weight")
    if(!is.null(E(gg)$edgeWeight)) gg <- remove.edge.attribute(gg, "edgeWeight")
    bt <- igraph::edge_betweenness(gg, directed = T)
    bt <- (1 - bt/max(bt))
    E(gg)$edgeWeight <- bt
    #--- set node and font sizes
    V(gg)$nodeLabel <- V(gg)$name
    V(gg)$nodeSize <- 5
    V(gg)$nodeSize[V(gg)$isLeaf] <- 30
    V(gg)$nodeLabelSize <- 15
    V(gg)$nodeLabelSize[!V(gg)$isLeaf] <- 1
    V(gg)$nodeColor <- "#ffcccc"
    V(gg)$nodeLineColor <- "#9999ff"
    V(gg)$nodeLineWidth <- 2
    E(gg)$edgeColor <- "#9999ff"
    E(gg)$edgeWidth <- 2
    gg$gtype <- "TreeAndLeaf"
    class(gg) <- c("tal","igraph")
    return(gg)
}

#-------------------------------------------------------------------------------
# This function gets coordinates from an unrooted 'phylo' objects,
# derived from a currectly not exported function from the ape package
.get.tree.coords <- function(phy, use.edge.length = FALSE){
    if(is.null(phy$edge.length)) use.edge.length <- FALSE
    phy <- ape::reorder.phylo(phy)
    Ntip <- length(phy$tip.label)
    Nedge <- dim(phy$edge)[1]
    Nnode <- phy$Nnode
    nb.sp <- ape::node.depth(phy)
    if(use.edge.length){
        XY <- .unrooted.xy(Ntip, Nnode, nb.sp, edge=phy$edge, 
                           edge.length=phy$edge.length)
    } else {
        XY <- .unrooted.xy(Ntip, Nnode, nb.sp, edge=phy$edge, 
                           edge.length=rep(1, Nedge))
    }
    xx <- XY[, 1] - min(XY[, 1])
    yy <- XY[, 2] - min(XY[, 2])
    layout <- cbind(x=xx,y=yy)
    if(is.null(phy$node.label)) phy <- ape::makeNodeLabel(phy)
    if(anyDuplicated(c(phy$tip.label, phy$node.label))) 
        stop("Duplicated labels!")
    rownames(layout) <- c(phy$tip.label,phy$node.label)
    centralVertex <- phy$node.label[1]
    return(list(phylo=phy, layout=layout, centralVertex=centralVertex))
}

#-------------------------------------------------------------------------------
# This function gets coordinates from an unrooted 'phylo' object; it is
# derived from a currectly not exported function from the ape package.
.unrooted.xy <- function(Ntip, Nnode, nb.sp, edge, edge.length){
    #--- Get nodes' parent/sons ordering
    # nodeseq <- function(node) {
    #     ind <- which(edge[, 1] == node)
    #     sons <- edge[ind, 2]
    #     nodes <<- c(nodes, node)
    #     for (i in sons) if (i > Ntip) nodeseq(i)
    # }
    # nodes <- NULL
    # nodeseq(node=Ntip + 1L)
    #--- Get nodes' parent/sons ordering (for large dendrograms)
    node <- Ntip + 1L
    nodes <- NULL
    nstack <- NULL
    flag <- TRUE
    while(flag){
        nodes <- c(nodes, node)
        ind <- which(edge[, 1] == node)
        sons <- edge[ind, 2]
        sons <- sons[!sons%in%nodes]
        sons <- sons[sons > Ntip]
        ns <- length(nstack)
        if(length(sons)>0){
            if(length(sons)>1)
                nstack <- c(nstack, rev(sons[-1]))
            node <- sons[1]
        } else if(ns > 0){
            node <- nstack[ns]
            nstack <- nstack[-ns]
        } else {
            flag <- FALSE
        }
    }
    #--- Gets coordinates
    yy <- xx <- numeric(Ntip + Nnode)
    axis <- angle <- numeric(Ntip + Nnode)
    angle[nodes[1]] <- 2*pi
    for(i in 1:length(nodes)){
        node <- nodes[i]
        ind <- which(edge[, 1] == node)
        sons <- edge[ind, 2]
        start <- axis[node] - angle[node]/2
        for (j in 1:length(sons)) {
            h <- edge.length[ind[j]]
            angle[sons[j]] <- alpha <- angle[node]*nb.sp[sons[j]]/nb.sp[node]
            axis[sons[j]] <- beta <- start + alpha/2
            start <- start + alpha
            xx[sons[j]] <- h*cos(beta) + xx[node]
            yy[sons[j]] <- h*sin(beta) + yy[node]
        }
    }
    M <- cbind(xx, yy)
    return(M)
}

