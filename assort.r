assortativity <- function(network){

       diag(network) = 0
       node_degree=rowSums(network)}
        
	# Weighted network
        	indices = which(!is.na(network), arr.ind=T )
                w = network[indices]/ sum(network[indices])
		x =  node_degree[indices[,1]] - sum( node_degree[indices[,1]] * w )
		y =  node_degree[indices[,2]] - sum( node_degree[indices[,2]] * w )
		vx = sum( w * x * x )
		vy = sum( w * y * y )
		vxy = sum(y * x * w)
                assort = vxy / sqrt(vx * vy)

	return(assort) 
}
