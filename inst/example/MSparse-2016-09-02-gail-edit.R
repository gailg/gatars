MSparse <- function(k)
{
#   for optimally assigning chromosome segments to sampling sets.  GI. 8/30/2016
#   input: matrix k(nc, m)
#   output: optimal solution vector x_opt with dimension nc*m + m + 1
#
#   generates the CA optimization problem for Gurobi in sparse matrix notation.
#   Only columns for which k_cm > 0 are generated.
#   Then solves the mixed integer lineare programming problem using Gurobi MIP.
#   x_opt[1:(nc*m)] optimal assignments, x_opt[(nc*m+1):(nc*m+m)] optimum number of markers for each sampling set m
#   x_opt[nc*m + m + 1] optimal objective
#
#   Uses the gurobi package and the Matrix package.

nc <- nrow(k)
m <- ncol(k)

npk <- sum(k>0)
nnz <- npk*2 + 3*m


A <- numeric(nnz)
irow <- numeric(nnz)+1
jcol <- numeric(nnz)+1
jind <- numeric(npk + m + 1)

ic=0
for (i in 1:nc){


    ip <- which(k[i,] > 0)
    icount = sum(ip>0)
    i1 = ic + 1
    ic = ic + icount
    jind[i1:ic] <- ip + (i-1)*m


    A[i1:ic] <- 1
    irow[i1:ic] <- i
    jcol[i1:ic] <- (1:icount) +  i1 - 1

    i2 <- i1 + npk
    i3 <- ic + npk
    A[i2:i3] <- k[i, ip]
    irow[i2:i3] <- nc + ip
    jcol[i2:i3] <- (1:icount) +  i1 - 1

}

n1 = npk*2
A[(n1+1):(n1+m)] <- -1
irow[(n1+1):(n1+m)] <- (1:m) + nc
jcol[(n1+1):(n1+m)] <- (1:m) + npk

n1 = npk*2 + m
A[(n1+1):(n1+m)] <- 1
irow[(n1+1):(n1+m)] <- (1:m) + nc + m
jcol[(n1+1):(n1+m)] <- (1:m) + npk

n1 = npk*2 + 2*m
A[(n1+1):(n1+m)] <- -1
irow[(n1+1):(n1+m)] <- (1:m) + nc + m
jcol[(n1+1):(n1+m)] <- npk + m + 1

jind[(npk+1):(npk+m+1)] <- (1:(m+1)) + nc*m

mrows <- nc + 2*m
ncols <- npk + m + 1

CA <- list()

CA$A          <- spMatrix(mrows, ncols, irow, jcol, A)
CA$obj        <- c(rep(0, ncols-1), 1)
CA$modelsense <- "max"
CA$rhs        <- c(rep(1, nc), rep(0, 2*m))
CA$sense      <- c(rep("=", nc+m), rep(">=", m))
CA$vtype      <- c(rep("B", npk), rep("C", m+1))
CA$jind       <- jind

result <- list()

result <- gurobi(CA)

result
#print(result$objval)

x_opt <- numeric(nc*m+m+1)
x_opt[jind] <- result$x

#print(result$x)
#print(x_opt)
list(CA = CA,
     x_opt = x_opt)
}
