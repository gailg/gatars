% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genome_resampling_fn.R
\name{genome_resampling_fn}
\alias{genome_resampling_fn}
\title{Simulate null observations of the optimized statistics}
\usage{
genome_resampling_fn(MMM, N_simulated_nulls_interval, optimized_names, Phi,
  sampling_set, theta, WWW, y_1, y_2, statistics)
}
\arguments{
\item{MMM}{A positive integer equal to the length of \code{target_markers}
in \code{gatars_sampling_set}}

\item{N_simulated_nulls_interval}{A positive integer equal to the number
of rows of simulated optimized statistics desired.}

\item{optimized_names}{A character vector naming the optimized statistcs,
e.g. \code{c("BS", "BT", "ST", "BST")}.}

\item{Phi}{A numerical matrix of dimension \code{2} by \code{2}.
\code{Phi_{k_1, k_2} = y_{k_1} Psi y_{k_2}}.  
This  matrix is a useful intermediate calculation for getting
\code{V_z}: \code{V_z = kronecker(Phi, W_VG_W)}. 
It is of dimension \code{2} by \code{2} because there are two entities 
\code{y_1} and \code{y_2}.}

\item{sampling_set}{#' A list of \code{MMM} matrices, one matrix for each target marker.
The \code{mmm}-th matrix is the sampling set for the \code{mmm}-th
target marker and has \code{NNN} rows and up to \code{1000} columns, each 
column containing a column from \code{genotype}.  These columns
do not intersect with any of the target markers or exclusion regions
and the minor allele frequencies of the columns in \code{mmm}-th sampling set match
the minor allele frequency of the \code{mmm}-th target marker. One of the objects
returned by \code{gatars_sampling_set}.}

\item{theta}{A vector of length \code{2} that is the initial value of the
reparametrization of alpha when I am finding
minimum p-value in the full triangle \code{(alpha_B, alpha_S, alpha_T)}}

\item{WWW}{A diagonal (numerical) matrix of dimension \code{MMM} by \code{MMM}
with the diagonals equal to the \code{weights}.  (The user will specify
\code{weights} in her call to \code{gatars_test_size}.)}

\item{y_1}{A numerical vector of length \code{NNN} equal to what is referred
to in the manuscript as \code{y}, the vector of subjects' coded trait 
phenotypes.}

\item{y_2}{A numerical vector of length \code{NNN} equal to what is referred
to in the manuscript as \code{mu}, the vector of user-specified phenotype
predictions.}

\item{statistics}{A character vector reflecting which optimized statistics you
would like \code{uno_experimento_fn} to compute.}
}
\value{
A \code{list} containing the two objects
\itemize{
\item{\code{simulated}: } {
A matrix with \code{N_simulated_nulls_interval} rows and a column for each of the
optimized statistics.  Each row contains the optimized statistcs for one simulated
genotype matrix.
}

\item{\code{so_far_so_good}: } {
A logical equal to \code{TRUE} if the function successfully obtained 
\code{N_simulated_nulls_interval} rows of simulated optimized statistics.
To construct each row, only 1000 bad tries are allowed, so if 
\code{genome_resampling_fn} ever exceeds 1000, the function aborts.
}
}
}
\description{
Accumulate a matrix with \code{N_simulated_nulls_interval}
rows and a column for each of the optimized statistics.  Each row
contains the optimized statistics for one simulated genotype matrix.
The reason this function looks so complicated is there are two 
\code{while} loops.  The inside loop generates a genotype matrix
called \code{genotype_sim} and then checks to see if it has full rank.
If not, it tries again; there is a limit of 1000 tries, and it counts
the number of bad tries before it successfully gets a matrix of full rank.
The outside loop accumulates \code{N_simulated_nulls_interval} good
tries.
}
\examples{
library(Matrix)
bim = gatars_example$bim
genotype = gatars_example$genotype
phenotype = gatars_example$phenotype
Psi = gatars_example$Psi
target_markers = gatars_example$target_markers[3:5]
g_target = genotype[, target_markers]
MMM = ncol(g_target)
NNN = nrow(g_target)
e_g_target_1 = colMeans(g_target)
p_target = e_g_target_1/2
e_g_target = matrix(rep(e_g_target_1, nrow(g_target)), nrow = nrow(g_target), byrow = TRUE)
y_1 = yyy = phenotype$y
y_2 = mu = phenotype$mu
Phi = Phi_fn(Psi, y_1, y_2)
www_num = rep(1, MMM)
www = www_num/sum(www_num) * MMM
WWW = diag(www)
zzz_etc = zzz_and_first_two_moments_fn(g_target, Phi, WWW, y_1, y_2)
zzz = zzz_etc$zzz
mu_z = zzz_etc$mu_z
V_z = zzz_etc$V_z
AAA = AAA_fn(1, 0, 0, MMM)
theta_init = rep(pi/3, 2)
statistics = c("BS", "BT", "ST", "BST")
bo = basic_and_optimized_lu_fn(g_target, Phi, theta_init, WWW, y_1, y_2, statistics)
bo$xxx
bo$theta
theta = bo$theta
x_observed = bo$xxx
calculate_optimized = TRUE
epsilon = 0.01
exclusion_region = NULL
sampling_set = gatars_sampling_set(
  bim, epsilon, exclusion_region,
  genotype, hotspot, target_markers)
print(sampling_set)
sampling_set = sampling_set$sampling_set
str(sampling_set)
set.seed(1)
N_simulated_nulls_interval = 7
optimized_names = names(x_observed)
sss = genome_resampling_fn(MMM, N_simulated_nulls_interval, optimized_names, Phi, 
                           sampling_set, theta, WWW, y_1, y_2, statistics)
sss
str(sss)

}
