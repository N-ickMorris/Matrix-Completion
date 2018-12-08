# clear the workspace

rm(list = ls())

# ----------------------------------------------------------------------------------- 
# ---- Load Packages ----------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

{

# data handling

require(data.table)
require(gdata)

# parallel processing

require(foreach)
require(doParallel)

# truncated svd

require(irlba)
require(Matrix)

# city to city data

require(TSP)

}

# ----------------------------------------------------------------------------------- 
# ---- The SVT Algorithm ------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

{

# parameter descriptions:
	# sample.set:
		# a matrix with three columns: row index (column 1), column index (column 2), value (column 3)
		# each row of sample.set corresponds to one entry of the unknown low rank matrix
		
	# matrix.size
		# indicates the number of rows and columns of the square low rank unknown matrix
		
	# step.size
		# a scalar used to decide how to converge onto the completion of the low rank unknown matrix
			# too large of a value may constantly jump back and forth over the point of convergence
			# too low of value may take too many iterations to reach the point of convergence
			
	# tolerance
		# a scalar used to decide if the algorithm has converged
		
	# threshold
		# the value that will be used to subtract from singular values at each iteration
		
	# increment
		# a scalar used to increase the number of required singular values to compute, until we find a singula value < threshold
		
	# max.iterations
		# a scalar used to stop the algorithm
			# if max.iterations have occured then this means we did not satisfy the tolerance condition for convergence

	# the.matrix
		# if the true low rank unknown matrix is given then we compute:
			# relative error (full.relative.error)
			# relative reconstruction error (relative.error)
			# best possible error (best.error)
		# otherwise we compute:
			# relative reconstruction error (relative.error)
	
	# noise.aware
		# if set to TRUE: tolerance = 1.05 * tolerance
		
svt = function(sample.set = NULL, matrix.size = NULL, step.size = NULL, tolerance = NULL, threshold = NULL, increment = NULL, max.iterations = NULL, the.matrix = NULL, noise.aware = FALSE)
{
	require(irlba)
	
	# check to make sure sample.set, matrix.size, and max.iterations have been given properly
	
	if(missing(sample.set))
	{
		return("You must give a sample.set")
		
	} else if(missing(matrix.size))
	{
		return("You must give a matrix.size")
		
	} else if(missing(max.iterations))
	{
		return("You must give a max.iterations")
		
	} else if(length(matrix.size) > 1)
	{
		return("matrix.size must be one number")
		
	} else if(!is.matrix(sample.set))
	{
		return("sample.set must be a matrix")
		
	} else if(ncol(sample.set) != 3)
	{
		return("sample.set must have three columns")
		
	} else
	{
		time.start = Sys.time()

		# lets compute values for step.size, tolerance, threshold, and increment if they were not already given
		# these computed values are given in the SVT paper
		
		if(missing(step.size))
		{
			step.size = 1.2 * ((matrix.size^2) / nrow(sample.set))
		}
		
		if(missing(tolerance))
		{
			tolerance = 10e-4
		}
		
		if(noise.aware == TRUE)
		{
			tolerance = 1.05 * tolerance
		}
		
		if(missing(threshold))
		{
			threshold = 5 * matrix.size
		}
		
		if(missing(increment))
		{
			increment = 5
		}
		
		# lets compute P.M where M represents the low rank unknown matrix
		# P.M will be a matrix of dimension (matrix.size by matrix.size) with values filled in according to sample.set and zeros everywhere else
		
		P.M = sparseMatrix(i = sample.set[,1], j = sample.set[,2], x = sample.set[,3], dims = c(matrix.size, matrix.size))
		
		# lets compute k.0
		# k.0 represents a kick start integer value for where the matrix completion iterations should start
		
		k.0 = round(ceiling(threshold / (step.size * base::norm(P.M, "2"))), 0)
		
		# initialize a value for r
		# r represents the rank of the current matrix X
		# X is the matrix that will iteratively be updated to converge onto the desired matrix M
		# X is initially a zero matrix so the rank is initially 0
		
		r = 0
		
		# initialize matrix Y
		# Y is the matrix with the singular values of interest 
		# These singular values will have the threshold operator applied to
		# then the updated singular values and the left and right singular vectors of Y will be used to compute X
		
		Y = k.0 * step.size * P.M
		
		# run the algorithm
		
		for(k in 1:max.iterations)
		{	
			# pick the number of singular values to compute
			
			s = r + 1
			
			# lets compute enough singular values until we found all singular values > threshold
			
			done = FALSE
			it.svd = 0
			
			while(done == FALSE)
			{	
				it.svd = it.svd + 1
				
				# if we have already computed a truncated svd then restart where we left off
				# otherwise, compute a truncated svd for the first time
				
				if("S" %in% ls())
				{
					S = irlba(Y, s, v = S)
					
				} else
				{
					S = irlba(Y, s)
				}
				
				# check if the smallest singular value is less than threshold
					# if so, then we are done
					# otherwise increase the number of singular values to compute
					
				if(S$d[s] <= threshold)
				{
					done = TRUE
					
				} else
				{
					s = s + increment
				}
			}
			
			# keep track of how many iterations it took to find the truncated svd
			
			if("it.svd.history" %in% ls())
			{
				it.svd.history = c(it.svd.history, it.svd)
				
			} else
			{
				it.svd.history = it.svd
			}
			
			# reset the value of r
			
			r = length(which(S$d > threshold))
			
			# keep track of the rank of X
			
			X.rank = r
			
			if("X.rank.history" %in% ls())
			{
				X.rank.history = c(X.rank.history, X.rank)
				
			} else
			{
				X.rank.history = X.rank
			}
			
			# compute X
			
			X = S$u[,1:r] %*% diag(S$d[1:r] - threshold, nrow = r, ncol = r) %*% t(S$v[,1:r])
			
			rm(S)

			# compute the full relative error of X if needed
			
			if(!missing(the.matrix))
			{
				full.relative.error = sqrt(sum((X - the.matrix)^2)) / sqrt(sum((the.matrix)^2))
				
				if("full.relative.error.history" %in% ls())
				{
					full.relative.error.history = c(full.relative.error.history, full.relative.error)
					
				} else
				{
					full.relative.error.history = full.relative.error
				}
			}
			
			# compute the best possible error of X if needed
			
			if(!missing(the.matrix))
			{
				S.M = irlba(the.matrix, r)
				
				M.r = S.M$u %*% diag(S.M$d, nrow = r, ncol = r) %*% t(S.M$v)
				
				rm(S.M)
				
				best.error = sqrt(sum((the.matrix - M.r)^2)) / sqrt(sum((the.matrix)^2))
				
				rm(M.r)
				
				if("best.error.history" %in% ls())
				{
					best.error.history = c(best.error.history, best.error)
					
				} else
				{
					best.error.history = best.error
				}
			}
			
			# compute the relative reconstruction error of X
			
			relative.error = sqrt(sum((X[sample.set[,1:2]] - P.M[sample.set[,1:2]])^2)) / sqrt(sum((P.M[sample.set[,1:2]])^2))
			
			# keep track of how the relative error converges
			
			if("relative.error.history" %in% ls())
			{
				relative.error.history = c(relative.error.history, relative.error)
				
			} else
			{
				relative.error.history = relative.error
			}
			
			# check for convergence
			
			if(relative.error <= tolerance) break
			
			# update Y
			
			Y[sample.set[,1:2]] = Y[sample.set[,1:2]] + (step.size * (P.M[sample.set[,1:2]] - X[sample.set[,1:2]]))
		}
		
		# return the output
		
		output = list(X = X,
						duration = difftime(time2 = time.start, time1 = Sys.time(), units = "secs"),
						iterations = k, 
						relative.errors = relative.error.history,
						svd.iterations = it.svd.history,
						X.ranks = X.rank.history)

		if(!missing(the.matrix))
		{
			output$full.relative.errors = full.relative.error.history
			output$best.errors = best.error.history
		}
		
		return(output)
	}
}

}

# ----------------------------------------------------------------------------------- 
# ---- Building a Low Rank Square Matrix --------------------------------------------
# ----------------------------------------------------------------------------------- 

{

# the following function computes a n by n matrix of rank r with entries sampled from a standard normal distribution

M = function(n, r, seeds = c(42, 5))
{
	set.seed(seeds[1])
	M.L = matrix(rnorm(n = n * r), nrow = n)
	
	set.seed(seeds[2])
	M.R = matrix(rnorm(n = n * r), nrow = n)
	
	return(M.L %*% t(M.R))
}

}

# ----------------------------------------------------------------------------------- 
# ---- Sampling a Random Matrix Subset ----------------------------------------------
# ----------------------------------------------------------------------------------- 

{

# the following function randomly samples m entries of a matrix M

omega = function(M, m, seed = 21)
{
	# randomly sample m entries from M

	set.seed(seed)
	ind = sample(1:length(M), m)
	
	# compute the row and column index of each entry
	
	n = nrow(M)
	
	col.id = round(ceiling(ind / n), 0)
	
	row.id = round(n * ((ind / n) - floor(ind / n)), 0)
	row.id = ifelse(row.id == 0, n, row.id)
	
	# build output
	# output is a matrix with three columns
		# column 1 contains row indicies
		# column 2 contains column indicies
		# column 3 contains the randomly sampled entries from M
	
	output = cbind(row.id, col.id)
	output = cbind(output, M[output])
	
	colnames(output) = c("row", "col", "value")
	
	return(output)
}

}

# ----------------------------------------------------------------------------------- 
# ---- Plot Figure 5.2 --------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

path = "\\\\labsvcsredirect\\ ... \\figure5-2.RData"

{

# ----  ----------------------------------------------------------------------- 

# lets 

data(USCA312)
dat = as.matrix(USCA312)

# lets 

tau = 10^7
delta = 2
m = 0.3 * length(dat)

# lets 

S3 = irlba(dat, 3)

M3 = S3$u %*% diag(S3$d, nrow = 3, ncol = 3) %*% t(S3$v)

# randomly sample entries

omega3 = omega(M = M3, m = m, seed = 21)

# complete the matrix

result1a = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 2000, threshold = tau, step.size = delta, the.matrix = dat)
result1b = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 2000, threshold = tau, step.size = delta, the.matrix = dat, noise.aware = TRUE)

save.image(path)

# lets plot figure 5.2

svt.fre = result1a$full.relative.errors
svt.re = result1a$relative.errors
svt.be = result1a$best.errors

nav.fre = result1b$full.relative.errors
nav.re = result1b$relative.errors
nav.be = result1b$best.errors

svt.rank = result1a$X.ranks
nav.rank = result1b$X.ranks

font = 3

par(mfrow = c(1, 3), mar = c(4, 8, 6, 4))

plot.svt = plot(svt.fre, 
				main = "SVT Convergence\n of Relative Error", 
				xlab = "Iteration Number", 
				ylab = "Error \n ", ylim = c(0, 0.82), xlim = c(1, 1569),
				las = 1,
				type = "l", 
				col = "blue",
				cex.main = font,
				cex.lab = font - 0.5,
				cex.axis = font - 0.5,
				cex.sub = font - 0.5,
				tcl = -0.25,
				lwd = 3,
				lty = 2)

plot.svt				
lines(svt.re, col = "red", lty = 2, lwd = 3)
lines(svt.be, col = "black", lwd = 3)
legend(x = 300, 
		y = 0.75, 
		legend = c("Relative Error", "Relative Residual Error", "Best Possible Relative Error"), 
		lty = c(2, 2, 1), 
		lwd = c(3, 3, 3),
		cex = 1.5,
		col = c("blue", "red", "black"))

plot.nav = plot(nav.fre, 
				main = "NAV Convergence\n of Relative Error", 
				xlab = "Iteration Number", 
				ylab = "Error \n ", ylim = c(0, 0.82), xlim = c(1, 1569),
				las = 1,
				type = "l", 
				col = "blue",
				cex.main = font,
				cex.lab = font - 0.5,
				cex.axis = font - 0.5,
				cex.sub = font - 0.5,
				tcl = -0.25,
				lwd = 3,
				lty = 2)
				
plot.nav
lines(nav.re, col = "red", lty = 2, lwd = 3)
lines(nav.be, col = "black", lwd = 3)
legend(x = 300, 
		y = 0.75, 
		legend = c("Relative Error", "Relative Residual Error", "Best Possible Relative Error"), 
		lty = c(2, 2, 1), 
		lwd = c(3, 3, 3),
		cex = 1.5,
		col = c("blue", "red", "black"))

plot.rank = plot(svt.rank, 
				main = "SVT & NAV Convergence\n of Matrix Rank", 
				xlab = "Iteration Number", 
				ylab = "Rank \n ", xlim = c(1, 1569), 
				las = 1,
				type = "l", 
				col = "blue",
				cex.main = font,
				cex.lab = font - 0.5,
				cex.axis = font - 0.5,
				cex.sub = font - 0.5,
				tcl = -0.25,
				lwd = 3,
				lty = 1)
				
plot.rank
lines(nav.rank, col = "red", lty = 2, lwd = 3)
legend(x = 1000, 
		y = 1.5, 
		legend = c("SVT", "NAV"), 
		lty = c(1, 2), 
		lwd = c(3, 3),
		cex = 1.5,
		col = c("blue", "red"))

keep(plot.svt, plot.nav, plot.rank, path, 
		svt.fre, svt.re, svt.be,
		nav.fre, nav.re, nav.be,
		svt.rank, nav.rank, sure = TRUE)

save.image(path)

}



