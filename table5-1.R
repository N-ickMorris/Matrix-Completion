# clear the workspace

rm(list = ls())

# ----------------------------------------------------------------------------------- 
# ---- Load Packages ----------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

{

# data handling

require(data.table)
require(gdata)

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
# ---- Build Table 5.1 --------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

{

# ----------------------------------------------------------------------------------- 
# ---- Build n1000r10 ---------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n1000r10.RData"

{

# 1000 x 1000 rank 10 matrix

M1 = M(n = 1000, r = 10, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.12 * (1000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 10 matrix

M2 = M(n = 1000, r = 10, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.12 * (1000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 10 matrix

M3 = M(n = 1000, r = 10, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.12 * (1000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 10 matrix

M4 = M(n = 1000, r = 10, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.12 * (1000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 10 matrix

M5 = M(n = 1000, r = 10, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.12 * (1000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n1000r50 ---------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n1000r50.RData"

{

# 1000 x 1000 rank 50 matrix

M1 = M(n = 1000, r = 50, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.39 * (1000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 50 matrix

M2 = M(n = 1000, r = 50, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.39 * (1000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 50 matrix

M3 = M(n = 1000, r = 50, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.39 * (1000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 50 matrix

M4 = M(n = 1000, r = 50, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.39 * (1000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 50 matrix

M5 = M(n = 1000, r = 50, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.39 * (1000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)
}

# ----------------------------------------------------------------------------------- 
# ---- Build n1000r100 --------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n1000r100.RData"

{

# 1000 x 1000 rank 100 matrix

M1 = M(n = 1000, r = 100, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.57 * (1000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 100 matrix

M2 = M(n = 1000, r = 100, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.57 * (1000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 100 matrix

M3 = M(n = 1000, r = 100, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.57 * (1000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 100 matrix

M4 = M(n = 1000, r = 100, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.57 * (1000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 1000 x 1000 rank 100 matrix

M5 = M(n = 1000, r = 100, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.57 * (1000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n5000r10 ---------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n5000r10.RData"

{

# 5000 x 5000 rank 10 matrix

M1 = M(n = 5000, r = 10, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.024 * (5000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 10 matrix

M2 = M(n = 5000, r = 10, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.024 * (5000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 10 matrix

M3 = M(n = 5000, r = 10, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.024 * (5000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 10 matrix

M4 = M(n = 5000, r = 10, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.024 * (5000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 10 matrix

M5 = M(n = 5000, r = 10, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.024 * (5000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n5000r50 ---------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n5000r50.RData"

{

# 5000 x 5000 rank 50 matrix

M1 = M(n = 5000, r = 50, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.1 * (5000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 50 matrix

M2 = M(n = 5000, r = 50, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.1 * (5000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 50 matrix

M3 = M(n = 5000, r = 50, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.1 * (5000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 50 matrix

M4 = M(n = 5000, r = 50, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.1 * (5000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 50 matrix

M5 = M(n = 5000, r = 50, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.1 * (5000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)
}

# ----------------------------------------------------------------------------------- 
# ---- Build n5000r100 --------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n5000r100.RData"

{

# 5000 x 5000 rank 100 matrix

M1 = M(n = 5000, r = 100, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.158 * (5000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 100 matrix

M2 = M(n = 5000, r = 100, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.158 * (5000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 100 matrix

M3 = M(n = 5000, r = 100, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.158 * (5000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 100 matrix

M4 = M(n = 5000, r = 100, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.158 * (5000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 5000 x 5000 rank 100 matrix

M5 = M(n = 5000, r = 100, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.158 * (5000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n10000r10 --------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n10000r10.RData"

{

# 10000 x 10000 rank 10 matrix

M1 = M(n = 10000, r = 10, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.012 * (10000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 10 matrix

M2 = M(n = 10000, r = 10, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.012 * (10000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 10 matrix

M3 = M(n = 10000, r = 10, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.012 * (10000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 10 matrix

M4 = M(n = 10000, r = 10, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.012 * (10000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 10 matrix

M5 = M(n = 10000, r = 10, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.012 * (10000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n10000r50 --------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n10000r50.RData"

{

# 10000 x 10000 rank 50 matrix

M1 = M(n = 10000, r = 50, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.05 * (10000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 50 matrix

M2 = M(n = 10000, r = 50, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.05 * (10000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 50 matrix

M3 = M(n = 10000, r = 50, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.05 * (10000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 50 matrix

M4 = M(n = 10000, r = 50, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.05 * (10000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 50 matrix

M5 = M(n = 10000, r = 50, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.05 * (10000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)
}

# ----------------------------------------------------------------------------------- 
# ---- Build n10000r100 -------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n10000r100.RData"

{

# 10000 x 10000 rank 100 matrix

M1 = M(n = 10000, r = 100, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.08 * (10000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 100 matrix

M2 = M(n = 10000, r = 100, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.08 * (10000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 100 matrix

M3 = M(n = 10000, r = 100, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.08 * (10000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 100 matrix

M4 = M(n = 10000, r = 100, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.08 * (10000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 10000 x 10000 rank 100 matrix

M5 = M(n = 10000, r = 100, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.08 * (10000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n20000r10 --------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n20000r10.RData"

{

# 20000 x 20000 rank 10 matrix

M1 = M(n = 20000, r = 10, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.006 * (20000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 20000 x 20000 rank 10 matrix

M2 = M(n = 20000, r = 10, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.006 * (20000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 20000 x 20000 rank 10 matrix

M3 = M(n = 20000, r = 10, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.006 * (20000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 20000 x 20000 rank 10 matrix

M4 = M(n = 20000, r = 10, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.006 * (20000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 20000 x 20000 rank 10 matrix

M5 = M(n = 20000, r = 10, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.006 * (20000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n20000r50 --------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n20000r50.RData"

{

# 20000 x 20000 rank 50 matrix

M1 = M(n = 20000, r = 50, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.025 * (20000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 20000 x 20000 rank 50 matrix

M2 = M(n = 20000, r = 50, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.025 * (20000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 20000 x 20000 rank 50 matrix

M3 = M(n = 20000, r = 50, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.025 * (20000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 20000 x 20000 rank 50 matrix

M4 = M(n = 20000, r = 50, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.025 * (20000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 20000 x 20000 rank 50 matrix

M5 = M(n = 20000, r = 50, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.025 * (20000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)
}


# ----------------------------------------------------------------------------------- 
# ---- Build n30000r10 --------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# pick a location and file name to save results

path = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n30000r10.RData"

{

# 30000 x 30000 rank 10 matrix

M1 = M(n = 30000, r = 10, seeds = c(42, 5))

# randomly sample entries

omega1 = omega(M = M1, m = 0.004 * (30000^2), seed = 21)

# complete the matrix

result1 = svt(sample.set = omega1, matrix.size = nrow(M1), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 30000 x 30000 rank 10 matrix

M2 = M(n = 30000, r = 10, seeds = c(43, 6))

# randomly sample entries

omega2 = omega(M = M2, m = 0.004 * (30000^2), seed = 22)

# complete the matrix

result2 = svt(sample.set = omega2, matrix.size = nrow(M2), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 30000 x 30000 rank 10 matrix

M3 = M(n = 30000, r = 10, seeds = c(44, 7))

# randomly sample entries

omega3 = omega(M = M3, m = 0.004 * (30000^2), seed = 23)

# complete the matrix

result3 = svt(sample.set = omega3, matrix.size = nrow(M3), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 30000 x 30000 rank 10 matrix

M4 = M(n = 30000, r = 10, seeds = c(45, 8))

# randomly sample entries

omega4 = omega(M = M4, m = 0.004 * (30000^2), seed = 24)

# complete the matrix

result4 = svt(sample.set = omega4, matrix.size = nrow(M4), max.iterations = 200)

save.image(path)

# ----------------------------------------------------------------------------------- 

# 30000 x 30000 rank 10 matrix

M5 = M(n = 30000, r = 10, seeds = c(46, 9))

# randomly sample entries

omega5 = omega(M = M5, m = 0.004 * (30000^2), seed = 25)

# complete the matrix

result5 = svt(sample.set = omega5, matrix.size = nrow(M5), max.iterations = 200)

save.image(path)

rm(M1, M2, M3, M4, M5,
	omega1, omega2, omega3, omega4, omega5,
	result1, result2, result3, result4, result5)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n1000 Rows -------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# 

rm(list = ls())

# 

path.table = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\table5-1.RData"
path10 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n1000r10.RData"
path50 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n1000r50.RData"
path100 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n1000r100.RData"

{

# ----------------------------------------------------------------------------------- 

# 

load(path10)

# 

tab5.1 = data.table(n = nrow(M1), 
					r = 10, 
					m.dr = 6,
					m.n2 = 0.12,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration), 
										as.numeric(result3$duration), 
										as.numeric(result4$duration), 
										as.numeric(result5$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations, 
										result3$iterations, 
										result4$iterations, 
										result5$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1), 
											tail(result3$relative.errors, 1), 
											tail(result4$relative.errors, 1), 
											tail(result5$relative.errors, 1))))

# 

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

# ----------------------------------------------------------------------------------- 

# 

load(path50)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 50, 
					m.dr = 4,
					m.n2 = 0.39,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration), 
										as.numeric(result3$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations, 
										result3$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1), 
											tail(result3$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

# ----------------------------------------------------------------------------------- 

# 

load(path100)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 100, 
					m.dr = 3,
					m.n2 = 0.57,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration), 
										as.numeric(result3$duration), 
										as.numeric(result4$duration), 
										as.numeric(result5$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations, 
										result3$iterations, 
										result4$iterations, 
										result5$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1), 
											tail(result3$relative.errors, 1), 
											tail(result4$relative.errors, 1), 
											tail(result5$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n5000 Rows -------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# 

path10 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n5000r10.RData"
path50 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n5000r50.RData"
path100 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n5000r100.RData"

{

# ----------------------------------------------------------------------------------- 

# 

load(path10)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 10, 
					m.dr = 6,
					m.n2 = 0.024,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration), 
										as.numeric(result3$duration), 
										as.numeric(result4$duration), 
										as.numeric(result5$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations, 
										result3$iterations, 
										result4$iterations, 
										result5$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1), 
											tail(result3$relative.errors, 1), 
											tail(result4$relative.errors, 1), 
											tail(result5$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

# ----------------------------------------------------------------------------------- 

# 

load(path50)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 50, 
					m.dr = 5,
					m.n2 = 0.10,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration), 
										as.numeric(result3$duration), 
										as.numeric(result4$duration), 
										as.numeric(result5$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations, 
										result3$iterations, 
										result4$iterations, 
										result5$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1), 
											tail(result3$relative.errors, 1), 
											tail(result4$relative.errors, 1), 
											tail(result5$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

# ----------------------------------------------------------------------------------- 

# 

load(path100)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 100, 
					m.dr = 4,
					m.n2 = 0.158,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration), 
										as.numeric(result3$duration), 
										as.numeric(result4$duration), 
										as.numeric(result5$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations, 
										result3$iterations, 
										result4$iterations, 
										result5$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1), 
											tail(result3$relative.errors, 1), 
											tail(result4$relative.errors, 1), 
											tail(result5$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n10000 Rows ------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# 

path10 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n10000r10.RData"
path50 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n10000r50.RData"
path100 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n10000r100.RData"

{

# ----------------------------------------------------------------------------------- 

# 

load(path10)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 10, 
					m.dr = 6,
					m.n2 = 0.012,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration), 
										as.numeric(result3$duration), 
										as.numeric(result4$duration), 
										as.numeric(result5$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations, 
										result3$iterations, 
										result4$iterations, 
										result5$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1), 
											tail(result3$relative.errors, 1), 
											tail(result4$relative.errors, 1), 
											tail(result5$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

# ----------------------------------------------------------------------------------- 

# 

load(path50)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 50, 
					m.dr = 5,
					m.n2 = 0.050,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration), 
										as.numeric(result3$duration), 
										as.numeric(result4$duration), 
										as.numeric(result5$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations, 
										result3$iterations, 
										result4$iterations, 
										result5$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1), 
											tail(result3$relative.errors, 1), 
											tail(result4$relative.errors, 1), 
											tail(result5$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

# ----------------------------------------------------------------------------------- 

# 

load(path100)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 100, 
					m.dr = 4,
					m.n2 = 0.080,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration), 
										as.numeric(result3$duration), 
										as.numeric(result4$duration), 
										as.numeric(result5$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations, 
										result3$iterations, 
										result4$iterations, 
										result5$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1), 
											tail(result3$relative.errors, 1), 
											tail(result4$relative.errors, 1), 
											tail(result5$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

}

# ----------------------------------------------------------------------------------- 
# ---- Build n20000 Rows ------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# 

path10 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n20000r10.RData"
path50 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n20000r50.RData"

{

# ----------------------------------------------------------------------------------- 

# 

load(path10)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 10, 
					m.dr = 6,
					m.n2 = 0.006,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

# ----------------------------------------------------------------------------------- 

# 

load(path50)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 50, 
					m.dr = 5,
					m.n2 = 0.025,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration), 
										as.numeric(result3$duration), 
										as.numeric(result4$duration), 
										as.numeric(result5$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations, 
										result3$iterations, 
										result4$iterations, 
										result5$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1), 
											tail(result3$relative.errors, 1), 
											tail(result4$relative.errors, 1), 
											tail(result5$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

tab5.1.update = data.table(n = 20000, 
					r = 50, 
					m.dr = 5,
					m.n2 = 0.025,
					time.sec = NA,
					iterations = NA,
					relative.error = NA)

tab5.1 = rbind(tab5.1, tab5.1.update)

keep(tab5.1, path.table, path10, path50, path100, sure = TRUE)

# 

}

# ----------------------------------------------------------------------------------- 
# ---- Build n30000 Row -------------------------------------------------------------
# ----------------------------------------------------------------------------------- 

# 

path10 = "\\\\labsvcsredirect\\users$\\njm2868\\Documents\\n30000r10.RData"

{

# 

load(path10)

# 

tab5.1.update = data.table(n = nrow(M1), 
					r = 10, 
					m.dr = 6, 
					m.n2 = 0.004,
					time.sec = mean(c(as.numeric(result1$duration), 
										as.numeric(result2$duration))),
										
					iterations = mean(c(result1$iterations, 
										result2$iterations)),
										
					relative.error = mean(c(tail(result1$relative.errors, 1), 
											tail(result2$relative.errors, 1))))

tab5.1 = rbind(tab5.1, tab5.1.update)

# 

tab5.1.update = data.table(n = 30000, 
					r = 10, 
					m.dr = 6, 
					m.n2 = 0.004,
					time.sec = NA,
					iterations = NA,
					relative.error = NA)

tab5.1 = rbind(tab5.1, tab5.1.update)

keep(tab5.1, path.table, sure = TRUE)

}

save.image(path.table)

}

tab5.1


