# 'fuzzy.match' and 'fuzzy.match2' pre-dated nearest, etc.

# For query and subject both width 1 ranges, fuzzy.match
# returned the indices of the nearest subject ranges, in
# case of tie choosing the subject *following* the query.

# The subject for fuzzy.match2 is taken without validation
# to specify a "set of exons", as either a 2-column matrix,
# hence the function name, or a sorted and reduced IRanges
# object.  There are two nearest exons exactly when the query
# lies precisely midway between them.  In that case, the
# second exon is deemed the winner.

# The typical use-case is by matchGenes where the query x
# is a single integer, although that might get vectorized
# some day.

# Returned: a length(x)-by-2 matrix consisting of signed
# distances, whereby "x > s" gets -distance(x,s), and the
# nearest subject indices.

# Note: this is strand unaware!
fuzzy.match2 <- function(x, M) {

	# x: integer vector or an IRanges of constant width 1

	# M: 2-column matrix representing a set of exons, i.e.
	# a sorted list of disjoint intervals, or a sorted and
	# reduced IRanges object

	X <- if (is(x, "IRanges")) {
		if (any(width(x) != 1))
			stop("'x' must have constant width 1")	
		x
	     }
	     else
		IRanges(start=x, end=x)

	Y <- if (is(M, "IRanges"))
		M
	     else
		IRanges(start=M[,1], end=M[,2])

	# ?distance says:
	# In Bioconductor >=2.12 the distance calculation has been
        # changed to accommodate zero-width ranges in a consistent and
        # intuitive manner.  Because of this change, a warning will be
        # emitted when 'distance' is called. This warning is temporary
        # and will be removed in Bioconductor 2.13. To suppress the
        # warning, code can be wrapped in 'suppressWarnings()'.
	D <- suppressWarnings(distanceToNearest(X, Y, select="all"))

	# There are two cases:
	# 1) x[i] is midway between 2 subject ranges <==>
	#  there are 2 nearest's in D, and we want the 2nd,
	# when select="arbitrary" would have given the 1st.
	# 2) x[i] has a unique nearest range

	if (nrow(D) > length(X))	# double hits
		D <- D[cumsum(rle(D$queryHits)$lengths),]

	n <- D$subjectHits	# match indices
	N <- Y[n]		# nearest ranges
	d <- D$distance		# unsigned distance (see above)

	following <- X>N
	preceding <- X<N

	# convert back to real distances, if necessary
	if (suppressWarnings(distance(IRanges(1,1), IRanges(2,2))) == 0) {
		d[d>0] <- d[d>0] + 1L
		d[d==0 & (following|preceding)] <- 1L
	}

	# queries past their nearest exon get negative distances
	sign <- ifelse(following, -1L, 1L)
	d <- d * sign
	cbind(distance=d, nearest=n, row.names=NULL)
}
