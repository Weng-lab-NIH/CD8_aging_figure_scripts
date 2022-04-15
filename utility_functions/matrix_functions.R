
R_removeAllZeroRows <- function(matrix) {
	print('This function can be used with Classes: matrix and sparseMatrix')
	matrix_trimmed <- matrix[which(rowSums(matrix) != 0),]
	return(matrix_trimmed)
}
