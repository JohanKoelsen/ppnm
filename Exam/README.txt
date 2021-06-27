EXAM QUESTION 8
One-sided Jacobi algorithm for Singular Value Decomposition

By Johan Kølsen de Wit, 201709130 - mod(30,22) = 8.

-----------OVERVIEW------------
main.c	the main function with int main() included. It contains matrix generation of A, V, U, D, runs the SVD and prints the results.

SVD.c	makes the one-sided Jacobi SVD algorithm.


---------The Project---------
In this exercise the one-sided Jacobi algorithm is used on a real, square matrix, A, which decomposes A into A-> UDV^T (with V^T being the transposed of matrix V).
Matrix D is diagonal with non-negative elements, and matrices U and V are orthogonal.

The method uses a cyclic sweep-method until the iterations have converged when the diagonal elements does not change under further rotations.
Then U, D, and V can be computed:
V = R, Dii=||a'i|| and ui=a'i/||a'i||.

The exercise is reminiscent of that of the "matrix diagonalization" homework exercise, where we used Jacobi algorithm for matrix diagonalization.
Only few elements had to altered in order for the SVD to work. The most important difference was the one-sided nature of SVD and that A was not constrained to a symmetric matrix.

The two conditions for A was that it had to be a real matrix with dimensions m, n with m = n.
The condition that m = n could be bypassed by first use QR decomposition of A -> QR, where Q is the n × m orthogonal matrix and R is a square triangular m × m
matrix.
Then we could use SVD on R, yielding R = U'DV^T. The SVD on A is given as A = UDV^T with U = QU'.


The output text of the program is written in out.txt.
