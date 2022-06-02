\documentclass{article}
%include polycode.fmt

\begin{document}

All the imported modules and Language Extensions which are used in Interactive are shown below. 
\begin{code}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ExplicitNamespaces #-}


module Interactive where

import Algebra
import Algebra (VectorSpace(..))

import Eigen (detNN, detGauss, newton, roots, showSol, evalMat)
import qualified Eigen as E
import Gauss (ElimOp(..), elimOpToMat, elimOpToFunc, foldElimOpsFunc, gauss, gaussTrace)
import Subspaces hiding (mEx, Quot, Sub)


import Matrix hiding (Matrix, pjMat)
import qualified Matrix as M
import QuadTree (QuadM)
import qualified QuadTree as Q
import SparseCSR (CSR)
import qualified SparseCSR as C
import ListVector (Vector, type (^), vec, dot, cross, Matrix, ToMat (..), toMatT, append', appendV, pureVec)
import qualified ListVector as L
import ListVecGauss (particularSol)
import qualified ListVecGauss as LG

import Prelude hiding ((+), (-), (*), (/), (^), (**), sum, product, recip, fromRational, span, elem)
import qualified Prelude as P
\end{code}

\section{Data types}

\begin{description}
\item \textbf{ElimOp n a}\\
    \begin{description}
        \item \textbf{Swap (Fin n) (Fin n)}\\
            Represents the swapping of rows
        \item \textbf{Mul (Fin n) a }\\
            Represents the scaling of a row
        \item \textbf{MulAdd (Fin n) (Fin n) a}\\
            Represents addition between rows where row 1 is scaled by \textbf{a}
    \end{description}
\item \textbf{Exp}\\
    \begin{description}
        \item \textbf{Const R}\\
            Represents a constant in an expression
        \item \textbf{X}\\
            Represents an unknown variable
        \item \textbf{Exp :+: Exp}\\
            Represents addition
        \item \textbf{Negate Exp}\\
            Represents negation, can be used in conjunction with addition
            to get subtraction
        \item \textbf{Exp :*: Exp}\\
            Represents multiplication
        \item \textbf{Recip Exp}\\
            Represents the multiplicative inverse, a.k.a. $\frac{1}{Exp}$
    \end{description}
\item \textbf{Matrix representations}\\
    \begin{description}
        \item \textbf{Matrix f m n}\\
            A list in list matrix representation for dense or small matrices,
            with \textbf{m} rows, \textbf{n} columns 
            and with elements of type \textbf{f}
        \item \textbf{CSR f m n}\\
            A compressed sparse row matrix representation for sparse matrices,
            with \textbf{m} rows, \textbf{n} columns 
            and with elements of type \textbf{f}
        \item \textbf{QuadM f m n}\\
            A quadtree matrix representation for sparse matrices,
            with \textbf{m} rows, \textbf{n} columns 
            and with elements of type \textbf{f}
    \end{description}
\end{description}


\section{Functions for vectors and matrices}

\subsection{Creation of matrices of vectors}
\begin{description}
\item \textbf{vec :: [f] $\to$ Vector f n}\\
    Creates a vector of size n with elements of type f
\item \textbf{toMat   :: x $\to$ Matrix (Under' x) m n}\\
    Converts input into a matrix
\item \textbf{toMatT  ::  x $\to$ Matrix (Under' x) m n}\\
    Row indexed version of 
\item \textbf{fromKeyValue :: [(Int,Int,f)] $\to$ mat f m n}\\
    Transforms a, (x,y, element), triple into a matrix
\end{description}

Some examples of these functions being used: 

\begin{code}

v1,v2 :: Vector 3 Double
v1 = vec [2,7,1] 
v2 = vec [8,2,8]

m1 :: Matrix 4 4 Double
m1 = toMat [[1,0,0,3],[2,0,8,0],[3,0,3,1],[4,0,2,0]]

m1T :: Matrix 4 4 Double
m1T = toMatT [
            [1,2,3,4],
            [0,0,0,0],
            [0,8,3,2],
            [3,0,1,0]]

m2 :: Matrix 4 4 Double
m2 = fromKeyValue [(1,1,1),(1,2,2),(1,3,3),(1,4,4),
                    (3,2,8),(3,3,3),(3,4,2),
                    (4,1,3),(4,3,1)]

m3 :: QuadM 4 4 Double
m3 = fromKeyValue [(1,1,1),(1,2,2),(1,3,3),(1,4,4),
                    (3,2,8),(3,3,3),(3,4,2),
                    (4,1,3),(4,3,1)]

m4 :: CSR 4 4 Double
m4 = fromKeyValue [(1,1,1),(1,2,2),(1,3,3),(1,4,4),
                    (3,2,8),(3,3,3),(3,4,2),
                    (4,1,3),(4,3,1)]
\end{code}





\subsection{Matrix class functions and matrix representations}

\begin{description}
\item \textbf{zero :: a}\\
    Creates the zero matrix/vector, also returns the additive identity
    for types that implement \textbf{AddGroup} which in most cases is 0
\item \textbf{one :: a}\\
    Creates the identity matrix, also returns the multiplicative identity
    for types that implement \textbf{Mul} which in most cases is 1
\item \textbf{fin :: Int $\to$ Fin n}\\
    Safe constructor of Fin, checks that the value is in range
\item \textbf{finToInt :: Fin n $\to$ Int}\\
    Returns the corresponding Int
\item \textbf{raiseBound :: Fin n $\to$ Fin m}\\
    Increases the max bound of the index
\item \textbf{values :: mat f m n $\to$ [((Fin m, Fin n), f)]}\\
    Returns all values with corresponding index 
\item \textbf{tabulate :: [((Fin m, Fin n), f)] $\to$ mat f m n}\\
    Builds a matrix from a list of positions and values
    Initializes with the zero matrix
\item \textbf{getDiagonal :: mat f n n $\to$ [(Fin n, f)]}\\
    Returns a list of elements and positions corresponding to the diagonal
\item \textbf{transpose :: mat f m n $\to$ mat f n m}\\
    Transposes a matrix
\item \textbf{changeRep :: mat1 f m n $\to$ mat2 f m n}\\
    Changes the underlying matrix type
\item \textbf{purgeToList :: mat f m n $\to$ [((Fin m, Fin n), f)]}\\
    Like values but also removes zeros
\item \textbf{purge :: mat f m n $\to$ mat f m n}\\
    Removes all zeroes from a matrix
\item \textbf{toSparse :: mat1 f m n $\to$ mat2 f m n}\\
    Similar to changeRep but also removes zeros
\item \textbf{fromMat :: Matrix (Under' x) m n $\to$ x}\\
    Converts a matrix into its underlying type
\item \textbf{extend :: mat f m n $\to$ [((Fin m, Fin n), f)] $\to$ mat f m n}\\
    Creates a matrix given a list of keyvalue pairs
\item \textbf{get :: mat f m n $\to$ (Fin m, Fin n) $\to$ f}\\
    Indexes into a matrix and gets a value
\item \textbf{getRow :: mat f m n $\to$ Fin m $\to$ [(Fin n, f)]}\\
    Returns a list of elements and positions corresponding to a row
\item \textbf{getCol :: mat f m n $\to$ Fin n $\to$ [(Fin m, f)]}\\
    Returns a list of elements and positions corresponding to a column
\item \textbf{size :: mat f m n $\to$ (Int, Int)}\\
    Returns the size of a matrix in the form of (#rows, #columns)
\item \textbf{append :: mat f m n1 $\to$ mat f m n2 $\to$ mat f m (n1 + n2)}\\
    Appends two matrices, analogous to (++)
\item \textbf{append' :: mat f m1 n $\to$ mat f m2 n $\to$ mat f (m1 + m2) n}\\
    Like appends but places the second matrix under the first
\item \textbf{set :: mat f m n $\to$ (Fin m, Fin n) $\to$ f $\to$ mat f m n }\\
    Sets the value at a given position
\item \textbf{setRow :: mat f m n $\to$ Fin m $\to$ [(Fin n, f)] $\to$ mat f m n}\\
    Sets a given row in a matrix into the given values
\item \textbf{identity :: mat f n n}\\
    A general implementation of the identity matrix
\end{description}

Some code examples:

\begin{code}

m5 :: CSR 4 4 Double
m5 = changeRep m1 

m6 :: CSR 4 4 Double
m6 = toSparse m1 

m7 :: QuadM 4 4 Double
m7 = toSparse m1 

\end{code}


\subsection{Matrix and vector operations}

\begin{description}
\item \textbf{neg :: a $\to$ a }\\
    Negates the matrix/vector
\item \textbf{(+)  :: a $\to$ a $\to$ a}\\
    Addition between matrices/vectors
\item \textbf{(-)  :: a $\to$ a $\to$ a}\\
    Subtraction between matrices/vectors 
\item \textbf{(*) :: a $\to$ a $\to$ a}\\
    Multiplication between square matrices
\item \textbf{(£) :: b $\to$ a $\to$ a}\\
    Scalar multiplication where b is a scalar
    and a is a matrix/vector
\item \textbf{(**) :: a $\to$ b $\to$ c}\\
    Matrix Vector multiplication and matrix matrix multiplication
    of all correct sizes
\item \textbf{dot :: Vector f n $\to$ Vector f n $\to$ f}\\
    dotproduct of two vectors
\item \textbf{cross :: Vector f 3 $\to$ Vector f 3 $\to$ Vector f 3}\\
    crossproduct of two vectors
\item \textbf{appendV :: Matrix f m n $\to$ Vector f m $\to$ Matrix f m (n+1)}\\
    Appends a vector to the right of a matrix
\item \textbf{gauss :: mat f m n $\to$ mat f m n}\\
    Transforms a given matrix into row echelon form
\item \textbf{gaussTrace :: mat f m n $\to$ [ElimOp m f]}\\
    Returns the required row operations to 
    transform a given matrix into row echelon form
\item \textbf{showSol :: Matrix f m n $\to$ IO()}\\
    Prints solution to a given solved list matrix
    for example, use showSol \$ gauss m1
\item \textbf{particularSol :: Matrix f m n $\to$ Vector f (n -1)}\\
    apply when solving systems of equations
    each element in list represents variable values 
\item \textbf{solveQ :: Matrix f m n $\to$ Vector f m $\to$ QuotientSpace (Vector f n)}\\
    Returns the set of solutions to Ax=v
\item \textbf{solveQ' :: Matrix f m (n+1)  $\to$ QuotientSpace (Vector f n)}\\
    Equivalent to solveQ but takes a matrix A `appendV` v representing Ax=v
\item \textbf{elimOpToMat :: ElimOp n f $\to$ mat f n n}\\
    Representation of an elementary row operation as a matrix 
\item \textbf{elimOpToFunc :: ElimOp m f $\to$ (mat f m n $\to$ mat f m n)}\\
    Representation of an elementary row operation as a function 
\item \textbf{foldElimOpsFunc :: [ElimOp m f] $\to$ (mat f m n $\to$ mat f m n)}\\
    Reduces a list/trace of elimOps to a single function
\item \textbf{detNN :: Matrix f n n $\to$ f}\\
    Determinant of a matrix using newtons method, bad complexity but works 
    for smaller matrices
\item \textbf{detGauss :: Matrix f n n $\to$ f }\\
    Determinant of a matrix using Gaussian elimination, fast but produces ugly
    expression in the end
\item \textbf{roots :: Exp $\to$ [a] $\to$ [a]}\\
    Finds the roots of an exp in a given span
\item \textbf{evalMat :: Matrix Exp m n $\to$ f $\to$ Matrix f m n}\\
    Evaluates the expressions in a given list matrix
\item \textbf{newton :: Exp $\to$ a $\to$ a $\to$ a}\\
    finds the zeros of the characteristic equation 
    which equals the eigen values using newtons method
\item \textbf{span :: [v] $\to$ Subspace v}\\
    Creates a subspace of a span of vectors
\item \textbf{elem :: Vector f n $\to$ Subspace (Vector f n) $\to$ Bool}\\
    Checks if a vector exists within a subspace 
\item \textbf{span' :: Matrix f m n $\to$ Vector f m $\to$ Bool}\\
    Checks if a vector is in the span of a list of vectors
    Normally span is defined as a set, but here we use it as a condition such that
    span [w1..wn] v = v `elem` span(w1..w2)
\item \textbf{spans :: Matrix f m n1 $\to$ Matrix f m n2 $\to$ Bool}\\
    Checks that m1 spans at least as much as m2 
\item \textbf{linIndep :: Matrix f m n $\to$ Bool}\\
    Checks if the vectors in a matrix are linearly independent 
\item \textbf{makeLinIndep :: [Vector f n] $\to$ [Vector f n]}\\
    Any list of vector can be made linearly independent by 
    iteratively removing vectors that are in the span of the remaining vectors
\item \textbf{isBasis :: Matrix f m n $\to$ Bool}\\
    Checks if the vectors in a matrix forms a basis of their vectorSpace, where
    a basis of V is a list of vectors in V that is linearly independent and spans V
\item \textbf{nullSpace :: Matrix f m n $\to$ Subspace (Vector f n)}\\
    The null space of a linear transformation is the set of vectors that maps to 0
    In terms of linear equations it is the solution to Ax=0
\item \textbf{range :: Matrix f m n $\to$ Subspace (Vector f m)}\\
    The range of a linear map is the set of possible outputs
\end{description}

\end{document}
