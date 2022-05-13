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
import Subspaces hiding (mEx, QuotientSpace(..), Subspace(..))


import Matrix hiding (Matrix, pjMat)
import qualified Matrix as M
import QuadTree (QuadM)
import qualified QuadTree as Q
import SparseCSR (CSR)
import qualified SparseCSR as C
import ListVector (Vector, type (^), vec, dot, cross, Matrix, ToMat (..))
import qualified ListVector as L
import ListVecGauss (particularSol)

import Prelude hiding ((+), (-), (*), (/), (^), (**), sum, product, recip, fromRational, span, elem)
import qualified Prelude as P
\end{code}

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

v1,v2 :: Vector Double 3
v1 = vec [2,7,1] 
v2 = vec [8,2,8]

m1 :: Matrix Double 4 4
m1 = toMat [[1,0,0,3],[2,0,8,0],[3,0,3,1],[4,0,2,0]]

m1T :: Matrix Double 4 4
m1T = toMatT [
            [1,2,3,4],
            [0,0,0,0],
            [0,8,3,2],
            [3,0,1,0]]

m2 :: Matrix Double 4 4
m2 = fromKeyValue [(1,1,1),(1,2,2),(1,3,3),(1,4,4),
                    (3,2,8),(3,3,3),(3,4,2),
                    (4,1,3),(4,3,1)]

\end{code}



\begin{code}

m3 :: CSR Double 4 4
m3 = changeRep m1 

m4 :: CSR Double 4 4
m4 = toSparse m1 


m5 :: CSR Double 4 4
m5 = fromKeyValue [(1,1,1),(1,2,2),(1,3,3),(1,4,4),
                    (3,2,8),(3,3,3),(3,4,2),
                    (4,1,3),(4,3,1)]


m6 :: QuadM Double 4 4
m6 = toSparse m1 

m7 :: QuadM Double 4 4
m7 = fromKeyValue [(1,1,1),(1,2,2),(1,3,3),(1,4,4),
                    (3,2,8),(3,3,3),(3,4,2),
                    (4,1,3),(4,3,1)]


\end{code}
\item \textbf{}\\


\subsection{Matrix class functions and matrix representations}

\begin{description}
\item \textbf{}\\

\item \textbf{}\\

\item \textbf{}\\

\item \textbf{}\\

\end{description}


\subsection{Matrix and vector operations}





Constants 

    Creates the zero matrix/vector
    zero :: a 

    Creates the identity matrix
    one :: a

    A general implementation of the identity matrix
    identity :: mat f n n

Unary operations

    Negates the matrix/vector
    neg :: a $\to$ a 


    Determinant of a list matrix
        Recommended for small matrixes
        detNN :: Matrix f n n $\to$ f

        Fast but has floating point errors
        detGauss :: Matrix f n n $\to$ f 
    
    Prints solution to a given solved list matrix
    for example, use showSol \$ gauss m1
    showSol :: Matrix f m n $\to$ IO()

    Representation of an elementary row operation as a matrix 
    elimOpToMat :: ElimOp n f $\to$ mat f n n

    Representation of an elementary row operation as a function 
    elimOpToFunc :: ElimOp m f $\to$ (mat f m n $\to$ mat f m n)

    Reduces a list/trace of elimOps to a single function
    foldElimOpsFunc :: [ElimOp m f] $\to$ (mat f m n $\to$ mat f m n)

    Transforms a given matrix into row echelon form
    gauss :: mat f m n $\to$ mat f m n

    Returns the required row operations to 
    transform a given matrix into row echelon form
    gaussTrace :: mat f m n $\to$ [ElimOp m f]

    Creates a subspace of a span of vectors
    span :: [v] $\to$ Subspace v

    Checks if the vectors in a matrix forms a basis of their vectorSpace, where
    a basis of V is a list of vectors in V that is linearly independent and spans V
    isBasis :: Matrix f m n $\to$ Bool

    The null space of a linear transformation is the set of vectors that maps to 0
    In terms of linear equations it is the solution to Ax=0
    nullSpace :: Matrix f m n $\to$ Subspace (Vector f n)

    The range of a linear map is the set of possible outputs
    range :: Matrix f m n $\to$ Subspace (Vector f m)

    Equivalent to solveQ but takes a matrix A `append` v representing Ax=v
    solveQ' :: Matrix f m (n+1)  $\to$ QuotientSpace (Vector f n)

    Safe constructor of Fin, checks that the value is in range
    fin :: KnownNat n => Int $\to$ Fin n

    Returns the corresponding Int
    finToInt :: Fin n $\to$ Int

    Increases the max bound of the index
    raiseBound :: Fin n $\to$ Fin m

    Returns all values with corresponding index 
    values :: mat f m n $\to$ [((Fin m, Fin n), f)]

    Builds a matrix from a list of positions and values
    Initializes with the zero matrix
    tabulate :: AddGroup (mat f m n) => [((Fin m, Fin n), f)] $\to$ mat f m n

    Returns a list of elements and positions corresponding to the diagonal
    getDiagonal :: mat f n n $\to$ [(Fin n, f)]

    Transposes a matrix
    transpose :: mat f m n $\to$ mat f n m

    Changes the underlying matrix type
    changeRep :: mat1 f m n $\to$ mat2 f m n

    Like values but also removes zeros
    purgeToList :: (Matrix mat, Eq f, AddGroup f) => mat f m n $\to$ [((Fin m, Fin n), f)]

    Removes all zeroes from a matrix
    purge :: (Matrix mat, Eq f, AddGroup f, AddGroup (mat f m n)) => mat f m n $\to$ mat f m n

    Similar to changeRep but also removes zeros
    toSparse :: mat1 f m n $\to$ mat2 f m n


    Converts a matrix into its underlying type
    fromMat :: Matrix (Under' x) m n $\to$ x

    apply when solving systems of equations
    each element in list represents variable values 
    particularSol :: Matrix f m n $\to$ Vector f (n -1)

Binary operations

    Addition between matrices/vectors
    (+)  :: a $\to$ a $\to$ a
    Subtraction between matrices/vectors 
    (-)  :: a $\to$ a $\to$ a
    
    Multiplication between square matrices
    (*) :: a $\to$ a $\to$ a

    Scalar multiplication where b is a scalar
    and a is a matrix/vector
    (£) :: b $\to$ a $\to$ a

    Matrix Vector multiplication and matrix matrix multiplication
    of all correct sizes
    (**) :: a $\to$ b $\to$ c

    Finds the roots of an exp in a given span
    roots :: Exp $\to$ [a] $\to$ [a]

    Evaluates the expressions in a given list matrix
    evalMat :: Matrix Exp m n $\to$ f $\to$ Matrix f m n

    Checks if a vector exists within a subspace 
    elem :: Vector f n $\to$ Subspace (Vector f n) $\to$ Bool

    Checks if a vector is in the span of a list of vectors
    Normally span is defined as a set, but here we use it as a condition such that
    span [w1..wn] v = v `elem` span(w1..w2)
    span' :: Matrix f m n $\to$ Vector f m $\to$ Bool

    Checks that m1 spans at least as much as m2 
    spans :: Matrix f m n1 $\to$ Matrix f m n2 $\to$ Bool
    
    Checks if the vectors in a matrix are linearly independent 
    linIndep :: (Eq f, Field f) => Matrix f m n $\to$ Bool

    Any list of vector can be made linearly independent by 
    iteratively removing vectors that are in the span of the remaining vectors
    makeLinIndep :: (Eq f, Field f) => [Vector f n] $\to$ [Vector f n]

    Returns the set of solutions to Ax=v
    solveQ :: Matrix f m n $\to$ Vector f m $\to$ QuotientSpace (Vector f n)

    Creates a matrix given a list of keyvalue pairs
    extend :: mat f m n $\to$ [((Fin m, Fin n), f)] $\to$ mat f m n


    Indexes into a matrix and gets a value
    get :: mat f m n $\to$ (Fin m, Fin n) $\to$ f
    
    Returns a list of elements and positions corresponding to a row
    getRow :: mat f m n $\to$ Fin m $\to$ [(Fin n, f)]

    Returns a list of elements and positions corresponding to a column
    getCol :: mat f m n $\to$ Fin n $\to$ [(Fin m, f)]

    Appends two matrices, analogous to (++)
    append :: mat f m n1 $\to$ mat f m n2 $\to$ mat f m (n1 + n2)

    dot :: Vector f n $\to$ Vector f n $\to$ f

    crossproduct of two vectors
    cross :: Vector f 3 $\to$ Vector f 3 $\to$ Vector f 3

Trinary operations
    Sets the value at a given position
    set :: mat f m n $\to$ (Fin m, Fin n) $\to$ f $\to$ mat f m n 

    finds the zeros of the characteristic equation 
    which equals the eigen values using newtons method
    newton :: Exp $\to$ a $\to$ a $\to$ a

    Sets a given row in a matrix into the given values
    setRow :: mat f m n $\to$ Fin m $\to$ [(Fin n, f)] $\to$ mat f m n

---------- Datatypes ----------

Represents elementary row operations
data ElimOp n a = Swap (Fin n) (Fin n) 
                | Mul (Fin n) a 
                | MulAdd (Fin n) (Fin n) a

\end{document}