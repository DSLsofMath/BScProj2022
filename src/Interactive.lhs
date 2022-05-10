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
import qualified ListVecGauss as LG

import Prelude hiding ((+), (-), (*), (/), (^), (**), sum, product, recip, fromRational, span, elem)
import qualified Prelude as P
\end{code}

\begin{code}

v1,v2 :: Vector Double 3
v1 = vec [2,7,1] 
v2 = vec [8,2,8]
\end{code}

example matrix from a list of lists

\begin{code}
m1 :: Matrix Double 4 4
m1 = toMat [[1,0,0,3],[2,0,8,0],[3,0,3,1],[4,0,2,0]]
\end{code}

using key value pairs to create the same matrix,
useful for creating larger sparse matrices

\begin{code}
m2 :: Matrix Double 4 4
m2 = fromKeyValue [(1,1,1),(1,2,2),(1,3,3),(1,4,4),
                    (3,2,8),(3,3,3),(3,4,2),
                    (4,1,3),(4,3,1)]
\end{code}

changeRep can be used to convert any matrix representation to another.
However this keeps zeroes

\begin{code}
m3 :: CSR Double 4 4
m3 = changeRep m1 
\end{code}


toSparse does the same thing as changeRep 
but removes zeros before reconstructing

\begin{code}
m4 :: CSR Double 4 4
m4 = toSparse m1 
\end{code}

fromKeyValue can also be used to 
make any matrix representation directly

\begin{code}
m5 :: CSR Double 4 4
m5 = fromKeyValue [(1,1,1),(1,2,2),(1,3,3),(1,4,4),
                    (3,2,8),(3,3,3),(3,4,2),
                    (4,1,3),(4,3,1)]
\end{code}

Likewise the same can be done
to create a QuadTree sparse matrix

\begin{code}
m6 :: QuadM Double 4 4
m6 = toSparse m1 

m7 :: QuadM Double 4 4
m7 = fromKeyValue [(1,1,1),(1,2,2),(1,3,3),(1,4,4),
                    (3,2,8),(3,3,3),(3,4,2),
                    (4,1,3),(4,3,1)]
\end{code}

---------- A list of functions for vectors and matrices ----------

Constants 

    Creates the zero matrix/vector
    zero :: a 

    Creates the identity matrix
    one :: a

    A general implementation of the identity matrix
    identity :: mat f n n

Unary operations

    Negates the matrix/vector
    neg :: a -> a 

    Creates a vector of size n with elements of type f
    vec :: [f] -> Vector f n

    Determinant of a list matrix
        Recommended for small matrixes
        detNN :: Matrix f n n -> f

        Fast but has floating point errors
        detGauss :: Matrix f n n -> f 
    
    Prints solution to a given solved list matrix
    for example, use showSol \$ gauss m1
    showSol :: Matrix f m n -> IO()

    Representation of an elementary row operation as a matrix 
    elimOpToMat :: ElimOp n f -> mat f n n

    Representation of an elementary row operation as a function 
    elimOpToFunc :: ElimOp m f -> (mat f m n -> mat f m n)

    Reduces a list/trace of elimOps to a single function
    foldElimOpsFunc :: [ElimOp m f] -> (mat f m n -> mat f m n)

    Transforms a given matrix into row echelon form
    gauss :: mat f m n -> mat f m n

    Returns the required row operations to 
    transform a given matrix into row echelon form
    gaussTrace :: mat f m n -> [ElimOp m f]

    Creates a subspace of a span of vectors
    span :: [v] -> Subspace v

    Checks if the vectors in a matrix forms a basis of their vectorSpace, where
    a basis of V is a list of vectors in V that is linearly independent and spans V
    isBasis :: Matrix f m n -> Bool

    The null space of a linear transformation is the set of vectors that maps to 0
    In terms of linear equations it is the solution to Ax=0
    nullSpace :: Matrix f m n -> Subspace (Vector f n)

    The range of a linear map is the set of possible outputs
    range :: Matrix f m n -> Subspace (Vector f m)

    Equivalent to solveQ but takes a matrix A `append` v representing Ax=v
    solveQ' :: Matrix f m (n+1)  -> QuotientSpace (Vector f n)

    Safe constructor of Fin, checks that the value is in range
    fin :: KnownNat n => Int -> Fin n

    Returns the corresponding Int
    finToInt :: Fin n -> Int

    Increases the max bound of the index
    raiseBound :: Fin n -> Fin m

    Returns all values with corresponding index 
    values :: mat f m n -> [((Fin m, Fin n), f)]

    Builds a matrix from a list of positions and values
    Initializes with the zero matrix
    tabulate :: AddGroup (mat f m n) => [((Fin m, Fin n), f)] -> mat f m n

    Returns a list of elements and positions corresponding to the diagonal
    getDiagonal :: mat f n n -> [(Fin n, f)]

    Transposes a matrix
    transpose :: mat f m n -> mat f n m

    Changes the underlying matrix type
    changeRep :: mat1 f m n -> mat2 f m n

    Like values but also removes zeros
    purgeToList :: (Matrix mat, Eq f, AddGroup f) => mat f m n -> [((Fin m, Fin n), f)]

    Removes all zeroes from a matrix
    purge :: (Matrix mat, Eq f, AddGroup f, AddGroup (mat f m n)) => mat f m n -> mat f m n

    Similar to changeRep but also removes zeros
    toSparse :: mat1 f m n -> mat2 f m n

    Converts input into a matrix
    toMat   :: x -> Matrix (Under' x) m n

    Converts a matrix into its underlying type
    fromMat :: Matrix (Under' x) m n -> x

    apply when solving systems of equations
    each element in list represents variable values 
    particularSol :: Matrix f m n -> Vector f (n -1)

Binary operations

    Addition between matrices/vectors
    (+)  :: a -> a -> a
    Subtraction between matrices/vectors 
    (-)  :: a -> a -> a
    
    Multiplication between square matrices
    (*) :: a -> a -> a

    Scalar multiplication where b is a scalar
    and a is a matrix/vector
    (£) :: b -> a -> a

    Matrix Vector multiplication and matrix matrix multiplication
    of all correct sizes
    (**) :: a -> b -> c

    Finds the roots of an exp in a given span
    roots :: Exp -> [a] -> [a]

    Evaluates the expressions in a given list matrix
    evalMat :: Matrix Exp m n -> f -> Matrix f m n

    Checks if a vector exists within a subspace 
    elem :: Vector f n -> Subspace (Vector f n) -> Bool

    Checks if a vector is in the span of a list of vectors
    Normally span is defined as a set, but here we use it as a condition such that
    span [w1..wn] v = v `elem` span(w1..w2)
    span' :: Matrix f m n -> Vector f m -> Bool

    Checks that m1 spans at least as much as m2 
    spans :: Matrix f m n1 -> Matrix f m n2 -> Bool
    
    Checks if the vectors in a matrix are linearly independent 
    linIndep :: (Eq f, Field f) => Matrix f m n -> Bool

    Any list of vector can be made linearly independent by 
    iteratively removing vectors that are in the span of the remaining vectors
    makeLinIndep :: (Eq f, Field f) => [Vector f n] -> [Vector f n]

    Returns the set of solutions to Ax=v
    solveQ :: Matrix f m n -> Vector f m -> QuotientSpace (Vector f n)

    Creates a matrix given a list of keyvalue pairs
    extend :: mat f m n -> [((Fin m, Fin n), f)] -> mat f m n

    Transforms a triple into a matrix, 
    recommended for easier usage
    fromKeyValue :: [(Int,Int,f)] -> mat f m n

    Indexes into a matrix and gets a value
    get :: mat f m n -> (Fin m, Fin n) -> f
    
    Returns a list of elements and positions corresponding to a row
    getRow :: mat f m n -> Fin m -> [(Fin n, f)]

    Returns a list of elements and positions corresponding to a column
    getCol :: mat f m n -> Fin n -> [(Fin m, f)]

    Appends two matrices, analogous to (++)
    append :: mat f m n1 -> mat f m n2 -> mat f m (n1 + n2)

    dot :: Vector f n -> Vector f n -> f

    crossproduct of two vectors
    cross :: Vector f 3 -> Vector f 3 -> Vector f 3

Trinary operations
    Sets the value at a given position
    set :: mat f m n -> (Fin m, Fin n) -> f -> mat f m n 

    finds the zeros of the characteristic equation 
    which equals the eigen values using newtons method
    newton :: Exp -> a -> a -> a

    Sets a given row in a matrix into the given values
    setRow :: mat f m n -> Fin m -> [(Fin n, f)] -> mat f m n

---------- Datatypes ----------

Represents elementary row operations
data ElimOp n a = Swap (Fin n) (Fin n) 
                | Mul (Fin n) a 
                | MulAdd (Fin n) (Fin n) a
