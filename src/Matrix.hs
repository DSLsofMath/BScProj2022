{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE NoStarIsType #-}
{-# LANGUAGE UndecidableInstances #-}

module Matrix where

import GHC.TypeLits
import Data.Kind (Type)
import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), (^), sum, product, recip, fromRational)

import Algebra
import FiniteIndex
import Finite
import Expression.Exp

-- | Generic class for a matrix with focus on sparse representations
--   The class is roughly based on the paper "APLicative Programming with Naperian Functors"
--   by Jeremy Gibbons. 
--   However we make some changes to better suite sparse matrices,
--   most notably we only allow implementations for two dimensions.
--   
class Matrix (mat :: Nat -> Nat -> Type -> Type) where
    {-# MINIMAL tabulate, values, mulMat, addMat #-}

    -- | Returns all values with corresponding index 
    values :: mat m n f -> [((Fin m, Fin n), f)]

    -- | Builds a matrix from a list of positions and values
    tabulate :: (KnownNats m n, AddGroup f) => [((Fin m, Fin n), f)] -> mat m n f

    

    -- | Matrix multiplication is needed in the class to 
    --   define Composable for matrices in a nice way.
    mulMat :: (KnownNat a, Ring f) => mat a b f -> mat b c f -> mat a c f

    -- | Addition of multiplication
    --   Needed to define Add- instances for all matrices
    addMat :: AddGroup f => mat m n f -> mat m n f -> mat m n f

    -- | Transforms a matrix into a linear map
    --   The vector are represented as (Fin n -> f) to be compatible with module Finite.
    --   TODO could be optimized a lot, and should be per instance
    linMap' :: (Ring f, Matrix mat) => mat m n f -> (Fin n -> f) -> (Fin m -> f)
    linMap' mat f = sum [ (\i -> if i == i' then f j * s else zero) | ((i', j), s) <- values mat]


-- | Transforms keyvalue pairs into a matrix
fromKeyValue :: (KnownNats m n, Matrix mat, AddGroup f) => [(Int,Int,f)] -> mat m n f
fromKeyValue as = tabulate [ ((fin a,fin b),f) | (a,b,f) <- as ]

-- | Indexes into a matrix and gets a value
get :: (Matrix mat, AddGroup f) => mat m n f -> (Fin m, Fin n) -> f
get m i = case lookup i (values m) of Just a  -> a
                                      Nothing -> zero
    
-- | Returns a list of elements and positions corresponding to a row
getRow :: Matrix mat => mat m n f -> Fin m -> [(Fin n, f)]
getRow m i = [ (j, a) | ((i', j), a) <- values m, i' == i]

-- | Returns a list of elements and positions corresponding to a column
getCol :: Matrix mat => mat m n f -> Fin n -> [(Fin m, f)]
getCol m j = [ (i, a) | ((i, j'), a) <- values m, j' == j]

-- | Returns a list of elements and positions corresponding to the diagonal
getDiagonal :: Matrix mat => mat n n f -> [(Fin n, f)]
getDiagonal m = [ (i, a) | ((i, j), a) <- values m, i == j ]

-- | Sets a given row in a matrix into the given values
setRow :: (KnownNats m n, Matrix mat, AddGroup f) => mat m n f -> Fin m -> [(Fin n, f)] -> mat m n f
setRow m i r = tabulate $ new ++ old
    where old = filter ((i /=) . fst . fst) $ values m
          new = [ ((i, j), a) | (j, a) <- r ]

-- | Returns the size of a matrix in the form of (#rows, #columns)
size :: forall m n mat f. KnownNats m n => mat m n f -> (Int, Int) 
size _ = (natInt' @m, natInt' @n)


-- | Appends two matrices, analogous to (++)
append :: (KnownNats m (n1 + n2), KnownNat n1, Matrix mat, AddGroup f) => mat m n1 f -> mat m n2 f -> mat m (n1 + n2) f
append m1 m2 = tabulate $ m1' ++ m2'
    where m1' = [ ((i, Fin j),       a) | ((i, Fin j), a) <- values m1 ]
          m2' = [ ((i, Fin (n + j)), a) | ((i, Fin j), a) <- values m2 ]
          (_,n) = size m1

-- | Like appends but places the second matrix under the first 
append' :: (KnownNats (m1 + m2) n, KnownNat m1, Matrix mat, AddGroup f) => mat m1 n f -> mat m2 n f -> mat (m1 + m2) n f
append' m1 m2 = tabulate $ m1' ++ m2'
    where m1' = [ ((Fin  i,      j), a) | ((Fin i, j), a) <- values m1 ]
          m2' = [ ((Fin (n + i), j), a) | ((Fin i, j), a) <- values m2 ]
          (n,_) = size m1

transpose :: (KnownNats m n, Matrix mat, AddGroup f) => mat m n f -> mat n m f
transpose = tabulate . map (\((i,j),a) -> ((j,i),a)) . values

-- | Changes the underlying matrix type
--   useful when converting from one representation to another
--   changeRep (identity :: CSR R 5 5) :: Matrix R 5 5
changeRep :: (Matrix mat2, KnownNats m n, Matrix mat1, AddGroup f) => mat1 m n f -> mat2 m n f
changeRep = tabulate . values

changeUnder :: (AddGroup f2, KnownNats m n, Matrix mat) => (f1 -> f2) -> mat m n f1 -> mat m n f2
changeUnder f = tabulate . map (\(i,a) -> (i, f a)) . values

diagMat :: (KnownNat n, Matrix mat, AddGroup f) => f -> mat n n f
diagMat f = tabulate [ ((i,i), f) | i <- [minBound .. maxBound]]

-- | A general implementation of the identity matrix
identity :: (KnownNat n, Matrix mat, Ring f) => mat n n f
identity = diagMat one

-- | A general implementation of the zero matrix
zeroMat :: (KnownNats m n, Matrix mat, AddGroup f) => mat m n f
zeroMat = tabulate []

-- | Like values but also removes zeros
purgeToList :: (Matrix mat, Eq f, AddGroup f) => mat m n f -> [((Fin m, Fin n), f)]
purgeToList = (filter ((zero /=) . snd)) . values

-- | Removes all zeroes from a matrix
purge :: (KnownNats m n, Matrix mat, Eq f, AddGroup f) => mat m n f -> mat m n f
purge = toSparse

-- | Like changeRep but also removes zeros
toSparse :: (Matrix mat2, KnownNats m n, Matrix mat1, Eq f, AddGroup f) => mat1 m n f -> mat2 m n f
toSparse = tabulate . purgeToList


toConst :: (KnownNats m n, Matrix mat, AddGroup f) => mat m n f -> mat m n (Exp f)
toConst = changeUnder Const


----------------------------------------------------------------------------------------
-- Composability with Finite 
--
-- Any matrix is a linear map between arbitrary finite vectorspaces of
-- compatible dimensions

linMap :: (Finite w, Finite v, Matrix mat, f ~ Under w, f ~ Under v) => mat (Dim w) (Dim v) f -> (v -> w)
linMap mat = fromIndexing . linMap' mat . scalar 


----------------------------------------------------------------------------------------
-- Instances for Matrix

instance (KnownNats m n, Matrix mat, AddGroup f) => AddGroup (mat m n f) where
    (+) = addMat
    neg = changeUnder neg
    zero = zeroMat

instance (KnownNats m n, Matrix mat, Ring f) => VectorSpace (mat m n f) where
    type Under (mat m n f) = f
    s £ m = changeUnder (s *) m

-- TODO optimize, 
instance (KnownNats m n, KnownNat (m * n), Matrix mat, Ring f) => Finite (mat m n f) where
    type Dim (mat m n f) = m * n
    basis i = tabulate [((splitFin i), one)]
    scalar mat i = get mat (splitFin i)
    

instance (Matrix mat, AddGroup f, Eq f ) => Eq (mat m n f) where
    m1 == m2 = purgeToList m1 == purgeToList m2

instance (KnownNat a, Matrix mat, Ring f, mat ~ mat', f ~ f', b ~ b') => 
            Composable (mat a b f) (mat' b' c f') (mat a c f) where
    (**) = mulMat

instance (KnownNat n, Matrix mat, Ring f) => Mul (mat n n f) where
    (*) = mulMat
    one = identity

--------------------------------------------------------------------------
-- examples

pjMat :: (KnownNat n, Matrix mat, Ring f) => mat n n f
pjMat = tabulate $ concat [ f i | i <- [minBound .. maxBound] ]
    where f i | i == minBound = tail $ stencil i
              | i == maxBound = init $ stencil i
              | otherwise     =        stencil i
          stencil i = [ ((pred i, i), one), ((i,i), neg (one+one)), ((succ i, i), one) ]

-- | Represents derivation on polynomials
deriveMat :: (KnownNat n, Matrix mat, Ring f) => mat n n f
deriveMat = tabulate $ inc one [ (i, succ i) | i <- init [minBound .. maxBound] ]
    where inc _ []     = []
          inc v (i:is) = (i, v) : inc (v + one) is

