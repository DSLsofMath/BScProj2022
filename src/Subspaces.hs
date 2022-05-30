{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE StandaloneDeriving #-}

module Subspaces where


import GHC.TypeLits hiding (type (^))
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), recip, sum, product, (**), span, elem)

import Algebra
import ListVector
import ListVecGauss

-- Mathematically, a subspace is a set and we would therefore like to represent a 
-- given subspace as its own type -- in the same way we do with normal vector spaces.
-- However, we need to create subspaces from terms such as the range or null space 
-- of a concrete linear map or from a given system of equation. This is a great 
-- example of dependent types, which unfortunately is not available in Haskell.
--
-- For our purposes we should therefore consider a subspace to be a term, 
-- more precisely we could define it as a list of vectors that spans the subspace. 
-- On the type level we would get something like "Subspace v" which corresponds to
-- the set of all subspaces of v.


-- | A subspace is defined by a list of vectors that spans the space.
--   We might want to add an invariant to the list such that it is linearly independent.
newtype Subspace v = Sub [v] deriving (Show, Eq)

-- | We define Addgroup over a subspace as "sum of subsets"
instance (AddGroup v) => AddGroup (Subspace v) where
    Sub v1 + Sub v2 = Sub $ v1 ++ v2 -- We might want to remove dependent elements
    Sub v1 - Sub v2 = undefined -- Unsure of how we should handle subtraction
    zero = Sub []
    

-------------------------------------------
-- Properties on Finite-Dimensional Vector Spaces
-- Span, linear independence and basis

-- | Creates a subspace of a span of vectors
span :: [v] -> Subspace v
span = Sub

-- | Checks if a vector exists within a subspace 
elem :: (Eq f, Field f) => Vector n f -> Subspace (Vector n f) -> Bool
elem v (Sub vs) = span' (M $ V vs) v

-- | Checks if a vector is in the span of a list of vectors
--   Normally span is defined as a set, but here we use it as a condition such that
--   span [w1..wn] v = v `elem` span(w1..w2)
span' :: (Eq f, Field f) => Matrix m n f -> Vector m f -> Bool
span' m v = all (\x -> pivotLen x /= 1) . unpack . transpose . gauss $ append m v'
    where v' = M (V @1 ([v])) -- Forces the matrix to size n 1
          pivotLen xs = length (dropWhile (==zero) xs)
-- M (Vector @_ @1 (L [v])) -- Forces the matrix to size n 1
-- | Checks that m1 spans at least as much as m2 
spans :: (KnownNat m, KnownNat n2, Eq f, Field f) => Matrix m n1 f -> Matrix m n2 f -> Bool
m1 `spans` m2 = all (span' m1) (matToList m2)

-- | Checks that m spans the whole vectorSpace
spansSpace :: (KnownNat m, Eq f, Field f) => Matrix m n f -> Bool
spansSpace m = m `spans` idm

-- | Checks if the vectors in a matrix are linearly independent
linIndep :: (Eq f, Field f) => Matrix m n f -> Bool
linIndep = not . any (\(v, m) -> m `span'` v ) . separateCols 

-- | Any list of vector can be made linearly independent by iteratively removing vectors 
-- that are in the span of the remaining vectors
makeLinIndep :: (Eq f, Field f) => [Vector n f] -> [Vector n f]
makeLinIndep = foldl (\vs v -> if span' (M $ V vs) v then vs else v:vs) [] 


-- 2.27
-- Definition: basis
-- A basis of V is a list of vectors in V that is linearly independent and spans V

-- | Checks if the vectors in a matrix forms a basis of their vectorSpace
isBasis :: (KnownNat m, KnownNat (n-1), Eq f, Field f) => Matrix m n f -> Bool
isBasis m = spansSpace m && linIndep m

-- | The null space of a linear transformation is the set of vectors that maps to 0
--   In terms of linear equations it is the solution to Ax=0
nullSpace :: (KnownNat n, Field f, Eq f) =>  Matrix m n f -> Subspace (Vector n f)
nullSpace m = Sub $ [ V b | (a,b) <- splitAt height <$> reduceCol m, all (== zero) a]
    where height = length (head $ unpack m)
          reduceCol m = unpack $ transpose $ gauss (transpose m `append` idm)

-- | The range of a linear map is the set of possible outputs
range :: (Eq f, Field f) => Matrix m n f -> Subspace (Vector m f)
range = Sub . makeLinIndep . matToList


-- | The dimension of a Subspace
dim :: Subspace v -> Int
dim (Sub vs) = length vs


          
-- | A quotient space is a subspace translated by a vector
data QuotientSpace v = Quot v (Subspace v) deriving (Show)

-- | Checks if a vector exists within a quotient space, 
--   i.e. is a solution to a Ax=v  
elemQ :: (KnownNat n, Eq f, Field f) => Vector n f -> QuotientSpace (Vector n f) -> Bool
elemQ v (Quot p sub) = (v - p) `elem` sub 

-- | Returns the set of solutions to Ax=v
solveQ :: (KnownNat n, Field f, Eq f, (n ~ (n+1-1)) ) => Matrix m n f -> Vector m f -> QuotientSpace (Vector n f)
solveQ m v = Quot (particularSol $ m `appendV` v) (nullSpace m)


-- | Equivalent to solveQ but takes a matrix A `appendV` v representing Ax=v
solveQ' :: (KnownNat n, Field f, Eq f, (n ~ (n+1-1)) ) =>  Matrix m (n+1) f  -> QuotientSpace (Vector n f)
solveQ' m = let (v', m') = last $ separateCols m in solveQ m' v'


mEx :: Matrix 3 3 R
mEx = toMat [[1,2,1],[0,1,1],[1,2,1]]


