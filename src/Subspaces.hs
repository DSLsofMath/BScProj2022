{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE FlexibleContexts #-}

module Subspaces where


import GHC.TypeLits hiding (type (^))
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), recip, sum, product, (**), span)

import Algebra
import ListVector


-- | A subspace is defined by its overlying vector space, its dimension
--   and a unique identifier. The identifier is needed only to avoid composition
--   of two subspaces with different bases.
--
--   Internally Subspace contains its basis and a vector
data Subspace x v (n :: Nat) = Sub { getBasis  :: List v n,
                                     getVector :: Vector (Under v) n
                                   }

instance (VectorSpace v, Show v) => Show (Subspace x v n) where
    show = show . getVec

getVec :: VectorSpace v => Subspace x v n -> v
getVec (Sub b v) = v `eval` b

instance (KnownNat n, VectorSpace v) => AddGroup (Subspace x v n) where
    -- We need a way to grantee that b1 == b2
    -- Might be solved with the introduction of type x
    Sub b1 v1 + Sub b2 v2 = Sub b1 (v1 + v2)
    Sub b1 v1 - Sub b2 v2 = Sub b1 (v1 - v2)
    zero = Sub undefined zero
    
instance (KnownNat n, VectorSpace v) => VectorSpace (Subspace x v n) where
    type Under (Subspace x v n) = Under v
    s £ (Sub b v) = Sub b (s£v)

instance (KnownNat n, VectorSpace v) => Finite (Subspace x v n) where
    type Dim (Subspace x v n) = n
    type BasisVec (Subspace x v n) = v
    basis' (Sub b _) = b

-------------------------------------------
-- Properties on Finite-Dimensional Vector Spaces
-- Span, linear independence and basis

-- | Checks if a vector is in the span of a list of vectors
--   Normaly span is defined as a set, but here we use it as a condition such that
--   span [w1..wn] v = v `elem` span(w1..w2)
span :: (Eq f, Field f) => Matrix f m n -> Vector f m -> Bool
span m v = all (\x -> pivotLen x /= 1) . unpack . transpose . utf $ append m v'
    where v' = M (Vector @_ @1 (L [v])) -- Forces the matrix to size n 1
          pivotLen xs = length (dropWhile (==zero) xs)

-- | Checks that m1 spans atleast as much as m2 
spans :: (KnownNat m, KnownNat n2, Eq f, Field f) => Matrix f m n1 -> Matrix f m n2 -> Bool
m1 `spans` m2 = all (span m1) (matToList m2)

-- | Checks that m spans the whole vectorSpace
spansSpace :: (KnownNat m, Eq f, Field f) => Matrix f m n -> Bool
spansSpace m = m `spans` idm

-- | Checks if the vectors in a matrix are linearly independant
linIndep :: (Eq f, Field f) => Matrix f m n -> Bool
linIndep = not . any (\(v, m) -> m `span` v ) . separateCols 


-- 2.27
-- Definition: basis
-- A basis of V is a list of vectors in V that is linearly independent and spans V

-- | Checks if the vectors in a matrix forms a basis of their vectorSpace
isBasis :: (KnownNat m, KnownNat (n-1), Eq f, Field f) => Matrix f m n -> Bool
isBasis m = spansSpace m && linIndep m


    
