{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
-- {-# LANGUAGE StandaloneKindSignatures #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE NoStarIsType #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}

module QuadTree where

import Prelude hiding ((+), (-), (*), (/), sum, )
import GHC.TypeLits
import Data.Kind

import ListVector
import Algebra

{- Taken from
 - Matrix Algorithms using Quadtrees Invited Talk, ATABLE-92 Technical Report 357

type Quadrants a = [Matrx a] --list of *four* submatrices.

data Matrx a = Zero | ScalarM a | Mtx (Quadrants a) | ID --used mostly in permutations.
    deriving Show

instance (Num a) => Num (Matrx a) where
    fromInteger 0 = Zero
    fromInteger 1 = ID
    fromInteger _ = error "fromInteger defined only on zero/1 Matrces."
    negate Zero = Zero
    negate (ScalarM x) = ScalarM (negate x)
    negate (Mtx quads) = Mtx (map negate quads)

    x + Zero = x
    Zero + y = y
    ScalarM x + ScalarM y = ScalarM (x + y)    
    Mtx x + Mtx y = case zipWith (+) x y of [Zero,Zero,Zero,Zero] -> Zero
                                            quads     -> Mtx quads
    _ + _ = error "Matrx addition undefined on IdentM"

    Zero * _ = Zero
    _ * Zero = Zero
    ID * y = y --NB: multiplication accepts IdentM
    x * ID = x
    ScalarM x * ScalarM y = ScalarM (x*y)
    --Except with infinitesimal floats: case x*y of zero->zero; z->ScalarM z
    Mtx x * Mtx y = case zipWith (+)
                            (zipWith (*) (colExchange x)(offDiagSqsh y))
                            (zipWith (*) x (prmDiagSqsh y))
                    of 
                        [Zero,Zero,Zero,Zero] -> Zero 
                        quads                 -> Mtx quads 
        where colExchange [nw,ne,sw,se] = [ne,nw,se,sw]
              prmDiagSqsh [nw,ne,sw,se] = [nw,se,nw,se]
              offDiagSqsh [nw,ne,sw,se] = [sw,ne,sw,ne]
              abs _ = error "abs not implemented yet on Matrces."
              signum _ = error "signum not implemented yet on Matrces."

-}

data Nat4 = One | Suc Nat4

type family Pred (n :: Nat4) :: Nat4 where
     Pred (Suc n)  = n

-- n+n instead of 2*n helps ghc to deduce constraints
type family ToNat (n :: Nat4) :: Nat where
     ToNat One = 1
     ToNat (Suc n)  = ToNat n + ToNat n 


-- Represents a matrix of size n*n where n = 2^Nat4
data Quad (n :: Nat4) a where 
    Zero :: Quad n a
    Scalar :: a -> Quad One a
    Mtx :: ToInt n => {
        nw :: Quad n a,
        ne :: Quad n a,
        sw :: Quad n a, 
        se :: Quad n a 
           } -> Quad (Suc n) a 
    -- Id :: Quad n a --used mostly in permutations.

class ToInt (n :: Nat4) where
    toInt :: forall proxy. proxy n -> Int

instance ToInt One where
    toInt _ = 1

instance forall n. ToInt n => ToInt (Suc n) where
    toInt _ = 2 * toInt (undefined :: undefined n)

-- | Returns the size of the Quad matrix, since a Quad
--   is always square we only return one value
quadSize :: forall n a. ToInt n => Quad n a -> Int
quadSize _ = toInt (undefined :: undefined n)

-- | Applies a function on the Quads scalars
mapQuad :: (a -> b) -> Quad n a -> Quad n b
mapQuad f Zero = Zero
mapQuad f (Scalar s) = Scalar (f s)
mapQuad f (Mtx nw ne sw se) = Mtx (mapQuad f nw) (mapQuad f ne) 
                                  (mapQuad f sw) (mapQuad f se)

instance Functor (Quad n) where
   fmap = mapQuad 

-- Quad is a vector
instance AddGroup a => AddGroup (Quad n a) where
    Zero + x    = x
    x    + Zero = x
    Scalar a + Scalar b = Scalar (a+b)
    Mtx nw1 ne1 sw1 se1 + Mtx nw2 ne2 sw2 se2 = Mtx (nw1+nw2) (nw1+nw2) 
                                                    (nw1+nw2) (nw1+nw2)
    neg = fmap neg
    zero = Zero 

instance Ring a => VectorSpace (Quad n a) where
    type Under (Quad n a) = a
    s £ q = fmap (s*) q


-- | Converts a Quad matrix into a list of lists matrix
--
--   TODO: zero requires a KnownNat constraint but we cannot deduce 
--   KnownNat in the recursive call even though ToNat is trivial.
toDense :: (ToInt n, Ring a) => Quad n a -> Matrix a (ToNat n) (ToNat n) 
toDense z@Zero = let n = quadSize z in pack $ replicate n (replicate n zero)
toDense (Scalar a) = a £ idm    -- Will always be of size 1x1
toDense (Mtx nw ne sw se) = (toDense nw `append` toDense ne) `append'` 
                            (toDense sw `append` toDense se)



-- Potentally use parallelism for big enough Quads
-- class Mul a where
-- 
-- instance (n <= 50) => Mul (Quad n a) where
-- 
-- instance (12 <= n) => Mul (Quad n a) where


