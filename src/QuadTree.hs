{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE NoStarIsType #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE MultiParamTypeClasses #-}

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

-- | Singleton type over the kind Nat4
data SNat4 (n :: Nat4) where
    SOne :: SNat4 One 
    SSuc :: Sized n => SNat4 n -> SNat4 (Suc n)

type family Pred (n :: Nat4) :: Nat4 where
     Pred (Suc n)  = n

-- n+n instead of 2*n helps ghc to deduce constraints
type family ToNat (n :: Nat4) :: Nat where
     ToNat One = 1
     ToNat (Suc n)  = ToNat n + ToNat n 


-- | A matrix representation based on QuadTrees
--   Represents a matrix of size n*n where n = 2^Nat4
data Quad (n :: Nat4) a where 
    Zero :: Quad n a
    Scalar :: a -> Quad One a
    Mtx :: Sized n => {
        nw :: Quad n a,
        ne :: Quad n a,
        sw :: Quad n a, 
        se :: Quad n a 
           } -> Quad (Suc n) a 
    -- Id :: Quad n a --used mostly in permutations.

instance (Sized n, Ring a, Show a) => Show (Quad n a) where
    show = show . toDense

-- | A class to get size information from type
class Sized (n :: Nat4) where
    toInt   :: forall proxy. proxy n -> Int
    toSNat4 :: forall proxy. proxy n -> SNat4 n

instance Sized One where
    toInt   _ = 1
    toSNat4 _ = SOne

instance forall n. Sized n => Sized (Suc n) where
    toInt   _ = 2 * toInt (undefined :: undefined n)
    toSNat4 _ = SSuc $ toSNat4 (undefined :: undefined n)

-- | Returns the size of the Quad matrix, since a Quad
--   is always square we only return one value
quadSize :: forall n a. Sized n => Quad n a -> Int
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
    Mtx nw1 ne1 sw1 se1 + Mtx nw2 ne2 sw2 se2 = Mtx (nw1+nw2) (ne1+ne2) 
                                                    (sw1+sw2) (se1+se2)
    neg = fmap neg
    zero = Zero 

instance Ring a => VectorSpace (Quad n a) where
    type Under (Quad n a) = a
    s £ q = fmap (s*) q


-- | Identity matrix in Quad representation
idQ :: forall n a. (Sized n, Ring a) => Quad n a
idQ = case toSNat4 (undefined :: undefined n) of 
    SOne   -> Scalar one
    SSuc _ -> Mtx idQ Zero Zero idQ

-- | Multiplication on Quad matrices
mulQ :: Ring a => Quad n a -> Quad n a -> Quad n a
Zero `mulQ` _ = Zero
_ `mulQ` Zero = Zero
Scalar a `mulQ` Scalar b = Scalar (a * b)
x@(Mtx _ _ _ _) `mulQ` y@(Mtx _ _ _ _) = case zipQuad (+) 
                        (zipQuad mulQ (colExchange x) (offDiagSqsh y))
                        (zipQuad mulQ x               (prmDiagSqsh y))
                of 
                    Mtx Zero Zero Zero Zero -> Zero 
                    quads                   -> quads 
    where colExchange :: Quad n a -> Quad n a
          colExchange (Mtx nw ne sw se) = Mtx ne nw se sw
          prmDiagSqsh :: Quad n a -> Quad n a
          prmDiagSqsh (Mtx nw ne sw se) = Mtx nw se nw se
          offDiagSqsh :: Quad n a -> Quad n a
          offDiagSqsh (Mtx nw ne sw se) = Mtx sw ne sw ne
          zipQuad :: (Quad n a -> Quad n a -> Quad n a) -> Quad (Suc n) a -> Quad (Suc n) a -> Quad (Suc n) a
          zipQuad (*) (Mtx a b c d) (Mtx e f g h) = Mtx (a*e) (b*f) (c*g) (d*h)

instance (Ring a, n ~ n', a ~ a') => Composable (Quad n a) (Quad n' a') where
    type Out (Quad n a) (Quad n' a')  = Quad n a
    comp = mulQ

instance (Sized n, Ring a) => Mul (Quad n a) where
    (*) = mulQ
    one = idQ


-- | Converts a Quad matrix into a list of lists matrix
--
--   TODO: zero requires a KnownNat constraint but we cannot deduce 
--   KnownNat in the recursive call even though ToNat is trivial.
toDense :: (Sized n, Ring a) => Quad n a -> Matrix a (ToNat n) (ToNat n) 
toDense z@Zero = let n = quadSize z in pack $ replicate n (replicate n zero)
toDense (Scalar a) = a £ idm    -- Will always be of size 1x1
toDense (Mtx nw ne sw se) = (toDense nw `append` toDense ne) `append'` 
                            (toDense sw `append` toDense se)

-- TODO a nice example would be to implement toCSN - from one sparse
-- format to another.

-- Potentally use parallelism for big enough Quads
-- class Mul a where
-- 
-- instance (n <= 50) => Mul (Quad n a) where
-- 
-- instance (12 <= n) => Mul (Quad n a) where


----------------
-- Just some example code
type Eight = Suc (Suc (Suc One))
pjQ8 :: Quad Eight Double
pjQ8 = stencil3 1 (-2) 1

-- pjMat8 :: MatR 8 8
-- pjMat8 = toMat pjQ8  -- TODO: I got some type error here. No time to debug now.

class Corner1 n where corner1 :: a -> Quad n a
class Corner2 n where corner2 :: a -> Quad n a
class (Corner1 n, Corner2 n) => Stencil n where stencil3 :: a -> a -> a -> Quad n a

instance Corner1 One where corner1 = Scalar
instance Corner2 One where corner2 = Scalar
instance Stencil One where stencil3 _ b _ = Scalar b

instance (Sized n, Corner1 n) => Corner1 (Suc n) where corner1 = corner1Suc
corner1Suc :: (Sized n, Corner1 n) => a -> Quad (Suc n) a
corner1Suc c = Mtx Zero Zero (corner1 c) Zero

instance (Sized n, Corner2 n) => Corner2 (Suc n) where corner2 = corner2Suc
corner2Suc :: (Sized n, Corner2 n) => a -> Quad (Suc n) a
corner2Suc c = Mtx Zero (corner2 c) Zero Zero

instance (Sized n, Stencil n) => Stencil (Suc n) where stencil3 = stencil3Suc
stencil3Suc :: (Sized n, Stencil n) => a -> a -> a -> Quad (Suc n) a
stencil3Suc a b c = Mtx  (stencil3 a b c)  (corner1 c)
                         (corner2 a)  (stencil3 a b c)
