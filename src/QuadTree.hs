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
{-# LANGUAGE ExplicitForAll #-}

module QuadTree where

import Prelude hiding ((+), (-), (*), (/), sum, )
import GHC.TypeLits
import Data.Type.Bool
import Data.Kind

import Data.List

import qualified ListVector as L
import Matrix
import Algebra


data Nat4 = One | Suc Nat4

-- | Singleton type over the kind Nat4
data SNat4 (n :: Nat4) where
    SOne :: SNat4 One 
    SSuc :: Sized n => SNat4 n -> SNat4 (Suc n)

-- n+n instead of 2*n helps ghc to deduce constraints
type family ToNat (n :: Nat4) :: Nat where
     ToNat One = 1
     ToNat (Suc n)  = ToNat n + ToNat n



type family ToNat4 (m :: Nat) (n :: Nat) :: Nat4 where
            ToNat4 m n = If ( Max m n <=? 2 ^ Log2 (Max m n) ) (ToNat4' (Max m n))
                                                          (Suc (ToNat4' (Max m n)))

type family ToNat4' (n :: Nat) :: Nat4 where
            -- ToNat4' 0 =  -- We should handle case 0
            ToNat4' 1 = One
            ToNat4' n = Suc (ToNat4' (n `Div` 2))

type family Max (m :: Nat) (n :: Nat) :: Nat where
            Max n n = n
            Max m n = If ( m <=? n ) n m



-- | A wrapper for Quad, needed for Matrix instance
data QuadM f (m :: Nat) (n :: Nat) = Sized (ToNat4 m n) => QuadM (Quad (ToNat4 m n) f)

instance (AddGroup (L.Matrix f m n), Show (L.Matrix f m n)) => Show (QuadM f m n) where
    show = show . toMat
        where toMat :: QuadM f m n -> L.Matrix f m n
              toMat = changeRep
    


instance (Sized (ToNat4 m n), AddGroup f) => AddGroup (QuadM f m n) where
    QuadM m1 + QuadM m2 = QuadM (m1 + m2)
    QuadM m1 - QuadM m2 = QuadM (m1 - m2)
    zero = QuadM zero

instance (Eq f, AddGroup f, AddGroup (QuadM f m n)) => Eq (QuadM f m n) where
    q1 == q2 = let (QuadM q1', QuadM q2') = (purge q1, purge q2) in q1' == q2'
    

instance Matrix QuadM where 
    -- get (QuadM q) (Fin i, Fin j) = getQ q (i-1, j-1)

    values (QuadM q) = map toFin $ valuesQ q 
        where toFin ((i,j),a) = ((Fin (i+1), Fin (j+1)), a)

    extend (QuadM m) xs = QuadM $ setMultipleQ m (map removeFin xs)
        where removeFin ((Fin i, Fin j), a) = ((i-1, j-1), a) 

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

instance (Sized n, Ring a, Show a) => Show (Quad n a) where
    show = show . toDense

instance (Eq a) => Eq (Quad n a) where
    Zero == Zero = True
    Scalar a == Scalar b = a == b
    Mtx nw1 ne1 sw1 se1 == Mtx nw2 ne2 sw2 se2 = nw1 == nw2 && ne1 == ne2 
                                              && sw1 == sw2 && se1 == se2 
    _ == _ = False

getQ :: AddGroup a => Quad n a -> (Int, Int) -> a
getQ Zero       _ = zero
getQ (Scalar a) _ = a
getQ q@(Mtx nw ne sw se) (i, j) = getQ quad (i `mod` m, j `mod` m)
    where m = quadSize q `div` 2
          quad = case (i < m, j < m) of (True,  True) -> nw; (True,  False) -> ne 
                                        (False, True) -> se; (False, False) -> sw 

valuesQ :: Quad n a -> [((Int, Int), a)]
valuesQ Zero = []
valuesQ (Scalar s) = [((0,0), s)]
valuesQ q@(Mtx nw ne sw se) = concat [s (0,0) nw, s (0,m) ne, s (m,0) sw, s (m,m) se] 
    where m = quadSize q `div` 2
          s (i,j) = map (\((i',j'),a) -> ((i'+i, j'+j),a)) . valuesQ


setQ :: Sized n => Quad n a -> (Int, Int) -> a -> Quad n a
setQ q i a = setMultipleQ q [(i,a)]

setMultipleQ :: Sized n => Quad n a -> [((Int, Int), a)] -> Quad n a
setMultipleQ q = setMultipleQ' (sNatQ q) q

setMultipleQ' :: SNat4 n -> Quad n a -> [((Int, Int), a)] -> Quad n a
setMultipleQ' _    a [] = a 
setMultipleQ' SOne _ ((_, a):_) = Scalar a 
setMultipleQ' (SSuc n) a xs = Mtx (setM nw nws) (setM ne nes) 
                                  (setM sw sws) (setM se ses)
    where ((nws, nes), (sws, ses)) = splitOn snd `onBoth` splitOn fst xs
          ( nw,  ne,    sw,  se)   = case a of Zero            -> (Zero, Zero, Zero,Zero)
                                               Mtx nw ne sw se -> (nw,   ne,   sw,   se)
          splitOn f = span (\(i,_) -> f i < m) . sortOn (f . fst) -- TODO: get rid of sort
          f `onBoth` (a,b) = (f a, f b)
          setM q = setMultipleQ' n q . map (\(i,a) -> ((`mod` m) `onBoth` i, a)) 
          m = toInt n


-- | The height of the quadtree
height :: Quad n a -> Int
height Zero = 0
height (Scalar _) = 0
height (Mtx nw ne sw se) = 1 + maximum [height nw, height ne,height sw,height se]

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

sNatQ :: forall n a. Sized n => Quad n a -> SNat4 n
sNatQ _ = toSNat4 (undefined :: undefined n)

-- | Applies a function on the Quads scalars
mapQuad :: (a -> b) -> Quad n a -> Quad n b
mapQuad f Zero = Zero
mapQuad f (Scalar s) = Scalar (f s)
mapQuad f (Mtx nw ne sw se) = Mtx (mapQuad f nw) (mapQuad f ne) 
                                  (mapQuad f sw) (mapQuad f se)

addQ :: AddGroup a => Quad n a -> Quad n a -> Quad n a
Zero     `addQ` x        = x
x        `addQ` Zero     = x
Scalar a `addQ` Scalar b = Scalar (a+b)
Mtx nw1 ne1 sw1 se1 `addQ` Mtx nw2 ne2 sw2 se2 = Mtx (nw1+nw2) (ne1+ne2) 
                                                     (sw1+sw2) (se1+se2)

instance Functor (Quad n) where
   fmap = mapQuad 

-- Quad is a vector
instance AddGroup a => AddGroup (Quad n a) where
    (+) = addQ
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

mulQM :: (Ring f, ToNat4 m n ~ ToNat4 n c, ToNat4 m c ~ ToNat4 n c) => QuadM f m n -> QuadM f n c -> QuadM f m c
mulQM (QuadM q1) (QuadM q2) = QuadM $ mulQ q1 q2

instance (Sized n, Ring a) => Mul (Quad n a) where
    (*) = mulQ
    one = idQ

instance (Ring f, f ~ f', b ~ b', ToNat4 a b ~ ToNat4 b' c, ToNat4 b' c ~ ToNat4 a c) => Composable (QuadM f a b) (QuadM f' b' c) (QuadM f a c) where
   (**) = (mulQM)

-- | Transposes a Quad matrix
transposeQ :: Quad n a -> Quad n a
transposeQ Zero = Zero
transposeQ (Scalar s) = Scalar s
transposeQ (Mtx nw ne sw se) = Mtx (transposeQ nw) (transposeQ sw)
                                   (transposeQ ne) (transposeQ se)

-- | Converts a Quad matrix into a list of lists matrix
--
--   TODO: zero requires a KnownNat constraint but we cannot deduce 
--   KnownNat in the recursive call even though ToNat is trivial.
toDense :: (Sized n, Ring a) => Quad n a -> L.Matrix a (ToNat n) (ToNat n) 
toDense z@Zero = let n = quadSize z in L.pack $ replicate n (replicate n zero)
toDense (Scalar a) = a £ L.idm    -- Will always be of size 1x1
toDense (Mtx nw ne sw se) = (toDense nw `L.append` toDense ne) `L.append'` 
                            (toDense sw `L.append` toDense se)


-- Potentally use parallelism for big enough Quads
-- class Mul a where
-- 
-- instance (n <= 50) => Mul (Quad n a) where
-- 
-- instance (12 <= n) => Mul (Quad n a) where


testQ :: Quad (Suc (Suc One)) R
testQ = Mtx (Mtx Zero (Scalar 3) (Scalar 5) (Scalar 7)) (Mtx Zero       (Scalar 6) (Scalar 5) Zero) 
            Zero                                        (Mtx (Scalar 2) Zero       (Scalar 7) Zero)


----------------
-- Just some example code
type Eight = Suc (Suc (Suc One))
pjQ8 :: Quad Eight Double
pjQ8 = stencil3 1 (-2) 1

pj10000 :: Quad (Suc (Suc (Suc (Suc (Suc (Suc (Suc (Suc (Suc (Suc (Suc (Suc (Suc One))))))))))))) Double
pj10000 = stencil3 1 (-2) 1

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
