{-# LANGUAGE TypeFamilies #-}

module FiniteIndex where

import GHC.TypeLits
import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/))

import Data.List (sortOn)

import Algebra


-- | Bounded integer
--   A Fin n should always be in the range of (1,n) 
newtype Fin (n :: Nat) = Fin Int deriving (Eq, Ord) 

-- | Safe constructor of Fin, checks that the value is in range
fin' :: forall n. KnownNat n => Int -> Maybe (Fin n)
fin' i | 1 <= i && i <= natInt' @n = Just $ Fin i
      | otherwise                 = Nothing

fin :: forall n. KnownNat n => Int -> Fin n
fin i = case fin' i of 
    Just fi -> fi
    Nothing -> error $ unwords ["Index is out of bounds, got: ", show i, 
                                " in constraint 0<", show i, "<=", show (natInt' @n)]

instance KnownNat n => Show (Fin n) where
    show n = show (finToInt n) ++ " of " ++ show (natVal n)

-- | Num (Fin n) is only used to get fromInteger
instance KnownNat n => Num (Fin n) where
    fromInteger = fin . fromInteger
    (+) = undefined
    (*) = undefined
    abs = undefined 
    signum = undefined 
    negate = undefined

-- Enum and Bounded allows for list sequencing 
-- i.e. [Fin 2 .. maxBound]
instance KnownNat n => Enum (Fin n) where
    toEnum   = fin
    fromEnum = finToInt
    
instance KnownNat n => Bounded (Fin n) where
    minBound = Fin 1
    maxBound = Fin (natInt' @n)

indices :: KnownNat n => [Fin n] 
indices = [minBound .. maxBound]

-- | Returns the corresponding Int
finToInt :: Fin n -> Int
finToInt (Fin n) = n

-- | Increases the max bound of the index
raiseBound :: n <= m => Fin n -> Fin m
raiseBound (Fin n) = Fin n

-- | Like raiseBound but also adds the difference to the index
--   pushBound (3 :: Fin 4) :: Fin 6 => 5 :: Fin 6
pushBound :: forall m n. (KnownNat (m - n), n <= m) => Fin n -> Fin m
pushBound (Fin i) = Fin (diff + i)
    where diff = natInt' @(m - n)

-- | Merges two indexed lists, adding their values if they have the same index
addL :: (Ord i, AddGroup a) => [(i, a)] -> [(i, a)] -> [(i, a)] 
addL xs ys = merge (sortOn fst xs) (sortOn fst ys)
    where merge [] bs = bs
          merge as [] = as
          merge ((i, a):as) ((j, b):bs) | i == j    = (i, a + b) : merge as bs
                                        | i <= j    = (i, a)     : merge as ((j, b):bs)
                                        | otherwise = (j, b)     : merge ((i, a):as) bs

scaleL :: (Ord i, Mul a) => a -> [(i, a)] -> [(i, a)] 
scaleL s l = [ (i, s*a) | (i,a) <- l ]

instance (Ord i, Ring a) => AddGroup [(i, a)] where
    (+) = addL
    neg = scaleL (neg one)
    zero = []

instance (Ord i, Ring a) => VectorSpace [(i,a)] where
    type Under [(i,a)] = a
    (£) = scaleL



