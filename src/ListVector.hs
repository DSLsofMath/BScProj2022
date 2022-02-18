{-# LANGUAGE DataKinds #-} 
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeOperators #-} 
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module ListVector where

import GHC.TypeLits
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), sum)

import qualified Data.List as L 
import Algebra

-- This file contains an example of using TypeLits to handle vectors with a given size. 
-- The implementation here is based on lists and should be replaced.
--
-- To try the code in ghci the DataKinds language extension needs to be enabled:
-- type ":set -XDataKinds" in the ghci prompt
--

-------------------------------------------
-- Vector definitions

-- | Dependent typed vector
newtype Vector f (n :: Nat) = V [f] deriving (Show, Eq)

type VecR = Vector Double

vecLen :: KnownNat n => Vector f n -> Int
vecLen v = fromInteger $ natVal v

-- | Vector that is all zero
zeroVec :: (KnownNat n, AddGroup f) => Vector f n
zeroVec = let v = V $ replicate (vecLen v) zero in v

-- | converts a list to a vector with a typed size
vec :: (KnownNat n, AddGroup f) => [f] -> Vector f n
vec ss = V (ss ++ repeat zero) + zero

-- | e i is the i:th basis vector
e :: (KnownNat n, Field f) => Int -> Vector f n
e i = V (replicate (i-1) zero ++ (one : repeat zero)) + zero


-- | Definition of a VectorSpace
class (Field s, AddGroup v) => VectorSpace v s where
    (£) :: s -> v -> v

-- Vector is a vector space over a field
instance (KnownNat n, AddGroup f) => AddGroup (Vector f n) where
    V as + V bs = V $ zipWith (+) as bs
    V as - V bs = V $ zipWith (-) as bs
    zero = zeroVec

instance (KnownNat n, Field f) => VectorSpace (Vector f n) f where
    s £ V ss = V $ map (s*) ss

-- | Dot product of two vectors 
dot :: (Field f) => Vector f n -> Vector f n -> f
V v1 `dot` V v2 = sum $ zipWith (+) v1 v2

-- | Cross product of two vectors of size 3
cross :: (Field f) => Vector f 3 -> Vector f 3 -> Vector f 3
V [a1,a2,a3] `cross` V [b1,b2,b3] = V [a2*b3-a3*b2, 
                                       a3*b1-a1*b3, 
                                       a1*b2-a2*b1]


-------------------------------------------
-- Matrix definitions

-- | Matrix as a vector of vectors
--   note that the outmost vector is not matimatically a vector 
newtype Matrix f m n = M (Vector (Vector f m) n)  deriving (Eq)

instance Show f => Show (Matrix f m n) where
    show m = let (M (V vs)) = transpose m in unlines $ map (\(V ss) -> show ss) vs


type MatR = Matrix Double

-- | Identity matrix
idm :: (KnownNat n, KnownNat m, Field f) => Matrix f m n
idm = let v = V [ e i | i <- [1 .. vecLen v] ] in M v


-- | Matrix vector multiplication
(££) :: (KnownNat m, Field f) => 
                    Matrix f m n -> Vector f n -> Vector f m 
M (V vs) ££ (V ss) = sum $ zipWith (£) ss vs

-- | Matrix matrix multiplication
(£££) :: (KnownNat a, Field f) => 
                    Matrix f a b -> Matrix f b c -> Matrix f a c
m £££ M (V vs) = M . V $ map (m££) vs


-- Matrices also forms a vector space over a field
instance (KnownNat m, KnownNat n, AddGroup f) => AddGroup (Matrix f m n) where
    M (V as) + M (V bs) = M . V $ zipWith (+) as bs
    M (V as) - M (V bs) = M . V $ zipWith (-) as bs
    zero = let v = V $ replicate (vecLen v) zero in M v

instance (KnownNat m, KnownNat n, Field f) => VectorSpace (Matrix f m n) f where
    s £ M (V vs) = M . V $ map (s£) vs



-- Convert objects to Matrices
class ToMat f m n x where
    toMat :: x -> Matrix f m n

instance ToMat Double 1 1 Double where
    toMat s = M (V [V [s]])

instance ToMat f n 1 (Vector f n) where
    toMat v = M (V [v])

instance ToMat f 1 n (Vector f n) where
    toMat (V ss) = M . V $ map (\s -> V [s]) ss

instance ToMat f m n (Matrix f m n) where
    toMat = id
    
instance (KnownNat m, KnownNat n, AddGroup f) => ToMat f m n [[f]] where
    toMat ls = M . vec $ map vec ls

-- Convert Matrices to objects 
matToVec :: Matrix f n 1 -> Vector f n
matToVec (M (V [v])) = v

matToField :: Matrix f 1 1 -> f
matToField (M (V (V (s:_):_))) = s


transpose :: Matrix f m n -> Matrix f n m
transpose (M (V a)) = M . V $ map V $ L.transpose $ map (\(V a) -> a) a

get :: Matrix f m n -> (Int, Int) -> f
get (M (V a)) (x,y) = map (\(V a) -> a) a !! x !! y

unPack :: Matrix f m n -> [[f]]
unPack (M (V a)) = map (\(V a) -> a) a 

pack :: [[f]] -> Matrix f m n 
pack = M . V . map V


utf :: (Eq f, Field f) => Matrix f m n -> Matrix f m n
utf = transpose . pack . sort . f . unPack . transpose 
    where 
          f []     = []
          f (x:[]) = let (x', n) = pivot x in [x']
          f (x:xs) = let (x', n) = pivot x in x' : (f $ map (reduce n x') xs)
          pivot [] = ([],-1) 
          pivot (x:xs) | x == zero = let (ys, n) = pivot xs in (x:ys, n + 1)
                       | otherwise = (map (/x) (x:xs), 0 :: Int)
          reduce n p x = zipWith (-) x (map ((x!!n)*) p)
          sort = L.sortOn (length . takeWhile (==zero))

