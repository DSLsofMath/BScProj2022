{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE TypeApplications #-}


module ListVector where

import GHC.TypeLits
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), sum, (**), span)

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

lenTest = V[1,2,3]::VecR 3
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


-- Vector is a vector space over a field
instance (KnownNat n, AddGroup f) => AddGroup (Vector f n) where
    V as + V bs = V $ zipWith (+) as bs
    V as - V bs = V $ zipWith (-) as bs
    zero = zeroVec

instance (KnownNat n, Field f) => VectorSpace (Vector f n) where
    s £ V ss = V $ map (s*) ss

-- | Dot product of two vectors 
v1, v2 :: Vector Double 3
v1 = V [2,7,1]::VecR 3
v2 = V [8,2,8]::VecR 3
-- är dot prod rätt def? ändra till zipWith(*)
-- test dot product
testdot = dot v1 v2 == 2*8 + 7*2 + 8*1

dot :: (Field f) => Vector f n -> Vector f n -> f
V v1 `dot` V v2 = sum $ zipWith (*) v1 v2

-- | Cross product of two vectors of size 3
-- test cross product
testross = cross v1 v2 == V[7*8-1*2, 1*8-2*8, 2*2-7*8]

cross :: (Field f) => Vector f 3 -> Vector f 3 -> Vector f 3
V [a1,a2,a3] `cross` V [b1,b2,b3] = V [a2*b3-a3*b2,
                                       a3*b1-a1*b3,
                                       a1*b2-a2*b1]


-------------------------------------------
-- Matrix definitions

-- | Matrix as a vector of vectors
--   note that the outmost vector is not matimatically a vector 
newtype Matrix f m n = M (Vector (Vector f m) n)  deriving (Eq)

type MatR = Matrix Double

instance Show f => Show (Matrix f m n) where
    show m = let (M (V vs)) = transpose m in unlines $ map (\(V ss) -> show ss) vs


-- | Identity matrix
--id3x3
idm3x3 :: Matrix Double 3 3
idm3x3 = idm :: MatR 3 3

idm :: (KnownNat n, Field f) => Matrix f n n
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

instance (KnownNat m, KnownNat n, Field f) => VectorSpace (Matrix f m n) where
    s £ M (V vs) = M . V $ map (s£) vs


-- | Defines a unified multiplication between scalars, vectors and matrices
class Mult a b c | a b -> c where
    (**) :: a -> b -> c

instance                          Field f  => Mult f              f              f              where (**) = (*)
instance (VectorSpace v, Under v ~ f)      => Mult f              v   v                         where (**) = (£)
instance (KnownNat n, KnownNat m, Field f) => Mult (Matrix f m n) (Vector f n)   (Vector f m)   where (**) = (££)
instance (KnownNat n, KnownNat m, Field f) => Mult (Matrix f m n) (Matrix f n o) (Matrix f m o) where (**) = (£££)



type instance Under (Vector f _)   = f
type instance Under (Matrix f _ _) = f
type instance Under [Vector f _]   = f

-- | Converts objects to and from Matrices.
--   Requires that the object is an instance of the type family Under.
class ToMat m n x where
    toMat   :: x -> Matrix (Under x) m n
    fromMat :: Matrix (Under x) m n -> x

instance (KnownNat n) => ToMat n n Double where
    toMat s = s £ idm
    fromMat (M (V (V (s:_):_))) = s

instance ToMat n 1 (Vector f n) where
    toMat v = M (V [v])
    fromMat (M (V [v])) = v

instance ToMat 1 n (Vector f n) where
    toMat   (V ss) = M . V $ map (\s -> V [s]) ss
    fromMat (M (V vs)) = V $ map (\(V (x:_)) -> x) vs

-- | Diagonal matrix
instance (KnownNat n, Field f) => ToMat n n (Vector f n) where
    toMat (V ss) = M . V $ zipWith (\s i-> s £ e i) ss [1..]
    fromMat m = vec $ zipWith (!!) (fromMat m) [0..]

instance ToMat m n (Matrix f m n) where
    toMat = id
    fromMat = id

instance (KnownNat m, KnownNat n, AddGroup f) => ToMat m n [Vector f m] where
    toMat vs = M . vec $ vs
    fromMat = matToList

instance (KnownNat m, KnownNat n, AddGroup f) => ToMat m n [[f]] where
    toMat = M . vec . map vec
    fromMat = unPack

-- test transpose
m1 :: Matrix Double 3 3
m1 = toMat [[1,1,1],[2,2,2],[3,3,3]]::MatR 3 3
testtran = transpose m1 == (toMat [[1,2,3],[1,2,3],[1,2,3]])

transpose :: Matrix f m n -> Matrix f n m
transpose = pack . L.transpose . unPack

-- test get
testget = get m1 (1,1) == 2      -- m1 row 1, col 1

get :: Matrix f m n -> (Int, Int) -> f
get m (x,y) = unPack m !! x !! y

-- matrix to array and back
testunpack = unPack m1 == [[1,1,1],[2,2,2],[3,3,3]]
testpack = pack(unPack m1) == m1

unPack :: Matrix f m n -> [[f]]
unPack (M (V vs)) = map (\(V a) -> a) vs

pack :: [[f]] -> Matrix f m n
pack = M . V . map V

-- | Analogous to (++)
--   useful for Ex. Guassian elimination
m2 :: Matrix Double 3 1
m2 = toMat[[4,4,4]] :: MatR 3 1

-- append adds m2 as the last col, new matrix = 3 rows, 4 cols
testappend = append m1 m2 == (toMat [[1,1,1],[2,2,2],[3,3,3],[4,4,4]]:: MatR 3 4)

append :: Matrix f m n1 -> Matrix f m n2 -> Matrix f m (n1+n2)
append m1 m2 = pack $ unPack m1 ++ unPack m2

-- | Converts a Matrix to a list of Vectors 
matToList :: Matrix f m n -> [Vector f m]
matToList (M (V vs)) = vs


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


-- | Takes a vector and a basis and returns the linear combination
--   For Example eval [a b c] [^0 ^1 ^2] returns the polinomial \x -> ax + bx + cx^2
eval :: (VectorSpace v, Under v ~ f) => Vector f n -> Vector v n -> v
eval (V fs) (V vs) = sum $ zipWith (£) fs vs


-- | Checks if a vector is in the span of a list of vectors
--   Normaly span is defined as a set, but here we use it as a condition such that
--   span [w1..wn] v = v `elem` span(w1..w2)
span :: (Eq f, Field f) => Matrix f m n -> Vector f m -> Bool
span m v = all (\x -> pivotLen x /= 1) . L.transpose . unPack . utf $ append m v'
    where v' = M (V @_ @1 [v]) -- Forces the matrix to size n 1
          pivotLen xs = length (dropWhile (==zero) xs)

-- | Checks that m1 spans atleast as much as m2 
spans :: (KnownNat m, KnownNat n2, Eq f, Field f) => Matrix f m n1 -> Matrix f m n2 -> Bool
m1 `spans` m2 = all (span m1) (matToList m2)

-- | Checks that m spans the whole vectorSpace
spansSpace :: (KnownNat m, Eq f, Field f) => Matrix f m n -> Bool
spansSpace m = m `spans` idm

-- | For each vector in a matrix, returns the vector and the remaining matrix
separate :: (KnownNat m, KnownNat n, AddGroup f) => Matrix f m n -> [(Matrix f m n, Vector f m)]
separate m = map (\(m',v) -> (toMat m',v)) $ separate' [] (matToList m)
    where separate' :: [Vector f m] -> [Vector f m] -> [([Vector f m], Vector f m)]
          separate' _   []     = []
          separate' acc (x:[]) = [(acc,     x)]
          separate' acc (x:xs) = (acc++xs, x) : separate' (acc++[x]) xs

-- | Checks if the vectors in a matrix are linearly independant
linIndep :: (KnownNat m, KnownNat n, Eq f, Field f) => Matrix f m n -> Bool
linIndep m = not . any (uncurry span) $ separate m


-- 2.27
-- Definition: basis
-- A basis of V is a list of vectors in V that is linearly independent and spans V

-- | Checks if the vectors in a matrix forms a basis of their vectorSpace
basis :: (KnownNat m, KnownNat n, Eq f, Field f) => Matrix f m n -> Bool
basis m = spansSpace m && linIndep m
