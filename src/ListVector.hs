{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE TypeApplications #-}


module ListVector where

import GHC.TypeLits hiding (type (^))
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), recip, sum, product, (**), span)

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

-- PROPOSAL: change the type of vector to be 
-- Vector f n = Matrix f n 1
-- to make the code more general

-- | Dependent typed vector
newtype Vector f (n :: Nat) = V [f] deriving (Show, Eq)

type VecR = Vector Double

-- | Nice type notation for vectors
--   v = vec [1,2,3,4] :: R^4
type f ^ n = Vector f n

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
    type Under (Vector f n) = f
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
    show = showMat

-- Show function for matrices 
showMat :: (Show f) => Matrix f m n -> String
showMat m = unlines $ map (\s -> "| " ++ unwords s ++ "|") $ L.transpose $ map padCol $ unpack m
    where
        getLongest = maximum . map (length)
        padString n s = (replicate (n-length s) ' ') ++ s ++ " "
        padCol l = let s = map show l in map (padString (getLongest s)) s



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
    type Under (Matrix f m n) = f
    s £ M (V vs) = M . V $ map (s£) vs


-- | Defines a unified multiplication between scalars, vectors and matrices
class Mult a b c | a b -> c where
    (**) :: a -> b -> c

instance                          Field f  => Mult f              f              f              where (**) = (*)
instance (VectorSpace v, Under v ~ f)      => Mult f              v   v                         where (**) = (£)
instance (KnownNat n, KnownNat m, Field f) => Mult (Matrix f m n) (Vector f n)   (Vector f m)   where (**) = (££)
instance (KnownNat n, KnownNat m, Field f) => Mult (Matrix f m n) (Matrix f n o) (Matrix f m o) where (**) = (£££)



-- | Converts objects to and from Matrices.
--   Requires that the object is an instance of the type family Under.
class ToMat m n x where
    type Under' x 
    toMat   :: x -> Matrix (Under' x) m n
    fromMat :: Matrix (Under' x) m n -> x

instance (KnownNat n) => ToMat n n Double where
    type Under' Double = Double
    toMat s = s £ idm
    fromMat (M (V (V (s:_):_))) = s

instance ToMat n 1 (Vector f n) where
    type Under' (Vector f n) = f
    toMat v = M (V [v])
    fromMat (M (V [v])) = v

instance ToMat 1 n (Vector f n) where
    type Under' (Vector f n) = f
    toMat   (V ss) = M . V $ map (\s -> V [s]) ss
    fromMat (M (V vs)) = V $ map (\(V (x:_)) -> x) vs

-- | Diagonal matrix
instance (KnownNat n, Field f) => ToMat n n (Vector f n) where
    type Under' (Vector f n) = f
    toMat (V ss) = M . V $ zipWith (\s i-> s £ e i) ss [1..]
    fromMat m = vec $ zipWith (!!) (fromMat m) [0..]

instance ToMat m n (Matrix f m n) where
    type Under' (Matrix f m n) = f
    toMat = id
    fromMat = id

instance (KnownNat m, KnownNat n, AddGroup f) => ToMat m n [Vector f m] where
    type Under' [Vector f m] = f
    toMat vs = M . vec $ vs
    fromMat = matToList

instance (KnownNat m, KnownNat n, AddGroup f) => ToMat m n [[f]] where
    type Under' [[f]] = f
    toMat = M . vec . map vec
    fromMat = unpack



-- test transpose
m1 :: Matrix Double 3 3
m1 = toMat [[1,1,1],[2,2,2],[3,3,3]]::MatR 3 3
testtran = transpose m1 == (toMat [[1,2,3],[1,2,3],[1,2,3]])

transpose :: Matrix f m n -> Matrix f n m
transpose = pack . L.transpose . unpack

-- test get
testget = get m1 (1,1) == 2      -- m1 row 1, col 1

get :: Matrix f m n -> (Index, Index) -> f
get m (x,y) = unpack m !! x !! y

-- matrix to array and back
testunpack = unpack m1 == [[1,1,1],[2,2,2],[3,3,3]]
testpack = pack(unpack m1) == m1

unpack :: Matrix f m n -> [[f]]
unpack (M (V vs)) = map (\(V a) -> a) vs

pack :: [[f]] -> Matrix f m n
pack = M . V . map V

-- | Analogous to (++)
--   useful for Ex. Guassian elimination
m2 :: Matrix Double 3 1
m2 = toMat[[4,4,4]] :: MatR 3 1

-- append adds m2 as the last col, new matrix = 3 rows, 4 cols
testappend = append m1 m2 == (toMat [[1,1,1],[2,2,2],[3,3,3],[4,4,4]]:: MatR 3 4)

append :: Matrix f m n1 -> Matrix f m n2 -> Matrix f m (n1+n2)
append m1 m2 = pack $ unpack m1 ++ unpack m2

-- | Converts a Matrix to a list of Vectors 
matToList :: Matrix f m n -> [Vector f m]
matToList (M (V vs)) = vs


-------------------------------------------
-- Elementary row operations definition 
-- Reduction and equation solver functions

type Index = Int 
-- | Represents elementary row operations
data ElimOp a = Swap Index Index 
              | Mul Index a 
              | MulAdd Index Index a
              deriving (Eq)

instance Show a => Show (ElimOp a) where show = showElimOp

-- | Prettier show function for ElimOp a
showElimOp :: Show a => ElimOp a -> String
showElimOp op = concat $ case op of 
                  Swap    i j   -> [         row i,              " <-> ", row j ]
                  Mul     i   s -> [ show s, row i,               " -> ", row i ]
                  MulAdd  i j s -> [ show s, row i, " + ", row j, " -> ", row j ]
    where row i = "R(" ++ show i ++ ")"


-- | Shows step by step how a matrix is transformed by a ElimOp trace
showElimOnMat :: (Field f, Show f) => [ElimOp f] -> Matrix f m n -> String
showElimOnMat t m0 = let matTrace = scanl (flip elimOpToFunc) m0 t 
                     in unlines [ show m ++ "\n" ++ show op | (m, op) <- zip matTrace t ]
                        ++ show (last matTrace)

-- | Representation of an elementary row operation as a matrix 
elimOpToMat :: (KnownNat n, Field f) => ElimOp f -> Matrix f n n
elimOpToMat (Swap x y) = let v = V [ e' i | i <- [1 .. vecLen v] ] in M v
    where e' i | i == x    = e y
               | i == y    = e x
               | otherwise = e i
elimOpToMat (Mul x n)  = let v = V [ e' i | i <- [1 .. vecLen v] ] in M v
    where e' i | i == x    = n £ e i
               | otherwise = e i
elimOpToMat (MulAdd x y n) = let v = V [ e' i | i <- [1 .. vecLen v] ] in M v
    where e' i | i == x    = e i + (n £ e y)
               | otherwise = e i

m :: Matrix Double 3 4
m = toMat [[2, -3, -2],[1, -1, 1],[-1, 2, 2],[8, -11, -3]]

es :: [ElimOp Double]
es = [MulAdd 1 2 (3/2), MulAdd 1 3 1, MulAdd 2 3 (-4), MulAdd 3 2 (1/2), MulAdd 3 1 (-1), Mul 2 2, Mul 3 (-1), MulAdd 2 1 (-1), Mul 1 (1/2)]

foldElemOps :: (Field f, KnownNat n) => [ElimOp f] -> Matrix f n n
foldElemOps = foldr (flip (£££)) idm . map elimOpToMat 


-- | Representation of an elementary row operation as a function 
--   The reason that we unpack all the way to list is to avoid a KnownNat constraint
elimOpToFunc :: Field f => ElimOp f -> (Matrix f m n -> Matrix f m n)
elimOpToFunc e = pack . L.transpose  . f . L.transpose . unpack
    where f m = case e of
            Swap i j     -> let 
                            (i', j') = (min i j - 1, max i j - 1)
                            (m1,x:xs) = splitAt i' m
                            (xs',y:m2) = splitAt (j' - i'-1) xs
                            in m1++y:xs'++x:m2
            Mul i s      -> let
                            (m1,x:m2) = splitAt (i-1) m
                            in m1++(map (s*) x): m2
            MulAdd i j s -> let 
                            (i', j') = (i - 1, j - 1)
                            (_,x:_) = splitAt i' m
                            (m1,y:m2) = splitAt j' m
                            y' = zipWith (+) y (map (s*) x)
                            in m1++y':m2

-- | Reduces a trace of elimOps to a single function
--   TODO: We should add a rewrite rule such that fold elemOpToFunc only packs and unpacks once
foldElemOpsFunc :: Field f => [ElimOp f] -> (Matrix f m n -> Matrix f m n)
foldElemOpsFunc = foldr (flip (.)) id . map elimOpToFunc


-- | Transform a matrix to upper triangular form
utf :: (Eq f, Field f) => Matrix f m n -> Matrix f m n
utf = transpose . pack . sort . f . unpack . transpose
    where
          f []     = []
          f (x:[]) = let (x', n) = pivot x in [x']
          f (x:xs) = let (x', n) = pivot x in x' : (f $ map (reduce n x') xs)
          pivot [] = ([],-1)
          pivot (x:xs) | x == zero = let (ys, n) = pivot xs in (x:ys, n + 1)
                       | otherwise = (map (/x) (x:xs), 0 :: Int)
          reduce n p x = zipWith (-) x (map ((x!!n)*) p)
          sort = L.sortOn (length . takeWhile (==zero))

-- | Generate a trace of ElimOps from reducing a matrix to upper triangular form
utfTrace :: (Field f, Eq f) => Matrix f m n -> [ElimOp f]
utfTrace m0 = case separateCols m0 of
      []               -> []
      (V (x:xs), m):_  -> let 
                      trace  = mul x ++ [ mulAdd s j | (s,j) <- zip xs [2..], s /= zero ]
                      augM = foldElemOpsFunc trace m
                      in trace ++ case (m, separateRows augM) of
                         (M (V (_:_)), (_,m'):_ ) -> map incIndex $ utfTrace m'
                         _                        -> []
    where
        mulAdd s j = MulAdd 1 j (neg s)
        mul s = if s /= one then [Mul 1 (recip s)] else []
        
        incIndex :: ElimOp f -> ElimOp f
        incIndex (Swap i j)     = Swap   (i+1) (j+1)
        incIndex (Mul i s)      = Mul    (i+1)       s
        incIndex (MulAdd i j s) = MulAdd (i+1) (j+1) s

-- !! Work in progress
-- Doesn't work for m n matricies where m is larger than n
ref :: (KnownNat m, Eq f, Field f) => Matrix f m n -> [ElimOp f]
ref mat = ref' mat 0

ref' :: (KnownNat m, Eq f, Field f) => Matrix f m n -> Index -> [ElimOp f]
ref' mat i | i >= length lm = []
           | otherwise     = do
                let q = pivot (lm!!i)
                let mulOp = Mul (i+1) (recip (lm!!i!!q))
                let addOp = [MulAdd (i+1) (j+1) (neg (lm!!j!!q)) | j <- [(i+1)..(length lm-1)]]
                let t = mulOp:addOp
                t ++ ref' (foldElemOps t £££ mat) (i+1)
                where
                    lm = unpack $ transpose mat

pivot :: (Eq f, Field f) => [f] -> Int
pivot [] = 0
pivot (x:xs) | x == zero = 1 + pivot xs
             | otherwise = 0


jordan :: (Eq f, Field f) => Matrix f m n -> [ElimOp f]
jordan = undefined


-- | Solve systems of equations
--   Check in Demo.hs how to use
solve :: (Field f) => [[f]] -> [f]
solve m = foldr next [last (last m)] (init m)
    where
        next row found = let
            subpart = init $ drop (length m - length found) row
            solved = last row - sum (zipWith (*) found subpart)
            in solved : found

-- | apply when solving systems of equations
--   each element in list represents variable values
solvesys :: (Eq f, Field f) => Matrix f n m -> [f]
solvesys = solve . unpack . transpose . utf


-------------------------------------------
-- Properties on Finite-Dimensional Vector Spaces
-- Span, linear independence and basis

-- | Takes a vector and a basis and returns the linear combination
--   For Example eval [a b c] [^0 ^1 ^2] returns the polinomial \x -> ax + bx + cx^2
eval :: (VectorSpace v, Under v ~ f) => Vector f n -> Vector v n -> v
eval (V fs) (V vs) = sum $ zipWith (£) fs vs


-- | Checks if a vector is in the span of a list of vectors
--   Normaly span is defined as a set, but here we use it as a condition such that
--   span [w1..wn] v = v `elem` span(w1..w2)
span :: (Eq f, Field f) => Matrix f m n -> Vector f m -> Bool
span m v = all (\x -> pivotLen x /= 1) . L.transpose . unpack . utf $ append m v'
    where v' = M (V @_ @1 [v]) -- Forces the matrix to size n 1
          pivotLen xs = length (dropWhile (==zero) xs)

-- | Checks that m1 spans atleast as much as m2 
spans :: (KnownNat m, KnownNat n2, Eq f, Field f) => Matrix f m n1 -> Matrix f m n2 -> Bool
m1 `spans` m2 = all (span m1) (matToList m2)

-- | Checks that m spans the whole vectorSpace
spansSpace :: (KnownNat m, Eq f, Field f) => Matrix f m n -> Bool
spansSpace m = m `spans` idm

-- | Seperates the first column from a matrix and returns it as a Vector along with the remaining Matrix
--   SeparateCol is safe since the given matrix has width >= 1
separateCol :: ( n ~ (n+1-1) ) => Matrix f m (n+1) -> (Vector f m, Matrix f m n)
separateCol = head . separateCols

-- | For each column vector in a matrix, returns the vector and the remaining matrix
separateCols :: Matrix f m n -> [(Vector f m, Matrix f m (n-1))]
separateCols m = map (\(v, m') -> (v, M (V m'))) $ separate' [] (matToList m)
    where separate' :: [Vector f m] -> [Vector f m] -> [(Vector f m, [Vector f m])]
          separate' _   []     = []
          separate' acc (v:[]) = [(v,     acc)]
          separate' acc (v:vs) =  (v, acc++vs) : separate' (acc++[v]) vs


-- | Seperates the first row from a matrix and returns it as a Vector along with the remaining Matrix
--   SeparateRow is safe since the given matrix has width >= 1
separateRow :: ( m ~ (m+1-1) ) => Matrix f (m+1) n -> (Vector f n, Matrix f m n)
separateRow = head . separateRows

-- | For each row vector in a matrix, returns the vector and the remaining matrix
separateRows :: Matrix f m n -> [(Vector f n, Matrix f (m-1) n)]
separateRows = map (\(v, m) -> (v, transpose m)) . separateCols . transpose

-- | Checks if the vectors in a matrix are linearly independant
linIndep :: (Eq f, Field f) => Matrix f m n -> Bool
linIndep = not . any (\(v, m) -> m `span` v ) . separateCols 


-- 2.27
-- Definition: basis
-- A basis of V is a list of vectors in V that is linearly independent and spans V

-- | Checks if the vectors in a matrix forms a basis of their vectorSpace
basis :: (KnownNat m, KnownNat (n-1), Eq f, Field f) => Matrix f m n -> Bool
basis m = spansSpace m && linIndep m

