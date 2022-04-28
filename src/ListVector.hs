{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}

module ListVector where

import GHC.TypeLits hiding (type (^))
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), (^), recip, sum, product, (**), span)
import Data.Coerce

import qualified Data.List as L
import Algebra
import qualified Matrix as M
import Matrix (Fin(..))
import qualified Gauss as G
import Data.List (sortOn, groupBy)
import Data.Function 

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
newtype Vector f (n :: Nat) = V [f] deriving (Eq)

mapV :: (a->b) -> Vector a n -> Vector b n
mapV f (V l) = V (map f l)

instance Show f => Show (Vector f n) where
    show v = show . M $ V [v]

-- | Allows us to patternmatch the vector elements
--   without explicitly unwrapping the underlying List
-- pattern V a = Vector (L a)

type VecR = Vector Double

-- | Nice type notation for vectors
--   v = vec [1,2,3,4] :: R^4
type f ^ n = Vector f n

-- | The length of the vector given by its type
vecLen :: KnownNat n => Vector f n -> Int
vecLen v = fromInteger $ natVal v

-- | Vector that is all zero
zeroVec :: (KnownNat n, AddGroup f) => Vector f n
zeroVec = let v = V $ replicate (vecLen v) zero in v

-- | converts a list to a vector with a typed size
vec :: (KnownNat n, AddGroup f) => [f] -> Vector f n
vec ss = V (ss ++ repeat zero) + zero

-- | e i is the i:th basis vector
e :: (KnownNat n, Ring f) => Int -> Vector f n
e i = V (replicate (i-1) zero ++ (one : repeat zero)) + zero

zipWithV :: (a -> b -> c) -> Vector a n -> Vector b n -> Vector c n
zipWithV op (V as) (V bs) = V $ zipWith op as bs

-- Vector is a vector space over a field
instance (KnownNat n, AddGroup f) => AddGroup (Vector f n) where
    (+) = zipWithV (+)
    (-) = zipWithV (-) 
    zero = zeroVec

instance (KnownNat n, Ring f) => VectorSpace (Vector f n) where
    type Under (Vector f n) = f
    s £ v = mapV (s*) v

-- | Vector is further a finite Vector Field
instance (KnownNat n, Ring f) => Finite (Vector f n) where
    type Dim (Vector f n) = n
    basis' _ = let M (V bs) = idm in L bs


-- | Dot product of two vectors 
v1, v2 :: Vector Double 3
v1 = V [2,7,1]::VecR 3
v2 = V [8,2,8]::VecR 3 
-- test dot product
testdot = dot v1 v2 == 2*8 + 7*2 + 8*1

dot :: Ring f => Vector f n -> Vector f n -> f
V v1 `dot` V v2 = sum $ zipWith (*) v1 v2

-- | Cross product of two vectors of size 3
-- test cross product
testross = cross v1 v2 == V[7*8-1*2, 1*8-2*8, 2*2-7*8]

cross :: Ring f => Vector f 3 -> Vector f 3 -> Vector f 3
V [a1,a2,a3] `cross` V [b1,b2,b3] = V [a2*b3-a3*b2,
                                       a3*b1-a1*b3,
                                       a1*b2-a2*b1]


-- | Takes a Vector and a List of vectors and returns the linear combination
linComb :: VectorSpace v => Vector v n -> Vector (Under v) n -> v
linComb (V vs) (V fs) = sum $ zipWith (£) fs vs

-- | Takes a Vector and a List of vectors and returns the linear combination
--   For Example eval [a b c] [^0 ^1 ^2] returns the polinomial \x -> ax + bx + cx^2
eval :: (VectorSpace v, Under v ~ f) => Vector f n -> List v n -> v
eval (V fs) (L vs) = sum $ zipWith (£) fs vs


-------------------------------------------
-- Matrix definitions

-- | Matrix as a vector of vectors
--   note that the outmost vector is not matimatically a vector 
newtype Matrix f (m :: Nat) (n :: Nat) = M (Vector (Vector f m) n)  deriving (Eq)

type MatR = Matrix Double

instance Show f => Show (Matrix f m n) where
    show = showMat

-- Show function for matrices 
showMat :: (Show f) => Matrix f m n -> String
showMat = ("\n"++) . unlines . map formatRow . L.transpose . map padCol . unpack 
    where
        getLongest = maximum . map length
        padString n s = replicate (n-length s) ' ' ++ s ++ " "
        padCol l = let s = map show l in map (padString (getLongest s)) s
        formatRow s = "| " ++ unwords s ++ "|"


-- | Identity matrix
idm :: (KnownNat n, Ring f) => Matrix f n n
idm = let v = V [ e i | i <- [1 .. vecLen v] ] in M v

-- | Matrix vector multiplication
(££) :: (KnownNat m, Ring f) => Matrix f m n -> Vector f n -> Vector f m
M vs ££ v = linComb vs v

-- | Matrix matrix multiplication
(£££) :: (KnownNat a, Ring f) => Matrix f a b -> Matrix f b c -> Matrix f a c
a £££ b = (a££) `onCols` b


-- | Applies a function on each column vector
--   Represents composition f M 
onCols :: (Vector f b -> Vector f a) -> Matrix f b c -> Matrix f a c
onCols f (M v) = M $ mapV f v

-- | Returns the matrix representation of a linear function
--   We should have that
--   f `onCols` m == funToMat f £££ m
--   UNSAFE: the function should be linear.
funToMat :: (KnownNat n, Ring f) => (Vector f n -> Vector f m) -> Matrix f m n
funToMat f = f `onCols` idm

-- | Puts all element in the diagonal in a vector
getDiagonal :: Matrix f n n -> Vector f n
getDiagonal m = V $ zipWith (!!) (unpack m) [0..]
  
-- Matrices also forms a vector space
instance (KnownNat m, KnownNat n, AddGroup f) => AddGroup (Matrix f m n) where
    M as + M bs = M $ zipWithV (+) as bs
    M as - M bs = M $ zipWithV (-) as bs
    zero = let v = V $ replicate (vecLen v) zero in M v

instance (KnownNat m, KnownNat n, Ring f) => VectorSpace (Matrix f m n) where
    type Under (Matrix f m n) = f
    s £ m = (s£) `onCols` m


instance M.Matrix Matrix where 
    set (M (V vs)) (Fin i, Fin j) a = M . V $ as ++ V (as' ++ a:bs') : bs
        where (as, V v:bs) = splitAt (j-1) vs
              (as', _:bs') = splitAt (i-1) v

   -- extend listM = tabulate' . merge' (M.values listM)
   --     where merge' as bs = map (\x -> last x : []) . groupBy ((==) `on` fst) $ sortOn fst (as ++ bs)
   --           tabulate' xs = M $ V [V [ b| ((_, _), b) <- a] | a <- xs] 


    -- TODO: Not sure how we should handle zeros since we do not want a Eq constraint in the class
    values (M (V vs)) = [ ((Fin i, Fin j), a) | (j, V v) <- zip [1..] vs
                                              , (i, a)   <- zip [1..] v ]

-- Composition on matrices
-- Note that we write n ~ n' instead of writing n on both places. 
-- This tells GHC that this is the only match for Matrix*Vector or Matrix*Matrix,
-- and allows it to infer the type of e.g. m44 `comp` idm
instance (KnownNat m, Ring f, f ~ f', n ~ n') => Composable (Matrix f m n) (Vector f' n') (Vector f m) where
    comp = (££)

instance (KnownNat a, Ring f, f ~ f', b ~ b' ) => Composable (Matrix f a b) (Matrix f' b' c) (Matrix f a c) where
    comp = (£££)

-- | Square matrices form a multiplicative group
instance (KnownNat n, Ring f) => Mul (Matrix f n n) where
    (*) = (£££)
    one = idm

-- | Converts objects to and from Matrices.
--   PROPOSAL: Should we get rid of this class and simply define functions instead?
class ToMat m n x where
    type Under' x 
    toMat   :: x -> Matrix (Under' x) m n
    fromMat :: Matrix (Under' x) m n -> x

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

instance (m' ~ m, n' ~ n) => ToMat m' n' (Matrix f m n) where
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


-- | Transposes the matrix
transpose :: Matrix f m n -> Matrix f n m
transpose = pack . L.transpose . unpack

-- | Gets the value at (column, row)
get :: Matrix f m n -> (Index, Index) -> f
get m (x,y) = unpack m !! x !! y

-- | Converts a Matrix to a list of lists 
unpack :: Matrix f m n -> [[f]]
unpack = coerce

-- | Converts a list of lists to a Matrix
--   UNSAFE: should only be used when the correct dimensions can be guaranteed.
--   For a safer alternative se toMat
pack :: [[f]] -> Matrix f m n
pack = coerce

-- | Appends the second matrix to the right of the first, analogous to (++)
--   useful for Ex. Guassian elimination
append :: Matrix f m n1 -> Matrix f m n2 -> Matrix f m (n1+n2)
append m1 m2 = pack $ unpack m1 ++ unpack m2

-- | Appends the second matrix below the first
append' :: Matrix f m1 n -> Matrix f m2 n -> Matrix f (m1+m2) n
append' m1 m2 = pack $ zipWith (++) (unpack m1) (unpack m2)

-- | Converts a Matrix to a list of Vectors 
matToList :: Matrix f m n -> [Vector f m]
matToList = coerce


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

-- | Elementary row functions.
--   It might be better to define the as (Vector f n -> Vector f n) e.g. a linear map. 
--   If we want to apply it to a matrix we can then use onCols. 
--   Doing so will also remove the need for transpose.
swap :: Index -> Index -> Matrix f m n -> Matrix f m n
swap i j = onCols $ \(V v) -> 
    let (i', j') = (min i j - 1, max i j - 1)
        (m1,x:xs) = splitAt i' v
        (xs',y:m2) = splitAt (j' - i'-1) xs
    in V $ m1++y:xs'++x:m2

mul :: Mul f => Index -> f -> Matrix f m n -> Matrix f m n
mul i s = onCols $ \(V v) -> 
    let (m1,x:m2) = splitAt (i-1) v in V $ m1++(s*x): m2

muladd :: Ring f => Index -> Index -> f -> Matrix f m n -> Matrix f m n
muladd i j s = onCols $ \(V v) -> 
    let (i', j') = (i - 1, j - 1)
        (_,x:_) = splitAt i' v
        (m1,y:m2) = splitAt j' v
        y' = (s*x) + y
    in V $ m1++y':m2

instance (Ring f, f ~ f') => Composable (ElimOp f') (Matrix f m n) (Matrix f m n) where
    (Swap i j) `comp` m = swap i j m
    (Mul i s) `comp` m = mul i s m
    (MulAdd i j s) `comp` m = muladd i j s m


instance (Ring f, n ~ n', f ~ f') => Composable (G.ElimOp n' f') (Matrix f m n) (Matrix f m n) where
    (G.Swap (Fin i) (Fin j)) `comp` m = swap i j m
    (G.Mul (Fin i) s) `comp` m = mul i s m
    (G.MulAdd (Fin i) (Fin j) s) `comp` m = muladd i j s m


-- | Representation of an elementary row operation as a matrix 
elimOpToMat :: (KnownNat n, Ring f) => ElimOp f -> Matrix f n n
elimOpToMat e = elimOpToFunc e idm

foldElemOps :: (KnownNat n, Ring f) => [ElimOp f] -> Matrix f n n
foldElemOps = product . map elimOpToMat . reverse

-- | Representation of an elementary row operation as a function 
elimOpToFunc :: Ring f => ElimOp f -> (Matrix f m n -> Matrix f m n)
elimOpToFunc (Swap   i j  ) = swap   i j
elimOpToFunc (Mul    i   s) = mul    i   s
elimOpToFunc (MulAdd i j s) = muladd i j s
                         
-- | Reduces a trace of elimOps to a single function
--   TODO: We should add a rewrite rule such that fold elemOpToFunc only packs and unpacks once
foldElemOpsFunc :: Ring f => [ElimOp f] -> (Matrix f m n -> Matrix f m n)
foldElemOpsFunc = foldr (.) id . map elimOpToFunc . reverse


-- | Applies a function on a unpacked and transposed matrix before transposing it back
--   UNSAFE: The function can not change the dimension of the matrix 
onUnpackedTrans :: ([[f]] -> [[f]]) -> Matrix f m n -> Matrix f m n
onUnpackedTrans f = pack . L.transpose . f . L.transpose . unpack

-- | Transform a matrix to upper triangular form
utf :: (Eq f, Field f) => Matrix f m n -> Matrix f m n
utf = onUnpackedTrans (sort . f) 
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


-- | Solve systems of equations
--   Check in Demo.hs how to use
solve :: Ring f => [[f]] -> [f]
solve m = foldr next [last (last m)] (init m)
    where
        next row found = let
            subpart = init $ drop (length m - length found) row
            solved = last row - sum (zipWith (*) found subpart)
            in solved : found

-- | apply when solving systems of equations
--   each element in list represents variable values
--   TODO: handle case when there is no solution, return Maybe 
solvesys :: (Eq f, Field f) => Matrix f m n -> Vector f (n -1)
solvesys = V . solve . unpack . transpose . utf


-- | Separates the first column from a matrix and returns it as a Vector along with the remaining Matrix
--   SeparateCol is safe since the given matrix has width >= 1
separateCol :: ( n ~ (n+1-1) ) => Matrix f m (n+1) -> (Vector f m, Matrix f m n)
separateCol = head . separateCols

-- | For each column vector in a matrix, returns the vector and the remaining matrix
separateCols :: Matrix f m n -> [(Vector f m, Matrix f m (n-1))]
separateCols = map (\(v, m) -> (v, M (V m)) ) . separate' [] . matToList 
    where separate' :: [Vector f m] -> [Vector f m] -> [(Vector f m, [Vector f m])]
          separate' _   []     = []
          separate' acc (v:[]) = [(v, acc    )]
          separate' acc (v:vs) =  (v, acc++vs) : separate' (acc++[v]) vs

-- | Separates the first row from a matrix and returns it as a Vector along with the remaining Matrix
--   SeparateRow is safe since the given matrix has width >= 1
separateRow :: ( m ~ (m+1-1) ) => Matrix f (m+1) n -> (Vector f n, Matrix f m n)
separateRow = head . separateRows

-- | For each row vector in a matrix, returns the vector and the remaining matrix
separateRows :: Matrix f m n -> [(Vector f n, Matrix f (m-1) n)]
separateRows = map (\(v, m) -> (v, transpose m)) . separateCols . transpose


scaleM :: Ring f => f -> Matrix f m n -> Matrix f m n
scaleM c (M vs) = M (mapV (scaleV c) vs)

scaleV :: Ring f => f -> Vector f m -> Vector f m
scaleV c = mapV (c*)

