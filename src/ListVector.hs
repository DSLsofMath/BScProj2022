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

type Index = Int 

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
-- and allows it to infer the type of e.g. m44 (**) idm
instance (KnownNat m, Ring f, f ~ f', n ~ n') => Composable (Matrix f m n) (Vector f' n') (Vector f m) where
    (**) = (££)

instance (KnownNat a, Ring f, f ~ f', b ~ b', mat ~ Matrix) => Composable (Matrix f a b) (mat f' b' c) (Matrix f a c) where
    (**) = (£££)

-- | Square matrices form a multiplicative group
instance (KnownNat n, Ring f) => Mul (Matrix f n n) where
    (*) = (£££)
    one = idm

-- | Row first alternative of toMat
--   toMatT [[1,2,3],
--           [4,5,6]]
toMatT :: ToMat m n x => x -> Matrix (Under' x) n m
toMatT = transpose . toMat

-- | Converts objects to and from Matrices.
--   PROPOSAL: Should we get rid of this class and simply define functions instead?
class ToMat m n x where
    type Under' x 
    toMat   :: x -> Matrix (Under' x) m n
    fromMat :: Matrix (Under' x) m n -> x

instance (M.Matrix mat, m ~ m', n ~ n', AddGroup (mat f m n), AddGroup (Matrix f m n)) => ToMat m' n' (mat f m n) where
    type Under' (mat f m n) = f
    toMat = M.changeRep
    fromMat = M.changeRep

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

appendV :: Matrix f m n -> Vector f m -> Matrix f m (n + 1)
appendV m (V vs) = pack $ unpack m ++ [vs]

-- | Converts a Matrix to a list of Vectors 
matToList :: Matrix f m n -> [Vector f m]
matToList = coerce



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

