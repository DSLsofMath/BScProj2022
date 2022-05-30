{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}

{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE DeriveFoldable #-}

module ListVector where

import GHC.TypeLits hiding (type (^))
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), (^), recip, sum, product, (**), span)
import Data.Coerce
import Data.Proxy

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
-- Vector n f = Matrix f n 1
-- to make the code more general

-- | Dependent typed vector
newtype Vector (n :: Nat) f = V [f] 
    deriving (Eq, Functor, Foldable)

instance Show f => Show (Vector n f) where
    show v = show . M $ V [v]

-- | Nice type notation for vectors
--   v = vec [1,2,3,4] :: R^4
type f ^ n = Vector n f

-- | The length of the vector given by its type
vecLen :: forall n f. KnownNat n => Vector n f -> Int
vecLen _ = fromInteger $ natVal (Proxy @n)

-- | Vector containing a single term
pureVec :: KnownNat n => f -> Vector n f
pureVec f = let v = V $ replicate (vecLen v) f in v

-- | Vector that is all zero
zeroVec :: (KnownNat n, AddGroup f) => Vector n f
zeroVec = pureVec zero

-- | converts a list to a vector with a typed size
vec :: KnownNat n => [f] -> Vector n f
vec ss = if (vecLen v == length ss) then v
                                    else error errorMsg
    where v = V ss
          errorMsg = "Vector is of dimension " ++ show (vecLen v) ++ 
                     " but was given a list of length " ++ show (length ss)

-- | e i is the i:th basis vector
e :: (KnownNat n, Ring f) => Int -> Vector n f
e i = V (replicate (i-1) zero ++ (one : repeat zero)) + zero

-- Applies a funcion pairwise on two vectors of same length
zipWithV :: (a -> b -> c) -> Vector n a -> Vector n b -> Vector n c
zipWithV op (V as) (V bs) = V $ zipWith op as bs

instance KnownNat n => Applicative (Vector n) where
    pure = pureVec
    (<*>) = zipWithV ($)

-- Vector is a vector space over a field
instance (KnownNat n, AddGroup f) => AddGroup (Vector n f) where
    (+) = zipWithV (+)
    (-) = zipWithV (-) 
    zero = zeroVec

instance (KnownNat n, Ring f) => VectorSpace (Vector n f) where
    type Under (Vector n f) = f
    s £ v = fmap (s*) v

-- | Vector is further a finite Vector Field
instance (KnownNat n, Ring f) => Finite (Vector n f) where
    type Dim (Vector n f) = n
    basis' _ = let M (V bs) = idm in L bs


-- | Dot product of two vectors 
v1, v2 :: Vector 3 Double
v1 = V [2,7,1]
v2 = V [8,2,8]
-- test dot product
testdot = dot v1 v2 == 2*8 + 7*2 + 8*1

dot :: Ring f => Vector n f -> Vector n f -> f
v1 `dot` v2 = sum $ zipWithV (*) v1 v2

-- | Cross product of two vectors of size 3
-- test cross product
testross = cross v1 v2 == V[7*8-1*2, 1*8-2*8, 2*2-7*8]

cross :: Ring f => Vector 3 f -> Vector 3 f -> Vector 3 f
V [a1,a2,a3] `cross` V [b1,b2,b3] = V [a2*b3-a3*b2,
                                       a3*b1-a1*b3,
                                       a1*b2-a2*b1]


-- | Takes a Vector and a List of vectors and returns the linear combination
linComb :: VectorSpace v => Vector n v -> Vector n (Under v) -> v
linComb vs fs = sum $ zipWithV (£) fs vs

-- | Takes a Vector and a List of vectors and returns the linear combination
--   For Example eval [a b c] [^0 ^1 ^2] returns the polinomial \x -> ax + bx + cx^2
eval :: (VectorSpace v, Under v ~ f) => Vector n f -> List v n -> v
eval (V fs) (L vs) = sum $ zipWith (£) fs vs


-------------------------------------------
-- Matrix definitions

-- | Matrix as a vector of vectors
--   note that the outmost vector is not matimatically a vector 
newtype Matrix (m :: Nat) (n :: Nat) f = M (Vector n (Vector m f)) 
    deriving (Functor)

instance Show f => Show (Matrix m n f) where
    show = showMat

-- Show function for matrices 
showMat :: (Show f) => Matrix m n f -> String
showMat = ("\n"++) . unlines . map formatRow . L.transpose . map padCol . unpack 
    where
        getLongest = maximum . map length
        padString n s = replicate (n-length s) ' ' ++ s ++ " "
        padCol l = let s = map show l in map (padString (getLongest s)) s
        formatRow s = "| " ++ unwords s ++ "|"

zipWithM :: (a -> b -> c) -> Matrix m n a -> Matrix m n b -> Matrix m n c
zipWithM op (M as) (M bs) = M $ zipWithV (zipWithV op) as bs

-- | Identity matrix
idm :: (KnownNat n, Ring f) => Matrix n n f
idm = let v = V [ e i | i <- [1 .. vecLen v] ] in M v

-- | Matrix vector multiplication
(££) :: (KnownNat m, Ring f) => Matrix m n f -> Vector n f -> Vector m f
M vs ££ v = linComb vs v

-- | Matrix matrix multiplication
(£££) :: (KnownNat a, Ring f) => Matrix a b f -> Matrix b c f -> Matrix a c f
a £££ b = (a ££) `onCols` b


-- | Applies a function on each column vector
--   Represents composition f M 
onCols :: (Vector b f -> Vector a f) -> Matrix b c f -> Matrix a c f
onCols f (M v) = M $ fmap f v

-- | Returns the matrix representation of a linear function
--   We should have that
--   f `onCols` m == funToMat f £££ m
--   UNSAFE: the function should be linear.
funToMat :: (KnownNat n, Ring f) => (Vector n f -> Vector m f) -> Matrix m n f
funToMat f = f `onCols` idm

-- | Puts all element in the diagonal in a vector
getDiagonal :: Matrix n n f -> Vector n f
getDiagonal m = V $ zipWith (!!) (unpack m) [0..]

zeroMat :: (KnownNats m n, AddGroup f) => Matrix m n f
zeroMat = M (pureVec zeroVec)

set :: Matrix m n f -> ((Fin m, Fin n), f) -> Matrix m n f
set (M (V vs)) ((Fin i, Fin j), a) = M . V $ as ++ V (as' ++ a:bs') : bs
    where (as, V v:bs) = splitAt (j-1) vs
          (as', _:bs') = splitAt (i-1) v

instance M.Matrix Matrix where 
    
    tabulate = foldl set zeroMat 

    -- extend listM = tabulate' . merge' (M.values listM)
    --    where merge' as bs = map (\x -> last x : []) . groupBy ((==) `on` fst) $ sortOn fst (as ++ bs)
    --          tabulate' xs = M $ V [V [ b| ((_, _), b) <- a] | a <- xs] 


    -- TODO: Not sure how we should handle zeros since we do not want a Eq constraint in the class
    values (M (V vs)) = [ ((Fin i, Fin j), a) | (j, V v) <- zip [1..] vs
                                              , (i, a)   <- zip [1..] v ]

    mulMat = (£££)

    addMat = zipWithM (+)

-- Composition on matrices
-- Note that we write n ~ n' instead of writing n on both places. 
-- This tells GHC that this is the only match for Matrix*Vector or Matrix*Matrix,
-- and allows it to infer the type of e.g. m44 (**) idm
instance (KnownNat m, Ring f, f ~ f', n ~ n') => Composable (Matrix m n f) (Vector n' f') (Vector m f) where
    (**) = (££)


-- | Row first alternative of toMat
--   toMatT [[1,2,3],
--           [4,5,6]]
toMatT :: (KnownNats m n, ToMat n m x) => x -> Matrix m n (Under' x) 
toMatT = transpose . toMat

-- | Converts objects to and from Matrices.
--   PROPOSAL: Should we get rid of this class and simply define functions instead?
class ToMat m n x where
    type Under' x 
    toMat   :: x -> Matrix m n (Under' x) 
    fromMat :: Matrix m n (Under' x) -> x

instance (KnownNats m n, M.Matrix mat, m ~ m', n ~ n', AddGroup f) => ToMat m' n' (mat m n f) where
    type Under' (mat m n f) = f
    toMat = M.changeRep
    fromMat = M.changeRep

instance ToMat n 1 (Vector n f) where
    type Under' (Vector n f) = f
    toMat v = M (V [v])
    fromMat (M (V [v])) = v

instance ToMat 1 n (Vector n f) where
    type Under' (Vector n f) = f
    toMat   (V ss) = M . V $ map (\s -> V [s]) ss
    fromMat (M (V vs)) = V $ map (\(V (x:_)) -> x) vs

-- | Diagonal matrix
instance (KnownNat n, Field f) => ToMat n n (Vector n f) where
    type Under' (Vector n f) = f
    toMat (V ss) = M . V $ zipWith (\s i-> s £ e i) ss [1..]
    fromMat m = vec $ zipWith (!!) (fromMat m) [0..]

instance (KnownNats m n) => ToMat m n [Vector m f] where
    type Under' [Vector m f] = f
    toMat vs = M . vec $ vs
    fromMat = matToList

instance (KnownNats m n) => ToMat m n [[f]] where
    type Under' [[f]] = f
    toMat = M . vec . map vec
    fromMat = unpack


-- | Transposes the matrix
transpose :: Matrix m n f -> Matrix n m f
transpose = pack . L.transpose . unpack

-- | Gets the value at (column, row)
get :: Matrix m n f -> (Int, Int) -> f
get m (x,y) = unpack m !! x !! y

-- | Converts a Matrix to a list of lists 
unpack :: Matrix m n f -> [[f]]
unpack = coerce

-- | Converts a list of lists to a Matrix
--   UNSAFE: should only be used when the correct dimensions can be guaranteed.
--   For a safer alternative se toMat
pack :: [[f]] -> Matrix m n f
pack = coerce

-- | Appends the second matrix to the right of the first, analogous to (++)
--   useful for Ex. Guassian elimination
append :: Matrix m n1 f -> Matrix m n2 f -> Matrix m (n1+n2) f
append m1 m2 = pack $ unpack m1 ++ unpack m2

-- | Appends the second matrix below the first
append' :: Matrix m1 n f -> Matrix m2 n f-> Matrix (m1+m2) n f
append' m1 m2 = pack $ zipWith (++) (unpack m1) (unpack m2)

-- | Appends a vector to the right of a matrix
appendV :: Matrix m n f -> Vector m f -> Matrix m (n + 1) f
appendV m (V vs) = pack $ unpack m ++ [vs]

-- | Converts a Matrix to a list of Vectors 
matToList :: Matrix m n f -> [Vector m f]
matToList = coerce



-- | Separates the first column from a matrix and returns it as a Vector along with the remaining Matrix
--   SeparateCol is safe since the given matrix has width >= 1
separateCol :: ( n ~ (n+1-1) ) => Matrix m (n+1) f -> (Vector m f, Matrix m n f)
separateCol = head . separateCols

-- | For each column vector in a matrix, returns the vector and the remaining matrix
separateCols :: Matrix m n f -> [(Vector m f, Matrix m (n-1) f)]
separateCols = map (\(v, m) -> (v, M (V m)) ) . separate' [] . matToList 
    where separate' :: [Vector f m] -> [Vector f m] -> [(Vector f m, [Vector f m])]
          separate' _   []     = []
          separate' acc (v:[]) = [(v, acc    )]
          separate' acc (v:vs) =  (v, acc++vs) : separate' (acc++[v]) vs

-- | Separates the first row from a matrix and returns it as a Vector along with the remaining Matrix
--   SeparateRow is safe since the given matrix has width >= 1
separateRow :: ( m ~ (m+1-1) ) => Matrix (m+1) n f -> (Vector n f, Matrix m n f)
separateRow = head . separateRows

-- | For each row vector in a matrix, returns the vector and the remaining matrix
separateRows :: Matrix m n f -> [(Vector n f, Matrix (m-1) n f)]
separateRows = map (\(v, m) -> (v, transpose m)) . separateCols . transpose


