{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleInstances #-}

module Matrix where

import GHC.TypeLits
import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), (^), sum, product, recip, fromRational)

import Data.List (sortOn, groupBy)
import Data.Function 

import Algebra

-- | Bounded integer
--   A Fin n should always be in the range of (1,n) 
newtype Fin (n :: Nat) = Fin Int deriving (Eq, Ord) 

-- | Safe constructor of Fin, checks that the value is in range
fin :: KnownNat n => Int -> Fin n
fin i = finite
    where finite | 1 <= i && i <= fromInteger (natVal finite) = Fin i
                 | otherwise = error $ "Index is out of bounds, got: " ++ 
                   show i ++ " in constraint 0<=" ++ show i ++ "<=" ++ show (natVal finite)

instance KnownNat n => Show (Fin n) where
    show n = show (finToInt n) ++ " of " ++ show (natVal n)

-- Enum and Bounded allows for list sequencing 
-- i.e. [Fin 2 .. maxBound]
instance KnownNat n => Enum (Fin n) where
    toEnum = fin
    fromEnum = finToInt
    
instance KnownNat n => Bounded (Fin n) where
    minBound = Fin 1
    maxBound = let f = Fin (fromInteger (natVal f)) in f

-- | Returns the corresponding Int
finToInt :: Fin n -> Int
finToInt (Fin n) = n

-- | Increases the max bound of the index
raiseBound :: n <= m => Fin n -> Fin m
raiseBound (Fin n) = Fin n


-- | Merges two indexed lists, adding their values if they have the same index
addL :: (Ord i, AddGroup a) => [(i, a)] -> [(i, a)] -> [(i, a)] 
addL xs ys = map (foldl1 add) . groupBy ((==) `on` fst) $ sortOn fst (xs ++ ys)
    where (i, a) `add` (_,b) = (i, a + b) 

scaleL :: (Ord i, Mul a) => a -> [(i, a)] -> [(i, a)] 
scaleL s l = [ (i, s*a) | (i,a) <- l ]

instance (Ord i, Ring a) => AddGroup [(i, a)] where
    (+) = addL
    neg = scaleL (neg one)
    zero = []

instance (Ord i, Ring a) => VectorSpace [(i,a)] where
    type Under [(i,a)] = a
    (£) = scaleL

-- | Generic class for a matrix with focus on sparse representations
--   The class is roughly based on the paper "APLicative Programming with Naperian Functors"
--   by Jeremy Gibbons. 
--   However we make some changes to better suite sparse matrices,
--   most notably we only allow implementations for two dimensions.
--   
--
--   When implementing an instance, the only functions that are needed are set and values.
--   Although, performance may greatly increase with a custom tabulate implementation.
--   
--   TODO/Notes: To use tabulate we currently need a AddGroup constraint
--   since we need a zero representation. This is bit awkward since it carries
--   over to seemingly unrelated functions, like transpose and changeRep.
--   It would be nice if the zero state was part of the Matrix class directly.
--   But it is unclear how one would do that in a reasonable way. 
--
--   Another related annoyance is that AddGroup (mat f m n)
--   does not imply AddGroup f, or vice verse. 
--   It would be nice, and it seams plausible, to have at least 
--   VectorSpace (mat f m n) => Ring f
--   As it is now we have to deal with some rather hefty type constraint, 
--   see for example identity or pjMat
--   
class Matrix (mat :: * -> Nat -> Nat -> *) where
    {-# MINIMAL (set | extend) , values #-}

    -- | Sets the value at a given position
    set :: mat f m n -> (Fin m, Fin n) -> f -> mat f m n 
    set m i f = extend m [(i,f)]

    -- | Returns all nonzero values with corresponding index 
    values :: mat f m n -> [((Fin m, Fin n), f)]

    extend :: mat f m n -> [((Fin m, Fin n), f)] -> mat f m n
    extend = foldl (\m (i, a) -> set m i a) 

-- | Indexes into a matrix and gets a value
get :: (Matrix mat, AddGroup f) => mat f m n -> (Fin m, Fin n) -> f
get m i = case lookup i (values m) of Just a  -> a
                                      Nothing -> zero
    
-- | Returns a list of elements and positions corresponding to a row
getRow :: Matrix mat => mat f m n -> Fin m -> [(Fin n, f)]
getRow m i = [ (j, a) | ((i', j), a) <- values m, i' == i]

-- | Returns a list of elements and positions corresponding to a column
getCol :: Matrix mat => mat f m n -> Fin n -> [(Fin m, f)]
getCol m j = [ (i, a) | ((i, j'), a) <- values m, j' == j]

-- | Returns a list of elements and positions corresponding to the diagonal
getDiagonal :: Matrix mat => mat f n n -> [(Fin n, f)]
getDiagonal m = [ (i, a) | ((i, j), a) <- values m, i == j ]

setRow :: (Matrix mat, AddGroup (mat f m n)) => mat f m n -> Fin m -> [(Fin n, f)] -> mat f m n
setRow m i r = extend m [ ((i, j), a) | (j, a) <- r ]
-- setRow m i r = tabulate $ new ++ old
--     where old = filter ((i /=) . fst . fst) $ values m
--           new = [ ((i, j), a) | (j, a) <- r ]

-- | Builds a matrix from a list of positions and values
--   Initializes with the zero matrix
tabulate :: (Matrix mat, AddGroup (mat f m n)) => [((Fin m, Fin n), f)] -> mat f m n
tabulate = extend zero 

-- | Appends two matrices, analogous to (++)
append :: (KnownNat n1, Matrix mat, AddGroup (mat f m (n1 + n2)) ) => mat f m n1 -> mat f m n2 -> mat f m (n1 + n2)
append m1 m2 = tabulate $ m1' ++ m2'
    where m1' = [ ((i, Fin j),       a) | ((i, Fin j), a) <- values m1 ]
          m2' = [ ((i, Fin (n + j)), a) | ((i, Fin j), a) <- values m2 ]
          n = fromInteger (natVal m1)

transpose :: (Matrix mat, AddGroup (mat f n m)) => mat f m n -> mat f n m
transpose = tabulate . map (\((i,j),a) -> ((j,i),a)) . values

-- | Changes the underlying matrix type
--   useful when converting from one representation to another
--   changeRep (identity :: CSR R 5 5) :: Matrix R 5 5
changeRep :: (Matrix mat1, Matrix mat2, AddGroup (mat2 f m n)) => mat1 f m n -> mat2 f m n
changeRep = tabulate . values

-- | A general implementation of the identity matrix
identity :: (KnownNat n, Matrix mat, AddGroup (mat f n n), Mul f) => mat f n n
identity = tabulate [ ((i,i), one) | i <- [minBound .. maxBound]]


-- | Like values but also removes zeros
purgeToList :: (Matrix mat, Eq f, AddGroup f) => mat f m n -> [((Fin m, Fin n), f)]
purgeToList = (filter ((zero /=) . snd)) . values

purge :: (Matrix mat, Eq f, AddGroup f, AddGroup (mat f m n)) => mat f m n -> mat f m n
purge = toSparse

-- | Like changeRep but also removes zeros
toSparse :: (Matrix mat1, Matrix mat2, Eq f, AddGroup f, AddGroup (mat2 f m n)) => mat1 f m n -> mat2 f m n
toSparse = tabulate . purgeToList

----------------------------------------------------------------------------------------
-- Instances of Matrix

 

--------------------------------------------------------------------------
-- examples

pjMat :: (KnownNat n, Matrix mat, AddGroup (mat f n n), Ring f) => mat f n n
pjMat = tabulate $ concat [ f i | i <- [minBound .. maxBound] ]
    where f i | i == minBound = tail $ stencil i
              | i == maxBound = init $ stencil i
              | otherwise     =        stencil i
          stencil i = [ ((pred i, i), one), ((i,i), neg (one+one)), ((succ i, i), one) ]


