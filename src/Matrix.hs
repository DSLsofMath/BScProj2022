{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}

module Matrix where

import GHC.TypeLits
import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), (^), sum, product, recip, fromRational)

import Data.List (sortOn, groupBy)
import Data.Function 

import Algebra
import qualified ListVector as L
import qualified SparseCSR as CSR

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
    {-# MINIMAL set , values #-}

    -- | Indexes into a matrix and gets a value
    get :: AddGroup f => mat f m n -> (Fin m, Fin n) -> f
    get m i = case lookup i (values m) of Just a  -> a
                                          Nothing -> zero

    -- | Sets the value at a given position
    set :: mat f m n -> (Fin m, Fin n) -> f -> mat f m n 

    -- | Returns all nonzero values with corresponding index 
    values :: mat f m n -> [((Fin m, Fin n), f)]

    -- | Builds a matrix from a list of positions and values
    --   Initializes with the zero matrix
    tabulate :: AddGroup (mat f m n) => [((Fin m, Fin n), f)] -> mat f m n
    tabulate = foldl (\m (i, a) -> set m i a) zero 

    
-- | Returns a list of elements and positions corresponding to a row
getRow :: Matrix mat => mat f m n -> Fin m -> [(Fin n, f)]
getRow m i = [ (j, a) | ((i', j), a) <- values m, i' == i]

-- | Returns a list of elements and positions corresponding to a column
getCol :: Matrix mat => mat f m n -> Fin n -> [(Fin m, f)]
getCol m j = [ (i, a) | ((i, j'), a) <- values m, j' == j]

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


----------------------------------------------------------------------------------------
-- Instances of Matrix

instance Matrix L.Matrix where 
    set (L.M (L.V vs)) (Fin i, Fin j) a = L.M . L.V $ as ++ L.V (as' ++ a:bs') : bs
        where (as, L.V v:bs) = splitAt (j-1) vs
              (as', _:bs') = splitAt (i-1) v

    -- TODO: Not sure how we should handle zeros since we do not want a Eq constraint in the class
    values (L.M (L.V vs)) = [ ((Fin i, Fin j), a) | (j, L.V v) <- zip [1..] vs
                                                  , (i, a)     <- zip [1..] v ]

instance Matrix CSR.CSR where

    -- TODO: There must be a better way to implement this
    set (CSR.CSR elems col row) (Fin i, Fin j) a = 
            let (col', elems') = unzip (as' ++ row' ++ bs') 
            in CSR.CSR elems' col' (as ++ n:map (+(length row' - length roww)) (n':bs))
        where (as, n:n':bs) = splitAt (i - 1) row
              (as', (roww,bs')) = splitAt (n'-n) <$> splitAt n (zip col elems)
              row' = insert (j - 1) a roww
              insert i x ys = sortOn fst $ (i,x) : filter ((i /=) . fst) ys

    tabulate xs = CSR.CSR elems col row
        where sorted = sortOn fst xs
              grouped = groupBy ((==) `on` fst . fst) sorted
              (col, elems) = unzip $ map (\((_,Fin i),a) -> (i - 1, a)) sorted
              row = scanl (+) 0 (map length grouped)

    values (CSR.CSR elems col row) = concat $ merge (zip col elems) perRow
        where perRow = zip [1..] $ zipWith (-) (tail row) row 
              merge _ [] = []
              merge xs ((i,n):ys) = let (cur, next) = splitAt n xs in
                    [ ((Fin (j+1),Fin i), a) | (j,a) <- cur ] : merge next ys


--------------------------------------------------------------------------
-- examples

pjMat :: (KnownNat n, Matrix mat, AddGroup (mat f n n), Ring f) => mat f n n
pjMat = tabulate $ concat [ f i | i <- [minBound .. maxBound] ]
    where f i | i == minBound = tail $ stencil i
              | i == maxBound = init $ stencil i
              | otherwise     =        stencil i
          stencil i = [ ((pred i, i), one), ((i,i), neg (one+one)), ((succ i, i), one) ]


