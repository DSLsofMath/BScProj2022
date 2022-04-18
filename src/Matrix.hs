{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}


import GHC.TypeLits
import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), (^), sum, product, recip, fromRational)

import Data.List (sortOn, groupBy)
import Data.Function 

import Algebra
import qualified ListVector as L
import qualified SparseCSR as CSR

newtype Fin (n :: Nat) = Fin Int deriving (Eq, Ord) 

fin :: KnownNat n => Int -> Fin n
fin i = finite
    where finite | 1 <= i && i <= fromInteger (natVal finite) = Fin i
                 | otherwise = error $ "Index is out of bounds, got: " ++ 
                   show i ++ " in constraint 0<=" ++ show i ++ "<=" ++ show (natVal finite)

instance KnownNat n => Show (Fin n) where
    show n = show (finToInt n) ++ " of " ++ show (natVal n)

instance KnownNat n => Enum (Fin n) where
    toEnum = fin
    fromEnum = finToInt
    
instance KnownNat n => Bounded (Fin n) where
    minBound = Fin 1
    maxBound = let f = Fin (fromInteger (natVal f)) in f

finToInt :: Fin n -> Int
finToInt (Fin n) = n

raiseBound :: n <= m => Fin n -> Fin m
raiseBound (Fin n) = Fin n


-- | Generic class for a matrix with focus on sparse representations
class Matrix (mat :: * -> Nat -> Nat -> *) where
    {-# MINIMAL set , values #-}
    -- type Test mat f m n :: Constraint

    get :: AddGroup f => mat f m n -> (Fin m, Fin n) -> f
    get m i = case lookup i (values m) of Just a  -> a
                                          Nothing -> zero

    set :: mat f m n -> (Fin m, Fin n) -> f -> mat f m n 

    -- | Returns all nonzero values in a matrix
    values :: mat f m n -> [((Fin m, Fin n), f)]

    -- | Builds a matrix from a list of positions and values
    --   Initializes with the zero matrix
    tabulate :: AddGroup (mat f m n) => [((Fin m, Fin n), f)] -> mat f m n
    tabulate = foldl (\m (i, a) -> set m i a) zero 

    transpose :: (AddGroup (mat f n m)) => mat f m n -> mat f n m
    transpose = tabulate . map (\((i,j),a) -> ((j,i),a)) . values

    
getRow :: Matrix mat => mat f m n -> Fin m -> [(Fin n, f)]
getRow m i = [ (j, a) | ((i', j), a) <- values m, i' == i]

getCol :: Matrix mat => mat f m n -> Fin n -> [(Fin m, f)]
getCol m j = [ (i, a) | ((i, j'), a) <- values m, j' == j]

append :: (KnownNat n1, Matrix mat, AddGroup (mat f m (n1 + n2)) ) => mat f m n1 -> mat f m n2 -> mat f m (n1 + n2)
append m1 m2 = tabulate $ m1' ++ m2'
    where m1' = [ ((i, Fin j),       a) | ((i, Fin j), a) <- values m1 ]
          m2' = [ ((i, Fin (n + j)), a) | ((i, Fin j), a) <- values m2 ]
          n = fromInteger (natVal m1)

changeRep :: (Matrix mat1, Matrix mat2, AddGroup (mat2 f m n)) => mat1 f m n -> mat2 f m n
changeRep = tabulate . values

identity :: (KnownNat n, Matrix mat, AddGroup (mat f n n), Mul f) => mat f n n
identity = idm
    where idm = tabulate [ ((fin i,fin i), one) | i <- [1 .. n]]
          n = fromInteger $ natVal idm

----------------------------------------------------------------------------------------
-- Instances of Matrix

instance Matrix L.Matrix where 
    set (L.M (L.V vs)) (Fin i, Fin j) a = L.M . L.V $ as ++ L.V (as' ++ a:bs') : bs
        where (as, L.V v:bs) = splitAt (j-1) vs
              (as', _:bs') = splitAt (i-1) v

    -- Not sure how we should handle zeros since we do not want a Eq constraint in the class
    values (L.M (L.V vs)) = [ ((Fin i, Fin j), a) | (j, L.V v) <- zip [1..] vs
                                                  , (i, a)     <- zip [1..] v ]

instance Matrix CSR.CSR where
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

          

