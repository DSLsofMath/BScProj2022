{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}

import GHC.TypeLits
import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), (^), sum, product, recip, fromRational)

import Algebra


newtype Fin (n :: Nat) = Fin Int deriving (Eq, Ord) 

fin :: KnownNat n => Int -> Fin n
fin i = finite
    where finite | 1 <= i && i <= fromInteger (natVal finite) = Fin i
                 | otherwise = error $ "Index is out of bounds, got: " ++ 
                   show i ++ " in constraint 0<=" ++ show i ++ "<=" ++ show (natVal finite)

instance KnownNat n => Show (Fin n) where
    show n = show (finToInt n) ++ " of " ++ show (natVal n)

finToInt :: Fin n -> Int
finToInt (Fin n) = n


-- | Generic class for a matrix with focus on sparse representations
class VectorSpace mat => Matrix mat (m :: Nat) (n :: Nat) | mat -> m n where
    {-# MINIMAL set , values #-}

    get :: mat -> (Fin m, Fin n) -> Under mat
    get m i = case lookup i (values m) of Just a  -> a
                                          Nothing -> zero

    set :: mat -> (Fin m, Fin n) -> Under mat -> mat

    -- | Returns all nonzero values in a matrix
    values :: mat -> [((Fin m, Fin n), Under mat)]

    -- | Builds a matrix from a list of positions and values
    --   Initializes with the zero matrix
    tabulate :: [((Fin m, Fin n), Under mat)] -> mat
    tabulate xs = foldl (\m (i, a) -> set m i a) zero xs

    
getRow :: Matrix mat m n => mat -> Fin m -> [(Fin n, Under mat)]
getRow m i = [ (j, a) | ((i', j), a) <- values m, i' == i]

getCol :: Matrix mat m n => mat -> Fin n -> [(Fin m, Under mat)]
getCol m j = [ (i, a) | ((i, j'), a) <- values m, j' == j]



