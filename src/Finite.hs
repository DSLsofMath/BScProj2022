{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE AllowAmbiguousTypes #-}


module Finite where

import GHC.TypeLits
import Prelude hiding ((+), (-), (*), (/), sum)
import Data.Foldable (toList)

import Algebra
import FiniteIndex


-- | All finite vector spaces has a dimension
--   and each element can be represented by a linear combination of its basis.
--   That is, a vector v of dimension n can be written as:
--        v = a1 £ e1  +  a2 £ e2  +  ...  +  an £ en 
--   where ei are basis vectors and ai are scalars.
--
--   Notably, we can index the scalars given a basis.
--
--   Although there are arbitrarily many basis, Finite restricts a type to a single 
--   basis to make implementation reasonable.
--   Nonetheless, the choice of a "standard" basis is often very natural.
--
--   For using arbitrary basis see Subspaces (TODO rewrite/create new Subspaces)
--
class (KnownNat (Dim v), VectorSpace v) => Finite v where 
    type Dim v :: Nat

    -- | The standard basis vectors of v
    basis :: Fin (Dim v) -> v
    
    -- | Indexing the scalars of v in its standard basis
    --   That is, the scalars in the linear combination corresponding to v 
    scalar :: v -> (Fin (Dim v) -> Under v)
    

-- | Given an indexing of scalars we can construct a vector
fromIndexing :: Finite v => (Fin (Dim v) -> Under v) -> v
fromIndexing f = sum (\i -> f i £ basis i)


-- | Two vector spaces -- over the same field -- are isomorphic iff they have the same dimension
--   That is, there exists an bijective map between them.
type Isomorphic v w = (Finite v, Finite w, Dim v ~ Dim w, Under v ~ Under w)

-- | Map a vector to an isomorphic vector space
isomorphMap :: Isomorphic v w => v -> w
isomorphMap = fromIndexing . scalar



showFinite :: (Show (Under v), Finite v) => v -> String
showFinite = showIndixing . scalar

showIndixing :: (KnownNat n, Show f) => (Fin n -> f) -> String
showIndixing = ('\n':) . unlines . map formatRow . padCol . map show . toList


padCol :: [String] -> [String]
padCol col = map padString col
    where longest = maximum $ map length col
          padString s = replicate (longest - length s) ' ' ++ s ++ " "

formatRow :: String -> String
formatRow s = "| " ++ s ++ "|"

-------------------------------------------------------------------
-- Instances of  Fin n -> a 
--
-- Since functions from  Fin n  play an important role for our
-- definition of finite vector spaces we add some helpful instances.

instance KnownNat n => Foldable ((->) (Fin n)) where
    foldr op b f = foldr (\i acc -> f i `op` acc) b indices

instance (KnownNat n, Show a) => Show (Fin n -> a) where
    show = showIndixing 


-------------------------------------------------------------------
-- Basic instances

instance Finite Int      where type Dim Int      = 1; basis _ = 1; scalar x _ = x
instance Finite Integer  where type Dim Integer  = 1; basis _ = 1; scalar x _ = x
instance Finite Double   where type Dim Double   = 1; basis _ = 1; scalar x _ = x
instance Finite Rational where type Dim Rational = 1; basis _ = 1; scalar x _ = x

instance (KnownNat n, Ring f) => Finite (Fin n -> f) where
    type Dim (Fin n -> f) = n
    basis i = \j -> if i == j then one else zero
    scalar = id



