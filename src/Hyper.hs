{-# LANGUAGE DataKinds #-} 
{-# LANGUAGE TypeOperators #-} 
{-# LANGUAGE TypeFamilyDependencies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-} -- Required by type family N
{-# LANGUAGE GADTs #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE PolyKinds #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE NoStarIsType #-}
{-# LANGUAGE AllowAmbiguousTypes #-}

{-# LANGUAGE UndecidableSuperClasses #-}

module TypedVector where

import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), sum, (**), lookup)

import GHC.TypeLits 
import Data.Kind
import Data.Proxy 
import Data.Type.Equality 

import qualified Data.List as L 

import Algebra
import ListVector
import Matrix

import Control.Applicative (liftA2)


-- Vector inspired by Jeremy Gibbons paper APLicative Programming with Naperian Functors
-- 
-- The focus of the paper is about lifting lower dimensional types to higher dimensional,
-- for example scalars to vectors and vectors to matrices, and gives the required 
-- structure in the form of type classes. It manages this by using fixed shapes 
-- i.e. where the structure is known at compile time, such as Vector a n or (a,a,a). 
-- With these fixed shapes they generalises matrices as vectors of vectors.
--
-- Problematically for us it discourages sparse representations. 
-- Partly since a sparse matrix's structure changes based on its value 
-- and partly since it would be seen as a single dimension.
-- However, it might be possible for us to adjust some of the core ideas.
--
-- Concretely, for sparse matrices we know that the lack of a value corresponds to a 0,
-- so instead of requiring replicate or pure we require a "zero state".
-- In the Naperian class, or in a new Matrix class, we replace (Index f -> a) 
-- with [(Index f, a)], a list containing all nonzero elements and their position.
-- Alternatively, since we only care about matrices, as [(Int, Int, a)] 
-- or more well typed as [(Fin m, Fin n, a)]
-- If all matrices shared this common interface then conversion between different 
-- representations would be trivial and it would also allow us to define functions 
-- that work on all matrices. 
-- Furthermore, since class functions are over writable, we still allow for 
-- efficient implementations on base functions such as multiplication and utf.
--

-------------------------------------------
-- Help definitions

data TList :: (Nat -> Type) -> [Nat] -> Type where
    Nil  :: TList e '[]
    (:>) :: (KnownNat n, Shapely ns) => e n -> TList e ns -> TList e (n:ns)

data SomeVec v a = forall n. KnownNat n => Some (v n a)

-- deriving instance (forall i. Show     (v i a)) => Show     (Some v a)
-- deriving instance (forall i. Functor  (v i))   => Functor  (Some v)
-- deriving instance (forall i. Foldable (v i))   => Foldable (Some v)


-----------------------------------------------------------------------------
-- Hyper 


type a ^ ns = Hyper ns a
data Hyper :: [Nat] -> Type -> Type where
    Scalar :: a -> Hyper '[] a
    Prism :: (KnownNat n, Shapely ns) => Hyper ns (Vector n a) -> Hyper (n:ns) a

unPrism :: Hyper (n:ns) a -> Hyper ns (Vector n a)
unPrism (Prism a) = a



instance Show a => Show (Hyper ns a) where
    show (Scalar a) = show a
    show (Prism (Prism a)) = show $ fmap M a
    show (Prism a) = show a

deriving instance Functor (Hyper ns)

class Shapely ns where
    hreplicate :: a -> Hyper ns a
    href       :: TList Proxy ns

instance Shapely '[] where 
    hreplicate = Scalar
    href = Nil

instance forall n ns. (KnownNat n, Shapely ns) => Shapely (n:ns) where 
    hreplicate a = Prism (hreplicate (pure a))
    href = Proxy @n :> href

hzipWith :: (a -> b -> c) -> Hyper ns a -> Hyper ns b -> Hyper ns c
hzipWith f (Scalar a) (Scalar b) = Scalar (f a b)
hzipWith f (Prism a)  (Prism b)  = Prism (hzipWith (zipWithV f) a b) 

instance Shapely ns => Applicative (Hyper ns) where
    pure  = hreplicate
    (<*>) = hzipWith (\f a -> f a)

instance (AddGroup a, Shapely ns) => AddGroup (Hyper ns a) where
    (+) = liftA2 (+)
    (-) = liftA2 (-)
    zero = pure zero

hreduce :: (a -> a -> a) -> a -> Hyper (n:ns) a -> Hyper ns a
hreduce f b (Prism as) = fmap (foldl f b ) as

hsum :: AddGroup a => Hyper (n:ns) a -> Hyper ns a
hsum = hreduce (+) zero


class (Shapely ms, Shapely ns) => Alignable ms ns where 
    align :: Hyper ms a -> Hyper ns a

instance Alignable '[] '[] where
    align = id

instance (KnownNat n, Alignable ns' ns) => Alignable (n:ns') (n:ns) where
    align (Prism x) = Prism (align x)

instance (KnownNat n, Shapely ns) => Alignable '[] (n:ns) where
    align (Scalar a) = hreplicate a

alignLeft :: (Alignable ns ms) => (a -> b -> c) -> Hyper ns a -> Hyper ms b -> Hyper ms c
alignLeft f a = hzipWith f (align a)

alignRight :: (Alignable ms ns) => (a -> b -> c) -> Hyper ns a -> Hyper ms b -> Hyper ns c
alignRight f a b = hzipWith f a (align b)


-- | Matrix multiplication where the left argument can be a vector, matrix or higher dimensional 
mul :: (Ring a, KnownNat m, Alignable '[] ns) => Hyper '[n, m] a -> Hyper (m:ns) a -> Hyper (n:ns) a
mul (Prism m) v = Prism . hsum $ alignRight (£) v m



----------------------------------------------------------------
-- Transforming Hypers of equal dimension but different shape

type family Product ns where
            Product '[]    = 1
            Product (n:ns) = n * (Product ns) 


flatten :: Hyper ns a -> Vector (Product ns) a
flatten (Scalar a) = pure a
flatten (Prism a)  = flatVec (flatten a)

inflate :: forall ms a. Shapely ms => Vector (Product ms) a -> Hyper ms a 
inflate v = case href @ms of
    Nil                  -> Scalar . head . unVec $ v -- head is safe since Product '[] = 1
    (_ :: Proxy n) :> ns -> Prism (inflate $ groupVec @n v) 

reShape :: (Shapely ns, Product ms ~ Product ns) => Hyper ms a -> Hyper ns a
reShape = inflate . flatten


----------------------------------------------------------------
-- Tests

vec3 :: Hyper '[3] Int
vec3 = Prism . Scalar $ vec [0,1,2]

mat32 :: Hyper '[3,2] Int
mat32 = Prism . Prism . Scalar $ vec [vec [1, 2, 3], vec [4, 5, 6]]


mat23 :: Hyper '[2, 3] Int
mat23 = Prism . Prism . Scalar $ vec [vec [1, 2], vec [4, 5], vec [7, 8]] 

mat222 :: Hyper '[2,2,2] Int
mat222 = Prism . Prism . Prism . Scalar $ vec [vec [vec [1, 2], vec [4, 5]], zero]


