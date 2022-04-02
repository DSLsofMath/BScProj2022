{-# LANGUAGE DataKinds #-} 
{-# LANGUAGE TypeOperators #-} 
{-# LANGUAGE TypeFamilies #-} 
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-} -- Required by type family N
{-# LANGUAGE GADTs #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE PolyKinds #-}

module TypedVector where

import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), sum, (**), lookup)

import qualified GHC.TypeLits as Lit

import qualified Data.List as L 
import Algebra


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
-- Vector definitions

-- | Using this implementation of Nats instead of GHC.TypeLits allows us 
--   to use induction on the types which helps ghc to deduce types. 
data Nat = Z | S Nat

-- | Allows us to convert literals, like 3, into our Nat 
--   N 3 = S (S (S Z))
type family N (n :: Lit.Nat) :: Nat where
    N 0 = Z
    N n = S (N (n Lit.- 1))

-- | Nicer type syntax for Vec
--   Syntactically the same as in ListVector
--   Allows us to encode the type as
--   zero :: Vector R 4
--   instead of
--   zero :: Vec (S (S (S (S Z)))) R
type Vector f n = Vec (N n) f


class Count (n :: Nat) where
    vreplicate :: a -> Vec n a

instance Count Z where
    vreplicate _ = Nil

instance Count n =>  Count (S n) where
    vreplicate a = a :> vreplicate a


infixr 5 :>
-- | Dependent typed vector
--   Since we encode the length in the constructor the vector is guaranteed 
--   to be of correct size.
data Vec (n :: Nat) f where
    Nil  :: Vec Z f
    (:>) :: f -> Vec n f -> Vec (S n) f  

instance Show f => Show (Vec n f) where show = show . toList

instance Eq a => Eq (Vec n a) where 
    v1 == v2 = and $ vzipWith (==) v1 v2

toList :: Vec n f -> [f]
toList = foldr (:) []

-- | Syntactically identical to vec from ListVector 
vec :: (Count n, AddGroup a) => [a] -> Vec n a
vec xs = vreplace xs (vreplicate zero)

vmap :: (a -> b) -> Vec n a -> Vec n b
vmap _ Nil     = Nil
vmap f (x:>xs) = f x :> vmap f xs

instance Functor (Vec n) where
    fmap = vmap


vzipWith :: (a -> b -> c) -> Vec n a -> Vec n b -> Vec n c
vzipWith _ Nil Nil = Nil
vzipWith op (a:>as) (b:>bs) = a `op` b :> vzipWith op as bs

-- | replaces elements in a vector with elements in a list
vreplace :: [a] -> Vec n a -> Vec n a
vreplace [] v           = v
vreplace _  Nil         = Nil
vreplace (a:as) (_:>bs) = a :> vreplace as bs


vfoldr :: (a -> b -> b) -> b -> Vec n a -> b
vfoldr _  b Nil     = b
vfoldr op b (x:>xs) = x `op` vfoldr op b xs

instance Foldable (Vec n) where
    foldr = vfoldr

azipWith :: Applicative f => (a -> b -> c) -> f a -> f b -> f c
azipWith f a b = f <$> a <*> b

instance Count n => Applicative (Vec n) where
    pure = vreplicate
    (<*>) = vzipWith (\f x -> f x)


class Functor f => Naperian f where
    type Index f
    lookup :: f a -> (Index f -> a)

    tabulate :: (Index f -> a) -> f a
    tabulate h = fmap h positions

    positions :: f (Index f )
    positions = tabulate id


transpose :: (Naperian f, Naperian g) => f (g a) -> g (f a)
transpose = tabulate . fmap tabulate . flip . fmap lookup . lookup


-- | Bounded naturals
--   Contains all Nats up to n
data Fin (n :: Nat) where
    FZ :: Fin (S n)
    FS :: Fin n -> Fin (S n)

deriving instance Show (Fin n)
deriving instance Eq (Fin n)

-- | Type safe (!!)
vlookup :: Vec n a -> Fin n -> a
vlookup (x:>_ ) FZ      = x
vlookup (_:>xs) (FS n) = vlookup xs n

-- | Returns a Vec containing its Index values
viota :: Count n => Vec n (Fin n)
viota = viota' $ vreplicate () 
    where viota' :: Vec n () -> Vec n (Fin n)
          viota' Nil = Nil
          viota' (_:>xs) = FZ :> fmap FS (viota' xs)

instance Count n => Naperian (Vec n) where
    type Index (Vec n) = Fin n
    lookup = vlookup
    positions = viota


---------------------------------------------------------------------
-- Linear Algebra definitions

-- Vec is a vector space over a ring 
instance (Count n, AddGroup f) => AddGroup (Vec n f) where
    (+) = vzipWith (+)
    (-) = vzipWith (-)
    zero = vreplicate zero

instance (Count n, Ring f) => VectorSpace (Vec n f) where
    type Under (Vec n f) = f
    s £ v = fmap (s*) v


dotProd :: (Applicative f, Foldable f, Ring a) => f a -> f a -> a
dotProd f1 f2 = sum $ azipWith (*) f1 f2

dot :: (Count n, Ring f) => Vec n f -> Vec n f -> f
dot = dotProd

cross :: Ring f => Vector f 3 -> Vector f 3 -> Vector f 3
(a1:>a2:>a3:>Nil) `cross` (b1:>b2:>b3:>Nil) = a2*b3-a3*b2
                                           :> a3*b1-a1*b3
                                           :> a1*b2-a2*b1 :>Nil


baseVec :: (Naperian f, Eq (Index f), Ring a) => Index f -> f a
baseVec x = tabulate (\n -> if n == x then one else zero)

e :: (Count n, Ring a) => Fin n -> Vec n a 
e = baseVec

identityMat :: (Naperian f, Eq (Index f), Ring a) => f (f a)
identityMat = tabulate baseVec

idm :: (Count n, Ring a) => Vec n (Vec n a)
idm = identityMat



