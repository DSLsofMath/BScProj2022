{-# LANGUAGE DataKinds #-} 
{-# LANGUAGE TypeOperators #-} 
{-# LANGUAGE TypeFamilies #-} 
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE UndecidableInstances #-} -- Required by type family N
{-# LANGUAGE GADTs #-}

module TypedVector where

import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), sum, (**), lookup)

import qualified GHC.TypeLits as Lit

import qualified Data.List as L 
import Algebra


-- Vector inspired by Jeremy Gibbons paper APLicative Programming with Naperian Functors

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
vfoldr _ b Nil = b
vfoldr op b (x:>xs) = x `op` vfoldr op b xs

instance Foldable (Vec n) where
    foldr = vfoldr


instance Count n => Applicative (Vec n) where
    pure = vreplicate
    (<*>) = vzipWith (\f x -> f x)


class Functor f => Naperian f where
    type Log f
    lookup :: f a -> (Log f -> a)

    tabulate :: (Log f -> a) -> f a
    tabulate h = fmap h positions

    positions :: f (Log f )
    positions = tabulate id


transpose :: (Naperian f, Naperian g) => f (g a) -> g (f a)
transpose = tabulate . fmap tabulate . flip . fmap lookup . lookup


-- | Bounded naturals
--   Contains all Nats up to n
data Fin (n :: Nat) where
    FZ :: Fin (S n)
    FS :: Fin n -> Fin (S n)

vlookup :: Vec n a -> Fin n -> a
vlookup (x:>_ ) FZ      = x
vlookup (_:>xs) (FS n) = vlookup xs n

viota :: Count n => Vec n (Fin n)
viota = viota' $ vreplicate () 
    where viota' :: Vec n () -> Vec n (Fin n)
          viota' Nil = Nil
          viota' (_:>xs) = FZ :> fmap FS (viota' xs)

instance Count n => Naperian (Vec n) where
    type Log (Vec n) = Fin n
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


dot :: Ring f => Vec n f -> Vec n f -> f
v1 `dot` v2 = sum $ vzipWith (*) v1 v2

cross :: Ring f => Vector f 3 -> Vector f 3 -> Vector f 3
(a1:>a2:>a3:>Nil) `cross` (b1:>b2:>b3:>Nil) = a2*b3-a3*b2
                                           :> a3*b1-a1*b3
                                           :> a1*b2-a2*b1 :>Nil


