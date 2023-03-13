{-# Language FlexibleInstances #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE AllowAmbiguousTypes #-}


module Algebra where 

import GHC.TypeLits hiding ( type (^) )
import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), (^), (**), sum, product, recip)
import Data.Proxy 

type R = Double

type KnownNats m n = (KnownNat m, KnownNat n)

-- | Returns the type level Nat as an Int
--   Int wrapper for GHC.TypeLits natVal
natInt :: forall proxy n. KnownNat n => proxy n -> Int
natInt = fromInteger . natVal

natInt' :: forall n. KnownNat n => Int
natInt' = natInt (Proxy :: Proxy n)

-- | A list with a given length
newtype List a (n :: Nat) = L [a] deriving (Show, Eq)

-- Type signature important to enforce the "unchanged length" invariant
mapL :: (a -> b) -> List a n -> List b n
mapL f (L xs) = L (map f xs)


-------------------------------------------------------------------------------
-- Algebraic definitions

infixl 6 -
infixl 6 +

infixl 7 *
infixl 7 /

infixl 7 **
infixl 9 £

-- Class definitions
class AddGroup a where
    (+)  :: a -> a -> a

    (-)  :: a -> a -> a
    a - b = a + neg b

    neg :: a -> a
    neg a = zero - a

    zero :: a

sum :: (Foldable t, AddGroup a) => t a -> a
sum = foldr (+) zero 

-- | Allows us to define Matrix-matrix-,matrix-vector-multiplication as a single operation
--   Potentially this should replace or be used to define Mul 
class Composable a b c | a b -> c where
    (**) :: a -> b -> c

instance (b ~ b') => Composable (b -> c) (a -> b') (a -> c) where
    (**) = (.)

class Mul a where
    (*) :: a -> a -> a
    one :: a

product :: (Foldable t, Mul a) => t a -> a
product = foldr (*) one

(^) :: Mul a => a -> Int -> a
a ^ n = product $ replicate n a

type Ring a = (AddGroup a, Mul a)

two :: Ring a => a
two = one+one

class (Ring a) => Field a where
    (/) :: a -> a -> a
    a / b = a * recip b

    recip :: a -> a
    recip a = one / a


-- | Definition of a VectorSpace
class (AddGroup v, Ring (Under v)) => VectorSpace v where
    type Under v -- The underlying type of the VectorSpace
    type Under v = v

    (£) :: Under v -> v -> v

-- | All finite vectorspaces has a dimension
--   and can be represented by a basis
--   To generate a basis for a vector space v we can use TypeApplications @:
--   basis @(v)
--   e.g. basis @(R^2) = V [V [1.0,0.0],V [0.0,1.0]]
-- **TODO: please use working example - "R^2" is not valid Haskell in this file (perhaps works in ListVector? - if so, instruct the reader)
class VectorSpace x => Finite x where 
    type Dim x :: Nat
    type BasisVec x :: *
    type BasisVec x = x 
    basis' :: x -> List (BasisVec x) (Dim x)

basis :: forall x. (Finite x, AddGroup x, x ~ BasisVec x) => List x (Dim x)
basis = basis' (zero :: x)


-- | Approximated equality
--   Mostly useful for Double 
--
--   It might be useful to have two kinds of equalities,
--   one for absolute and one for relative.
--   For example 100001 and 100004 are "relatively" but not "absolute" equal.
--   Meanwhile 0.0001 and 0.000001 are "absolute" but not "relatively" equal
class Approx a where
    (~=) :: a -> a -> Bool

(~/=) :: Approx a => a -> a -> Bool
a ~/= b = not $ a ~= b


-------------------------------------------------------------------------------
-- Instance definitions
--

instance AddGroup Int       where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Integer   where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Double    where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Rational  where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup b => AddGroup (a -> b) where 
    f + g = \x -> f x + g x; neg f = \x -> neg (f x); zero = const zero

instance Mul Int          where (*) = (P.*); one = 1
instance Mul Integer      where (*) = (P.*); one = 1
instance Mul Double       where (*) = (P.*); one = 1
instance Mul Rational     where (*) = (P.*); one = 1
instance Mul b => Mul (a -> b)  where f * g = \x -> f x * g x; one = const one

instance Field Double       where (/) = (P./) 
instance Field Rational     where (/) = (P./) 
instance Field b => Field (a -> b) where recip f = \x -> recip (f x)  

instance VectorSpace Int      where (£) = (*)
instance VectorSpace Integer  where (£) = (*)
instance VectorSpace Double   where (£) = (*)
instance VectorSpace Rational where (£) = (*)
instance Ring b => VectorSpace (a -> b) where
    type Under (a -> b) = b
    s £ f = \x -> s * f x

instance Approx Double   where a ~= b = abs (a - b) < 0.00001
instance Approx Rational where a ~= b = abs (a - b) < 0.00001


