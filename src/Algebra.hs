{-# Language FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DataKinds #-}


module Algebra where 

import GHC.TypeLits hiding ( type (^) )
import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), sum, product, recip)

type R = Double


-- | A list with a given length
newtype List a (n :: Nat) = L [a] deriving (Show, Eq)


data Exp =  Const R
             |  X
             |  Exp :+: Exp
             |  Exp :*: Exp
             |  Recip Exp
             |  Negate Exp
    deriving (Eq, Ord, Show)

infixl 6 -
infixl 6 +

infixl 7 *
infixl 7 /

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


class Mul a where
    (*) :: a -> a -> a
    one :: a

product :: (Foldable t, Mul a) => t a -> a
product = foldr (*) one


class (AddGroup a, Mul a) => Field a where
    (/) :: a -> a -> a
    a / b = a * recip b

    recip :: a -> a
    recip a = one / a


-- | Definition of a VectorSpace
class (AddGroup v, AddGroup (Under v)) => VectorSpace v where
    type Under v -- The underlying type of the VectorSpace
    (£) :: Under v -> v -> v

-- | All finite vectorspaces has a dimension
--   and can be represented by a basis
--   To generate a basis for a vector space v we can use TypeApplications @:
--   basis @(v)
--   e.g. basis @(R^2) = V [V [1.0,0.0],V [0.0,1.0]]
--
class VectorSpace x => Finite x where 
    type Dim x :: Nat
    basis :: List x (Dim x)


-- Instance definitions
instance AddGroup Int       where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Integer   where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Double    where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Rational  where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Exp       where (+) = (:+:); neg = Negate; zero = Const 0
instance AddGroup b => AddGroup (a -> b)  where 
    f + g = \x -> f x + g x; neg f = \x -> neg (f x); zero = const zero

instance Mul Int          where (*) = (P.*); one = 1
instance Mul Integer      where (*) = (P.*); one = 1
instance Mul Double       where (*) = (P.*); one = 1
instance Mul Rational     where (*) = (P.*); one = 1
instance Mul Exp          where (*) = (:*:); one = Const 1
instance Mul b => Mul (a -> b)  where f * g = \x -> f x * f x; one = const one

instance Field Double       where (/) = (P./); 
instance Field Rational     where (/) = (P./); 
instance Field Exp          where recip = Recip; 
instance Field b => Field (a -> b)  where recip f = \x -> recip (f x); 

instance VectorSpace b => VectorSpace (a -> b) where
    type Under (a -> b) = Under b
    s £ f = \x -> s £ f x



