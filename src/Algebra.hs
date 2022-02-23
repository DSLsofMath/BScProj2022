{-# Language FlexibleInstances #-}
{-# LANGUAGE TypeFamilies #-}


module Algebra where 

import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), sum, recip)
type R = Double

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

class AddGroup a => Field a where
    (*) :: a -> a -> a

    (/) :: a -> a -> a
    a / b = a * recip b

    recip :: a -> a
    recip a = one / a

    one :: a

sum :: (Foldable t, AddGroup a) => t a -> a
sum = foldr (+) zero 


-- | The underlying type (often a field) of x  
type family   Under x 
type instance Under Double         = Double
type instance Under (_ -> f)       = f
type instance Under [[f]]          = f

-- | Definition of a VectorSpace
class (AddGroup v) => VectorSpace v where
    (£) :: Under v -> v -> v



-- Instance definitions
instance AddGroup Int       where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Integer   where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Double    where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Rational  where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Exp       where (+) = (:+:); neg = Negate; zero = Const 0
instance AddGroup b => AddGroup (a -> b)  where 
    f + g = \x -> f x + g x; neg f = \x -> neg (f x); zero = const zero

instance Field Double       where (*) = (P.*); (/) = (P./); one = 1
instance Field Rational     where (*) = (P.*); (/) = (P./); one = 1
instance Field Exp          where (*) = (:*:); recip = Recip; one = Const 1
instance Field b => Field (a -> b)  where 
    f * g = \x -> f x * f x; recip f = \x -> recip (f x); one = const one

