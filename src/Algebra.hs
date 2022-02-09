{-# Language TypeSynonymInstances #-}
{-# Language FlexibleInstances #-}

module Algebra where 

import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), sum)


infixl 6 -
infixl 6 +

infixl 7 *
infixl 7 /

-- Class definitions
class AddGroup a where
    (+)  :: a -> a -> a
    (-)  :: a -> a -> a
    zero :: a

class AddGroup a => Field a where
    (*) :: a -> a -> a
    (/) :: a -> a -> a
    one :: a

sum :: (Foldable t, AddGroup a) => t a -> a
sum = foldr (+) zero 


-- Instance definitions
instance AddGroup Int       where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Integer   where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Double    where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Rational  where (+) = (P.+); (-) = (P.-); zero = 0

instance Field Double       where (*) = (P.*); (/) = (P./); one = 1
instance Field Rational     where (*) = (P.*); (/) = (P./); one = 1


