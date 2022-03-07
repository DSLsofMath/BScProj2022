{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies, TypeOperators #-}
{-# LANGUAGE FlexibleContexts, FlexibleInstances #-}


-- | A module for creating generic linear maps
module LinearMap where

import GHC.TypeLits hiding ( type (^) )
import Data.Coerce

import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), recip, sum, product, (**), span)

import ListVector hiding (m)
import Algebra


-- | Nice notation for a linear map
--   f :: R^3 --> R^4
infixr 1 -->
type b --> a = LinearMapType b a


-- | A general interface for all linear maps
--   By using an existential type capturing a linear map we can hidde the implementation 
--   and treat all representation of linear maps in a single type.
data LinearMapType b a where
    Wrap :: (AddGroup x, LinearMap x) => x -> From x --> To x
    -- Comp :: b --> a -> c --> b -> c --> a
    Zero :: b --> a

    
-- | Contains the basic properties of a linear map
--   Needed to make us of LinearMapType
class (AddGroup x, VectorSpace (From x), VectorSpace (To x), Under (From x) ~ Under (To x)) => LinearMap x where 
    type To   x :: * 
    type From x :: *

    -- | Transforms the representation to a function
    toLinMapFun :: x -> From x -> To x

    -- | Optional function for more efficent compositions if the types match
    compOwn :: (LinearMap y, Coercible x y, From x ~ To y) =>  x -> y -> From y --> To x
    compOwn x y = Wrap x `comp` Wrap y


-- A linear map is a vector space
instance AddGroup (b --> a) where
    Wrap a + Wrap b = Wrap (toLinMapFun a + toLinMapFun b)
    Zero   + a      = a
    a      + Zero   = a
    neg (Wrap a) = Wrap (neg a)
    neg Zero     = Zero
    zero = Zero

type instance Under (b --> a) = Under a

instance (VectorSpace a) => VectorSpace (b --> a) where 
    s £ Wrap f = Wrap (\x -> s £ (toLinMapFun f x))


-- An example of composing linear maps with different implementations
f :: R^3 -> R^2 
f (V [x,y,z]) = V [2*x, 3*y] 

m :: MatR 2 2 
m = toMat [[0,1], [-1,0]]

mf :: R^3 --> R^2
mf = Wrap m `comp` Wrap f

-- | Coposition of linear maps
--   Alternativly we add a Comp constructor in LinearMap
comp :: (b --> a) -> (c --> b) -> (c --> a)
Wrap f `comp` Wrap g = Wrap $ toLinMapFun f . toLinMapFun g

-- | Computes the linear map on a vector
apply :: (b --> a) -> b -> a
Wrap f `apply` a = toLinMapFun f a



-- Instances of LinearMap
instance (KnownNat m, KnownNat n, Field f) => LinearMap (Matrix f m n ) where
    type To   (Matrix f m n) = Vector f m 
    type From (Matrix f m n) = Vector f n 
    toLinMapFun m = (m££)
    -- compOwn a b = Wrap $ (a) £££ (b)


instance (VectorSpace a, VectorSpace b, Under a ~ Under b) => LinearMap (b -> a) where
    type To   (b -> a) = a
    type From (b -> a) = b
    toLinMapFun f = f 







