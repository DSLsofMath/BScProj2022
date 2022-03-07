{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies, TypeOperators #-}
{-# LANGUAGE FlexibleContexts, FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}


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
    Wrap :: (LinearMap x) => x -> From x --> To x
    -- We might add a Comp constructor in LinearMap
    -- Comp :: b --> a -> c --> b -> c --> a
    Zero :: b --> a

    
-- | Contains the basic properties of a linear map
--   Needed to make us of LinearMapType
class (VectorSpace x, VectorSpace (From x), VectorSpace (To x), Under (From x) ~ Under (To x), Under (To x) ~ Under x) => LinearMap x where 
    type To   x :: * 
    type From x :: *

    -- | Transforms the representation to a function
    toLinMapFun :: x -> From x -> To x


-- | Coposition of linear maps and addition
--   Since this is an overlappable class we can add new instances with 
--   more efficent compositions for certaint types
class (LinearMap x, LinearMap y) => Composable x y where 
    comp :: (To y ~ From x) => x -> y -> (From y --> To x)
    f `comp` g = Wrap $ toLinMapFun f . toLinMapFun g

    add :: (From x ~ From y, To x ~ To y) => x -> y -> From x --> To x
    x `add` y = Wrap $ toLinMapFun x + toLinMapFun y

instance {-# OVERLAPPABLE #-} (LinearMap x, LinearMap y) => Composable x y where 
    

-- A linear map is a vector space
instance AddGroup (b --> a) where
    Wrap a + Wrap b = a `add` b
    Zero   + a      = a
    a      + Zero   = a
    neg (Wrap a) = Wrap (neg a)
    neg Zero     = Zero
    zero = Zero

instance (VectorSpace a) => VectorSpace (b --> a) where 
    type Under (b --> a) = Under a
    s £ Wrap f = Wrap (s £ f)


-- An example of composing linear maps with different implementations
f :: R^3 -> R^2 
f (V [x,y,z]) = V [2*x, 3*y] 

m :: MatR 2 2 
m = toMat [[0,1], [-1,0]]

mf :: R^3 --> R^2
mf = m `comp` f


-- | Computes the linear map on a vector
apply :: (b --> a) -> b -> a
Wrap f `apply` a = toLinMapFun f a


--------------------------------------------
---- Instances of LinearMap

instance (KnownNat m, KnownNat n, Field f) => LinearMap (Matrix f m n ) where
    type To   (Matrix f m n) = Vector f m 
    type From (Matrix f m n) = Vector f n 
    toLinMapFun m = (m££)

-- | Since we have a better way to compose matrices than transforming them to funcion 
--   we use our own implementation
instance (KnownNat a, KnownNat b, KnownNat c, Field f) => Composable (Matrix f a b) (Matrix f b c) where
    a `comp` b = Wrap $ a £££ b
    a `add`  b = Wrap $ a + b


instance (VectorSpace a, VectorSpace b, Under a ~ Under b) => LinearMap (b -> a) where
    type To   (b -> a) = a
    type From (b -> a) = b
    toLinMapFun f = f 




