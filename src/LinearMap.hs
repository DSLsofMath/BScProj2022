{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies, TypeOperators #-}
{-# LANGUAGE FlexibleContexts, FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}


-- | A module for creating generic linear maps
module LinearMap where

import GHC.TypeLits hiding ( type (^) )

import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), recip, sum, product, (**), span)

import ListVector hiding (m)
import Algebra


-- | Nice notation for a linear map
--   Especially useful for the interactive experince
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

-- | Show instance for a linear map if it can be represented as a matrix
instance (Finite b, Show f) => Show (b --> (Vector f m)) where
    show = show . toListMat

    
-- | Contains the basic properties of a linear map
--   Needed to make us of LinearMapType
class (VectorSpace x, VectorSpace (From x), VectorSpace (To x), Under (From x) ~ Under (To x), Under (To x) ~ Under x) => LinearMap x where 
    type To   x :: * 
    type From x :: *

    -- Save a list of prefered representations 
    -- Such that efficent metods can be used when 
    -- combining linear maps
    -- preferedRep :: x -> [Reps]

    -- | Transforms the representation to a function
    toLinMapFun :: x -> From x -> To x



-- | Coposition of linear maps and addition
--   Since this is an overlappable class we can add new instances with 
--   more efficent compositions for certaint types
--   NOTE: The use of this class does not currently work when extracting 
--   values from Wrap, which was the main motivation for its existance
class (LinearMap x, LinearMap y) => Composable x y where 
    comp :: (To y ~ From x) => x -> y -> (From y --> To x)
    f `comp` g = Wrap $ toLinMapFun f . toLinMapFun g

    add :: (From x ~ From y, To x ~ To y) => x -> y -> From x --> To x
    x `add` y = Wrap $ toLinMapFun x + toLinMapFun y

instance {-# INCOHERENT #-} (LinearMap x, LinearMap y) => Composable x y where 
    

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
f :: R^2 -> R^2 
f (V [x,y]) = V [2*x, 3*y] 

m :: MatR 2 2 
m = toMat [[0,1], [-1,0]]

mf :: R^2 --> R^2
mf = m `comp` f


-- | Computes the linear map on a vector
apply :: (b --> a) -> b -> a
Wrap f `apply` a = toLinMapFun f a


-- | If a linear map goes from a finite space to a space with our Vector type
--   it can be represented as a matrix
toListMat :: (Finite b) => (b --> Vector f m) -> Matrix f m (Dim b)
toListMat (Wrap x) = let L bs = basis in M . V $ map (toLinMapFun x) bs


prop_MatToFunToMat :: (KnownNat n, LinearMap (Matrix f m n), Field f, Eq f) => Matrix f m n -> Bool
prop_MatToFunToMat m = m == toListMat (Wrap (toLinMapFun m))



-- | We might be able to generalize this further such that any linear map 
--   from and to a finite vector space can be represented as a matrix
instance (KnownNat m, KnownNat n, Field f) => ToMat m n (Vector f n --> Vector f m) where
    type Under' (Vector f n --> Vector f m) = f
    toMat = toListMat 
    fromMat = Wrap


--------------------------------------------
---- Instances of LinearMap

instance (KnownNat m, KnownNat n, Field f) => LinearMap (Matrix f m n ) where
    type To   (Matrix f m n) = Vector f m 
    type From (Matrix f m n) = Vector f n 
    toLinMapFun m = (m££)

-- | Since we have a better way to compose matrices than transforming them to funcion 
--   we use our own implementation
--   NOTE: Curently LinearMapType will always use the most general instance 
instance (KnownNat a, KnownNat b, KnownNat c, Field f) => Composable (Matrix f a b) (Matrix f b c) where
    a `comp` b = Wrap $ a £££ b
    a `add`  b = Wrap $ a + b


instance (VectorSpace a, VectorSpace b, Under a ~ Under b) => LinearMap (b -> a) where
    type To   (b -> a) = a
    type From (b -> a) = b
    toLinMapFun f = f 



