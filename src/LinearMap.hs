{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies, TypeOperators #-}
{-# LANGUAGE FlexibleContexts, FlexibleInstances #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}


-- | A module for creating generic linear maps
module LinearMap where

import GHC.TypeLits hiding ( type (^) )

import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), recip, sum, product, (**), span)

import ListVector hiding (m)
import qualified SparseCSR
import qualified QuadTree 
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
    Wrap :: (LinearMap x) => Rep x -> x -> From x --> To x
    -- We might add a Comp constructor in LinearMap
    -- Comp :: b --> a -> c --> b -> c --> a
    Zero :: b --> a


-- | To be able to extract the implementation from a linear map
--   we wrap it with a semi-singelton representation of its type.
--   
--   Note: we might want to expand upon this structure to also capture invariants
--   such as invertability or upper/lower echelon
data Rep x where
    Mat   :: (KnownNat m, KnownNat n, Field f) => Rep (Matrix f m n)
    Fun   :: Rep (b -> a)
    CSR   :: Rep (SparseCSR.CSR f m n)
    Quad  :: Rep (QuadTree.Quad n f)
    Other :: Rep a
    ZeroR :: Rep ()

deriving instance Show (Rep x)

-- | Returns a string of the underlying type of the linear map
--   Useful for debuging
getRep :: (b --> a) -> String 
getRep Zero       = "Zero"
getRep (Wrap r _) = show r


-- | Show instance for a linear map if it can be represented as a matrix
instance (Finite b, Show f) => Show (b --> (Vector f m)) where
    show = show . toListMat
    

-- | Contains the basic properties of a linear map
--   Needed to make us of LinearMapType
class (VectorSpace x, VectorSpace (From x), VectorSpace (To x), 
       Under (From x) ~ Under (To x), Under (To x) ~ Under x) => LinearMap x where 
    type To   x :: * 
    type From x :: *

    wrap :: x -> From x --> To x

    -- | Transforms the representation to a function
    toLinMapFun :: x -> From x -> To x


-- A linear map is a vector space
instance AddGroup (b --> a) where
    Wrap rA a + Wrap rB b = case (rA, rB) of 
        (Mat, Mat) -> wrap $ a + b
        _          -> wrap $ toLinMapFun a + toLinMapFun b
    Zero   + a      = a
    a      + Zero   = a
    neg (Wrap _ a) = wrap (neg a)
    neg Zero     = Zero
    zero = Zero

instance (VectorSpace a) => VectorSpace (b --> a) where 
    type Under (b --> a) = Under a
    s £ Wrap _ f = wrap (s £ f)


-- An example of composing linear maps with different implementations
f :: R^2 -> R^2 
f (V [x,y]) = V [2*x, 3*y] 

m :: MatR 2 2 
m = toMat [[0,1], [-1,0]]

mf :: R^2 --> R^2
mf = wrap m `comp` wrap f

-- | Composition of linear maps
--   Here we may pattern match on the representations to use a more efficent way of composing
--
--   Note: We would like to use a class for composition, but its hard to implement it such
--   that it does not default to the most generall instance
comp :: (b --> a) -> (c --> b) -> (c --> a)
Zero      `comp` _         = Zero
_         `comp` Zero      = Zero
Wrap rA a `comp` Wrap rB b = case (rA, rB) of
    (Mat, Mat) -> wrap $ a £££ b 
    (_, _) -> wrap $ toLinMapFun a . toLinMapFun b


-- | Computes the linear map on a vector
apply :: (b --> a) -> b -> a
Wrap _ f `apply` a = toLinMapFun f a


-- | If a linear map goes from a finite space to a space with our Vector type
--   it can be represented as a matrix
toListMat :: (Finite b) => (b --> Vector f m) -> Matrix f m (Dim b)
toListMat (Wrap Mat x) = x
toListMat (Wrap _   x) = let L bs = basis in M . V $ map (toLinMapFun x) bs


prop_MatToFunToMat :: (KnownNat n, LinearMap (Matrix f m n), Field f, Eq f) => Matrix f m n -> Bool
prop_MatToFunToMat m = m == toListMat (wrap (toLinMapFun m))


-- | We might be able to generalize this further such that any linear map 
--   from and to a finite vector space can be represented as a matrix
instance (KnownNat m, KnownNat n, Field f) => ToMat m n (Vector f n --> Vector f m) where
    type Under' (Vector f n --> Vector f m) = f
    toMat = toListMat 
    fromMat = wrap


--------------------------------------------
-- Instances of LinearMap
-- Note that these also needs to be represented internally if we want to use them efficently.
-- As of now they will simply be transformed to functions when composing

instance (KnownNat m, KnownNat n, Field f) => LinearMap (Matrix f m n ) where
    type To   (Matrix f m n) = Vector f m 
    type From (Matrix f m n) = Vector f n 
    wrap x = Wrap Mat x
    toLinMapFun m = (m££)

instance (VectorSpace a, VectorSpace b, Under a ~ Under b) => LinearMap (b -> a) where
    type To   (b -> a) = a
    type From (b -> a) = b
    wrap x = Wrap Fun x
    toLinMapFun f = f 



