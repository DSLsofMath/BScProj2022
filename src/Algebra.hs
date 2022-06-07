{-# Language FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}

{-# LANGUAGE GADTs #-}
{-# LANGUAGE StandaloneDeriving #-}


module Algebra where 

import GHC.TypeLits hiding ( type (^) )
import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), (^), (**), sum, product, recip, fromRational)
import qualified Data.Ratio

type R = Double

type KnownNats m n = (KnownNat m, KnownNat n)

-- | A list with a given length
newtype List a (n :: Nat) = L [a] deriving (Show, Eq)

-- Type signature important to enforce the "unchanged length" invariant
mapL :: (a -> b) -> List a n -> List b n
mapL f (L xs) = L (map f xs)

-- | Datatype for field expressions of one variable
data Exp a where
    Const  :: a -> Exp a
    X      :: Exp a
    (:+:)  :: AddGroup a => Exp a -> Exp a -> Exp a
    Negate :: AddGroup a => Exp a -> Exp a
    (:*:)  :: Mul a      => Exp a -> Exp a -> Exp a
    Recip  :: Field a    => Exp a -> Exp a
    -- (:£:)  :: VectorSpace a => Exp (Under a) -> Exp a -> Exp a

deriving instance Eq a => Eq (Exp a)

-- Remove + const 0 from Exp matrices
-------------------------------------Only to simplify how to view matrix expressions and determinants--------------------------------

instance (Show a, Field a, Num a, Eq a) => Show (Exp a) where
   show = showE . simplifyE

showE :: Show a => Exp a -> String
showE X                 = "X"
showE (e1 :+: e2)       = showE e1 ++ " + " ++ showE e2
showE (e1 :*: Recip e2) = showFactor e1 ++ " / " ++ showFactor e2 
showE (e1 :*: e2)       = showFactor e1 ++ " * " ++ showFactor e2 
showE (Recip e)         = "1 / " ++ showFactor e 
showE (Negate e)        = "-" ++ showFactor e
showE (Const r)         = show r

-- | Handles cases when parenthesis need to be shown
showFactor :: Show a => Exp a -> String
showFactor e@(_ :+: _)  = "(" ++ showE e ++ ")"
showFactor e@(Negate _) = "(" ++ showE e ++ ")"
showFactor e            = showE e 


------------------------------------------------------------------------------------------------------------------------------------

-- Eval functions for expressions

instance Expr Exp where
    eval = evalExp
    derive = deriveExp

-- Eval for expressions, apply a value for X
evalExp :: Exp a -> a -> a
evalExp (Const a)   = const a
evalExp X           = id
evalExp (e1 :+: e2) = evalExp e1 + evalExp e2
evalExp (e1 :*: e2) = evalExp e1 * evalExp e2
evalExp (Negate e)  = neg (evalExp e)
evalExp (Recip e)   = recip (evalExp e)

-- derive  ::  FunExp        ->  FunExp
deriveExp :: Ring a => Exp a -> Exp a
deriveExp (Const _)   = Const zero
deriveExp X           = Const one
deriveExp (e1 :+: e2) = deriveExp e1 :+: deriveExp e2
deriveExp (e1 :*: e2) = (deriveExp e1 :*: e2) :+: (e1 :*: deriveExp e2)
deriveExp (Negate e)  = neg (deriveExp e)
deriveExp (Recip e)   = neg (deriveExp e :*: (recip (e^2)))

evalExp' :: Ring a => Exp a -> a -> a
evalExp' =  evalExp . deriveExp

evalExp'' :: Ring a => Exp a -> a -> a
evalExp'' = evalExp' . deriveExp


-- Added this implementation from DSLsofMath to simplify algebraic expressions
simplifyE :: (Field a, Num a, Eq a) => Exp a -> Exp a
simplifyE (x :+: y)   = simplifyAdd (simplifyE x) (simplifyE y)
simplifyE (x :*: y)   = simplifyMul (simplifyE x) (simplifyE y)
simplifyE (Negate x)  = simplifyNeg (simplifyE x)
simplifyE (Recip x)   = simplifyRec (simplifyE x)
simplifyE X           = X
simplifyE (Const c)   = Const c

simplifyAdd :: (Ring a, Num a, Eq a) => Exp a -> Exp a -> Exp a
simplifyAdd (Const a)   (Const b)   = Const (a+b)
simplifyAdd (Const 0)   x           = x
simplifyAdd x           (Const 0)   = x
simplifyAdd (x:+:y)     z           = simplifyAdd x (y:+:z)
simplifyAdd (Const a)   (Const b:+:x) = simplifyAdd (Const (a+b)) x
simplifyAdd x           y             = case scaledEq x y of
    Left    (a,b,x) -> simplifyMul (Const (a+b)) x
    Right   (x,y)   -> x :+: y

simplifyMul :: (Ring a, Num a, Eq a) => Exp a -> Exp a -> Exp a
simplifyMul (Const a)  (Const b)  = Const (a*b)
simplifyMul x          (Const c)  = simplifyMul (Const c) x
simplifyMul (Const 0)  _          = Const 0
simplifyMul (Const 1)  x          = x
simplifyMul (x:*:y)    z          = simplifyMul x (y:*:z)
simplifyMul (Const a)  (Const b :*: x) = simplifyMul (Const (a*b)) x
simplifyMul x          (Const c :*: y) = simplifyMul (Const c) (x :*: y)
simplifyMul (Negate a)  b         = simplifyNeg $ simplifyMul a b
simplifyMul a          (Negate b) = simplifyMul (simplifyNeg a) b
simplifyMul x          y          = x :*: y

scaledEq :: (Num a, Eq a) => Exp a -> Exp a -> Either (a, a, Exp a) (Exp a, Exp a)
scaledEq x y                         | x==y = Left (1,1,x)
scaledEq            x  (Const b:*:y) | x==y = Left (1,b,x)
scaledEq (Const a:*:x)            y  | x==y = Left (a,1,x)
scaledEq (Const a:*:x) (Const b:*:y) | x==y = Left (a,b,x)
scaledEq x y = Right (x,y)

simplifyNeg :: AddGroup a => Exp a -> Exp a
simplifyNeg (Const c)   = Const (neg c)
simplifyNeg (Negate x)  = x
simplifyNeg x           = Negate x

simplifyRec :: Field a => Exp a -> Exp a
simplifyRec (Const c)   = Const (recip c)
simplifyRec (Recip x)   = x
simplifyRec (Negate x)  = simplifyNeg $ simplifyRec x
simplifyRec x           = Recip x

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


-- | General class for Ring expression
class Expr (exp :: * -> *) where
    eval :: Ring a => exp a -> a -> a

    derive :: Ring a => exp a -> exp a

eval' :: (Expr exp, Ring a) => exp a -> a -> a
eval' = eval . derive


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


-- Instance definitions
instance AddGroup Int       where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Integer   where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Double    where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup Rational  where (+) = (P.+); (-) = (P.-); zero = 0
instance AddGroup a => AddGroup (Exp a)  where (+) = (:+:); neg = Negate; zero = Const zero
instance AddGroup b => AddGroup (a -> b) where 
    f + g = \x -> f x + g x; neg f = \x -> neg (f x); zero = const zero

instance Mul Int          where (*) = (P.*); one = 1
instance Mul Integer      where (*) = (P.*); one = 1
instance Mul Double       where (*) = (P.*); one = 1
instance Mul Rational     where (*) = (P.*); one = 1
instance Mul a => Mul (Exp a)   where (*) = (:*:); one = Const one
instance Mul b => Mul (a -> b)  where f * g = \x -> f x * g x; one = const one

instance Field Double       where (/) = (P./) 
instance Field Rational     where (/) = (P./) 
instance Field a => Field (Exp a)  where recip = Recip 
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


