{-# Language FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}



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

data Exp  =  Const R
          |  X
          |  Exp :+: Exp
          |  Negate Exp
          |  Exp :*: Exp
          |  Recip Exp
    deriving (Eq, Ord)


-- Remove + const 0 from Exp matrices
-------------------------------------Only to simplify how to view matrix expressions and determinants--------------------------------
instance Show Exp where
   show = showE . simplifyE

showE :: Exp -> String
showE X = "X"
showE (e :+: Const 0) = showE e
showE (Const 0 :+: e) = showE e
showE (Const x :+: Const y) = showE (Const(x+y))
showE (e1 :+: e2)     = showE e1 ++ " + " ++ showE e2

showE (_ :*: Const 0) = showE zero
showE (Const 0 :*: _) = showE zero
showE (e :*: Const 1) = showE e
showE (Const 1 :*: e) = showE e
showE (Const x :*: Const y) = showE (Const(x*y))
showE (e1 :*: e2)     = "(" ++ showE e1 ++ ") * (" ++ showE e2 ++ ")"

showE (Recip e)       = "(" ++ "1 / " ++ showE e ++ ")"
showE (Negate e)      = "-" ++ "("  ++ showE e ++ ")"
showE (Const r)       = show r
------------------------------------------------------------------------------------------------------------------------------------

-- Eval functions for expressions
fromRational :: (Field a, Num a) => Data.Ratio.Ratio Integer -> a
fromRational x  =  fromInteger (Data.Ratio.numerator x) / fromInteger (Data.Ratio.denominator x)


-- Eval for expressions, apply a value for X
evalExp :: (Field a, Num a) => Exp -> a -> a
evalExp (Const alpha)  =  const (fromRational (toRational alpha))
evalExp X              =  id
evalExp (e1 :+: e2)    =  evalExp e1 + evalExp e2
evalExp (e1 :*: e2)    =  evalExp e1 * evalExp e2
evalExp (Negate e)     =  neg (evalExp e)
evalExp (Recip e)      =  recip (evalExp e)

-- derive  ::  FunExp        ->  FunExp
derive :: Exp -> Exp
derive      (Const alpha)  =  Const 0
derive      X              =  Const 1
derive      (e1 :+: e2)    =  derive e1 :+: derive e2
derive      (e1 :*: e2)    =  (derive e1 :*: e2) :+: (e1 :*: derive e2)
derive      (Negate e)     =  neg (derive e)
derive      (Recip e)      =  neg (derive e :*: (recip (e^2)))

evalExp' :: (Field a, Num a) => Exp -> a -> a
evalExp' =  evalExp . derive

evalExp'' :: (Field a, Num a) => Exp -> a -> a
evalExp'' = evalExp'.derive



-- Added this implementation from DSLsofMath to simplify algebraic expressions
type SimplifiedExp = Exp
simplifyE  ::  Exp -> SimplifiedExp
simplifyE (x :+: y)   = simplifyAdd (simplifyE x) (simplifyE y)
simplifyE (x :*: y)   = simplifyMul (simplifyE x) (simplifyE y)
simplifyE (Negate x)  = simplifyNeg (simplifyE x)
simplifyE (Recip x)   = simplifyRec (simplifyE x)
simplifyE X           = X
simplifyE (Const c)   = Const c

simplifyAdd :: Exp -> Exp -> Exp
simplifyAdd (Const a)   (Const b)   = Const (a+b)
simplifyAdd (Const 0)   x           = x
simplifyAdd x           (Const 0)   = x
simplifyAdd (x:+:y)     z           = simplifyAdd x (y:+:z)
simplifyAdd (Const a)   (Const b:+:x) = simplifyAdd (Const (a+b)) x
simplifyAdd x           y             = case scaledEq x y of
    Left    (a,b,x) -> simplifyMul (Const (a+b)) x
    Right   (x,y)   -> x :+: y

simplifyMul :: Exp -> Exp -> Exp
simplifyMul (Const a)  (Const b)  = Const (a*b)
simplifyMul (Const 0)  _x         = Const 0
simplifyMul _x         (Const 0)  = Const 0
simplifyMul (Const 1)  x          = x
simplifyMul x          (Const 1)  = x
simplifyMul (x:*:y)    z          = simplifyMul x (y:*:z)
simplifyMul x          (Const c)  = simplifyMul (Const c) x
simplifyMul (Const a)  (Const b :*: x) = simplifyMul (Const (a*b)) x
simplifyMul x          (Const c :*: y) = simplifyMul (Const c) (x :*: y)
simplifyMul x          y          = x :*: y

scaledEq x y                         | x==y = Left (1,1,x)
scaledEq            x  (Const b:*:y) | x==y = Left (1,b,x)
scaledEq (Const a:*:x)            y  | x==y = Left (a,1,x)
scaledEq (Const a:*:x) (Const b:*:y) | x==y = Left (a,b,x)
scaledEq x y = Right (x,y)

simplifyNeg :: Exp -> Exp
simplifyNeg (Const c)   = Const (negate c)
simplifyNeg (Negate x)  = x
simplifyNeg x           = Negate x

simplifyRec :: Exp -> Exp
simplifyRec (Const c)   = Const (recip c)
simplifyRec (Recip x)   = x
simplifyRec x           = Recip x

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

class (Ring a) => Field a where
    (/) :: a -> a -> a
    a / b = a * recip b

    recip :: a -> a
    recip a = one / a


-- | Definition of a VectorSpace
class (AddGroup v, Ring (Under v)) => VectorSpace v where
    type Under v -- The underlying type of the VectorSpace
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
instance Mul b => Mul (a -> b)  where f * g = \x -> f x * g x; one = const one

instance Field Double       where (/) = (P./); 
instance Field Rational     where (/) = (P./); 
instance Field Exp          where recip = Recip; 
instance Field b => Field (a -> b)  where recip f = \x -> recip (f x); 

instance Ring b => VectorSpace (a -> b) where
    type Under (a -> b) = b
    s £ f = \x -> s * f x



