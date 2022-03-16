{-# Language FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ConstraintKinds #-}



module Algebra where 

import GHC.TypeLits hiding ( type (^) )
import qualified Prelude as P
import Prelude hiding ((+), (-), (*), (/), sum, product, recip, fromRational)
import qualified Data.Ratio



type R = Double


-- | A list with a given length
newtype List a (n :: Nat) = L [a] deriving (Show, Eq)


data Exp =  Const R
             |  X
             |  Exp :+: Exp
             |  Exp :*: Exp
             |  Recip Exp
             |  Negate Exp
    deriving (Eq, Ord)


-- Remove + const 0 from Exp matrices
instance Show Exp where
   show s = replaced (showE s) zerosToReplace  


-- Remove all substrings in zerosToReplace
replaceS :: Eq a => [a] -> [a] -> [a] -> [a]
replaceS [] _ _ = []
replaceS s find repl =
    if take (length find) s == find
        then repl ++ (replaceS (drop (length find) s) find repl)
        else [head s] ++ (replaceS (tail s) find repl)

zerosToReplace :: [String]
zerosToReplace = [ " + (-0.0)", " + 0.0", " - (-0.0)", " - 0.0", "0.0 + ", "(-0,0) + ","0.0 - ","(-0.0) - " ]

replaced :: String -> [String] -> String
replaced s []      = s
replaced s (x:xs)  = replaced sub xs
    where
        sub = replaceS s x ""



showE :: Exp -> String
showE X = "X"
showE (e :+: Const 0) = showE e
showE (Const 0 :+: e) = showE e
showE (e1 :+: e2) = "(" ++ showE e1 ++ " + " ++ showE e2 ++ ")"
showE (e1 :*: e2) = "(" ++ showE e1 ++ " * " ++ showE e2 ++ ")"
showE (Recip e)   = "(" ++ "Recip" ++ showE e ++ ")"
showE (Negate e)  = "(" ++ "-" ++ showE e ++ ")"
showE (Const r)   = show r


-- Eval functions for expressions
fromRational :: (Field a, Num a) => Data.Ratio.Ratio Integer -> a
fromRational x  =  fromInteger (Data.Ratio.numerator x) / fromInteger (Data.Ratio.denominator x)


-- Vad är det som är fel med *, kör ^, inte * ?
evalExp :: (AddGroup a, Mul a, Field a, Num a) => Exp -> a -> a
evalExp (Const alpha)  =  const (fromRational (toRational alpha))
evalExp X              =  id
evalExp (e1 :+: e2)    =  evalExp e1 + evalExp e2
evalExp (e1 :*: e2)    =  evalExp e1 * evalExp e2
evalExp (Negate e)     =  neg (evalExp e)
evalExp (Recip e)      =  recip (evalExp e)

{- 
simplifyE :: Exp -> Exp
simplifyE X               = X
simplifyE (Const x)       = Const x
simplifyE (e :+: Const 0) = simplifyE e
simplifyE (Const 0 :+: e) = simplifyE e
simplifyE (e1 :+: e2)     = simplifyE e1 :+: simplifyE e2
simplifyE (e1 :*: e2)     = e1 :*: e2
simplifyE (Negate e)      = Negate e
simplifyE (Recip e)       = Recip e

-}


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

type Ring a = (AddGroup a, Field a)

class (AddGroup a, Mul a) => Field a where
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
--
class VectorSpace x => Finite x where 
    type Dim x :: Nat
    type BasisVec x :: *
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

instance VectorSpace b => VectorSpace (a -> b) where
    type Under (a -> b) = Under b
    s £ f = \x -> s £ f x



