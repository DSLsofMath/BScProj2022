{-# LANGUAGE GADTs #-}
{-# LANGUAGE StandaloneDeriving #-}


module Expression.Exp where

import Prelude hiding ((+), (-), (*), (/), (^), sum, product, recip, fromRational)

import Algebra
import Expression

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
evalExp' = evalExp . deriveExp

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


instance AddGroup a => AddGroup (Exp a) where 
    (+) = (:+:)
    neg = Negate
    zero = Const zero

instance Mul a => Mul (Exp a) where 
    (*) = (:*:)
    one = Const one

instance Field a => Field (Exp a) where 
    recip = Recip 


