{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE TypeFamilies #-}

module Polynomial where

import Prelude hiding ((+), (-), (*), (/), (^), (**), sum, product, recip, fromRational)
import Algebra


-- | Poly is a polynomial representation 
--   P [1,2,3,0,5] = 1 + 2x + 3x^2 + 5x^4
newtype Poly a = P [a] deriving (Functor, Eq)


showP :: (Show a, AddGroup a, Eq a) => Poly a -> String
showP (P xs) = case [ format x i | (x, i) <- zip xs [0 ..], x /= zero ] of
        []   -> "0"
        x:xs -> foldl (\s t -> s ++ " + " ++ t) x xs
    where format x 0 = show x 
          format x 1 = show x ++ "x"
          format x i = show x ++ "x^" ++ show i

instance (Show a, AddGroup a, Eq a) => Show (Poly a) where
    show = showP

-----------------------------------------------------------
-- Eval functions on Poly

evalP :: Ring a => Poly a -> a -> a 
evalP (P cs) x = sum $ zipWith (*) cs exps 
    where exps = iterate (* x) one

deriveP :: Ring a => Poly a -> Poly a
deriveP (P xs) = P . drop 1 $ zipWith (*) xs powers
    where powers = iterate (+ one) zero

instance Expr Poly where
    eval = evalP
    derive = deriveP

-----------------------------------------------------------
-- Operations and values on Poly

-- | polyX = x^1 
polyX :: Ring a => Poly a
polyX = P [zero, one]

addP :: AddGroup a => Poly a -> Poly a -> Poly a
P as `addP` P bs = P $ as `add` bs
    where []     `add` b      = b
          a      `add` []     = a
          (a:as) `add` (b:bs) = a + b : as `add` bs

mulP :: Ring a => Poly a -> Poly a -> Poly a
P as `mulP` bs = sumP $ map (£ bs) as
    where sumP = foldr (\p acc -> p + raise acc) zero 
          raise (P xs) = P (zero : xs)


instance AddGroup a => AddGroup (Poly a) where
    (+) = addP
    neg = fmap neg
    zero = P []

instance Ring a => VectorSpace (Poly a) where
    type Under (Poly a) = a
    s £ p = (s *) <$> p

instance Ring a => Mul (Poly a) where
    (*) = mulP 
    one = P [one]


