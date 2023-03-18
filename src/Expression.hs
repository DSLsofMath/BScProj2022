
module Expression where

import Prelude hiding ((+), (-), (*), (/))
import Data.Kind (Type)
import Data.List (nubBy)
import Data.Maybe (mapMaybe)

import Algebra 


-- | General class for Ring expression
class Expr (exp :: Type -> Type) where
    eval :: Ring a => exp a -> a -> a

    derive :: Ring a => exp a -> exp a

eval' :: (Expr exp, Ring a) => exp a -> a -> a
eval' = eval . derive

--------------------------------------------------------------------------
-- Root solving algorithms 
--
-- Use newton to find the zeros of the characteristic equation = the eigen values
-- 

-- | Returns a root of the expression if found
newton :: (Expr exp, Field a, Approx a) => exp a -> a -> Maybe a
newton exp = converge . take 200 . newtonList exp

-- | Returns, if found, the converging value of a sequence
--   Note: there is probably a better way to check for convergence,
--   but this is sufficient for newtonList.
--   For performance, we could exit early if the sequence begins to "jump".
converge :: Approx a => [a] -> Maybe a 
converge (x:x':xs) | x ~= x'   = Just x'
                   | otherwise = converge (x':xs)
converge _ = Nothing

-- | An infinite sequence approaching a root of the expression
newtonList :: (Expr exp, Field a) => exp a -> a  -> [a]
newtonList exp = iterate (\x -> x - f x / f' x)
    where (f, f') = (eval exp, eval' exp)

-- | Returns a list of found roots based on a list of guesses
roots :: (Expr exp, Field a, Approx a) => exp a -> [a] -> [a]
roots exp = nubBy (~=) . mapMaybe (newton exp) 



