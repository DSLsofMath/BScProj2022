{-# LANGUAGE KindSignatures #-}

module Expression where

import Algebra (Ring)


-- | General class for Ring expression
class Expr (exp :: * -> *) where
    eval :: Ring a => exp a -> a -> a

    derive :: Ring a => exp a -> exp a

eval' :: (Expr exp, Ring a) => exp a -> a -> a
eval' = eval . derive

