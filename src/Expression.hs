
module Expression where


import Data.Kind (Type)

import Algebra (Ring)


-- | General class for Ring expression
class Expr (exp :: Type -> Type) where
    eval :: Ring a => exp a -> a -> a

    derive :: Ring a => exp a -> exp a

eval' :: (Expr exp, Ring a) => exp a -> a -> a
eval' = eval . derive

