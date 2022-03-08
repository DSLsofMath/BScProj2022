{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
-- {-# LANGUAGE StandaloneKindSignatures #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE NoStarIsType #-}
{-# LANGUAGE UndecidableInstances #-}

module QuadTree where

import Prelude hiding ((+), (-), (*), (/), sum, )
import GHC.TypeLits
import Data.Kind

import ListVector
import Algebra

{- Taken from
 - Matrix Algorithms using Quadtrees Invited Talk, ATABLE-92 Technical Report 357

type Quadrants a = [Matrx a] --list of *four* submatrices.

data Matrx a = Zero | ScalarM a | Mtx (Quadrants a) | ID --used mostly in permutations.
    deriving Show

instance (Num a) => Num (Matrx a) where
    fromInteger 0 = Zero
    fromInteger 1 = ID
    fromInteger _ = error "fromInteger defined only on zero/1 Matrces."
    negate Zero = Zero
    negate (ScalarM x) = ScalarM (negate x)
    negate (Mtx quads) = Mtx (map negate quads)

    x + Zero = x
    Zero + y = y
    ScalarM x + ScalarM y = ScalarM (x + y)    
    Mtx x + Mtx y = case zipWith (+) x y of [Zero,Zero,Zero,Zero] -> Zero
                                            quads     -> Mtx quads
    _ + _ = error "Matrx addition undefined on IdentM"

    Zero * _ = Zero
    _ * Zero = Zero
    ID * y = y --NB: multiplication accepts IdentM
    x * ID = x
    ScalarM x * ScalarM y = ScalarM (x*y)
    --Except with infinitesimal floats: case x*y of zero->zero; z->ScalarM z
    Mtx x * Mtx y = case zipWith (+)
                            (zipWith (*) (colExchange x)(offDiagSqsh y))
                            (zipWith (*) x (prmDiagSqsh y))
                    of 
                        [Zero,Zero,Zero,Zero] -> Zero 
                        quads                 -> Mtx quads 
        where colExchange [nw,ne,sw,se] = [ne,nw,se,sw]
              prmDiagSqsh [nw,ne,sw,se] = [nw,se,nw,se]
              offDiagSqsh [nw,ne,sw,se] = [sw,ne,sw,ne]
              abs _ = error "abs not implemented yet on Matrces."
              signum _ = error "signum not implemented yet on Matrces."

-}

-- Reprsents a matrix size n*n where n = 4^Nat4
data Nat4 = Nil | Succ Nat4

type family Pred (n :: Nat4) :: Nat4 where
     Pred (Succ n)  = n

type family ToNat (n :: Nat4) :: Nat where
     ToNat Nil = 1
     ToNat (Succ n)  = 4 * ToNat n


-- type Quad :: Nat4 -> Type -> Type
data Quad (n :: Nat4) a where 
    Zero :: Quad n a
    Scalar :: a -> Quad Nil a
    Mtx :: {
        nw :: Quad (Pred n) a,
        ne :: Quad (Pred n) a,
        sw :: Quad (Pred n) a, 
        se :: Quad (Pred n) a 
           } -> Quad n a 
    Id :: Quad n a --used mostly in permutations.

instance AddGroup a => AddGroup (Quad n a) where
    Zero + x    = x
    x    + Zero = x
    Scalar a    + Scalar b = Scalar (a+b)
    Mtx nw1 ne1 sw1 se1 + Mtx nw2 ne2 sw2 se2 = Mtx (nw1+nw2) (nw1+nw2) (nw1+nw2) (nw1+nw2)
    a - b = undefined
    zero = Zero 


-- Potentally use parallelism for big enough Quads
-- class Mul a where
-- 
-- instance (n <= 50) => Mul (Quad n a) where
-- 
-- instance (12 <= n) => Mul (Quad n a) where


