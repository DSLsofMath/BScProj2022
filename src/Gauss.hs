{-# language DataKinds #-}

module Gauss where

import GHC.TypeLits 
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), (^), recip, sum, product, (**), span)

import Algebra 
import Matrix 


-------------------------------------------
-- Elementary row operations definition 
-- Reduction and equation solver functions

-- | Represents elementary row operations
data ElimOp n a = Swap (Fin n) (Fin n) 
                | Mul (Fin n) a 
                | MulAdd (Fin n) (Fin n) a
              deriving (Eq)

instance Show a => Show (ElimOp n a) where show = showElimOp

-- | Prettier show function for ElimOp a
showElimOp :: Show a => ElimOp n a -> String
showElimOp op = concat $ case op of 
        Swap    (Fin i) (Fin j)   -> [         row i,              " <-> ", row j ]
        Mul     (Fin i)         s -> [ show s, row i,               " -> ", row i ]
        MulAdd  (Fin i) (Fin j) s -> [ show s, row i, " + ", row j, " -> ", row j ]
    where row i = "R(" ++ show i ++ ")"


-- | Shows step by step how a matrix is transformed by a ElimOp trace
-- showElimOnMat :: (Field f, Show f) => [ElimOp m f] -> Matrix f m n -> String
-- showElimOnMat t m0 = let matTrace = scanl (flip elimOpToFunc) m0 t 
--                      in unlines [ show m ++ "\n" ++ show op | (m, op) <- zip matTrace t ]
--                         ++ show (last matTrace)


-- | Representation of an elementary row operation as a matrix 
elimOpToMat :: (KnownNat n, Matrix mat, AddGroup (mat f n n), Ring f) => ElimOp n f -> mat f n n
elimOpToMat e = elimOpToFunc e identity

-- | Representation of an elementary row operation as a function 
elimOpToFunc :: (Matrix mat, Ring f, AddGroup (mat f m n)) => ElimOp m f -> (mat f m n -> mat f m n)
elimOpToFunc e m = case e of Swap   i j   -> setRow (setRow m j (row i)) i (row j)
                             Mul    i   s -> setRow m i (s `scale` row i)
                             MulAdd i j s -> setRow m j (s `scale` row i `merge` row j)
    where row = getRow m 
          scale s l = [ (i, s*a) | (i,a) <- l ]

-- | Reduces a trace of elimOps to a single function
foldElemOpsFunc :: (Matrix mat, Ring f, AddGroup (mat f m n)) => [ElimOp m f] -> (mat f m n -> mat f m n)
foldElemOpsFunc = foldr (.) id . map elimOpToFunc . reverse



