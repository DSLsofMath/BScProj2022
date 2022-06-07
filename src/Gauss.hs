{-# language DataKinds #-}
{-# language FlexibleContexts #-}

module Gauss where

import GHC.TypeLits 
import qualified Prelude
import Algebra 
import Prelude hiding ((+), (-), (*), (/), (^), recip, sum, product, (**), span)
import Data.List 

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

instance Show a => Show (ElimOp n a) where 
    show = showElimOp

-- | Prettier show function for ElimOp a
showElimOp :: Show a => ElimOp n a -> String
showElimOp op = concat $ case op of 
        Swap    (Fin i) (Fin j)   -> [ " ",              row i,              " <-> ", row j ]
        Mul     (Fin i)         s -> [ " ", show s, "*", row i,              " -> ", row i ]
        MulAdd  (Fin i) (Fin j) s -> [ " ", show s, "*", row i, " + ", row j, " -> ", row j ]
    where row i = "R(" ++ show i ++ ")"


-- -- | Shows step by step how a matrix is transformed by a ElimOp trace
-- showElimOnMat :: (Field f, Show f) => [ElimOp m f] -> Matrix f m n -> String
-- showElimOnMat t m0 = let matTrace = scanl (flip elimOpToFunc) m0 t 
--                      in unlines [ show m ++ "\n" ++ show op | (m, op) <- zip matTrace t ]
--                         ++ show (last matTrace)


-- | Representation of an elementary row operation as a matrix 
elimOpToMat :: (KnownNat n, Matrix mat, Ring f) => ElimOp n f -> mat n n f
elimOpToMat e = elimOpToFunc e identity

swap :: (KnownNats m n, Matrix mat, AddGroup f) => 
          Fin m -> Fin m -> mat m n f -> mat m n f
swap i j m = setRow (setRow m i (getRow m j)) j (getRow m i) 

mul :: (KnownNats m n, Matrix mat, Ring f) => 
        Fin m -> f -> mat m n f -> mat m n f
mul i s m = setRow m i (s £ getRow m i) 

mulAdd :: (KnownNats m n, Matrix mat, Ring f) => 
            Fin m -> Fin m -> f -> mat m n f -> mat m n f
mulAdd i j s m = setRow m j (s £ getRow m i + getRow m j)

-- | Representation of an elementary row operation as a function 
elimOpToFunc :: (KnownNats m n, Matrix mat, Ring f) => 
                  ElimOp m f -> (mat m n f -> mat m n f)
elimOpToFunc e = case e of Swap   i j   -> swap i j 
                           Mul    i   s -> mul i s 
                           MulAdd i j s -> mulAdd i j s 

-- | Reduces a trace of elimOps to a single function
foldElimOpsFunc :: (KnownNats m n,Matrix mat, Ring f) => 
                    [ElimOp m f] -> (mat m n f -> mat m n f)
foldElimOpsFunc = foldr (.) id . map elimOpToFunc . reverse

-- | Transforms a given matrix into row echelon form
gauss :: (KnownNats m n, Matrix mat, Field f, Approx f) =>
           mat m n f -> mat m n f
gauss = snd . gauss' 

-- | Returns the required row operations to 
-- transform a given matrix into row echelon form
gaussTrace :: (KnownNats m n, Matrix mat, Field f, Approx f) =>
               mat m n f -> [ElimOp m f]
gaussTrace = fst . gauss' 


gauss' :: (KnownNats m n, Matrix mat, Field f, Approx f) => 
                mat m n f -> ([ElimOp m f], mat m n f)
gauss' m = gauss'' m [] [minBound .. maxBound] [minBound .. maxBound] 

gauss'' :: (KnownNats m n, Matrix mat, Field f, Approx f) => 
            mat m n f -> [ElimOp m f] -> [Fin m] -> [Fin n] -> ([ElimOp m f], mat m n f)
gauss'' m t _      []     = (t, m)
gauss'' m t []     _      = (t, m)
gauss'' m t (i:is) (j:js) = case getCol' j of 
              []           -> gauss'' m t (i:is) js
              (i', a) : rs -> let xs = swap' i' ++ [Mul i (recip a)] ++ map mulAdd' rs 
                              in  gauss'' (foldElimOpsFunc xs m) (t ++ xs) is js
    where getCol' j = sortOn fst [ (i', a) | (i', a) <- getCol m j, a ~/= zero, i' >= i ]
          mulAdd' (j,b) = MulAdd i j (neg (b))
          swap' j = if j /= i then [Swap i j] else []


jordan :: (KnownNat m, Matrix mat, Composable (ElimOp  m f) (mat m n f) (mat m n f), Ord f, Field f, Fractional f) => 
            mat m n f -> [ElimOp m f] -> [Fin m] -> [Fin n] -> ([ElimOp m f], mat m n f)
jordan m t (_) [] = (t, m)
jordan m t [] (_) = (t, m)
jordan m t (i:is) (j:js) = case getCol' j of 
              [] -> jordan m t (i:is) js
              (i', a) : rs -> let xs = (if i' /= i then [Swap i i'] else []) ++ map (mulAdd' i a) rs 
                              in  jordan (foldr (**) m xs) (t ++ xs) is js
    where getCol' = dropWhile ((<i) . fst) . sortOn fst . filter filterZ . getCol m
          mulAdd' i a (j,b) = MulAdd i j (neg (b/a))
          filterZ (_,s) = s > 0.0001 || s < -0.0001

-- PJ: some example(s), please. Or at least a comment pointing the
-- user/reader to where examples can be found.
