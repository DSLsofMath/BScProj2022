{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FunctionalDependencies #-}

-- Old Gauss implementation, hard coded for ListVector
-- See Gauss for general version
module ListVecGauss where

import GHC.TypeLits
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), (^), recip, sum, product, (**), span)

import qualified Data.List as L

import Algebra
import ListVector
import FiniteIndex
import qualified Gauss as G



-------------------------------------------
-- Elementary row operations definition 
-- Reduction and equation solver functions

type Index = Int

-- | Represents elementary row operations
data ElimOp a = Swap Index Index 
              | Mul Index a 
              | MulAdd Index Index a
              deriving (Eq, Show)

-- instance Show a => Show (ElimOp a) where show = showElimOp

-- | Prettier show function for ElimOp a
showElimOp :: Show a => ElimOp a -> String
showElimOp op = concat $ case op of 
                  Swap    i j   -> [         row i,              " <-> ", row j ]
                  Mul     i   s -> [ show s, row i,               " -> ", row i ]
                  MulAdd  i j s -> [ show s, row i, " + ", row j, " -> ", row j ]
    where row i = "R(" ++ show i ++ ")"


-- | Shows step by step how a matrix is transformed by a ElimOp trace
showElimOnMat :: (Field f, Show f) => [ElimOp f] -> Matrix m n f -> String
showElimOnMat t m0 = let matTrace = scanl (flip elimOpToFunc) m0 t 
                     in unlines [ show m ++ "\n" ++ show op | (m, op) <- zip matTrace t ]
                        ++ show (last matTrace)

-- | Elementary row functions.
--   It might be better to define the as (Vector f n -> Vector f n) e.g. a linear map. 
--   If we want to apply it to a matrix we can then use onCols. 
--   Doing so will also remove the need for transpose.
swap :: Index -> Index -> Matrix m n f -> Matrix m n f
swap i j = onCols $ \(V v) -> 
    let (i', j') = (min i j - 1, max i j - 1)
        (m1,x:xs) = splitAt i' v
        (xs',y:m2) = splitAt (j' - i'-1) xs
    in V $ m1++y:xs'++x:m2

mul :: Mul f => Index -> f -> Matrix m n f -> Matrix m n f
mul i s = onCols $ \(V v) -> 
    let (m1,x:m2) = splitAt (i-1) v in V $ m1++(s*x): m2

muladd :: Ring f => Index -> Index -> f -> Matrix m n f -> Matrix m n f
muladd i j s = onCols $ \(V v) -> 
    let (i', j') = (i - 1, j - 1)
        (_,x:_) = splitAt i' v
        (m1,y:m2) = splitAt j' v
        y' = (s*x) + y
    in V $ m1++y':m2

instance (Ring f, f ~ f') => Composable (ElimOp f') (Matrix m n f) (Matrix m n f) where
    (Swap i j) ** m = swap i j m
    (Mul i s) ** m = mul i s m
    (MulAdd i j s) ** m = muladd i j s m


instance (Ring f, n ~ n', f ~ f') => Composable (G.ElimOp n' f') (Matrix m n f) (Matrix m n f) where
    (G.Swap (Fin i) (Fin j)) ** m = swap i j m
    (G.Mul s (Fin i)) ** m = mul i s m
    (G.MulAdd s (Fin i) (Fin j)) ** m = muladd i j s m


-- | Representation of an elementary row operation as a matrix 
elimOpToMat :: (KnownNat n, Ring f) => ElimOp f -> Matrix n n f
elimOpToMat e = elimOpToFunc e idm

foldElemOps :: (KnownNat n, Ring f) => [ElimOp f] -> Matrix n n f
foldElemOps = product . map elimOpToMat . reverse

-- | Representation of an elementary row operation as a function 
elimOpToFunc :: Ring f => ElimOp f -> (Matrix m n f -> Matrix m n f)
elimOpToFunc (Swap   i j  ) = swap   i j
elimOpToFunc (Mul    i   s) = mul    i   s
elimOpToFunc (MulAdd i j s) = muladd i j s
                         
-- | Reduces a trace of elimOps to a single function
--   TODO: We should add a rewrite rule such that fold elemOpToFunc only packs and unpacks once
foldElemOpsFunc :: Ring f => [ElimOp f] -> (Matrix m n f -> Matrix m n f)
foldElemOpsFunc = foldr (.) id . map elimOpToFunc . reverse


-- | Applies a function on a unpacked and transposed matrix before transposing it back
--   UNSAFE: The function can not change the dimension of the matrix 
onUnpackedTrans :: ([[f]] -> [[f]]) -> Matrix m n f -> Matrix m n f
onUnpackedTrans f = pack . L.transpose . f . L.transpose . unpack

-- | Transform a matrix to upper echelon form
gauss :: (Eq f, Field f) => Matrix m n f -> Matrix m n f
gauss = onUnpackedTrans (sort . f) 
    where
          f []     = []
          f (x:[]) = let (x', n) = pivot x in [x']
          f (x:xs) = let (x', n) = pivot x in x' : (f $ map (reduce n x') xs)
          pivot [] = ([],-1)
          pivot (x:xs) | x == zero = let (ys, n) = pivot xs in (x:ys, n + 1)
                       | otherwise = (map (/x) (x:xs), 0 :: Int)
          reduce n p x = zipWith (-) x (map ((x!!n)*) p)
          sort = L.sortOn (length . takeWhile (==zero))

-- | Generate a trace of ElimOps from reducing a matrix to upper echelon form
gaussTrace :: (Field f, Eq f) => Matrix m n f -> [ElimOp f]
gaussTrace m0 = case separateCols m0 of
      []               -> []
      (V (x:xs), m):_  -> let 
                      trace  = mul x ++ [ mulAdd s j | (s,j) <- zip xs [2..], s /= zero ]
                      augM = foldElemOpsFunc trace m
                      in trace ++ case (m, separateRows augM) of
                         (M (V (_:_)), (_,m'):_ ) -> map incIndex $ gaussTrace m'
                         _                        -> []
    where
        mulAdd s j = MulAdd 1 j (neg s)
        mul s = if s /= one then [Mul 1 (recip s)] else []
        incIndex :: ElimOp f -> ElimOp f
        incIndex (Swap i j)     = Swap   (i+1) (j+1)
        incIndex (Mul i s)      = Mul    (i+1)       s
        incIndex (MulAdd i j s) = MulAdd (i+1) (j+1) s


-- | Solve systems of equations
--   Check in Demo.hs how to use
solve :: (Eq f, Field f) => [[f]] -> [f]
solve m = foldr next [last (last m)] (init m)
    where
        next row found = let
            subpart = init $ drop (length m - length found) row
            solved = last row - sum (zipWith (*) found subpart)
            in solved : found

-- | apply when solving systems of equations
--   each element in list represents variable values
--   Recommend using showSol from eigen to verify result, or v `elem` range m
particularSol :: (Eq f, Field f) => Matrix m n f -> Vector (n -1) f
particularSol = V . solve . unpack . transpose . gauss
