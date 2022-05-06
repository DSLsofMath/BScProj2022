{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE TypeApplications #-}


module Eigen where

import GHC.TypeLits
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), (^), recip, sum, product, (**), span)

import qualified Data.List as L
import Algebra
import ListVector hiding (v1,v2,m1,m2)
import ListVecGauss

-- This file contains an example of using TypeLits to handle vectors with a given size. 
-- The implementation here is based on lists and should be replaced.
--
-- To try the code in ghci the DataKinds language extension needs to be enabled:
-- type ":set -XDataKinds" in the ghci prompt
--

m44 = toMat [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]] :: MatR 4 4
mtest = transpose $ toMat [[-5,4,1,7],[-9,3,2,-5],[-2,0,-1,1],[1,14,0,3]] :: MatR 4 4

-- | Determinant for any nxn matrix 
--   The algorithm is based on Laplace expansion
--   Note that the complexity is exponantial since it calls itself n times per iteration
detNN :: Field f => Matrix f n n -> f
detNN (M (V [V [x]])) = x
detNN m = sum $ zipWith (*) (cycle [one, neg one]) $ do 
    (V (s:_), subM) <- separateCols m      -- All combinations of columns and their remaining matrix
    let (_, m') = head $ separateRows subM -- Removes the upper row, creating a square matrix
    return $ s * detNN m'


-- | Uses Guassian elimination to compute the determinant
--   This algorithm is only O(n^3)
detGauss :: (Field f, Eq f) => Matrix f n n -> f
detGauss m = let V diag = getDiagonal utfM in product diag / traceProduct
    where trace = utfTrace m 
          utfM = foldElemOpsFunc trace m
          traceProduct = product [ s | Mul _ s <- trace ]


--tests, all give correct results
idm33 = idm :: MatR 3 3 -- det33 = 1
idm22 = idm :: MatR 2 2 -- det 22 = 1

m22   = toMat [[10,1],[1,10]] :: MatR 2 2 -- det22 = 99
m33   = toMat [[1,-3,2],[3,-1,3],[2,-3,1]] :: MatR 3 3 -- det33 = -15


exp22 :: Matrix Exp 2 2
exp22 = toMat [[X:*:Const 1,Const 2],[Const 2, zero]] --- (X £ idm :: Matrix Exp 2 2)


-- Use newton to find the zeros of the characteristic equation = the eigen values
newton :: (Field a, Num a, Eq a, Ord a) => Exp -> a -> a -> a
newton f eps x = if abs fx < eps then x
                              else if fx' /= 0 then newton f eps next
                                                       else  newton f eps (x+eps)
      where fx  = evalExp f x
            fx' = evalExp' f x
            next = x - (fx/fx')

newtonLast :: (Field a, Num a, Eq a, Ord a) => Exp -> a -> a -> a
newtonLast f eps x = last (newtonList f eps x)

newtonList :: (Field a, Num a, Eq a, Ord a) => Exp -> a -> a -> [a]
newtonList f eps x = x : if abs fx < eps then [ ]
                                  else if fx' /= 0 then newtonList f eps next
                                                           else newtonList f eps (x+eps)
      where fx  = evalExp f x
            fx' = evalExp' f x
            next = x - (fx/fx')

-- Todo, let user choose span
roots :: (Field a, Num a, Eq a, Ord a, Enum a, Fractional a) => Exp -> [a]
roots f = map (newton f 0.001) [0.0,0.5..2.0]

-- Eigen values = 1, 1/2 
matrix1 :: Matrix Exp 2 2
matrix1 = toMat[[Const(3/4), Const(1/4)], [Const(1/4), Const(3/4)]] - X £ idm

-- Using detNN, seems to get some error with detGauss, probably because of derive recip
eigenM1 :: (Field a, Num a, Eq a, Ord a, Enum a, Fractional a) => [a]
eigenM1 = roots(detNN matrix1)

-- Eigen values = -2, 2, 0
matrix2 :: Matrix Exp 3 3
matrix2 = toMat[[one, Const 2, one], [one , zero, neg one],[neg one, neg Const 2, neg one]] - X £ idm

eigenM2 :: (Field a, Num a, Eq a, Ord a, Enum a, Fractional a) => [a]
eigenM2 = roots(detNN matrix2)


-- Eigenvectors are solutions to (A − λI)x = 0 for eigen values
-- 

evalCol :: (Field a, Num a, Eq a, Ord a, Enum a, Fractional a) => [Exp] -> a -> [a]
evalCol [] _ = []
evalCol (x:xs) val = [evalExp x val] ++ evalCol xs val

evalColnRow :: (Field a, Num a, Eq a, Ord a, Enum a, Fractional a) => [[Exp]] -> a -> [[a]]
evalColnRow [] _       = []
evalColnRow (x:xs) val = evalCol x val : evalColnRow xs val

evalMat :: (Field f, Num f, Eq f, Ord f, Enum f, Fractional f) => Matrix Exp m n -> f -> Matrix f m n
evalMat m val = pack $ evalColnRow (unpack m) val

{-
For matrix1 => Eigenvalues = 0.5 & 1

Gauss $ evalMat matrix 0.5 `append` zeroVec => Solution = eigenvector

Gauss $ evalMat matrix 1 `append` zeroVec => Solution = eigenvector

-}

eigen05 :: Matrix Double 2 3
eigen05 = (evalMat matrix1 0.5) `append` toMat [[0,0,0]] -- utf eigen05 => x = -y 

eigen1 :: Matrix Double 2 3
eigen1 = (evalMat matrix1 1) `append` toMat [[0,0,0]]    -- utf eigen1  => x = y

----------------

two :: Ring a => a
two = one+one

{-
pjM :: Ring f => Matrix f 2 2
pjM = M (Vector (L [Vector (L [one,two]),
                    Vector (L [two,one])]))
      

pjExp :: Exp
pjExp = detNN (pjM - scaleM X idm)

pjFun :: Field f => f -> f
pjFun x = detNN (pjM - scaleM x idm)
-- zeroes in (-1) and 3

-- TODO pjPoly - using the polynomial instance from DSLsofMath
-}




-- | Idea for visualizing the answers of systems of equations as strings
--   Try running showSol on an row echelon form matrix (utf matrix)

leadingZeros :: (Field f, Eq f, Num f) => [f] -> Int
leadingZeros = length . takeWhile (== 0)

showVariableValues :: (Field f, Eq f, Num f, Ord f, Show f) => [f] -> [String] -> String
showVariableValues r var_names
  | not (null other_coefficients) = var_str ++ other_vars_str
  | otherwise = var_str
  where
    index = leadingZeros r
    coefficient = r !! index
    value = last r
    raw_row = reverse . drop 1 . reverse $ r -- row coefficients, except the free member
    elements_count = length raw_row
    other_coefficients = filter (\(k, k_idx) -> k /= 0 && k_idx /= index) (zip raw_row [0 .. elements_count])
    subtract_coefficient k = if k < 0 then " + " ++ show (- k) else " - " ++ show k
    other_vars_str = concatMap (\(k, k_idx) -> subtract_coefficient k ++ " * " ++ (var_names !! k_idx)) other_coefficients
    var_str = (var_names !! index) ++ " = " ++ show (value / coefficient)




variables :: [String]
variables = [ x ++ show(i) | x <- ["x"] , i <- [1..] ] -- ["x1", "x2", "x3" ...]

showColnRow :: (Field f, Eq f, Num f, Ord f, Show f) => [[f]] -> [String] -> String
showColnRow [] _        = []
showColnRow (x:xs) vars | skipRow == length x = showColnRow xs vars
                        | otherwise           = showVariableValues x vars ++ "\n" ++ showColnRow xs vars
    where
        skipRow = leadingZeros x


-- showSol $ utf eigen1
-- showSOl $ utf eigen05

showSol :: (Field f, Eq f, Num f, Ord f, Show f, KnownNat m, KnownNat n) => Matrix f m n -> IO()
showSol m = putStr $ showColnRow (unpack $ transpose m) vars
    where
        rows     = unpack $ transpose m
        nrOfVars = length(unpack m) - 1
        vars     = take (nrOfVars) variables 