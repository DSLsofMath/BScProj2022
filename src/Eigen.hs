{-# LANGUAGE DataKinds #-}

module Eigen where

import GHC.TypeLits
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), (^), recip, sum, product, (**), span)

import Data.List (nubBy)
import Data.Maybe (mapMaybe)
import Data.Function (on)

import Algebra
import ListVector hiding (v1,v2,m1,m2)
import ListVecGauss
import Polynomial
import Subspaces


m44, mtest :: Matrix 4 4 R
m44 = toMat [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]] 
mtest = toMatT [[-5,4,1,7],[-9,3,2,-5],[-2,0,-1,1],[1,14,0,3]] 

--------------------------------------------------------------------------
-- Determinants of matrices 
--

-- | Determinant for any nxn matrix 
--   The algorithm is based on Laplace expansion
--   Note that the complexity is exponantial since it calls itself n times per iteration
detNN :: Ring f => Matrix n n f -> f
detNN (M (V [V [x]])) = x
detNN m = sum $ zipWith (*) (cycle [one, neg one]) $ do 
    (V (s:_), subM) <- separateCols m      -- All combinations of columns and their remaining matrix
    let (_, m') = head $ separateRows subM -- Removes the upper row, creating a square matrix
    return $ s * detNN m'


-- | Uses Guassian elimination to compute the determinant
--   This algorithm is only O(n^3)
detGauss :: (Field f, Eq f) => Matrix n n f -> f
detGauss m = let diag = getDiagonal gaussM in product diag / traceProduct
    where trace = gaussTrace m 
          gaussM = foldElemOpsFunc trace m
          traceProduct = product [ s | Mul _ s <- trace ]


--tests, all give correct results
m22 :: Matrix 2 2 R
m22   = toMat [[10,1],[1,10]] -- detNN = 99

m33 :: Matrix 3 3 R
m33   = toMat [[1,-3,2],[3,-1,3],[2,-3,1]] -- detNN = -15


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


--------------------------------------------------------------------------
-- Eigen-related functions
--
-- Note: In general our root finder, newton, struggles with roots  
-- located on stationary points. To compensate for this we explicitly test 
-- if 0 and +- 1 are eigenvalues, as these are common. A better solution 
-- might be to check all stationary points. That however is costly and we
-- might still encounter stationary points. 
-- As a side note it would also be nice to not give instal guesses 
-- when finding eigenvalues/vectors


-- | The characteristic polynomial
--   The roots of this polynomial is the eigenvalues of the matrix
--   It also has the property that:  
--   prop_charPoly m = evalOver (charPoly m) m == zero
charPoly :: (KnownNat n, Ring f) => Matrix n n f -> Poly f
charPoly m = detNN (m' - polyX £ idm) 
    where m' = fmap (\x -> P [x]) m

-- | Computes the eigenvalues of a matrix
eigenValues :: (KnownNat n, Field f, Approx f) => Matrix n n f -> [f] -> [f]
eigenValues m = roots (charPoly m) 

-- | Computes the eigenvectors of a matrix
eigenVectors :: (KnownNat n, Field f, Approx f) => Matrix n n f -> [f] -> [Vector n f]
eigenVectors m = map snd . eigenValVec m

-- | Computes the eigen-pairs of a matrix
eigenValVec :: (KnownNat n, Field f, Approx f) => Matrix n n f -> [f] -> [(f, Vector n f)]
eigenValVec m = removeDuplicates . concatMap getVecs . (guesses++) . eigenValues m 
    where getVecs x = (\v -> (x, v)) <$> unSub ( nullSpace' (m - x £ idm) )
          removeDuplicates = nubBy ((~=) `on` snd)
          guesses = [zero, one, neg one] -- Gives nicer results on common eigenvalues, 
                                         -- see Note above



-- Eigenvalues = 1, 1/2 
matrix1 :: Matrix 2 2 R
matrix1 = toMat [[3/4, 1/4], 
                 [1/4, 3/4]] 

-- Eigenpairs (1, [1,1]), (1/2, [-1, 1])
eigenM1 :: [(R, Vector 2 R)]
eigenM1 = eigenValVec matrix1 [-3, -2.5 .. 3.0]

--
-- Eigen values = -2, 2, 0
matrix2 :: Matrix 3 3 R
matrix2 = toMat[[1, 2, 1], [1 , 0, -1],[-1, -2, -1]] 

-- Eigenpairs (0, [1,0,1]), (-2, [0,1,1]), (2 [1,1,0])
eigenM2 :: [(R, Vector 3 R)]
eigenM2 = eigenValVec matrix2 [-3, -2.5 .. 3.0]


-- Eigenvectors are solutions to (A − λI)x = 0 for eigen values
-- 

evalMat :: Matrix m n (Exp f) -> f -> Matrix m n f
evalMat m val = fmap (\x -> evalExp x val) m



-- | Idea for visualizing the answers of systems of equations as strings
--   Try running showSol on an row echelon form matrix (gauss matrix)

leadingZeros :: (Eq f, Num f) => [f] -> Int
leadingZeros = length . takeWhile (== 0)

showVariableValues :: (Field f, Num f, Ord f, Show f) => [f] -> [String] -> String
showVariableValues r var_names
  | not (null other_coefficients) = var_str ++ other_vars_str
  | otherwise = var_str
  where
    noSol = "This system has no solution"
    index = leadingZeros r
    coefficient = (r !! index)                       
    value = last r
    raw_row = reverse . drop 1 . reverse $ r -- row coefficients, except the free member
    elements_count = length raw_row
    other_coefficients = filter (\(k, k_idx) -> k /= 0 && k_idx /= index) (zip raw_row [0 .. elements_count])
    subtract_coefficient k = if k < 0 then " + " ++ show (- k) else " - " ++ show k
    other_vars_str = concatMap (\(k, k_idx) -> subtract_coefficient k ++ " * " ++ (var_names !! k_idx)) other_coefficients
    var_str = if(index < length var_names) then (var_names !! index) ++ " = " ++ show (value / coefficient)
                                        else error noSol
        
         




variables :: [String]
variables = [ x ++ show(i) | x <- ["x"] , i <- [1..] ] -- ["x1", "x2", "x3" ...]

showColnRow :: (Field f, Num f, Ord f, Show f) => [[f]] -> [String] -> String
showColnRow [] _        = []
showColnRow (x:xs) vars | skipRow == length x = showColnRow xs vars
                        | otherwise           = showVariableValues x vars ++ "\n" ++ showColnRow xs vars
    where
        skipRow = leadingZeros x


-- showSol $ gauss eigen1
-- showSOl $ gauss eigen05

-- if no solution exists for the system, error message thrown
showSol :: (Field f, Num f, Ord f, Show f, KnownNats m n) => Matrix m n f -> IO()
showSol m = putStr $ showColnRow (unpack $ transpose m) vars
    where
        rows     = unpack $ transpose m
        nrOfVars = length(unpack m) - 1
        vars     = take (nrOfVars) variables


