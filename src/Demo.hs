{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE TypeApplications #-}

module Main where


import Algebra
import ListVector hiding (v1,v2,m1,m2)
import qualified Matrix as M

import qualified Prelude
import Prelude hiding ((*),(+),(/),(-),(**))
import Eigen
import Subspaces
import QuadTree
import SparseCSR

--import Test.QuickCheck
-- To run QuickCheck:
-- type ":set -package QuickCheck" in the ghci prompt



-- To try the code in ghci the DataKinds language extension needs to be enabled:
-- type ":set -XDataKinds" in the ghci prompt



{-
Below we initialize some vector variables of different lengths and underlying types and perform
some operations.

A Vector varaible takes to the following parameters: 
- List of elements, ex : [1,2,3]
- Underlying element type, ex: Double, Int, Rational 
- Length of the vector

newtype Vector f (n :: Nat) = V [f] deriving (Show, Eq)

-}

v1 :: Vector Int 4
v1 = vec [1,2,3,4]


v2,v3 :: Vector Double 3
v2 = V [2,7,1]::VecR 3
v3 = V [8,2,8]::VecR 3

-- additon of vectors
addTest = v2 + v3 == vec [2+8,7+2,1+8]

-- Subtraction of vectors
subTest = v2 - v3 == vec [2-8,7-2,1-8]

-- | Dot product of two vectors 
-- test dot product
dotProd = dot v2 v3 == 2*8 + 7*2 + 1*8



-- | Cross product of two vectors of size 3
-- test cross product
{-
A good test of the cross
product is to check that it is orthogonal to both incoming vectors:
dot v (cross v w) == dot w (cross v w) == 0
-}
crossProd = cross v2 v3 == vec [7*8-1*2, 1*8-2*8, 2*2-7*8]


orth = dot v2 (cross v2 v3) == dot v3 (cross v2 v3)
                            && dot v2 (cross v2 v3) == 0


--Functions and example matrixes for performance tests.

-- use :set +s and :set +r when running tests

--Sparse matrixes, identity matrix

testM51 = M.identity :: Matrix Double 5 5
testM5001 = M.identity :: Matrix Double 500 500

testC51 = M.identity :: CSR Double 5 5
testC5001 = M.identity :: CSR Double 500 500

testQ51 =  M.identity :: QuadM Double 5 5
testQ5001 = M.identity :: QuadM Double 500 500

--dense sparse matrixes, tridiagonal matrix

testC5 = denseCSR5
testC500 = denseCSR500

-- to slow, takes long time constructing

testQ5 =  M.toSparse testC5 :: QuadM Double 5 5
testQ500 =  M.toSparse testC500 :: QuadM Double 500 500

-- Performance test functions
-- Designed to avoid GHCi prints when running m + m,
-- will add the equal time to solution which could be
-- deducted manually or ignored as it still gives a 
-- general idea of which is fast. 

performaddtest :: ((AddGroup (mat Double m n)), (Eq (mat Double m n)), 
                    M.Matrix mat) => mat Double m n -> Bool
performaddtest m = a == a
    where a = m + m

performmultest :: ((Mul (mat Double m n)), (Eq (mat Double m n)), 
                    M.Matrix mat) => mat Double m n -> Bool
performmultest m = a == a
    where a = m * m

-- note, add mulgroup for QuadM

testquadperf = a == a
    where a = testQ5001 ** testQ5001 

-- test transpose
m1 :: Matrix Double 3 3
m1 = toMat [[1,1,1],[2,2,2],[3,3,3]]::MatR 3 3
testtran = transpose m1 == (toMat [[1,2,3],[1,2,3],[1,2,3]])

m2 :: Matrix Double 3 1
m2 = toMat[[4,4,4]] :: MatR 3 1

-- append adds m2 as the last col, new matrix = 3 rows, 4 cols
testappend = append m1 m2 == (toMat [[1,1,1],[2,2,2],[3,3,3],[4,4,4]]:: MatR 3 4)

{- 
A system of linear equations can be represented in matrix form using a coefficient matrix, a variable matrix, and a constant matrix
Ax = B

System of equations
 2x + y - z  =  8 
-3x - y + 2z = -11
-2x + y + 2z = -3

-}

-- Column matrices of coefficients
mx, my, mz :: MatR 3 1
mx = toMat [[2,-3,-2]]
my = toMat [[1,-1,1]]
mz = toMat [[-1,2,2]] 


-- A matrix with infinite solutions
vx2, vy2, vz2, vb2 :: VecR 3
vx2 = vec [1,2,1]  :: VecR 3
vy2 = vec [0,1,1]  :: VecR 3
vz2 = vec [1,2,1]  :: VecR 3

vb2 = vec [1,3,2] :: VecR 3

aM2 :: MatR 3 3
aM2 = toMat [vx2,vy2,vz2]

aM3 :: MatR 2 2
aM3 = toMat [[1,0],[0,1]]

matx :: MatR 3 4
matx = toMat [vx2,vy2,vz2,vb2]

smatx = solveQ' matx

-- A matrix
aM :: MatR 3 3
aM =  (mx `append` my) `append` mz

-- B matrix
bM :: MatR 3 1
bM = toMat [[8,-11,-3]]

s1 = solveQ aM (V(head $ unpack bM))
s2 = solveQ' augmentedM
{-      A      b
    |2  1 -1 | 8 |
    |3 -1  2 |-11|
    |2  1  2 |-3 |
-}
-- Augmented matrix
augmentedM :: MatR 3 4
augmentedM = aM `append` bM


--solvedSystem :: MatR 3 1
--solvedSystem = toMat [solvesys augmentedM]
{- 
-- A matrix
aM :: Matrix Double 3 3
aM = toMat [[2,-3,-2],[1,-1,1],[-1,2,2]] :: MatR 3 3

-- B matrix
bM :: Matrix Double 3 1
bM = toMat [[8,-11,-3]] :: MatR 3 1

-- Augmented matrix
augM :: Matrix Double 3 4
augM = append aM bM

-- Triangular form
trMiM = utf augM
-}

main :: IO()
main = putStrLn "" -- $ show solvedSystem
