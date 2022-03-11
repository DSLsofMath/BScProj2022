{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE TypeApplications #-}



import Algebra
import ListVector hiding (v1,v2,m1,m2)

import qualified Prelude
import Prelude hiding ((*),(+),(/),(-))
import Eigen
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

{- 
A system of linear equations can be represented in matrix form using a coefficient matrix, a variable matrix, and a constant matrix
Ax = B

System of equations
 2x + y - z  =  8 
-3x - y + 2z = -11
-2x + y + 2z = -3

-}

-- Column vectors of coefficients
vx, vy, vz, vb :: VecR 3
vx = vec [2,-3,-2] :: VecR 3
vy = vec [1,-1,1]  :: VecR 3
vz = vec [-1,2,2]  :: VecR 3

vb = vec [8,-11,-3] :: VecR 3

-- A matrix with infinite solutions
vx2, vy2, vz2, vb2 :: VecR 3
vx2 = vec [1,2,1] :: VecR 3
vy2 = vec [0,1,1]  :: VecR 3
vz2 = vec [1,2,1]  :: VecR 3

vb2 = vec [1,3,2] :: VecR 3

matx :: MatR 3 4
matx = toMat [vx2,vy2,vz2,vb2]

-- | Should give [1-x,1,x], currently very wrong
solveMatx :: MatR 3 1
solveMatx = toMat [solvesys matx]

-- A matrix
aM :: MatR 3 3
aM =  (mx `append` my) `append` mz
    where
        mx = toMat [vx] :: MatR 3 1
        my = toMat [vy] :: MatR 3 1
        mz = toMat [vz] :: MatR 3 1

-- B matrix
bM :: MatR 3 1
bM = toMat [vb] :: MatR 3 1

{-      A      b
    |2  1 -1 | 8 |
    |3 -1  2 |-11|
    |2  1  2 |-3 |
-}
-- Augmented matrix
augM :: MatR 3 4
augM = aM `append` bM

-- Triangular form
triM :: MatR 3 4
triM = utf augM

solvedSystem :: MatR 3 1
solvedSystem = toMat [solvesys augM]
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