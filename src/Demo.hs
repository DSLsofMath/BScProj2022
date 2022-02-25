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

crossProd = cross v2 v3 == vec [7*8-1*2, 1*8-2*8, 2*2-7*8]



{- 
Ax = B

System of equations
 2x + y - z  =  8 
-3x - y + 2z = -11
-2x + y + 2z = -3

-}

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
triM = utf augM
