{-# Language DataKinds #-}

import Algebra

import Eigen (detNN, detGauss, newton, roots, showSol, evalMat)
import qualified Eigen as E
import Gauss (ElimOp(..), elimOpToMat, elimOpToFunc, foldElimOpsFunc, gauss, gaussTrace)
import Subspaces hiding (mEx, QuotientSpace(..))


import Matrix hiding (Matrix, pjMat)
import qualified Matrix as M
import QuadTree (QuadM)
import qualified QuadTree as Q
import SparseCSR (CSR)
import qualified SparseCSR as C
import ListVector (Vector, vec, dot, cross, Matrix, ToMat (..), toMatT)
import qualified ListVector as L
import ListVecGauss (particularSol)

import Prelude hiding ((+), (-), (*), (/), (^), (**), sum, product, recip, fromRational, span, elem)
import qualified Prelude as P

-- 1 

a1 :: Matrix Rational 3 3
a1 = toMatT [[ 1, -1,  5],
             [ 2,  0,  7],
             [-3, -5, -3]]

u1 :: Vector Rational 3
u1 = vec [1, 2, -5]


-- (a) Tillhör u1 nollrummet till matrisen a1? 
f1a = u1 `elem` nullSpace a1

-- (b) Tillhör u1 kolumnrummet till matrisen a1? 
f1b = u1 `elem` range a1

-- (c) Vilken rang har a1?
f1c = dim (range a1)




-- 2 

a2 ::  Matrix R 3 3
a2 = toMatT [[1, 0, 1], 
             [2, 1, 2], 
             [1, 1, 1]]

u2 :: (Field a, Num a) => Vector a 3
u2 = vec [1, 2, -5]

-- (a) Använd Gaussalgoritmen för att omforma matrisen A till trappstegsform.
f2a = gauss a2

-- (b) Visa alla elementära radoperationer använda i uppgift (a). 
f2b = gaussTrace a2

-- (c) Hitta lösningen till ekvationssystemet Ax = u
f2c = solveQ a2 u2


-- 3 

a3 :: (Fractional f, Field f) => Matrix f 2 2
a3 = toMatT [[  1, -1 ], 
             [ -2,  0 ]] 

-- (a) Bestäm As egenvärden och egenvektorer.

-- Karakteriska polynomet
eqv = detNN ( toConst a3 - X £ identity )

-- Hittar ekvationens rötter, vi får ut egenvärdena -1 och 2
eqvRoots :: [R]
eqvRoots = roots eqv [-1.0 .. 2]

-- Får ut egenvektorerna
eigen :: Subspace (Vector R 2)
eigen = nullSpace (a3 - (-1) £ identity ) + nullSpace (a3 - 2 £ identity )
          

-- (b) Visa för funna egenvärden x och egenvektorer v att Ax = v


