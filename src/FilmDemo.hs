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
import ListVector (Vector, vec, dot, cross, Matrix, ToMat (..), toMatT, appendV)
import qualified ListVector as L
import ListVecGauss (particularSol)

import Prelude hiding ((+), (-), (*), (/), (^), (**), sum, product, recip, fromRational, span, elem)
import qualified Prelude as P

-- Starta med DataKinds, TypeApplications och :set +t



-- Vektorer

v3, w3 :: Vector R 3
v3 = vec [1,2,3]
w3 = vec [2,2,3]

-- Stödjer addition, skalär- och kryss-produkt

-- Kan också skapas interaktivt
-- > vec @2 [3,1]

-- Notera att vektorer av olika längd inte kan adderas






-- LIL Matriser

m32 :: Matrix R 3 2
m32 = toMatT [[2,0],
              [1,1],
              [0,3]]

m23 :: Matrix R 2 3
m23 = toMatT [[2,0,0],
              [0,2,1]]

-- Matris vektor multiplication
-- > m23 ** v3

-- Multiplication av matriser med olika storlekar

m33 :: Matrix R 3 3
m33 = m32 ** m23

-- Kan också skapas interaktivt
-- > toMatT @2 @2 [[1,2],
--                [3,4]]


-- identity, zero och typsignaturer
-- > identity :: Matrix R 3 3
-- > zero :: Matrix R 2 6

-- > identity `appendV` v3








-- Delrum 

f1a = nullSpace m33

f1b = range m33

-- Visa igen att m33*x=w3 kan lösas
-- > w3 `elem` range m33

f1c = dim (range m33)







-- Ekvationssystem m33*x = w3

-- Använd Gaussalgoritmen för att omforma matrisen A till trappstegsform.
f2a = gauss m33

-- Visa alla elementära radoperationer använda i uppgift (a). 
f2b = gaussTrace m33

-- Hitta lösningen till ekvationssystemet Ax = u
f2c = solveQ m33 w3

-- Visa att det är lösning









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


