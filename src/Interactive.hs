module Interactive where

import Algebra 

import Eigen (detNN, detGauss, newton, roots, showSol, evalMat)
import qualified Eigen as E
import Gauss (ElimOp(..), elimOpToMat, elimOpToFunc, foldElimOpsFunc, utf, utfTrace)
import Subspaces hiding (mEx, QuotientSpace(..), Subspace(..))


import Matrix hiding (Matrix, pjMat)
import qualified Matrix as M
import QuadTree (QuadM)
import qualified QuadTree as Q
import SparseCSR (CSR)
import qualified SparseCSR as C
import ListVector (Vector, vec, dot, cross, Matrix, ToMat, append, )
import ListVector as L
import ListVecGauss (particularSol)

import Prelude hiding ((+), (-), (*), (/), (^), (**), sum, product, recip, fromRational, span, elem)
import qualified Prelude as P
