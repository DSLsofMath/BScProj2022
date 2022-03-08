{-# LANGUAGE DataKinds #-} 
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeOperators #-} 
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeFamilies #-}


module SparseCSR where

import GHC.TypeLits
import Algebra
import Prelude hiding ((+), (-), (*), (/), sum)
import ListVector 

-- A sparse matrix using the compressed sparse row format
data CSR f (m :: Nat) (n :: Nat) = CSR { elems :: [f],
                            col :: [Int],
                            row :: [Int]}

deriving instance Show f => Show (CSR f m n)

instance (KnownNat m, KnownNat n, AddGroup f) => ToMat m n (CSR f m n) where
    type Under' (CSR f m n) = f
    toMat csr = M ( V [ V [ getElem csr (x,y) | y <- [0..] ] | x <- [0..] ]) + zero

-- Test Matrix
test :: CSR Int 4 4
test = CSR {
    elems = [ 5, 8, 3, 6 ],
    col = [ 0, 1, 2, 1 ],
    row = [ 0, 1, 2, 3, 4 ]}

colVecTest :: CSR Int 4 1
colVecTest = CSR {
    elems = [4,1,3],
    col = [0, 0, 0],
    row = [0, 0, 1, 2, 3]
}


-- Matrix ?

--Returns the element of a given coordinate in the sparse matrix
getElem :: AddGroup f => CSR f m n -> (Int,Int) -> f
getElem csr (x,y) = maybe zero id $ lookup x $ uncurry zip $ getRow csr y

getRow :: CSR f m n -> Int -> ([Int], [f])
getRow (CSR elems col row) i = case td i 2 row of
        [a,b] -> unzip $ td a (b-a) (zip col elems) 
        _     -> ([],[])
    where td i j = take j . drop i

--
-- Matrix operations
--

-- Matrix Vector Multiplication                                                                                                                                                                                                                                                                                                                                                                       
