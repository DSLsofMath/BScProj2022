{-# LANGUAGE BlockArguments #-}
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

--Returns the element of a given coordinate in the sparse matrix
getElem :: AddGroup f => CSR f m n -> (Int,Int) -> f
getElem csr (x,y) = maybe zero id $ lookup x $ uncurry zip $ getRow csr y

getRow :: CSR f m n -> Int -> ([Int], [f])
getRow (CSR elems col row) i = case td i 2 row of
        [a,b] -> unzip $ td a (b-a) (zip col elems) 
        _     -> ([],[])
        where td i j = take j . drop i

getColumn :: CSR f m n -> Int -> ([Int], [f])
getColumn (CSR elems col row) i = undefined 

--getColumnIndex i col elems = unzip $ filter zip col elems
--    \(i',_) -> i' = i

dotL :: (AddGroup f, Mul f) => [f] -> [f] -> f
dotL v1 v2 = sum $ zipWith (*) v1 v2

--
-- Matrix operations
--  

-- Matrix Vector Multiplication                                                                                                                                                                                                                                                                                                                                                                       

-- Multiplies a CSR matrix with a Vector
-- Might be some problems with precise vector definitions 
smV :: (AddGroup f, Mul f) => CSR f a b -> Vector f b  -> Vector f a
smV m (V v) = V (smv m v)

-- Multiplies a CSR matrix with a list
-- getRow could be used instead of calculating difference to take j from elems&col
smv :: (AddGroup f, Mul f) => CSR f a b -> [f] -> [f]
smv (CSR _ _ (r:[])) v = []
smv (CSR [] _ _) v = []
smv (CSR elems col (r:row)) v = dotL (take j elems) (map (v!!) (take j col)) : 
                                     smv (CSR (drop j elems) (drop j col) row) v
            where j = head row - r 

-- Test values/functions

test :: CSR Double 4 4
test = CSR {
    elems = [ 5, 8, 3, 6 ],
    col = [ 0, 1, 2, 1 ],
    row = [ 0, 1, 2, 3, 4 ]}

test2 :: CSR Double 4 4
test2 = CSR {
    elems = [ 5, 4, 8, 3, 6 ],
    col = [ 0, 1, 1, 2, 1 ],
    row = [ 0, 2, 3, 4, 5 ]}

test3 :: CSR Double 2 2
test3 = CSR {
    elems = [ 5, 6 ],
    col = [ 0, 1],
    row = [ 0, 1, 2]}

colVecTest :: CSR Double 4 1
colVecTest = CSR {
    elems = [4,1,3],
    col = [0, 0, 0],
    row = [0, 0, 1, 2, 3]
}

testsmv = smv test [2,7,8,0]  == [10.0,56.0,24.0,42.0]

testsmv2 = smv test2 [1,1,1,1] -- == [9,8,3,6]

v11, v22 :: Vector Double 4
v11 = V [5,3,1,8]::VecR 4
v22 = V [8,2,8,3]::VecR 4

testsmv3 = smV test v11  