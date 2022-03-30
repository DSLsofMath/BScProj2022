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

import Data.List hiding (sum)

-- A sparse matrix using the compressed sparse row format
data CSR f (m :: Nat) (n :: Nat) = CSR { elems :: [f],
                            col :: [Int],
                            row :: [Int]}

deriving instance Show f => Show (CSR f m n)

instance AddGroup f => AddGroup (CSR f m n) where
    (+) = cSRAdd
    (-) = cSRSub
    neg (CSR e c r) = CSR (map neg e) c r 
    zero = CSR [] [] [0,0] -- should be [0,size+1, info from type] 

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

getColumn :: (AddGroup f, Eq f) => CSR f m n -> Int -> ([Int], [f])
getColumn (CSR elems col row) i = loop (CSR elems col row) i (length row - 2)
    where loop m i (-1) = ([],[])
          loop m i r | e == zero = loop m i (r-1) `addLT` ([],[]) 
                     | otherwise = loop m i (r-1) `addLT` ([r],[e])
                            where e = getElem m (i,r)

dotL :: (AddGroup f, Mul f) => [f] -> [f] -> f
dotL v1 v2 = sum $ zipWith (*) v1 v2

--
-- Matrix operations
--  

cSRSub :: AddGroup f => CSR f a b -> CSR f a b  -> CSR f a b
cSRSub m1 (CSR e c r) = cSRAdd  m1 (CSR (map neg e) c r)

cSRAdd :: AddGroup f => CSR f a b -> CSR f a b  -> CSR f a b
cSRAdd (CSR e1 c1 [r]) _ = CSR [] [] [r]
cSRAdd (CSR e1 c1 (r1:rw1)) 
       (CSR e2 c2 (r2:rw2)) = CSR e c [r1] `cSRcombine` cSRAdd (CSR (drop j1 e1) (drop j1 c1) (map (+(length e - j1)) rw1)) (CSR (drop j2 e2) (drop j2 c2) rw2)
       where (c,e) = cSRAddRow (take j1 c1,take j1 e1) (take j2 c2,take j2 e2)
             j1 = head rw1 - r1
             j2 = head rw2 - r2

-- Adds two rows of a csr matrix
cSRAddRow :: AddGroup f => ([Int], [f]) -> ([Int], [f])  -> ([Int], [f])
cSRAddRow (_,[]) (_,[]) = ([],[])
cSRAddRow (_,[]) (cs,es) = (cs,es)
cSRAddRow (cs, es) (_,[]) = (cs, es)
cSRAddRow (c1:cs1,e1:es1) 
           (c2:cs2,e2:es2) | c1 == c2 = ([c1],[e1 + e2]) `addLT` cSRAddRow (cs1,es1) (cs2,es2)
                           | c1 > c2  = ([c2],[e2]) `addLT` cSRAddRow (c1:cs1,e1:es1) (cs2,es2)
                           | c1 < c2  = ([c1],[e1]) `addLT` cSRAddRow (cs1,es1) (c2:cs2,e2:es2)

cSRcombine :: AddGroup f => CSR f a b -> CSR f a b  -> CSR f a b
cSRcombine (CSR e1 c1 r1) (CSR e2 c2 r2) = CSR (e1++e2) (c1++c2) (r1++r2)

addLT :: AddGroup f => ([Int], [f]) -> ([Int], [f])  -> ([Int], [f])
addLT (cs1, es1) (cs2, es2) = (cs1++cs2, es1++es2)                                   

-- Matrix Vector Multiplication                                                                                                                                                                                                                                                                                                                                                                       

-- Multiplies a CSR matrix with a Vector 
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

test1 :: CSR Double 4 4
test1 = CSR {
    elems = [ 2, 7],
    col = [ 3,3 ],
    row = [ 0, 1, 1, 2, 2 ]}

-- Large
bigBoi :: CSR Double 10000 10000
bigBoi = CSR {
    elems = [1,2..10000],
    col = [1,2..10000],
    row = [0,1..10000]}

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

v11, v22 :: Vector Double 4
v11 = V [5,3,1,8]::VecR 4
v22 = V [8,2,8,3]::VecR 4
