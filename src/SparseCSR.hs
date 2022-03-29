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

-- Test Matrix
test :: CSR Int 4 4
test = CSR {
    elems = [ 5, 8, 3, 6 ],
    col = [ 0, 1, 2, 1 ],
    row = [ 0, 1, 2, 3, 4 ]}

-- Test Matrix
test2 :: CSR Int 4 4
test2 = CSR {
    elems = [ 5, 4, 8, 3, 6 ],
    col = [ 0, 1, 1, 2, 1 ],
    row = [ 0, 2, 3, 4, 5 ]}

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

smv :: (AddGroup f, Mul f) => CSR f a b -> [f] -> [f]
smv (CSR _ _ (r:[])) v = []
smv (CSR [] _ _) v = []
smv (CSR elems col (r:row)) v = dotL (take j elems) (listIndex v (take j col)) : 
                                     smv (CSR (drop j elems) (drop j col) row) v
            where j = head row - r 

testsmv = smv test [1,1,1,1]  == [5,8,3,6]

testsmv2 = smv test2 [1,1,1,1]  == [9,8,3,6]

-- | Take a list uses a list of indexes to print index values
listIndex :: [f] -> [Int] -> [f]
listIndex as [] = []
listIndex as indexes = map (as!!) indexes -- i : listIndex as indexes
