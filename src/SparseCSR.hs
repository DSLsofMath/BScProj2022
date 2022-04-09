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

import Data.Function
import Data.List hiding (sum)
import qualified Data.List as L

-- A sparse matrix using the compressed sparse row format
data CSR f (m :: Nat) (n :: Nat) = CSR { elems :: [f],
                            col :: [Int],
                            row :: [Int]}
    deriving (Eq)

deriving instance Show f => Show (CSR f m n)

instance (AddGroup f, Eq f) => AddGroup (CSR f m n) where
    (+) = cSRAdd
    (-) = cSRSub
    neg (CSR e c r) = CSR (map neg e) c r 
    zero = CSR [] [] [0,0] -- row = [replicate (size+1) zero], size from type 

instance (AddGroup f, Mul f, Eq f) => Mul (CSR f m m) where
    (*) = cSRMM
    one = csrIdm

instance (KnownNat m, KnownNat n, AddGroup f, m ~ m', n ~ n') => ToMat m n (CSR f m' n') where
    type Under' (CSR f m' n') = f
    toMat csr = M ( V [ V [ getElem csr (x,y) | y <- [0..] ] | x <- [0..] ]) + zero
    fromMat m = undefined


csrIdm :: CSR f m n
csrIdm = undefined

--Returns the element of a given coordinate in the sparse matrix
getElem :: AddGroup f => CSR f m n -> (Int,Int) -> f
getElem csr (x,y) = maybe zero id $ lookup x $ uncurry zip $ getRow2 csr y

getRow2 :: CSR f m n -> Int -> ([Int], [f])
getRow2 (CSR elems col row) i = case td i 2 row of
        [a,b] -> unzip $ td a (b-a) (zip col elems) 
        _     -> ([],[])
        where td i j = take j . drop i

getElem2 :: AddGroup f => CSR f m n -> (Int,Int) -> f
getElem2 csr (x,y) = snd $ head $ filter ((==x).fst) (getRow csr y)

getRow :: CSR f m n -> Int -> [(Int, f)]
getRow (CSR elems col row) i = case td i 2 row of
        [a,b] -> td a (b-a) (zip col elems) 
        _     -> []
        where td i j = take j . drop i

getColumn :: (AddGroup f) => CSR f m n -> Int -> [(Int, f)]
getColumn (CSR elems col row) i = [ x | (x, y) <- zip cs col, y==i ]
    where
        cs = couple elems row
        couple [] _ = []
        couple e rw@(r:r':rs) = let (xs, ys) = splitAt (r'-r) e in 
                                zip (repeat (length row - length rw)) xs ++ couple ys (r':rs)

dotL :: (AddGroup f, Mul f) => [f] -> [f] -> f
dotL v1 v2 = sum $ zipWith (*) v1 v2

-- dot product between two "csr vectors".
dotCsr :: (AddGroup f, Mul f) => [(Int,f)] -> [(Int,f)]  -> f
dotCsr v [] = zero
dotCsr [] v = zero
dotCsr (v1:vs1) (v2:vs2) | c1 == c2 =  e1 * e2 + dotCsr vs1 vs2
                         | c1 > c2  =  dotCsr (v1:vs1) vs2
                         | c1 < c2  =  dotCsr vs1 (v2:vs2)
        where (c1,c2) = (fst v1, fst v2)
              (e1,e2) = (snd v1, snd v2)

--
-- Matrix operations
--  

-- currently slightly off, transposing answer corrects it
-- as it places column vectors as row vectors in the answer. 
cSRMM :: (AddGroup f, Mul f, Eq f) => CSR f a b -> CSR f b c  -> CSR f a c
cSRMM m1@(CSR e1 c1 r1) 
      m2@(CSR e2 c2 r2) = csrTranspose2 $ foldl comb emptyCSR bs
        where 
            bs = [ filter ((/=zero).snd) (csrMV m1 (getColumn m2 b)) | b <- [0..maximum c1]]
            emptyCSR = CSR [] [] (scanl (+) 0 (map length bs)) :: CSR f a c

cSRSub :: (AddGroup f, Eq f) => CSR f a b -> CSR f a b  -> CSR f a b
cSRSub m1 m2 = cSRAdd m1 (neg m2)

cSRAdd :: (AddGroup f, Eq f) => CSR f a b -> CSR f a b  -> CSR f a b
cSRAdd (CSR e1 c1 [r]) _ = CSR [] [] [r]
cSRAdd m1@(CSR e1 c1 r1) 
       m2@(CSR e2 c2 r2) = foldl comb emptyCSR bs
       where 
            bs = opRows m1 m2 cSRAddRow
            emptyCSR = CSR [] [] (scanl (+) 0 (map length bs)) :: CSR f a b

-- Could move this to be apart of cSRAdd
opRows :: (AddGroup f, Eq f) => CSR f a b -> CSR f a b -> ([(Int,f)] -> [(Int,f)]  -> [(Int,f)]) -> [[(Int, f)]]
opRows m1@(CSR e1 c1 r1) m2 op =  [ filter ((/=zero).snd) (op (getRow m1 a) (getRow m2 a)) | a <- [0..length r1 - 2]]

cSRAddRow :: AddGroup f => [(Int,f)] -> [(Int,f)]  -> [(Int,f)]
cSRAddRow [] [] = []
cSRAddRow [] as = as
cSRAddRow as [] = as
cSRAddRow v1@((c1,e1):as) v2@((c2,e2):bs) | c1 == c2 = (c1,e1+e2) : cSRAddRow as bs
                                          | c1 > c2  = (c2,e2) : cSRAddRow v1 bs
                                          | c1 < c2  = (c1,e1) : cSRAddRow as v2

csrTranspose :: (AddGroup f, Eq f) => CSR f a b -> CSR f a b
csrTranspose m1@(CSR e1 c1 r1) =  foldl comb emptyCSR bs
             where
                bs = [getColumn m1 a | a <- [0..maximum c1]]
                emptyCSR = CSR [] [] (scanl (+) 0 (map length bs))

csrTranspose2 :: (AddGroup f, Eq f) => CSR f a b -> CSR f a b
csrTranspose2 m1@(CSR e1 c1 r1) =  foldl comb emptyCSR $ bs
             where
                bs = map (map snd) $ groupBy ((==) `on` fst) $ sortOn fst $ concat qs
                qs = [ zip (map fst as) (zip (repeat rs) (map snd as))| (rs, as) <- zip [0..] [getRow m1 a | a <- [0..length r1 - 2]] ]
                emptyCSR = CSR [] [] (scanl (+) 0 (map length bs))

-- toList :: CSR a m n -> [(Int, Int a)]
-- toCSR  :: [(Int, Int a)] -> CSR a m n

-- Test values/functions

comb :: AddGroup f => CSR f x y -> [(Int, f)] -> CSR f x y
comb csr [] = csr
comb (CSR e1 c1 r1) as = CSR (e1++es) (c1++cs) r1
    where
        (cs, es) = unzip as                                

-- Matrix Vector Multiplication                                                                                                                                                                                                                                                                                                                                                                       

-- Multiplies a CSR matrix with a Vector 
smV :: (AddGroup f, Mul f) => CSR f a b -> Vector f b  -> Vector f a
smV m (V v) = V (smv m v)

-- Multiplies a CSR matrix with a list vector
-- getRow could be used instead of calculating difference to take j from elems&col
smv :: (AddGroup f, Mul f) => CSR f a b -> [f] -> [f]
smv (CSR _ _ (r:[])) v = []
smv (CSR [] _ _) v = []
smv (CSR elems col (r:row)) v = dotL (take j elems) (map (v!!) (take j col)) : 
                                     smv (CSR (drop j elems) (drop j col) row) v
            where j = head row - r 

-- Multiplies a CSR matrix with a CSR vector
csrMV :: (AddGroup f, Mul f, Eq f) => CSR f a b -> [(Int,f)]  -> [(Int,f)]
csrMV m1@(CSR e1 c1 r1) v1 = filter ((/=zero).snd) [(a,dotCsr (getRow m1 a) v1)| a <- [0..length r1 - 2]]

--- Test values/functions

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
    elems = [1..10000],
    col = [0,1..9999],
    row = [0,1..10000]}

-- Large
denseBigBoi :: CSR Double 10000 10000
denseBigBoi = CSR {
    elems = [1..50000],
    col = concat $ replicate 5 [0,1..9999],
    row = [0,5..50000]}

-- Large in one row
bigBoi2 :: CSR Double 10000 10000
bigBoi2 = CSR {
    elems = [1..10000],
    col = [0,1..9999],
    row = 0 : replicate 10000 10000}

mediumBoi :: CSR Double 1000 1000
mediumBoi = CSR {
    elems   = [1..1000],
    col     = [0..999 ],
    row     = [0..1000]
}

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

rowVecTest :: CSR Double 3 3
rowVecTest = CSR {
    elems = [4,1,9],
    col = [0, 1, 2],
    row = [0, 3, 3, 3]
}

v11, v22 :: Vector Double 4
v11 = V [5,4,0,0]::VecR 4
v22 = V [8,2,8,3]::VecR 4

bigVec :: Vector Double 10000
bigVec = V [1,2..10000]

-- tridiagonal with -2 on diagonal and +1 above and below diagonal
-- TODO make it parameterised on the size n
pjCSR :: CSR Double 10 10
pjCSR = CSR {
    elems = concat stencils,
    col   = concat scols,
    row = scanl (+) 0 slens
  }
  where  n = 10
         stencil = [1,-2,1]
         stencils = [drop 1 stencil] ++ replicate (n-2) stencil ++ [take 2 stencil]
         colfun c = [c-1,c,c+1]
         scols    = [drop 1 (colfun 0)] ++ map colfun [1..n-2]  ++ [take 2 (colfun (n-1))]
         slens = map length stencils

pjMat :: MatR 10 10
pjMat = toMat pjCSR
