{-# LANGUAGE BlockArguments #-}
{-# LANGUAGE DataKinds #-} 
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeOperators #-} 
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE UndecidableInstances #-}



module SparseCSR where

import GHC.TypeLits
import Algebra
import Prelude hiding ((+), (-), (*), (/), sum)
import ListVector 
import qualified Matrix as M
import Matrix (values, Fin(..))

import Data.Function
import Data.List hiding (sum)
import qualified Data.List as L

-- A sparse matrix using the compressed sparse row format
data CSR f (m :: Nat) (n :: Nat) = CSR { elems :: [f],
                            col :: [Int],
                            row :: [Int]}

deriving instance Show f => Show (CSR f m n)

instance (Eq f, AddGroup f, AddGroup (CSR f m n)) => Eq (CSR f m n) where
    m1 == m2 = and bs
        where
            (m1', m2') = (M.purge m1, M.purge m2)
            bs = [
                elems m1' == elems m2',
                col m1' == col m2',
                row m1' == row m2']

instance forall m n f. (KnownNat n, AddGroup f) => AddGroup (CSR f m n) where
    (+) = cSRAdd
    (-) = cSRSub
    neg (CSR e c r) = CSR (map neg e) c r 
    zero = CSR [] [] (replicate (size+1) 0) -- size from type 
        where size = fromInteger $ natVal (undefined :: undefined n)

instance (KnownNat m, AddGroup f, Mul f) => Mul (CSR f m m) where
    (*) = cSRMM
    one = csrIdm

instance (Ring f, f ~ f', n ~ n') => Composable (CSR f m n) (Vector f' n') (Vector f m) where
    (**) = (smV)

instance (Ring f, f ~ f', b ~ b') => Composable (CSR f a b) (CSR f' b' c) (CSR f a c) where
    (**) = (cSRMM)

instance (KnownNat m, KnownNat n, Ring f) => VectorSpace (CSR f m n) where
    type Under (CSR f m n) = f
    s £ m = s `scaleCSR` m

instance M.Matrix CSR where

    extend csr = tabulate' . merge' (values csr)
        where merge' as bs = map last . groupBy ((==) `on` fst) $ sortOn fst (as ++ bs)
              count [] i | i >= length (row csr) = [] 
              count [] i = 0 : count [] (i+1)
              count cs i = let (t, d) = span ((Fin i==) . fst . fst) cs in
                  length t : count d (i+1)
              tabulate' xs = CSR elems col row
                  where sorted = sortOn fst xs
                        (col, elems) = unzip $ map (\((_,Fin i),a) -> (i - 1, a)) sorted
                        row = scanl (+) 0 (count sorted 1)

                        
    values (CSR elems col row) = concat $ merge (zip col elems) perRow
        where perRow = zip [1..] $ zipWith (-) (tail row) row 
              merge _ [] = []
              merge xs ((i,n):ys) = let (cur, next) = splitAt n xs in
                    [ ((Fin i,Fin (j+1)), a) | (j,a) <- cur ] : merge next ys

-- | returns size of the sparse matrix
csrSize :: (KnownNat m, KnownNat n) => CSR f m n -> (Int,Int)
csrSize csr@(CSR e c r) = (length r - 1, fromInteger (natVal csr))

-- | returns only column size, useful for square sparse matrixes
csrLen :: KnownNat m => CSR f m m -> Int
csrLen csr = fromInteger (natVal csr)

-- | identity matrix for CSR representation
csrIdm :: (KnownNat m, Mul f) => CSR f m m
csrIdm = idm
  where  idm = CSR ones [0..mN-1] [0..mN]
         ones = replicate mN one
         mN = csrLen idm

-- | scale function for CSR
scaleCSR :: Mul f => f -> CSR f m n -> CSR f m n
scaleCSR c (CSR elems col row) = CSR (map (c*) elems) col row

-- | returns the element of a given coordinate in the sparse matrix
getElem :: AddGroup f => CSR f m n -> (Int,Int) -> f
getElem csr (x,y) = maybe zero id $ lookup x $ getRow csr y

-- | returns a given row in the sparse matrix
getRow :: CSR f m n -> Int -> [(Int, f)]
getRow (CSR elems col row) i = case td i 2 row of
        [a,b] -> td a (b-a) (zip col elems) 
        _     -> []
        where td i j = take j . drop i

-- | returns a given column in the sparse matrix
getColumn :: CSR f m n -> Int -> [(Int, f)]
getColumn (CSR elems col row) i = [ x | (x, y) <- zip cs col, y==i ]
    where
        cs = couple elems row
        couple [] _ = []
        couple e rw@(r:r':rs) = let (xs, ys) = splitAt (r'-r) e in 
                                zip (repeat (length row - length rw)) xs ++ couple ys (r':rs)

-- | dot product between two lists
dotL :: (AddGroup f, Mul f) => [f] -> [f] -> f
dotL v1 v2 = sum $ zipWith (*) v1 v2

-- | dot product between two "csr vectors".
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

-- | matrix multiplication
--  currently slightly off, transposing answer corrects it
-- as it places column vectors as row vectors in the answer. 
cSRMM :: (Ring f) => CSR f a b -> CSR f b c  -> CSR f a c
cSRMM m1@(CSR e1 c1 r1) 
      m2@(CSR e2 c2 r2) = csrTranspose $ foldl comb emptyCSR bs
        where 
            bs = [(csrMV m1 (getColumn m2 b)) | b <- [0..maximum c1]]
            emptyCSR = CSR [] [] (scanl (+) 0 (map length bs)) :: CSR f a c

cSRSub :: (KnownNat b, AddGroup f) => CSR f a b -> CSR f a b  -> CSR f a b
cSRSub m1 m2 = cSRAdd m1 (neg m2)

cSRAdd :: (AddGroup f) => CSR f a b -> CSR f a b  -> CSR f a b
cSRAdd (CSR e1 c1 [r]) _ = CSR [] [] [r]
cSRAdd m1@(CSR e1 c1 r1) 
       m2@(CSR e2 c2 r2) = foldl comb emptyCSR bs
       where 
            bs = opRows m1 m2 cSRAddRow
            emptyCSR = CSR [] [] (scanl (+) 0 (map length bs)) :: CSR f a b

-- Could move this to be apart of cSRAdd
opRows :: (AddGroup f) => CSR f a b -> CSR f a b -> ([(Int,f)] -> [(Int,f)]  -> [(Int,f)]) -> [[(Int, f)]]
opRows m1@(CSR e1 c1 r1) m2 op =  [ (op (getRow m1 a) (getRow m2 a)) | a <- [0..length r1 - 2]]

cSRAddRow :: AddGroup f => [(Int,f)] -> [(Int,f)]  -> [(Int,f)]
cSRAddRow [] [] = []
cSRAddRow [] as = as
cSRAddRow as [] = as
cSRAddRow v1@((c1,e1):as) v2@((c2,e2):bs) | c1 == c2 = (c1,e1+e2) : cSRAddRow as bs
                                          | c1 > c2  = (c2,e2) : cSRAddRow v1 bs
                                          | c1 < c2  = (c1,e1) : cSRAddRow as v2

csrTranspose :: CSR f a b -> CSR f a b
csrTranspose m1@(CSR e1 c1 r1) =  foldl comb emptyCSR $ bs
             where
                bs = map (map snd) $ groupBy ((==) `on` fst) $ sortOn fst $ concat qs
                qs = [ zip (map fst as) (zip (repeat rs) (map snd as))| (rs, as) <- zip [0..] [getRow m1 a | a <- [0..length r1 - 2]] ]
                emptyCSR = CSR [] [] (scanl (+) 0 (map length bs))

-- Slow/bad transpose
csrTranspose2 :: CSR f a b -> CSR f a b
csrTranspose2 m1@(CSR e1 c1 r1) =  foldl comb emptyCSR bs
             where
                bs = [getColumn m1 a | a <- [0..maximum c1]]
                emptyCSR = CSR [] [] (scanl (+) 0 (map length bs))

-- toList :: CSR a m n -> [(Int, Int a)]
-- toCSR  :: [(Int, Int a)] -> CSR a m n

-- Test values/functions

comb :: CSR f x y -> [(Int, f)] -> CSR f x y
comb csr [] = csr
comb (CSR e1 c1 r1) as = CSR (e1++es) (c1++cs) r1
    where
        (cs, es) = unzip as

-- Matrix Vector Multiplication

-- Multiplies a CSR matrix with a Vector 
smV :: (Ring f) => CSR f a b -> Vector f b  -> Vector f a
smV m (V v) = V (smv m v)

-- Multiplies a CSR matrix with a list vector
-- getRow could be used instead of calculating difference to take j from elems&col
smv :: (Ring f) => CSR f a b -> [f] -> [f]
smv (CSR _ _ (r:[])) v = []
smv (CSR [] _ _) v = []
smv (CSR elems col (r:row)) v = dotL (take j elems) (map (v!!) (take j col)) : 
                                     smv (CSR (drop j elems) (drop j col) row) v
            where j = head row - r 

-- Multiplies a CSR matrix with a CSR row/column
csrMV :: (Ring f) => CSR f a b -> [(Int,f)]  -> [(Int,f)]
csrMV m1@(CSR e1 c1 r1) v1 = [(a,dotCsr (getRow m1 a) v1)| a <- [0..length r1 - 2]]

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
denseCSR :: CSR Double 500 500
denseCSR = CSR {
    elems = [1..500*500],
    col = [0,1..(500*500 - 1)],
    row = [0,500..500*500]}

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
triCSR :: (Ring f, Num f, KnownNat m) => CSR f m m
triCSR = csr
  where  csr = CSR (concat stencils) (concat scols) (scanl (+) 0 slens)
         n = csrLen csr
         stencil = [one,negate (one + one),one]
         stencils = [drop 1 stencil] ++ replicate (n-2) stencil ++ [take 2 stencil]
         colfun c = [c-1,c,c+1]
         scols    = [drop 1 (colfun 0)] ++ map colfun [1..n-2]  ++ [take 2 (colfun (n-1))]
         slens = map length stencils

-- tomat for triCSR
-- "triMat :: Matrix Double 5 5" creates a 5 X 5 matrix with type Double
triMat :: (Ring f, Num f, KnownNat m) => Matrix f m m
triMat = toMat triCSR

triPJ :: Matrix Double 10 10
triPJ = triMat