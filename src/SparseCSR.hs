{-# LANGUAGE BlockArguments #-}
{-# LANGUAGE DataKinds #-} 
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeOperators #-} 
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE UndecidableInstances #-}



module SparseCSR where

import GHC.TypeLits
import Algebra
import Prelude hiding ((+), (-), (*), (/), sum)
import ListVector 
import qualified Matrix as M
import Matrix (values, Fin(..))

import Data.Proxy
import Data.Function
import Data.List hiding (sum)
import qualified Data.List as L

-- A sparse matrix using the compressed sparse row format
data CSR (m :: Nat) (n :: Nat) f = CSR { elems :: [f],
                            col :: [Int],
                            row :: [Int]}

deriving instance Show f => Show (CSR m n f)


instance forall m n f. (KnownNat n, AddGroup f) => AddGroup (CSR m n f) where
    (+) = cSRAdd
    (-) = cSRSub
    neg (CSR e c r) = CSR (map neg e) c r 
    zero = CSR [] [] (replicate (n + 1) 0) -- size from type
        where n = fromInteger $ natVal (Proxy @n)

instance (Ring f, f ~ f', n ~ n') => Composable (CSR m n f) (Vector n' f') (Vector m f) where
    (**) = (smV)

instance (KnownNat m, KnownNat n, Ring f) => VectorSpace (CSR m n f) where
    type Under (CSR m n f) = f
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

    mulMat = cSRMM

-- | returns size of the sparse matrix
csrSize :: forall m n f. KnownNats m n => CSR m n f -> (Int,Int)
csrSize _ = (m, n)
    where m = fromInteger $ natVal (Proxy @m)
          n = fromInteger $ natVal (Proxy @n)

-- | returns only column size, useful for square sparse matrixes
csrLen :: KnownNat m => CSR m m f -> Int
csrLen = snd . csrSize

-- | identity matrix for CSR representation
csrIdm :: (KnownNat m, Mul f) => CSR m m f
csrIdm = idm
  where  idm = CSR ones [0..mN-1] [0..mN]
         ones = replicate mN one
         mN = csrLen idm

-- | scale function for CSR
scaleCSR :: Mul f => f -> CSR m n f -> CSR m n f
scaleCSR c (CSR elems col row) = CSR (map (c*) elems) col row

-- | returns the element of a given coordinate in the sparse matrix
getElem :: AddGroup f => CSR m n f -> (Int,Int) -> f
getElem csr (x,y) = maybe zero id $ lookup x $ getRow csr y

-- | returns a given row in the sparse matrix
getRow :: CSR m n f -> Int -> [(Int, f)]
getRow (CSR elems col row) i = case td i 2 row of
        [a,b] -> td a (b-a) (zip col elems) 
        _     -> []
        where td i j = take j . drop i

-- | returns a given column in the sparse matrix
getColumn :: CSR m n f -> Int -> [(Int, f)]
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
cSRMM :: (Ring f) => CSR a b f -> CSR b c f  -> CSR a c f
cSRMM m1@(CSR e1 c1 r1) 
      m2@(CSR e2 c2 r2) = csrTranspose $ foldl comb emptyCSR bs
        where 
            bs = [(csrMV m1 (getColumn m2 b)) | b <- [0..maximum c1]]
            emptyCSR = CSR [] [] (scanl (+) 0 (map length bs)) :: CSR f a c

cSRSub :: (KnownNat a, AddGroup f) => CSR b a f -> CSR b a f  -> CSR b a f
cSRSub m1 m2 = cSRAdd m1 (neg m2)

cSRAdd :: (AddGroup f) => CSR b a f -> CSR b a f  -> CSR b a f
cSRAdd (CSR e1 c1 [r]) _ = CSR [] [] [r]
cSRAdd m1@(CSR e1 c1 r1) 
       m2@(CSR e2 c2 r2) = foldl comb emptyCSR bs
       where 
            bs = opRows m1 m2 cSRAddRow
            emptyCSR = CSR [] [] (scanl (+) 0 (map length bs)) :: CSR b a f

-- applies a given operation row by row,
-- between rows from two given matrices.
opRows :: (AddGroup f) => CSR b a f -> CSR b a f -> ([(Int,f)] -> [(Int,f)]  -> [(Int,f)]) -> [[(Int, f)]]
opRows m1@(CSR e1 c1 r1) m2 op =  [ (op (getRow m1 a) (getRow m2 a)) | a <- [0..length r1 - 2]]

-- Adds one csr vector with another, useful for sparse vectors.
cSRAddRow :: AddGroup f => [(Int,f)] -> [(Int,f)]  -> [(Int,f)]
cSRAddRow [] [] = []
cSRAddRow [] as = as
cSRAddRow as [] = as
cSRAddRow v1@((c1,e1):as) v2@((c2,e2):bs) | c1 == c2 = (c1,e1+e2) : cSRAddRow as bs
                                          | c1 > c2  = (c2,e2) : cSRAddRow v1 bs
                                          | c1 < c2  = (c1,e1) : cSRAddRow as v2

csrTranspose :: CSR b a f -> CSR b a f
csrTranspose m1@(CSR e1 c1 r1) =  foldl comb emptyCSR $ bs
             where
                bs = map (map snd) $ groupBy ((==) `on` fst) $ sortOn fst $ concat qs
                qs = [ zip (map fst as) (zip (repeat rs) (map snd as))| (rs, as) <- zip [0..] [getRow m1 a | a <- [0..length r1 - 2]] ]
                emptyCSR = CSR [] [] (scanl (+) 0 (map length bs))

-- Appends two csr matrices.
comb :: CSR x y f -> [(Int, f)] -> CSR x y f
comb csr [] = csr
comb (CSR e1 c1 r1) as = CSR (e1++es) (c1++cs) r1
    where
        (cs, es) = unzip as

--
-- Matrix Vector Multiplication
--

-- Multiplies a CSR matrix with a Vector 
smV :: (Ring f) => CSR a b f -> Vector b f  -> Vector a f
smV m (V v) = V (smv m v)

-- Multiplies a CSR matrix with a list vector
-- getRow could be used instead of calculating difference to take j from elems&col
smv :: (Ring f) => CSR b a f -> [f] -> [f]
smv (CSR _ _ (r:[])) v = []
smv (CSR [] _ _) v = []
smv (CSR elems col (r:row)) v = dotL (take j elems) (map (v!!) (take j col)) : 
                                     smv (CSR (drop j elems) (drop j col) row) v
            where j = head row - r 

-- Multiplies a CSR matrix with a CSR row/column
csrMV :: (Ring f) => CSR a b f -> [(Int,f)]  -> [(Int,f)]
csrMV m1@(CSR e1 c1 r1) v1 = [(a,dotCsr (getRow m1 a) v1)| a <- [0..length r1 - 2]]

--
-- Test values/example matrices
--

test :: CSR 4 4 R
test = CSR {
    elems = [ 5, 8, 3, 6 ],
    col = [ 0, 1, 2, 1 ],
    row = [ 0, 1, 2, 3, 4 ]}

test1 :: CSR 4 4 R
test1 = CSR {
    elems = [ 5, 4, 8, 3, 6 ],
    col = [ 0, 1, 1, 2, 1 ],
    row = [ 0, 2, 3, 4, 5 ]}

-- Large
denseCSR5 :: CSR 5 5 R
denseCSR5 = CSR {
    elems = [1..25],
    col = concat $ replicate 5 [0..4],
    row = [0,5..5*5]}

-- Large
denseCSR500 :: CSR 500 500 R
denseCSR500 = CSR {
    elems = [1..500*500],
    col = concat $ replicate 500 [0..499],
    row = [0,500..500*500]}

colVecTest :: CSR 4 1 R
colVecTest = CSR {
    elems = [4,1,3],
    col = [0, 0, 0],
    row = [0, 0, 1, 2, 3]
}

rowVecTest :: CSR 3 3 R
rowVecTest = CSR {
    elems = [4,1,9],
    col = [0, 1, 2],
    row = [0, 3, 3, 3]
}

-- tridiagonal with -2 on diagonal and +1 above and below diagonal
triCSR :: (Ring f, Num f, KnownNat m) => CSR m m f
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
triMat :: (Ring f, Num f, KnownNat m) => Matrix m m f
triMat = toMat triCSR

triPJ :: Matrix 10 10 Double
triPJ = triMat
