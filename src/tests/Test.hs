{-# Language DataKinds #-} 
{-# Language GADTs #-} 
{-# Language StandaloneDeriving #-} 
{-# Language ScopedTypeVariables #-} 
{-# Language RankNTypes #-} 
{-# Language TypeOperators #-} 
{-# Language TypeApplications #-} 
{-# Language KindSignatures #-} 

module Main where

import GHC.TypeLits hiding (type (^))
import qualified Prelude
import Prelude hiding ((+),(*), (-), (/), recip, sum, product, (**), span)

import Algebra
import ListVector
import Eigen
import HiddenNat hiding (Suc)
import QuadTree hiding (toInt)
import SparseCSR

import Test.QuickCheck


main :: IO ()
main = putStrLn "Test loaded"


-- | Generator for arbitrary vectors of a given size.
--   To see some example vectors of length 5 
--   > sample (arbitrary @(Vector R 5) ) 
instance forall n f. (KnownNat n, Arbitrary f) => 
                        Arbitrary (Vector f n) where
    arbitrary = V <$> vector ( vecLen (undefined :: Vector () n) )

-- | Generator for matrices of a given size
instance forall m n f. (KnownNat m, KnownNat n, Arbitrary f) => 
                        Arbitrary (Matrix f m n) where
    arbitrary = do 
        let v = arbitrary @(Vector f m)
        M . V <$> vectorOf (vecLen (undefined :: Vector () n)) v


instance Arbitrary f => Arbitrary (Quad One f) where
    arbitrary = do 
            s <- arbitrary
            elements [Zero, Scalar s]

instance (Sized n, Arbitrary f, Arbitrary (Quad n f)) => Arbitrary (Quad (Suc n) f) where
    arbitrary = do
            (nw, ne, sw, se) <- arbitrary
            frequency [(1, zero'), (3, mtx' nw ne sw se) ]
        where zero' = return Zero
              mtx' a b c d = return $ Mtx a b c d


-----------------------------------------------------------------------------------
-- Test for vectors 
--
-- As of now we can only run these test if we provide a given length.
-- See next section for tests on arbitrary sizes.
--
-- For example:
-- > quickCheck (prop_vectorAddAssoc @5)
-- tests that addition is associative for R^5
--

--  Tests vector addition properties for a given vector lenght
prop_vectorAddZero :: KnownNat n => Vector R n -> Bool
prop_vectorAddZero v = v + zero == v

-- Vector addition is commutative
prop_vectorAddComm :: KnownNat n => Vector R n -> Vector R n -> Bool
prop_vectorAddComm v1 v2 = v1 + v2 == v2 + v1

-- Vector addition is associative
prop_vectorAddAssoc :: KnownNat n => Vector Rational n -> 
                        Vector Rational n -> Vector Rational n -> Bool
prop_vectorAddAssoc v1 v2 v3 = (v1 + v2) + v3 == v1 + (v2 + v3)


-- | Tests cross product properties
prop_crossDot :: Vector Rational 3 -> Vector Rational 3 -> Bool
prop_crossDot v1 v2 = dot v1 (cross v1 v2) == dot v2 (cross v1 v2)
                            && dot v1 (cross v1 v2) == 0

prop_crossProduct :: Vector Rational 3 -> Vector Rational 3 -> Bool
prop_crossProduct v1 v2 = cross v1 v2 == neg (cross v2 v1)                            

--  Tests vector dot product properties for a given vector lenght
prop_dotProduct :: KnownNat n => Vector R n -> Vector R n -> Bool
prop_dotProduct v1 v2 = dot v1 v2 == dot v2 v1


--  Tests LIL matrix multipliciation properties
--  for a given matrix/vector lenght
prop_matMatmulIdm :: (KnownNat a, KnownNat b, Field f,Eq f) =>
                    Matrix f a b -> Bool
prop_matMatmulIdm m1 = m1 £££ idm == idm £££ m1 && m1 £££ idm == m1

prop_matMatmulAssoc :: (KnownNat a,KnownNat b, KnownNat c, Field f,Eq f) =>
                    Matrix f a b -> Matrix f b c -> Matrix f c d -> Bool
prop_matMatmulAssoc m1 m2 m3 = m1 £££ (m2 £££ m3) == (m1 £££ m2) £££ m3 


-- Test on Quad matrices for a given size
prop_quadAddZero :: (Sized n, AddGroup f,Eq f) => Quad n f -> Bool
prop_quadAddZero m1 = m1 + zero == zero + m1 && m1 + zero == m1

prop_quadAddComm :: Sized n => Quad n R -> Quad n R -> Bool
prop_quadAddComm m1 m2 = m1 + m2 == m2 + m1

prop_quadAddAssoc :: Sized n => Quad n Rational -> Quad n Rational -> Quad n Rational -> Bool
prop_quadAddAssoc m1 m2 m3 = (m1 + m2) + m3 == m1 + (m2 + m3)


prop_quadMulIdQ :: (Sized n, Field f,Eq f) => Quad n f -> Bool
prop_quadMulIdQ m1 = m1 `mulQ` idQ == idQ `mulQ` m1 && m1 `mulQ` idQ == m1

prop_quadMulAssoc :: (Sized n, Field f,Eq f) => Quad n f -> Quad n f -> Quad n f -> Bool
prop_quadMulAssoc m1 m2 m3 = m1 `mulQ` (m2 `mulQ` m3) == (m1 `mulQ` m2) `mulQ` m3 


-- Test on determinant
prop_detHomomorphism :: (KnownNat n, Field f, Eq f) => Matrix f n n -> Matrix f n n -> Bool
prop_detHomomorphism m1 m2 = detNN(m1 £££ m2) == detNN(m1) * detNN(m2) 

-- ElimOp and determinant


-- Test on eigenvalues/vectors

-- Matrix in characteristic polynomial



-- Test on utility functions

prop_transpose :: Eq f => Matrix f m n -> Bool
prop_transpose m = transpose (transpose m) == m

--prop_onUnpackedTrans :: Eq f => Matrix f m n -> Bool
--prop_onUnpackedTrans m = onUnpackedTrans id m == m

prop_funMatComp :: (KnownNat b, KnownNat a, Field f, Eq f) => 
    (Vector f b -> Vector f a) -> Matrix f b c -> Vector f c -> Bool
prop_funMatComp f m v = f (m ££ v) == (f `onCols` m) ££ v

-- matrix to array and back
prop_packUnpack :: Eq f => Matrix f m n -> Bool
prop_packUnpack m = pack(unpack m) == m


-- Test on generators for arbitrary sizes

-- | Tests that the length of the vectors list is equal to the value of HiddenNat
--   Since HiddenNat is an instance of Arbitrary we can use quickCheck 
--   > quickCheck prop_correctLengthVec
prop_genVecWithLen :: HiddenNat -> Property
prop_genVecWithLen (Hidden n) = forAll (genVecWithLen @R n) $ -- arbitrary vector in R^n
    \(V l) -> length l == toInt n

-- | Tests that genMatWithSize creates a matrix of the right size
prop_genMatWithSize :: HiddenNat -> HiddenNat -> Property
prop_genMatWithSize (Hidden m) (Hidden n) = forAll (genMatWithSize @R m n) $ 
    \(M (V vs)) -> let rect = allEqual $ map vecLen vs
                       n' = length vs
                       m' = vecLen (head vs)
                   in rect && n' == toInt n && if n' > 0 then m' == toInt m else True  
  where vecLen (V xs) = length xs
        allEqual []     = True
        allEqual (x:xs) = all (==x) xs

