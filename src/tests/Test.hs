{-# Language DataKinds #-} 
{-# Language GADTs #-} 
{-# Language StandaloneDeriving #-} 
{-# Language ScopedTypeVariables #-} 
{-# Language RankNTypes #-} 
{-# Language TypeOperators #-} 
{-# Language TypeApplications #-} 
{-# Language KindSignatures #-} 

module Test(tests) where

import GHC.TypeLits -- hiding (type (^))
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), recip, sum, product, (**), span)

import Algebra
import ListVector

import Distribution.TestSuite.QuickCheck
import Test.QuickCheck

-- The test suite runs a list of different properties
-- To run the test suite use: "cabal test test --allow-newer"
tests :: IO [ Test ]
tests = return [testProperty "Always true" True] 


-- | Generator for arbitrary vectors of a given size.
--   To see some example vectors of length 5 
--   > sample (arbitrary @(Vector R 5) ) 
instance forall n f. (KnownNat n, Arbitrary f) => Arbitrary (Vector f n) where
    arbitrary = V <$> vector ( vecLen (undefined :: Vector () n) )

-- | Generator for matrices of a given size
instance forall m n f. (KnownNat m, KnownNat n, Arbitrary f) => Arbitrary (Matrix f m n) where
    arbitrary = do 
        let v = arbitrary @(Vector f m)
        M . V <$> vectorOf (vecLen (undefined :: Vector () n)) v



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

-- | Tests vector addition properties for a given vector lenght

prop_vectorAddZero :: KnownNat n => Vector R n -> Bool
prop_vectorAddZero v = v + zero == v

prop_vectorAddAssoc :: KnownNat n => Vector R n -> Vector R n -> Bool
prop_vectorAddAssoc v w = v + w == w + v

-- | Tests cross product properties for a given vector lenght

prop_crossDot :: Vector Rational 3 -> Vector Rational 3 -> Bool
prop_crossDot v1 v2 = dot v1 (cross v1 v2) == dot v2 (cross v1 v2)
                            && dot v1 (cross v1 v2) == 0

prop_crossProduct :: Vector Rational 3 -> Vector Rational 3 -> Bool
prop_crossProduct v1 v2 = cross v1 v2 == neg (cross v2 v1)                            

-- | Tests vector dot product properties for a given vector lenght

prop_dotProduct :: KnownNat n => Vector R n -> Vector R n -> Bool
prop_dotProduct v1 v2 = dot v1 v2 == dot v2 v1

-- | Tests matrix multipliciation properties for a given matrix/vector lenght

prop_matVecmul :: (KnownNat m, Field f,Eq f) =>
                    Matrix f m n -> Vector f n -> Bool
prop_matVecmul m v = (m ££ v) == (m ££ v)

prop_matMatmul :: (KnownNat a, KnownNat b, Field f,Eq f) =>
                    Matrix f a b -> Bool
prop_matMatmul m1 = m1 £££ idm == idm £££ m1 && m1 £££ idm == m1

prop_matMatmul2 :: (KnownNat a,KnownNat b, KnownNat c, Field f,Eq f) =>
                    Matrix f a b -> Matrix f b c -> Matrix f c d -> Bool
prop_matMatmul2 m1 m2 m3 = m1 £££ (m2 £££ m3) == (m1 £££ m2) £££ m3 

-----------------------------------------------------------------------------------
-- Generator for arbitrary sizes 
-- 
-- Since the size of Vector and Matrix is part of its type we can not generate 
-- them with arbitrary sizes directly.
-- As a workaround we can generate hidden singletons corresponding to the size
-- and then use these to generate vectors and matrices of arbitrary size.
--
-- As an example of how we would like this to work see prop_arbitraryAddZero
-- @
-- prop_arbitraryAddZero :: HiddenNat -> Property
-- prop_arbitraryAddZero = liftVecProp1 prop_vectorAddZero
-- @
-- 
-- prop_arbitraryAddZero does not work because it can not deduce KnownNat n, 
-- which is required of the instance VectorSpace (Vector f n).
-- Since we use singletons we can deduce it by pattern matching but 
-- there should be a better way to do this automatically.
--
-- For a working example see have prop_genVecWithLen.


-- | Singleton type over the kind Nat
data SNat (n :: Nat) where
    Nil :: SNat 0 
    Suc :: SNat n -> SNat (n + 1)
deriving instance Show (SNat n) 

-- | Existential for SNat
data HiddenNat = forall n. Hidden (SNat n)

instance Show HiddenNat where
    show = ("Hidden " ++) . show . hiddenToInt

-- | Converts the type n to a corresponding value n
toInt :: SNat n -> Int
toInt Nil     = 0
toInt (Suc n) = 1 + toInt n

-- | Lifts a function on SNat to a function on HiddenNat
--   Restricted to not return any type information of n.
--   for example 
--   > onHidden id 
--   results in an error "...'n' would escape its scope..."
onHidden :: (forall n. SNat n -> a ) -> HiddenNat -> a
onHidden f (Hidden n) = f n

-- | Gets the value of the hidden SNat
hiddenToInt :: HiddenNat -> Int
hiddenToInt = onHidden toInt

-- | Converts a Int to a SNat and wraps it.
toHidden :: Int -> HiddenNat
toHidden 0 = Hidden Nil
toHidden n = go n Nil
    where go :: Int -> SNat n -> HiddenNat
          go 0 m = Hidden m
          go n m = go (n-1) (Suc m) 


-- | Generator for SNat
--   Since we cannot expose the type we wrap the result in a HiddenNat
instance Arbitrary HiddenNat where
    arbitrary = toHidden . abs <$> arbitrary

-- | Generates a arbitrary Vector with a given size
genVecWithLen :: Arbitrary f => SNat n -> Gen (Vector f n)
genVecWithLen n = V <$> vector (toInt n) -- vector comes from QuickCheck

-- | Generates a arbitrary Matrix with a given size
genMatWithSize :: Arbitrary f => SNat m -> SNat n -> Gen (Matrix f m n)
genMatWithSize m n = do 
    let v = genVecWithLen m  
    M . V <$> vectorOf (toInt n) v


-- | Takes a prop for vectors of arbitrary length and converts it to a prop
--   on HiddenNat. This enables quickCheck to actually test the property against 
--   vectors of arbitrary size. For example:
--   > quickCheck dummy_prop
--   is not valid since quickCheck needs to know the size of vectors to test. 
--   > quickCheck $ liftVecProp1 dummy_prop
--   will run and test vectors of different sizes.
liftVecProp1 :: (Testable prop, Arbitrary f, Show f) => (forall n. Vector f n -> prop) -> HiddenNat -> Property
liftVecProp1 p (Hidden n) = forAll (genVecWithLen n) p

-- | dummy_prop to test liftVecProp1
dummy_prop :: Vector R n -> Bool
dummy_prop (V x) = length (zipWith (+) x x) == length x 


-- | Tests that the length of the vectors list is equal to the value of HiddenNat
--   Since HiddenNat is an instance of Arbitrary we can use quickCheck 
--   > quickCheck prop_correctLengthVec
prop_genVecWithLen :: HiddenNat -> Property
prop_genVecWithLen (Hidden n) = forAll (genVecWithLen @R n) $ -- arbitrary vector in R^2
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





