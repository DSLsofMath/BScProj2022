{-# Language DataKinds #-} 
{-# Language GADTs #-} 
{-# Language StandaloneDeriving #-} 
{-# Language ScopedTypeVariables #-} 
{-# Language RankNTypes #-} 
{-# Language TypeOperators #-} 
{-# Language TypeApplications #-} 
{-# Language KindSignatures #-} 


module HiddenNat where

import GHC.TypeLits 
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), recip, sum, product, (**), span)

import Test.QuickCheck

import Algebra
import ListVector

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



