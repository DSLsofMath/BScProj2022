{-# Language DataKinds #-} 
{-# Language GADTs #-} 
{-# Language StandaloneDeriving #-} 
{-# Language ScopedTypeVariables #-} 
{-# Language RankNTypes #-} 
{-# Language TypeOperators #-} 
{-# Language TypeApplications #-} 
{-# Language KindSignatures #-} 
{-# Language FlexibleInstances #-} 

module Main where

import GHC.TypeLits -- hiding (type (^))
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), recip, sum, product, (**), span)

import Algebra
import ListVector


import Test.QuickCheck


main :: IO ()
main = putStrLn "Test loaded"



-- | Generator for arbitrary vectors of a given size.
--   To see some example vectors of length 5 
--   > sample (arbitrary @(Vector R 5) ) 
instance forall n f. (KnownNat n, Arbitrary f) => Arbitrary (Vector f n) where
    arbitrary = V <$> vector ( vecLen (undefined :: Vector () n) )



-----------------------------------------------------------------------------------
-- Test for vectors 
--
-- As of now we can only run a test if we provide a given length.
-- For example:
-- > quickCheck (prop_vectorAddAssoc @5)
-- tests that addition is associative for R^5


prop_vectorAddZero :: KnownNat n => Vector R n -> Bool
prop_vectorAddZero v = v + zero == v

prop_vectorAddAssoc :: KnownNat n => Vector R n -> Vector R n -> Bool
prop_vectorAddAssoc v w = v + w == w + v





-----------------------------------------------------------------------------------
-- Generator for arbitrary sizes 
-- 
-- Since the size of Vector and Matrix is part of its type we can not generate 
-- them with arbitrary sizes directly.
-- As a workaround we can generate hidden singletons corresponding to the size
-- and then use these to generate vectors and matrices of arbitrary size.
--
-- As an example of how we would like this to work see prop_arbitraryAddZero
--
-- @
-- prop_arbitraryAddZero :: HiddenNat -> Property
-- prop_arbitraryAddZero (Hidden n) = property $ do 
--         v <- genVecWithLen n
--         return $ prop_vectorAddZero v
-- @
-- 
-- prop_arbitraryAddZero does not work because it can not deduce KnownNat n, 
-- which is required of the instance VectorSpace (Vector f n).
-- Since we use singletons we can deduce it by pattern matching but 
-- there should be a better way to do this automatically.
--
-- For a working example we have prop_correctLength.


-- | Singleton type over the kind Nat
data SNat (n :: Nat) where
    Nil :: SNat 0 
    Suc :: SNat n -> SNat (n + 1)
deriving instance Show (SNat n) 
    
-- | Existential for SNat
data HiddenNat = forall n. Hidden (SNat n)
deriving instance Show HiddenNat

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
    arbitrary = sized $ return . toHidden 


-- | Generates a arbitrary Vector with a given size
genVecWithLen :: Arbitrary f => SNat n -> Gen (Vector f n)
genVecWithLen n = V <$> vector (toInt n) -- vector comes from QuickCheck



-- | Tests that the length of the vectors list is equal to the value of HiddenNat
prop_correctLength :: HiddenNat -> Property
prop_correctLength (Hidden n) = property $ do 
    v <- genVecWithLen @R n  -- arbitrary vector in R^n
    let V l = v 
    return $ length l == toInt n 






