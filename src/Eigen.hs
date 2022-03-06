{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE TypeApplications #-}


module Eigen where

import GHC.TypeLits
import qualified Prelude
import Prelude hiding ((+), (-), (*), (/), recip, sum, product, (**), span)

import qualified Data.List as L
import Algebra
import ListVector hiding (v1,v2,m1,m2)

-- This file contains an example of using TypeLits to handle vectors with a given size. 
-- The implementation here is based on lists and should be replaced.
--
-- To try the code in ghci the DataKinds language extension needs to be enabled:
-- type ":set -XDataKinds" in the ghci prompt
--

-- Implementation for 2X2 and 3x3 determinants
-- seems to work, tested several inputs and compared with wolfram


detPlus :: Field f => Matrix f m n -> [f]
detPlus m = zipWith (*) (zipWith (*) r1 r2) r3
    where
        dubM = unpack $ (transpose m) `append` (transpose m)
        r1   = head dubM
        r2   = take 2 (drop 1 (head $ drop 1 dubM)) ++ [head $ head $ drop 4 dubM]
        r3   =  [last $ head $ drop 2 dubM] ++ take 2 (last dubM)


detMinus :: Field f => Matrix f m n -> [f]
detMinus m = zipWith (*) (zipWith (*) r1 r2) r3
    where
        dubM = unpack $ (transpose m) `append` (transpose m)
        r1   = reverse $ head $ drop 3 dubM
        r2   = reverse $ [last $ head $ drop 1 dubM] ++ take 2 (head $ drop 4 dubM)
        r3   = reverse $ take 2 (drop 1 (head $ drop 2 dubM)) ++ [head $ last dubM]


det22 ::  Field f => Matrix f m n -> f
det22 m = p1 - p2
    where
        cols = unpack m
        p1 = (head $ head cols) * (last $ last cols)
        p2 = (head $ last cols) * (last $ head cols)

det33 ::  Field f => Matrix f m n -> f
det33 m = sum (detPlus m) - sum (detMinus m)


m44 = toMat [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]] :: MatR 4 4
mtest = transpose $ toMat [[-5,4,1,7],[-9,3,2,-5],[-2,0,-1,1],[1,14,0,3]] :: MatR 4 4


testdet = det mtest == 1693

-- Determinant 4x4 for matrices
det44 :: Field f => Matrix f m n -> f
det44 m = a*det33 m1 - b*det33 m2 + c*det33 m3 - d*det33 m4
    where
        m1 = pack $ L.transpose $ drop 1 (L.transpose $ drop 1 (unpack m))
        m2 = pack $ drop 1 (head $ unpack m) : drop 1 (unpack m1)
        m3 = pack $ head (unpack m2) : head (unpack m1) : drop 2 (unpack m1)
        m4 = pack $ head (unpack m2) : take 2 (unpack m1)
        
        row1 = head $ unpack $ transpose m  -- row1 = [a,b,c,d]
        a  = head row1
        b  = head $ drop 1 row1
        c  = head $ drop 2 row1
        d  = last row1


-- Determinants
det :: Field f => Matrix f m n -> f
det m | cols /= rows = error "Matrix not squared"
      | cols == 1    = head $ head $ unpack m    -- det [[x]] == x
      | cols == 2    = det22 m
      | cols == 3    = det33 m
      | cols == 4    = det44 m
      | otherwise    = error "Not defined yet"
    where
        cols = length $ unpack m
        rows = length $ head $ unpack m


--tests, all give correct results
idm33 = idm :: MatR 3 3 -- det33 = 1
idm22 = idm :: MatR 2 2 -- det 22 = 1

m22   = toMat [[10,1],[1,10]] :: MatR 2 2 -- det22 = 99
m33   = toMat [[1,-3,2],[3,-1,3],[2,-3,1]] :: MatR 3 3 -- det33 = -15