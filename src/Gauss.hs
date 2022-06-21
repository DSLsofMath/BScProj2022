{-# language DataKinds #-}

module Gauss where

import GHC.TypeLits 
import Prelude hiding ((+), (-), (*), (/), (^), recip, sum, product, (**), span)
import Data.List (sortOn)

import Control.Monad.State.Lazy

import Algebra 
import Matrix 


-------------------------------------------
-- Elementary row operations definition 
--

-- | Represents elementary row operations 
data ElimOp n a = Swap     (Fin n) (Fin n)
                | Mul    a (Fin n)
                | MulAdd a (Fin n) (Fin n)
              deriving (Eq)

instance Show a => Show (ElimOp n a) where 
    show = showElimOp

-- | Prettier show function for ElimOp a
showElimOp :: Show a => ElimOp n a -> String
showElimOp op = concat $ case op of 
        Swap   (Fin i) (Fin j)   -> [ " ",              row i,              " <-> ", row j ]
        Mul    s (Fin i)         -> [ " ", show s, "*", row i,               " -> ", row i ]
        MulAdd s (Fin i) (Fin j) -> [ " ", show s, "*", row i, " + ", row j, " -> ", row j ]
    where row i = "R(" ++ show i ++ ")"


-- -- | Shows step by step how a matrix is transformed by a ElimOp trace
-- showElimOnMat :: (Field f, Show f) => [ElimOp m f] -> Matrix f m n -> String
-- showElimOnMat t m0 = let matTrace = scanl (flip elimOpToFunc) m0 t 
--                      in unlines [ show m ++ "\n" ++ show op | (m, op) <- zip matTrace t ]
--                         ++ show (last matTrace)


-- | Representation of an elementary row operation as a matrix 
elimOpToMat :: (KnownNat n, Matrix mat, Ring f) => ElimOp n f -> mat n n f
elimOpToMat e = elimOpToFunc e identity

swap :: (KnownNats m n, Matrix mat, AddGroup f) => 
          Fin m -> Fin m -> mat m n f -> mat m n f
swap i j m = setRow (setRow m i (getRow m j)) j (getRow m i) 

mul :: (KnownNats m n, Matrix mat, Ring f) => 
        f -> Fin m -> mat m n f -> mat m n f
mul s i m = setRow m i (s £ getRow m i)

mulAdd :: (KnownNats m n, Matrix mat, Ring f) => 
            f -> Fin m -> Fin m -> mat m n f -> mat m n f
mulAdd s i j m = setRow m j (s £ getRow m i + getRow m j)

-- | Representation of an elementary row operation as a function 
elimOpToFunc :: (KnownNats m n, Matrix mat, Ring f) => ElimOp m f -> mat m n f -> mat m n f
elimOpToFunc e m = case e of Swap     i j -> swap i j m
                             Mul    s i   -> mul s i m
                             MulAdd s i j -> mulAdd s i j m

-- | Reduces a trace of elimOps to a single function
foldElimOpsFunc :: (KnownNats m n,Matrix mat, Ring f) => [ElimOp m f] -> (mat m n f -> mat m n f)
foldElimOpsFunc = foldr (.) id . map elimOpToFunc . reverse


--------------------------------------------------------------------------
-- Gauss elimination 
--


data GaussState mat m n f = GS {
            matrix :: mat m n f,
            trace  :: [ElimOp m f],
            pivotIndex :: (Fin m, Fin n)
            }

getTrace :: GaussState mar m n f -> [ElimOp m f]
getTrace = reverse . trace

type Gauss mat m n f a = State (GaussState mat m n f) a


-- | initial state of the gauss algorithm
initGauss :: (KnownNats m n, Matrix mat) => mat m n f -> GaussState mat m n f
initGauss m = GS m [] (minBound, minBound)

-- | Returns all nonzero elements of the column in sorted order
getColumn :: (Matrix mat, Approx f, AddGroup f) => Fin n -> Gauss mat m n f [(Fin m, f)]
getColumn j = sortOn fst . filter ((~/= zero) . snd) . flip getCol j <$> gets matrix

-- | Performs an elementary row operation and adds it to the trace
doRowOp :: (KnownNats m n, Matrix mat, Ring f) => ElimOp m f -> Gauss mat m n f ()
doRowOp rowOp = modify $ \(GS mat trace p) -> 
                    GS (elimOpToFunc rowOp mat) (rowOp : trace) p

-- Sets the position of current pivot element
setPivot :: (Fin m, Fin n) -> Gauss mat m n f ()
setPivot i = modify $ \gs -> gs { pivotIndex = i }

-- Gets the index of the current pivot element
getPivot :: Gauss mat m n f (Fin m, Fin n)
getPivot = gets pivotIndex

succMaybe :: (Enum a, Bounded a, Eq a) => a -> Maybe a
succMaybe i = if i /= maxBound then Just (succ i) else Nothing

predMaybe :: (Enum a, Bounded a, Eq a) => a -> Maybe a
predMaybe i = if i /= minBound then Just (pred i) else Nothing


-- | Transforms a given matrix into row echelon form
gauss :: (KnownNats m n, Matrix mat, Field f, Approx f) => mat m n f -> mat m n f
gauss = snd . gauss' 

-- | Returns the row operations that transforms a given matrix into row echelon form
gaussTrace :: (KnownNats m n, Matrix mat, Field f, Approx f) => mat m n f -> [ElimOp m f]
gaussTrace = fst . gauss' 

gauss' :: (KnownNats m n, Matrix mat, Field f, Approx f) => mat m n f -> ([ElimOp m f], mat m n f)
gauss' m = (getTrace gaussed, matrix gaussed)
    where gaussed = execState gauss'' (initGauss m) 

gauss'' :: (KnownNats m n, Matrix mat, Field f, Approx f) => Gauss mat m n f () 
gauss'' = do (i, j) <- getPivot
             col <- filter ((>= i) . fst) <$> getColumn j
             nextRowIndex <- case col of 
                 []               -> return (Just i)
                 ((i', s) : col') -> do 
                         unless (i' == i  ) $ doRowOp (Swap i' i)
                         unless (s  ~= one) $ doRowOp (Mul (recip s) i)
                         mapM_ doRowOp [MulAdd (neg a) i j | (j, a) <- col']
                         return (succMaybe i)
             case (nextRowIndex, succMaybe j) of
                 (Just i', Just j') -> setPivot (i', j') >> gauss''
                 _                  -> return ()

