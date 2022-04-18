open import Data.List
open import Data.Nat using (ℕ; zero; suc)
open import Data.Product
open import Data.Bool
open import Algebra
open import Level using (Level)


module test where
   
    variable
      a b c ℓ : Level
      A : Set a 
      B : Set b
      C : Set c
      m n : ℕ
      
    module Matrix where
        -- Inductive definition of a vector
        data Vector (A : Set a) : (n : ℕ) → Set a where
          [] : Vector A zero
          _::_ :  A → Vector A n → Vector A (suc n)

        infixr 5 _::_

        Veclength : Vector A n -> ℕ 
        Veclength {A} {a} {n} v = n
        
        -- Matrices are defined as vector of vectors 
        Matrix : (A : Set a) (m n : ℕ) → Set a
        Matrix A m n = Vector (Vector A n) m

        Matlength : Matrix A m n -> ℕ × ℕ 
        Matlength {A} {a} {m} {n} mat  = m , n

        v1 : Vector ℕ 4
        v1 = 1 :: 3 :: 4 :: 5 :: []

        m1 : Matrix ℕ 2 2 
        m1 =  (1 :: 2 :: []) :: (1 :: 2 :: []) :: []  

        mat1 = Matrix Bool 2 2

        zipV : (A → B → C) → (Vector A n → Vector B n → Vector C n)
        zipV f [] [] = []
        zipV f (x :: xs) (y :: ys) = f x y :: zipV f xs ys


-- Some operations

    module Operations (S : Semiring c ℓ) where
        open Semiring S   
        open Matrix             

        infixr 6 _+v_
        
        -- Vector addition
        _+v_ : (u v : Vector Carrier n) → Vector Carrier n 
        _+v_ = zipV _+_    
