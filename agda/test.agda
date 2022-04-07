open import Data.List
open import Data.Nat 
open import Data.Product
open import Data.Bool

variable A : Set
         m n : ℕ

data Vector (A : Set) : ℕ → Set where
  [] : Vector A zero
  _::_ : {n : ℕ} → A → Vector A n → Vector A (suc n)

Veclength : Vector A n -> ℕ 
Veclength {A} {n} v = n

infixr 5 _::_

Matrix : (A : Set ) (m n : ℕ) → Set
Matrix A m n = Vector (Vector A n) m

Matlength : Matrix A m n -> ℕ × ℕ 
Matlength {A} {m} {n} m₁ = m , n

v1 : Vector ℕ 4
v1 = 1 :: 3 :: 4 :: 5 :: []

m1 : Matrix ℕ 2 2 
m1 =  (1 :: 2 :: []) :: (1 :: 2 :: []) :: []  

mat1 = Matrix Bool 2 2





