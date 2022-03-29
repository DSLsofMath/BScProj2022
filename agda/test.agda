open import Data.List
open import Data.Nat
open import Data.Bool

variable A : Set
         m n : ℕ

data Vector (A : Set) : ℕ → Set where
  [] : Vector A zero
  _::_ : {n : ℕ} → A → Vector A n → Vector A (suc n)
  
infixr 5 _::_

Matrix : (A : Set ) (m n : ℕ) → Set
Matrix A m n = Vector (Vector A n) m


v1 : Vector ℕ 2
v1 = 1 :: 3 :: []

mat1 = Matrix Bool 2 2



