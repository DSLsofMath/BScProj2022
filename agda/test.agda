open import Data.List
open import Data.Nat using (ℕ; zero; suc)
open import Data.Product
open import Data.Bool
open import Function using (id; _∘_)
open import Algebra
open import Level using (Level; _⊔_)

module test where

    variable
      a b c ℓ : Level
      A : Set a
      B : Set b
      C : Set c
      m n o : ℕ

-- Vector representation

    module Vector where
        -- Inductive definition of a vector
        data Vector (A : Set a) : (n : ℕ) → Set a where
          [] : Vector A zero
          _::_ :  A → Vector A n → Vector A (suc n)

        infixr 5 _::_

        vecLength : Vector A n → ℕ
        vecLength {n = n} v = n

        headV : Vector A (suc n) → A
        headV (x :: _) = x

        tailV : Vector A (suc n) → Vector A n
        tailV (_ :: xs) = xs

        -- Matrices are defined as vector of vectors
        Matrix : (A : Set a) → (m n : ℕ) → Set a
        Matrix A m n = Vector (Vector A n) m

        matLength : Matrix A m n → ℕ × ℕ
        matLength {m = m} {n = n} mat  = m , n

        -- Some examples
        v1 : Vector ℕ 4
        v1 = 1 :: 3 :: 4 :: 5 :: []

        m1 : Matrix ℕ 2 2
        m1 =  (1 :: 2 :: []) :: (3 :: 4 :: []) :: []

        -- Some standard functions for working with vectors
        zipV : (A → B → C) → (Vector A n → Vector B n → Vector C n)
        zipV f [] [] = []
        zipV f (x :: xs) (y :: ys) = f x y :: zipV f xs ys

        mapV : (A → B) → Vector A n → Vector B n
        mapV f [] = []
        mapV f (x :: v) = f x :: mapV f v

        replicateV : A → Vector A n
        replicateV {n = zero} x = []
        replicateV {n = suc n} x = x :: replicateV x

        transpose : ∀{A m n} → Matrix A m n → Matrix A n m
        transpose [] = replicateV []
        transpose (x :: xs) = ((zipV _::_) x) (transpose xs)

        -- Pointwise equality on vectors (lifting _∼_ from elements to vectors)
        data EqV {A : Set a} (_∼_ : A → A → Set c) :
                       ∀ {m n} (xs : Vector A m) (ys : Vector A n) → Set (a ⊔ c)
                       where
          eq-[]  : EqV _∼_ [] []
          eq-::  : ∀ {m n x y} {xs : Vector A m} {ys : Vector A n}
                  (x∼y : x ∼ y) (xs∼ys : EqV _∼_ xs ys) →
                  EqV _∼_ (x :: xs) (y :: ys)

-- Operations

    module Operations (R : Ring c ℓ) where
        open Ring R
        open Vector

        infixr 6 _+v_
        infixr 7 _◁_
        infixr 7 _◁ₘ_
        infixr 7 _x_
        infixr 7 _⊗_
        infixr 7 _*m_ 

        sumV : Vector Carrier n → Carrier
        sumV {n = zero} v = 0#
        sumV {n = suc n} (x :: xs) = x + sumV {n} xs

        0v : Vector Carrier n
        0v = replicateV 0#

        -- Vector addition
        _+v_ : Vector Carrier n → Vector Carrier n  → Vector Carrier n
        _+v_ = zipV _+_

        -- Scale vector
        _◁_ : Carrier → Vector Carrier n → Vector Carrier n
        c ◁ v = mapV (c *_) v

        -- Dot product
        _•_ : Vector Carrier n → Vector Carrier n → Carrier
        u • v = sumV (zipV _*_ u v)

        -- Cross
        _x_ : Vector Carrier 3 → Vector Carrier 3 → Vector Carrier 3
        (v1 :: v2 :: v3 :: []) x (u1 :: u2 :: u3 :: []) = (v2 * u3 + -(v3 * u2) ::
                                                           v3 * u1 + -(v1 * u3) ::
                                                           v1 * u2 + -(v2 * u1) :: [])

        -- Outer product
        _⊗_ : Vector Carrier m → Vector Carrier n → Matrix Carrier m n
        _⊗_ [] ys = []
        _⊗_ (x :: xs) ys = mapV (x *_) ys :: xs ⊗ ys

        -- Scale matrix
        _◁ₘ_ : Carrier → Matrix Carrier m n → Matrix Carrier m n
        c ◁ₘ m = mapV (c ◁_) m

        -- Add matrix/matrix
        _+m_ : Matrix Carrier m n → Matrix Carrier m n → Matrix Carrier m n
        _+m_ = zipV _+v_

        -- Mul matrix/vector
        _*mv_ : Matrix Carrier m n → Vector Carrier n → Vector Carrier m
        m *mv v = mapV (_• v) m

        -- Mul matrix/matrix
        _*m_ : Matrix Carrier m n → Matrix Carrier n o → Matrix Carrier m o
        --m1 *m m2 = mapV (_*mv (transpose m2)) m1

        Sign = Carrier
        altSumVHelp : Sign → Vector Carrier n → Carrier
        altSumVHelp s [] = 0#
        altSumVHelp s (x :: xs) = s * x + altSumVHelp (- s) xs

        altSumV : Vector Carrier n → Carrier
        altSumV = altSumVHelp 1# -- alternating sum: multiply every second term by minus one

        submatricesStep : Matrix Carrier m (suc (suc n)) → Vector (Matrix Carrier m (suc n)) (suc (suc n))

        submatrices : Matrix Carrier m (suc n) → Vector (Matrix Carrier m n) (suc n)
        submatrices {m} {ℕ.zero} ma = replicateV [] :: []
        submatrices {m} {suc n}   ma = submatricesStep ma

        submatricesStep ma with mapV headV ma | mapV tailV ma 
        submatricesStep ma | heads | tails with submatrices tails
        submatricesStep ma | heads | tails | rec = tails :: mapV (zipV _::_ heads) rec

        -- Determinant
        det : Matrix Carrier m m → Carrier
        det [] = 1#
        det (v :: m) = altSumV (zipV _*_ v (mapV det (submatrices m)))

        module TestM (a11 a12 a21 a22 : Carrier) where
          m22 : Matrix Carrier 2 2
          m22 = (a11 :: a12 :: []) :: (a21 :: a22 :: []) :: []
          test : Carrier
          test = det m22

        -- Equality on our vectors is a lifted version of the
        -- underlying equality of the ring of components.
        _≈v_ : Vector Carrier n → Vector Carrier m → Set (c ⊔ ℓ)
        _≈v_ = EqV _≈_
        -- The equality type is basically just a vector of equality
        -- proofs between pairs of corresponding elements.

        -- Left and right proof of identity with vector addition
        vectorAddIdentityˡ : ∀ {n} (v1 : Vector Carrier n) → (0v +v v1) ≈v v1
        vectorAddIdentityˡ [] = eq-[]
        vectorAddIdentityˡ (x1 :: v1) = eq-::(+-identityˡ x1) (vectorAddIdentityˡ v1)

        vectorAddIdentityʳ : ∀ {n} (v1 : Vector Carrier n) → (v1 +v 0v) ≈v v1
        vectorAddIdentityʳ [] = eq-[]
        vectorAddIdentityʳ (x1 :: v1) = eq-::(+-identityʳ x1) (vectorAddIdentityʳ v1)
        
        -- Vector addition is commutative (statement, and inductive proof)
        vectorAddComm : ∀ {n} (v1 v2 : Vector Carrier n) →
                        (v1 +v v2) ≈v (v2 +v v1)
        vectorAddComm [] [] = eq-[]
        vectorAddComm (x1 :: v1) (x2 :: v2) =
            eq-:: (+-comm x1 x2) (vectorAddComm v1 v2)

        -- Vector addition is associative (statement, and inductive proof)
        vectorAddAssoc : ∀ {n} (v1 v2 v3 : Vector Carrier n) →
                         ((v1 +v v2) +v v3) ≈v (v1 +v (v2 +v v3))
        vectorAddAssoc [] [] [] = eq-[]
        vectorAddAssoc (x1 :: v1) (x2 :: v2) (x3 :: v3) =
            eq-:: (+-assoc x1 x2 x3) (vectorAddAssoc v1 v2 v3)
