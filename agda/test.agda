open import Data.List
open import Data.Nat using (ℕ; zero; suc)
open import Data.Product
open import Data.Vec.Base using (foldr)
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
      m n : ℕ

-- Vector representation

    module Vector where
        -- Inductive definition of a vector
        data Vector (A : Set a) : (n : ℕ) → Set a where
          [] : Vector A zero
          _::_ :  A → Vector A n → Vector A (suc n)

        infixr 5 _::_

        vecLength : Vector A n -> ℕ
        vecLength {n = n} v = n

        -- Matrices are defined as vector of vectors
        Matrix : (A : Set a) → (m n : ℕ) → Set a
        Matrix A m n = Vector (Vector A n) m

        matLength : Matrix A m n -> ℕ × ℕ
        matLength {m = m} {n = n} mat  = m , n

        -- Some examples
        v1 : Vector ℕ 4
        v1 = 1 :: 3 :: 4 :: 5 :: []

        m1 : Matrix ℕ 2 2
        m1 =  (1 :: 2 :: []) :: (1 :: 2 :: []) :: []

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

        zeroVec : Vector ℕ n
        zeroVec = replicateV 0

        -- Pointwise equality on vectors (lifting _∼_ from elements to vectors)
        data Pointwise {A : Set a} (_∼_ : A -> A -> Set c) :
                       ∀ {m n} (xs : Vector A m) (ys : Vector A n) → Set (a ⊔ c)
                       where
          eq-[]  : Pointwise _∼_ [] []
          eq-::   : ∀ {m n x y} {xs : Vector A m} {ys : Vector A n}
                  (x∼y : x ∼ y) (xs∼ys : Pointwise _∼_ xs ys) →
                  Pointwise _∼_ (x :: xs) (y :: ys)

-- Operations

    module Operations (R : Ring c ℓ) where
        open Ring R
        open Vector

        infixr 6 _+v_
        infixr 7 _◁_
        infixr 7 _◁ₘ_
        infixr 7 _x_
        infixr 7 _⊗_

        sumV : Vector Carrier n → Carrier
        sumV {n = zero} v = 0#
        sumV {n = suc n} (x :: xs) = x + sumV {n} xs

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

        -- Add matrix
        _+m_ : Matrix Carrier m n → Matrix Carrier m n → Matrix Carrier m n
        _+m_ = zipV _+v_

        altSumV : Vector Carrier n → Carrier
        altSumV v = {!!} -- alternating sum: multiply every second term by minus one

        submatrices : Matrix Carrier m (suc n) -> Vector (Matrix Carrier m n) (suc n)
        submatrices = {!!}

        -- Determinant
        det : Matrix Carrier m m → Carrier
        det [] = 1#
        det (v :: m) = altSumV (zipV _*_ v (mapV det (submatrices m)))

        -- Equality on our vectors is a lifted version of the
        -- underlying equality of the ring of components.
        eqVec : Vector Carrier n -> Vector Carrier m -> Set (c ⊔ ℓ)
        eqVec = Pointwise _≈_

        vectorAddComm :  (v1 v2 : Vector Carrier n) ->
                         eqVec (v1 +v v2) (v2 +v v1)
        vectorAddComm v1 v2 = {!!}

        -- Vector addition is associative (statement, and inductive proof)
        vectorAddAssoc :  (v1 v2 v3 : Vector Carrier n) ->
                          eqVec ((v1 +v v2) +v v3) (v1 +v (v2 +v v3))
        vectorAddAssoc [] [] [] = eq-[]
        vectorAddAssoc (x1 :: v1) (x2 :: v2) (x3 :: v3)
          = eq-:: (+-assoc x1 x2 x3) (vectorAddAssoc v1 v2 v3)
