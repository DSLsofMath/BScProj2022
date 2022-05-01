open import Data.Nat hiding ( _+_; _⊔_; zero ) 
open import Algebra hiding (Zero) 
open import Level using (Level; _⊔_)
open import Data.Product



module quad { c ℓ } ( ring : Ring c ℓ ) where
open Ring ring


data Nat4 : Set where
  One : Nat4
  Suc : Nat4 -> Nat4 


variable
  a d : Level
  A : Set a
  m n : Nat4


-- Indu:ctive definition of a Quad matrix
data Quad ( A : Set a ) : Nat4 → Set a where
  Zero : Quad A n 
  Scalar : A → Quad A One
  Mtx : ( nw ne sw se : Quad A n ) → Quad A (Suc n)


-- Zero with explicit size
zeroQ : ( n : Nat4) → Quad A n 
zeroQ n = Zero


-- Equality

-- eq-Zero is only valid when m ≡ n
data EqQ {A : Set a} (_∼_ : A -> A -> Set d) :
                       (xs : Quad A m) (ys : Quad A n) → Set (a ⊔ d)
                       where
                       
     eq-Zero : ∀ { n } → EqQ _∼_ (zeroQ n) (zeroQ n)

     eq-Scalar : { x y : A } ( x∼y : x ∼ y ) → EqQ _∼_ (Scalar x) (Scalar y)
          
     eq-Mtx  : { nw1 ne1 sw1 se1 : Quad A m } { nw2 ne2 sw2 se2 : Quad A n }
               ( nw : EqQ _∼_ nw1 nw2 ) ( ne : EqQ _∼_ ne1 ne2 )
               ( sw : EqQ _∼_ sw1 sw2 ) ( se : EqQ _∼_ se1 se2 ) →
               EqQ _∼_ (Mtx nw1 ne1 sw1 se1) (Mtx nw2 ne2 sw2 se2) 



eqQuad : Quad Carrier n -> Quad Carrier n -> Set ( c ⊔ ℓ)
eqQuad = EqQ _≈_

-- Functions

mapq : ∀ { A B : Set } → ( A → B ) → Quad A n -> Quad B n
mapq f Zero = Zero
mapq f (Scalar x) = Scalar (f x) 
mapq f (Mtx nw ne sw se) = Mtx (mapq f nw) (mapq f ne) (mapq f sw) (mapq f se)


-- Addition on Quads

_+q_ : ( x y : Quad Carrier n) → Quad Carrier n
Zero +q y = y
x +q Zero = x
Scalar x +q Scalar x₁ = Scalar (x + x₁)
Mtx x x₁ x₂ x₃ +q Mtx y y₁ y₂ y₃ = Mtx (x +q y) (x₁ +q y₁) (x₂ +q y₂) (x₃ +q y₃)


-- Example Quad

q1 : Quad ℕ (Suc One)
q1 = Mtx (Scalar 2) (Scalar 7) Zero Zero 

q2 : Quad ℕ (Suc (Suc One))
q2 = Mtx q1 q1 Zero (Mtx (Scalar 1) Zero Zero (Scalar 1))



-- Proofs

quad-refl : ( a : Quad Carrier n ) → eqQuad a a
quad-refl Zero = eq-Zero
quad-refl (Scalar x) = eq-Scalar refl
quad-refl (Mtx a a₁ a₂ a₃) = eq-Mtx (quad-refl a) (quad-refl a₁) (quad-refl a₂) (quad-refl a₃)

quad-sym : { a b : Quad Carrier n } → eqQuad a b → eqQuad b a
quad-sym eq-Zero = eq-Zero
quad-sym (eq-Scalar x∼y) = eq-Scalar (sym x∼y)
quad-sym (eq-Mtx x x₁ x₂ x₃) = eq-Mtx (quad-sym x) (quad-sym x₁) (quad-sym x₂) (quad-sym x₃)



quad-zerol : ∀ { n } ( x : Quad Carrier n ) → eqQuad ( Zero +q x ) x
quad-zerol Zero = eq-Zero
quad-zerol (Scalar x) = eq-Scalar refl
quad-zerol (Mtx x x₁ x₂ x₃) = eq-Mtx (quad-refl x) (quad-refl x₁) (quad-refl x₂) (quad-refl x₃)

quad-zeror : ∀ { n } ( x : Quad Carrier n ) → eqQuad ( x +q Zero ) x
quad-zeror Zero = eq-Zero
quad-zeror (Scalar x) = eq-Scalar refl
quad-zeror (Mtx x x₁ x₂ x₃) = eq-Mtx (quad-refl x) (quad-refl x₁) (quad-refl x₂) (quad-refl x₃)

quad-zero : ( x : Quad Carrier n ) → eqQuad ( Zero +q x ) ( x +q Zero )
quad-zero x = quad-sym (quad-zeror x)


quad-comm : ( a b : Quad Carrier n ) → eqQuad ( a +q b) ( b +q a )
quad-comm Zero b = quad-zero b
quad-comm a Zero = quad-zeror a
quad-comm (Scalar x) (Scalar x₁) = eq-Scalar (+-comm x x₁)
quad-comm (Mtx a a₁ a₂ a₃) (Mtx b b₁ b₂ b₃) = eq-Mtx (quad-comm a b) (quad-comm a₁ b₁) (quad-comm a₂ b₂) (quad-comm a₃ b₃)

