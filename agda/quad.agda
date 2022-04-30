open import Data.Nat
open import Algebra hiding (Zero)


data Nat4 : Set where
  One : Nat4
  Suc : Nat4 -> Nat4

-- Inductive definition of a Quad matrix
data Quad ( A : Set ) : Nat4 → Set where
  Zero : ∀ { n } → Quad A n 
  Scalar : A → Quad A One
  Mtx : ∀ { n } ( nw ne sw se : Quad A n ) → Quad A (Suc n)


-- Functions
R = ℕ -- Should be changed to a ring

variable
  m n : Nat4

mapq : ∀ { A B : Set } → ( A → B ) → Quad A n -> Quad B n
mapq f Zero = Zero
mapq f (Scalar x) = Scalar (f x) 
mapq f (Mtx nw ne sw se) = Mtx (mapq f nw) (mapq f ne) (mapq f sw) (mapq f se)

-- Addition on Quads
_+q_ : ( x y : Quad R n) → Quad R n
Zero +q y = y
x +q Zero = x
Scalar x +q Scalar y = Scalar (x + y)
Mtx nw1 ne1 sw1 se1 +q Mtx nw2 ne2 sw2 se2 = Mtx (nw1 +q nw2) (ne1 +q ne2)
                                                 (sw1 +q sw2) (se1 +q se2)




--Example Quad

q1 : Quad ℕ (Suc One)
q1 = Mtx (Scalar 2) (Scalar 7) Zero Zero 

q2 : Quad ℕ (Suc (Suc One))
q2 = Mtx q1 q1 Zero (Mtx (Scalar 1) Zero Zero (Scalar 1))
