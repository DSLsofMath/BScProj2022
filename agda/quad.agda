open import Data.Nat hiding ( _+_; _*_; _⊔_; zero ) 
open import Algebra hiding (Zero) 
open import Level using (Level; _⊔_)
open import Data.Product



module quad { c ℓ } ( ring : Ring c ℓ ) where
open Ring ring


data Nat4 : Set where
  One : Nat4
  Suc : Nat4 -> Nat4 


variable
  a b d : Level
  A : Set a
  B : Set b
  C : Set c
  m n : Nat4


-- Indu:ctive definition of a Quad matrix
data Quad ( A : Set a ) : Nat4 → Set a where
  Zero : Quad A n 
  Scalar : A → Quad A One
  Mtx : ( nw ne sw se : Quad A n ) → Quad A (Suc n)


-- Zero with explicit size
zeroQ : ( n : Nat4) → Quad A n 
zeroQ n = Zero

oneQ : ∀ { n } → Quad Carrier n
oneQ {One} = Scalar 1#
oneQ {Suc n} = Mtx oneQ Zero
                   Zero oneQ

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

_=q_ : Quad Carrier n -> Quad Carrier n -> Set ( c ⊔ ℓ)
_=q_ = EqQ _≈_

-- Functions

mapq : ∀ { A B : Set } → ( A → B ) → Quad A n -> Quad B n
mapq f Zero = Zero
mapq f (Scalar x) = Scalar (f x) 
mapq f (Mtx nw ne sw se) = Mtx (mapq f nw) (mapq f ne) (mapq f sw) (mapq f se)


-- Addition on Quads
_+q_ : ( x y : Quad Carrier n) → Quad Carrier n
Zero +q y = y
Scalar x +q Zero = Scalar x
Scalar x +q Scalar x₁ = Scalar (x + x₁)
Mtx x x₁ x₂ x₃ +q Zero = Mtx x x₁ x₂ x₃
Mtx x x₁ x₂ x₃ +q Mtx y y₁ y₂ y₃ = Mtx (x +q y) (x₁ +q y₁) (x₂ +q y₂) (x₃ +q y₃)


-- Multiplication

data Q' ( A : Set a) : Set a where
  Mtx' : ( nw ne sw se : A ) → Q' A

Q'toQuad : Q' (Quad A n) → Quad A (Suc n)
-- Q'toQuad (Mtx' Zero Zero Zero Zero) = Zero
Q'toQuad (Mtx' nw ne sw se) = Mtx nw ne sw se

zipQ' : ( A → B → C ) → Q' A → Q' B → Q' C
zipQ' _op_ (Mtx' nw ne sw se) (Mtx' nw₁ ne₁ sw₁ se₁) = Mtx' (nw op nw₁) (ne op ne₁)
                                                            (sw op sw₁) (se op se₁)

colExchange : Q' A → Q' A
colExchange (Mtx' nw ne sw se) = Mtx' ne nw se sw

prmDiagSqsh : Q' A → Q' A
prmDiagSqsh (Mtx' nw ne sw se) = Mtx' nw se nw se

offDiagSqsh : Q' A → Q' A
offDiagSqsh (Mtx' nw ne sw se) = Mtx' sw ne sw se


_*q_ : ( x y : Quad Carrier n) → Quad Carrier n
Zero *q y = Zero
Scalar x *q Zero = Zero
Scalar x *q Scalar x₁ = Scalar (x * x₁)
Mtx x x₁ x₂ x₃ *q Zero = Zero
Mtx x x₁ x₂ x₃ *q Mtx y y₁ y₂ y₃ =
      let X = Mtx' x x₁ x₂ x₃
          Y = Mtx' y y₁ y₂ y₃
      in Q'toQuad (zipQ' _+q_ (zipQ' _*q_ (colExchange X) (offDiagSqsh Y))
                              (zipQ' _*q_              X  (prmDiagSqsh Y)))



-- Example Quad

q1 : Quad ℕ (Suc One)
q1 = Mtx (Scalar 2) (Scalar 7) Zero Zero 

q2 : Quad ℕ (Suc (Suc One))
q2 = Mtx q1 q1 Zero (Mtx (Scalar 1) Zero Zero (Scalar 1))



-- Proofs

quad-refl : { a : Quad Carrier n } → a =q a
quad-refl { n } { Zero } = eq-Zero
quad-refl { n } { (Scalar x) } = eq-Scalar refl
quad-refl { n } { (Mtx a a₁ a₂ a₃) } = eq-Mtx quad-refl quad-refl quad-refl quad-refl

quad-sym : { a b : Quad Carrier n } → a =q b → b =q a
quad-sym eq-Zero = eq-Zero
quad-sym (eq-Scalar x∼y) = eq-Scalar (sym x∼y)
quad-sym (eq-Mtx x x₁ x₂ x₃) = eq-Mtx (quad-sym x) (quad-sym x₁) (quad-sym x₂) (quad-sym x₃)


quad-zerol : ∀ { n } ( x : Quad Carrier n ) → ( Zero +q x ) =q x
quad-zerol x = quad-refl

quad-zeror : ∀ { n } ( x : Quad Carrier n ) → ( x +q Zero ) =q x
quad-zeror Zero = quad-refl
quad-zeror (Scalar x) = quad-refl
quad-zeror (Mtx x x₁ x₂ x₃) = quad-refl

-- Need a proof that Zero is left and right +q-identity
-- Should be some kind of (quad-zerol, quad-zeror)



quad-comm : ( a b : Quad Carrier n ) → ( a +q b) =q ( b +q a )
quad-comm Zero b = quad-sym (quad-zeror b)
quad-comm (Scalar x) Zero = quad-refl
quad-comm (Scalar x) (Scalar x₁) = eq-Scalar (+-comm x x₁)
quad-comm (Mtx a a₁ a₂ a₃) Zero = quad-refl 
quad-comm (Mtx a a₁ a₂ a₃) (Mtx b b₁ b₂ b₃) = eq-Mtx (quad-comm a b) (quad-comm a₁ b₁)
                                                     (quad-comm a₂ b₂) (quad-comm a₃ b₃)


quad-assoc : ( a b c : Quad Carrier n ) → ((a +q b ) +q c ) =q ( a +q ( b +q c )) 
quad-assoc Zero b c = quad-refl 
quad-assoc (Scalar x) Zero c = quad-refl 
quad-assoc (Scalar x) (Scalar x₁) Zero = quad-refl
quad-assoc (Scalar x) (Scalar x₁) (Scalar x₂) = eq-Scalar (+-assoc x x₁ x₂)
quad-assoc (Mtx a a₁ a₂ a₃) Zero c = quad-refl
quad-assoc (Mtx a a₁ a₂ a₃) (Mtx b b₁ b₂ b₃) Zero = quad-refl
quad-assoc (Mtx a a₁ a₂ a₃) (Mtx b b₁ b₂ b₃) (Mtx c c₁ c₂ c₃) = eq-Mtx (quad-assoc a  b  c ) (quad-assoc a₁ b₁ c₁)
                                                                       (quad-assoc a₂ b₂ c₂) (quad-assoc a₃ b₃ c₃)



quad-*-zerol : ( x : Quad Carrier n ) → ( Zero *q x ) =q Zero
quad-*-zerol x = eq-Zero

quad-*-zeror : ( x : Quad Carrier n ) → ( x *q Zero ) =q Zero
quad-*-zeror Zero = eq-Zero
quad-*-zeror (Scalar x) = eq-Zero
quad-*-zeror (Mtx x x₁ x₂ x₃) = eq-Zero


quad-onel : ( x : Quad Carrier n ) → ( oneQ *q x ) =q x

quad-onel-lemma : ( x : Quad Carrier n ) → ((oneQ *q x) +q Zero) =q x
quad-onel-lemma {One} Zero = eq-Zero
quad-onel-lemma {Suc n} Zero = eq-Zero
quad-onel-lemma (Scalar x) = eq-Scalar (*-identityˡ x)
quad-onel-lemma (Mtx x x₁ x₂ x₃) = eq-Mtx (quad-onel x) (quad-onel-lemma x₁) (quad-onel-lemma x₂) (quad-onel x₃)

quad-onel Zero = quad-*-zeror oneQ
quad-onel (Scalar x) = eq-Scalar (*-identityˡ x)
quad-onel (Mtx x x₁ x₂ x₃) = eq-Mtx (quad-onel x) (quad-onel-lemma x₁) (quad-onel-lemma x₂) (quad-onel x₃)

