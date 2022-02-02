\documentclass{article}


\usepackage{fancyvrb} 
\DefineVerbatimEnvironment{code}{Verbatim}{fontsize=\small}


% Adds a blank row between paragraphs
\setlength{\parindent}{0pt}
\setlength{\parskip}{\baselineskip}

\begin{document}

\begin{code}

{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ConstraintKinds #-}


import qualified Prelude
import Prelude hiding ((*),(+),(/),(-))

\end{code}


\section{Definition of a vector Space}

The following quote is taken from Linear Algebra Done Right by Sheldon Axler.
In our text $F$ denotes a arbitrary field, however, 
in the book it is further specified to either be a real or a complex number.

``
The motivation for the definition of a vector space comes from properties of
addition and scalar multiplication in $F^n$: 
Addition is commutative, associative, and has an identity. 
Every element has an additive inverse. 
Scalar multiplication is associative. 
Scalar multiplication by 1 acts as expected. 
Addition and scalar multiplication are connected by distributive properties.

We will define a vector space to be a set V with an addition and a scalar
multiplication on V that satisfy the properties in the paragraph above.

\textbf{Definition: addition, scalar multiplication}.
An addition on a set V is a function that assigns an element $u+v \in V$
to each pair of elements $u, v \in V$.

A scalar multiplication on a set V is a function that assigns an 
element $\lambda v \in V$ to each $\lambda \in F$ and each $v \in V$.

Now we are ready to give the formal definition of a vector space.

\textbf{Definition: vector space}.
A vector space is a set V along with an addition on V and a scalar multi-
plication on $V$ such that the following properties hold:

\textbf{commutativity}
$u + v = v + u$ for all $u, v \in V$;

\textbf{associativity}
$(u + v) + w = u + (v + w)$ and $(ab)v = a(bv)$ for all $u, v, w \in V$
and all $a, b \in F;$

\textbf{additive identity}
there exists an element $0 \in V$ such that $v + 0 = v$ for all $v \in V ;$

\textbf{additive inverse}
for every $v \in V$, there exists $w \in V$ such that $v + w = 0$;

\textbf{multiplicative identity}
$1v = v$ for all $v \in V ;$

\textbf{distributive properties}
$a(u + v) = au + av$ and $(a + b)v = av + bv$ 
for all $a, b \in F$ and all $u, v \in V$.
''

We will later prove these properties for our datatypes in Agda, 
but first we will simply define a vector space as a class in Haskell.
To begin with we need a couple of help classes.
From above we know that $V$ forms a group over addition,
so we need a way to describe that in our code.

\begin{code}
class AddGroup a where
    (+)  :: a -> a -> a
    (-)  :: a -> a -> a
    zero :: a
\end{code}

We also need a way to represent a field.

\begin{code}
class AddGroup a => Field a where
    (*) :: a -> a -> a
    (/) :: a -> a -> a
    one :: a
\end{code}

From this we can define a vector space class.

\begin{code}
class (Field s, AddGroup v) => VectorSpace v s where
    (£) :: s -> v -> v
\end{code}

The interpretation of this is that if we have
a field s, an additive group v and can define a operation £ --
that follows the above rules -- we have a vector space.

Note: Often in mathematical texts the scalar multiplication operation is implicit.
Since we cannot do that in Haskell, 
we have arbitrarily chosen £ as our scalar multiplication operator.


Now that we have a class to represent vector spaces we may start creating instances.
As a simple example we can define a vector space over a field onto itself.

\begin{code}

instance Field a => VectorSpace a a where
    (£) = (*)

\end{code}

If we now make Double an instance of field we get our first concrete vector space.

\begin{code}

instance AddGroup Double where
    (+)  = (Prelude.+)
    (-)  = (Prelude.-)
    zero = 0

instance Field Double where
    (*) = (Prelude.*)
    (/) = (Prelude./)
    one = 1

\end{code}

Since this file is written in literate Haskell we can load it in GHCI 
and try out our vector space interactively. 


Another example that may come naturally is to define a vector 
as a list containing some Field a.

\begin{code}

instance Field a => AddGroup [a] where
    (+)  = zipWith (+)
    (-)  = zipWith (-)
    zero = error "Cannot make zero element"

\end{code}


We will notice however that there is no clear way to define zero 
such that the additive inverse holds. 
In essence we need a way to make sure the list are of the same length.
So instead we will define a vector space over functions 
$(Field\\ a) \Rightarrow g \rightarrow a $.
This works since the size is bound by its domain $g$.

\begin{code}

instance Field a => AddGroup (g -> a) where
    f + g = \x -> f x + g x
    f - g = \x -> f x - g x
    zero  = \_ -> zero

instance Field a => VectorSpace (g -> a) a where
    s £ f = \x -> s £ f x

\end{code}



\section{Definition of linear map}

Taken from Linear Algebra Done Right by Sheldon Axler

$V$ and $W$ are vector spaces over a field $F$.

A linear map from V to W is a function $T : V \to W$ 
with the following properties:

\textbf{additivity}
$T(u + v) = T u + T v$ for all $u, v \in V $;

\textbf{homogeneity}
$T(\lambda v) = \lambda T v$ for all $\lambda \in F$ and all $v \in V$.


\end{document}

