# Titel

Matematikens domänspecifika språk (DSLsofMath) för andra kurser

## Bakgrund:

DSLsofMath [1,2] är namnet på ett pedagogiskt projekt som lett till
en ny valfri kurs i årskurs 2-3 riktad till datavetare och matematiker
på Chalmers och GU. Kursen presenterar klassiska matematiska ämnen
från ett datavetenskapligt perspektiv: genom att specificera de
introducerade begreppen, vara uppmärksam på syntax och typer, och
slutligen genom att bygga domänspecifika språk for vissa matematiska
områden. (Exempelvis linjär algebra, Laplace-transform, potensserier,
derivator.)

Inspirerat av detta har flera studentgrupper genomfört kandidatarbetesprojekt under de senaste åren med följande resultat:
+ 2016: Programmering som undervisningsverktyg för Transformer, signaler och system - Utvecklingen av läromaterialet TSS med DSL
+ 2018: Ett komplementerande läromaterial för datastudenter som lär sig fysik - Läromaterialet Learn You a Physics for Great Good!
+ 2020: A Computer Science Approach to Teaching Control Theory - Developing Learning Material Using Domain-Specific Languages

## Projektbeskriving:

Det här kandidatprojektet går ut på att ta fram DSLsofMath-inspirerat
kompletterande material för andra närliggande kurser som exempelvis
* Matematisk statistik och diskret matematik, eller
* Grundläggande datorteknik, eller
* Linjär algebra, eller
* Datastrukturer och algoritmer, eller
* andra kurser som ni känner skulle må bra av mer fokus på syntax, typer och funktioner.

Implementationsspråk är Haskell och Agda och målet är dels att
förbättra förståelsen hos projektmedlemmarna av de kurser och ämnen
som väljs och dels att ge framtida studenter mer material att arbeta
med. Materialet som utvecklas skall finnas öppet tillgängligt på
github.

Efter tre omgångar med huvudfokus på lärmaterial är fokus i år mer inriktat mot korrekthet: DSL, typer, specifikation, test, bevis.

Att göra ("produkt"):
* Designa och implementera (ett par) DSL för det valda området
* Specificera lagar som bör gälla i Haskell eller Agda
* Testa de lagar som kan testas med QuickCheck
* Bevisa någon eller några lagar i Agda

Samt dokumentation i form av kandidatarbetesrapport mm.

## Litteraturförslag:

* [1] https://github.com/DSLsofMath/DSLsofMath
* [2] http://www.cse.chalmers.se/~patrikj/papers/Ionescu_Jansson_DSLsofMath_TFPIE_2015_paper_preprint.pdf
* [3] http://www.cse.chalmers.se/~patrikj/papers/Janssonetal_DSLsofMathCourseExamplesResults_preprint_2018-08-17.pdf

## Målgrupp:

DV, D, IT, TM

## Särskilda förkunskaper:

Funktionell programmering (Haskell) och kursen DSLsofMath eller gott om matematik (TM-programmet eller liknande).

(Det kan gå att ta kursen DSLsofMath parallellt med projektet, men det blir svårare.)

## Förslagslämnare:

Patrik Jansson

## Handledare:

Patrik Jansson
