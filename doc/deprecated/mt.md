
# Model theory

WORK IN PROGRESS.


Understanding a system requires observation and reasoning;
collaborative understanding involves collaborative observation and
reasoning, and collaboration requires communication.  When
communication itself seems to be inefficient or an impediment, model
theory can be a useful tool.  Model theory posits that the system of
communication and reasoning is separate from the system of study and
that it can be worthwhile to try to understand the relationship
between the two.

It is sentences (also called 'formulae') that are transmitted and that
are the outcome of reasoning.  We start with sentences capturing
observations that someone has made (or other kinds of proposition),
and reason about them, yielding new sentences.  Sometimes a conclusion
affects our predictions about the future of the system and can be
tested by observation.

  - The machinery of sentences and reasoning is called a 'theory'.
    When a sentence S follows logically (theoretically) from a
    collection of others sentences Γ, we say Γ _yields_ S, and write Γ
    ⊢ S.  (Γ is typically a set of premises or observations relevant
    to a system.)

  - A correspondence between the sentences of a theory and a given
    system is called an 'interpretation', because we interpret the
    sentences as expressing propositions that are about the system.
    Use the letter I to stand for an interpretation when the theory
    and system are understood.

  - An interpretation 'satisfies' the theory if logical apparatus of
    the formal language is interpreted sensibly, e.g. a sentence 'p
    and q' is interpreted to be the proposition that's true just when
    I(p) and I(q) are.

  - (Definition of 'satisfying interpretation' or 'model' ...?
    compositional?  If Γ ⊢ S then I(Γ) logically implies I(S)?  I
    respects the language's rules of deduction?)

Reasoning is connected to truth via the entailment relationship.
Write Γ _entails_ S, or write Γ ⊨ S, if I(S) holds in any model of the
system in which I(Γ) is the case.  That is, I(S) is true in all models
of Γ.

If Γ ⊨ S implies that Γ ⊢ S, then the deductive system is said to be
_complete_ - anything that is true in all models can be inferred in
the logical system.

I is faithful/satisfying: Γ ⊢ S implies that I(Γ) implies I(S) ?]

Model theory gives us a structure for diagnosing and correcting
communication failures.  No two collaborators will interpret abstract
terminology in the same way; they will have their own interpretations
of the words and their own private, perhaps unspoken, knowledge of
systems under study.  But if they can stick to a formal language, they
can at least agree on what a _possible_ interpretation might be, and
agree on what conclusions are _entailed_ by agreed upon premises.

## Things to notice.

_Propositions_: It's important to distinguish sentences from
propositions.  Propositions are the interpretations or meanings of
sentences, not the sentences themselves, but they can be _expressed_
as sentences.  For example, the same proposition could be expressed in
different languages.  Saying that a sentence is true or false is
shorthand for saying that its interpretation as a proposition is true
or false.

## Application: the word 'taxon'

...

## Application: RCC-5

...

