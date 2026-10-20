# Study Q8, amendment 2: a chain-aware lining, and the look that has already happened

Dated 2026-10-20, after the ground-truth arm has been run and read.

## The disclosure, first

This amendment is written **after** seeing the ground-truth arm's result, so it
begins with what was seen, in full:

> Of 137 annotated IHP copies, 28 produced a usable site lining and 5 cleared
> the sibling floor. The topology half rested on 4 resampling groups with a
> mean quantile of 0.279; the register half on 3 groups with 0.328. Both sit
> **below 0.5**, the direction the hypothesis predicts. Both are **no
> evidence** by the project's five-group floor, with minimum detectable
> differences of 0.37 and 0.47.

Everything below is therefore **exploratory**, and a repowered arm that agrees
with that direction is not a confirmation of it. The confirmatory set is named
at the end and has not been touched.

## What went wrong, and why it is a method fault rather than a threshold

The arm lost 109 of 137 copies before the statistics began, and the dominant
cause is a rule this plan's amendment 1 introduced:

> keep the ligand-bearing chain, and **exclude the copy** when another polymer
> chain in the same structure shares any residue number with it, because the
> pocket's residue list cannot then be assigned to a chain.

That rule exists because the detector's pocket residue lists are **chain-blind
numbers**. It is conservative and it was the right call to make in advance
rather than to invent a lining definition under time pressure. But it discards
homodimers, and deposited inositol-phosphate complexes are largely multimeric,
so the rule removed most of the dataset.

**It is also unnecessary.** fpocket itself writes, for every pocket, a file of
the atoms lining it — `pocket<N>_atm.pdb` — and those records carry the chain
identifier in column 22. The project's parser reads the residue number from
columns 23 to 26 of exactly those lines and drops the chain. Reading one more
column gives a **chain-aware lining from fpocket's own definition**: no
centroid heuristic, no second notion of what a lining is, no change to any
statistic.

So the fix is to read the column that was already there, and the chain rule
then has nothing to exclude: a pocket's lining in the ligand-bearing chain is
simply the residues of that chain among the pocket's own atoms.

## What changes, and what does not

1. **Changes.** The ground-truth arm's linings — site pocket and siblings
   alike — are taken from `pocket<N>_atm.pdb` restricted to the ligand-bearing
   chain. The chain-ambiguity exclusion is removed because the ambiguity is
   gone. A pocket that contributes no residue in that chain is not a pocket of
   that chain and is skipped, with a counted reason.
2. **Does not change.** The contact-shell definition (4.5 Å), the recall floor
   (0.3), the size tolerance (2 residues), the sibling floor (8), the cap (40),
   the mid-p quantile, the unresolved rule, the group bootstrap over
   `homology_group_strict`, the five-group floor, and Holm across the two
   halves. None of those is touched, and the candidate arm is not touched at
   all.

## The confirmatory set, named now and untouched

`results/redocking/census.csv` carries **222** eligible, complete,
non-symmetry IHP copies. Amendment 1 used the **137** the census marks
`selected`, one representative per site. The remaining **85** have never been
fetched, measured or looked at in this study.

If the repowered exploratory arm shows anything, the confirmatory run is those
85 copies, with every rule above fixed and no further changes. Its result is
the one that may be reported as a test.

**What that confirmatory set is and is not.** Those copies come from the same
deposited entries as the exploratory ones — they are other copies of the same
sites, not new proteins. So it tests whether the measurement **reproduces**
across copies, grouped as always by strict homology group; it does not test
generalisation to proteins the study has not seen, and it will not be
described as if it did.

## What this still cannot conclude

* Everything amendment 1 could not conclude, unchanged.
* An exploratory arm reported after a look is weaker evidence than a
  pre-registered one, whatever it shows, and this file exists so that weakness
  is visible rather than inferred.
