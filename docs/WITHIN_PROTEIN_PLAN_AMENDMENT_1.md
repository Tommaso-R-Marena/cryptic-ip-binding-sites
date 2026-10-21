# Study Q8, amendment 1: the same design on ground truth

Dated 2026-10-20, **before any crystal structure has been fetched for this arm
and before any annotated site has been passed to `describe_architecture`**.
Written while the candidate arm's pocket detection is still running in CI and
before any of its output has been read.

## Why the candidate arm needs this one

Q8 asks whether a **candidate** pocket's lining architecture differs from its
own protein's other pockets of the same size. Whatever it returns, one
objection survives:

> If the answer is null, is that because candidates are ordinary — or because
> the instrument cannot tell a real inositol-phosphate site from an ordinary
> pocket either?

The candidate arm cannot answer that about itself. Nothing in the architecture
line has ever been run on an **annotated** inositol-hexakisphosphate site: the
real-lining study used candidate pockets, the paired study used subsets of
candidate linings, and every power number comes from synthetic barrels.

So this arm runs the identical design on ground truth, and it is a **positive
control in the strict sense**: it can fail, and its failure would reinterpret
the candidate arm rather than excuse it.

* If annotated sites separate from their proteins' other pockets, the
  instrument registers real sites, and a null in the candidate arm is then a
  statement about candidates.
* If annotated sites do **not** separate, the instrument cannot see this class
  of site at all, and the candidate arm's result — null or not — says nothing
  about candidates. That conclusion would be the most consequential thing in
  the architecture line, and it would be recorded as such.

## The copies

`results/redocking/census.csv`, filtered exactly as the census itself defines:
`comp_id == IHP`, `status == eligible`, `complete`, not a symmetry contact, and
`selected` — the census's own representative copy per site. That is **137
copies across 135 PDB entries and 17 strict homology groups**, which clears the
standing five-group floor. The counts are read from the committed census, not
from anything fetched.

## What will be computed

1. **Structures.** Fetch the 135 entries with the project's existing fetcher,
   the same one the burial survey uses.
2. **The annotated site.** Protein residues with any heavy atom within
   `CONTACT_ANGSTROMS = 4.5` of any heavy atom of that IHP copy, in the copy's
   own chain and in every other chain. The distance is the one the census
   already uses for protein contacts and is not tuned here.
3. **Pockets.** `ProteinAnalyzer.detect_pockets()` on the same structure, with
   every pocket's lining residues kept — identical to the candidate arm.
4. **The site pocket.** The pocket whose lining best **covers the annotated
   site**, by recall `|pocket ∩ site| / |site|`. Recall rather than Jaccard,
   because a pocket may legitimately be larger than a ligand's contact shell,
   and a copy whose best recall is below `MIN_RECALL = 0.3` is **excluded with
   a counted reason**: no detected pocket contains the site, so there is
   nothing to compare. The recall distribution is reported whatever it shows.
5. **One chain, or nothing.** A crystal structure is often several chains, and
   the detector's pocket residue lists are **chain-blind residue numbers** — a
   homodimer's two residue 50s are indistinguishable in them. The architecture
   statistic reads positions along one chain, so a lining pooled across chains
   is not a lining at all. The rule, fixed here before any structure is
   fetched: keep the **ligand-bearing chain**, and **exclude the copy** when
   another polymer chain in the same structure shares any residue number with
   it, because the pocket's residue list cannot then be assigned to a chain.
   The number of copies this removes is reported. It is a conservative rule —
   it discards homodimers, which are common — and it is preferred to a second
   lining definition, because the arms must compute the same statistic on the
   same kind of input for their comparison to mean anything.

   *(Added to this file before any structure for this arm was fetched and
   before any of its results existed; the arm had not run.)*

6. **Siblings and statistics.** Identical to the candidate arm: other pockets
   of the same structure matched on lining size to within two residues, at
   least eight and at most forty, the candidate's mid-p quantile among them per
   half, an unresolved half left unresolved, a group bootstrap over
   `homology_group_strict`, and Holm across the two halves.

Both arms are therefore the same statistic on the same kind of input, and the
only difference is whether the pocket was **selected by the screen** or
**occupied by an inositol phosphate in a crystal**.

## Decision rules, fixed now

* **Primary:** the two halves' mean quantiles against 0.5, Holm-corrected
  within this arm's own family of two. The arms are not pooled and their
  families are kept separate, so neither arm's p-values depend on the other's.
* **Secondary:** the recall distribution; the fraction of sites called a
  helical or strand face against their siblings; the same split by the census's
  `burial_class`, since a surface site and a cryptic one are not the same
  object — reported, and **not** promoted to a finding, because a stratum with
  fewer than five groups gives no evidence.
* **The comparison between arms is descriptive.** A difference between the
  candidate arm's quantile and this arm's is reported with both intervals and
  is never itself tested: the two arms differ in structure source, in detection
  conditions and in what selected the pocket, and no single test can carry
  that.
* **Power is stated.** The minimum detectable difference comes from the same
  bootstrap, and the transfer ceiling stands over this arm too.

## What this still cannot conclude

* A crystal structure is a crystal structure. These sites are occupied by an
  inositol phosphate in one deposited conformation, which is a stronger label
  than anything else in this project and still not a binding constant.
* Many of these sites are at the surface (99 of 137 by the census's own class).
  The project's target is a **buried** site, so a positive result here is not
  automatically a result about cryptic sites, and the burial split is reported
  for that reason rather than as a finding.
* A sibling pocket is not a non-site, in this arm either.
