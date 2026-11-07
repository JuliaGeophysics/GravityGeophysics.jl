# Choosing between them

Both inversions read the same files and return the same kind of result, so the cheapest way
to decide is to run both and compare. That is exactly what the shipped demos do. Still, the
two have different characters, and it is worth knowing which one your problem suits.

## The short version

| Your situation | Try first |
|---|---|
| blocky, compact targets — a dyke, an ore body, a shallow block | **Gauss–Newton** |
| gradational structure — a basin, a weathering profile, an alteration halo | **INR** |
| dense, regular station grid | **Gauss–Newton** |
| sparse or irregular stations, gaps, no-access ground | **INR** |
| you need an answer in seconds | **Gauss–Newton** |
| depth of the target is the thing you care about | **INR** |
| you want parameters you can explain to a reviewer | **Gauss–Newton** |

## Why they differ

Gauss–Newton treats every cell as an independent unknown and then has to *argue* them back
into agreement with a smoothing penalty. The INR never lets them disagree in the first place
— a single small network produces all of them, so nearby cells get similar values because
nearby inputs give similar outputs.

That has two consequences you can see in the results.

**Edges.** Gauss–Newton can put a sharp contact exactly where the data demand one. The INR
rounds it off, because a smooth function is what it is good at.

**Empty ground.** Where there are no stations, Gauss–Newton falls back on its smoothing term,
which mostly pulls towards the prior. The INR fills the gap by continuing the function it has
learned elsewhere — which is usually closer to what a geologist would have drawn.

There is also a bias worth naming: with no term that prices the *volume* of a cell,
Gauss–Newton's sensitivity weighting makes it comfortable putting mass deep and wide, while
the INR has no such pull and tends to keep structure where the data put it. In practice that
shows up as depth.

## What the demos show

Two synthetic surveys ship with the package, one built for each case. Both are scored over
the survey core, against the true model.

**`demo1`** — two uniform ±200 kg/m³ blocks under a dense 100 m grid. Blocky truth, good
coverage, which is the classical setting:

| | RMS | iterations | correlation |
|---|---|---|---|
| Gauss–Newton | 0.99 | 1 | **0.66** |
| INR | 1.00 | 281 | 0.63 |

**`demo2`** — a basin whose fill grades continuously with depth, beside a block, under
150 scattered stations with a no-access gap:

| | RMS | iterations | correlation | deficit | its centroid |
|---|---|---|---|---|---|
| true | | | | −3.05 Mt | 52 m |
| Gauss–Newton | 0.95 | 3 | 0.61 | −6.34 Mt | 189 m |
| INR | 0.99 | 257 | **0.81** | **−2.89 Mt** | **54 m** |

The second case is the interesting one. Gauss–Newton reaches the target misfit perfectly
well — it simply does so with twice the mass it needs, buried more than three times too deep.
The INR lands the basin within 2 m of the right depth. And this is not a matter of tuning:
across a sweep of stronger smallness, heavier smoothing, depth weighting and no weighting,
Gauss–Newton's correlation stays between 0.41 and 0.65 and its basin never comes shallower
than 149 m.

## Practical advice

1. **Run Gauss–Newton first.** It is fast, and if your target is blocky and your coverage
   good, you may not need anything else.
2. **Look at the depths.** If the result sits suspiciously deep, or mass is smeared below
   where you expect it, try the INR.
3. **Run both and compare.** Where they agree, you can believe it. Where they disagree, you
   have found something the data do not actually constrain — which is worth knowing before
   you put it in a report.
4. **Change the INR's seed and re-run.** If the answer moves much, you are asking for more
   detail than the data support.

Remember that neither method can overcome gravity's fundamental ambiguity. A deeper, denser,
larger body and a shallower, lighter, smaller one can produce the same anomaly. Both demos
overshoot the true mass, because the data genuinely do not pin it down. Bounds, prior
information and geological sense do more for that than the choice of algorithm.
