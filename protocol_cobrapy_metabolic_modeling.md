# Protocol: Metabolic Modeling with COBRApy

*Prepared 2026-10-24*

## Overview

This protocol runs constraint-based metabolic modeling in Python with COBRApy, from loading a genome-scale model to designing a production strain. Every code block was tested with COBRApy 0.32.1 on the bundled *E. coli* core model (95 reactions, 72 metabolites, 137 genes); the outputs shown are from that run.

**What COBRApy does:** it represents metabolism as a stoichiometric matrix S and finds steady-state fluxes v that satisfy S·v = 0 within reaction bounds, while optimizing an objective such as biomass.

$$\max_{v} \; c^{T} v \quad \text{subject to} \quad S v = 0, \quad v_{lb} \le v \le v_{ub}$$

**Workflow:**

1. Install and load a model.
2. Inspect reactions, metabolites, genes and the medium.
3. Run flux balance analysis (FBA).
4. Change growth conditions.
5. Run parsimonious FBA and flux variability analysis.
6. Simulate gene and reaction knockouts.
7. Map production envelopes and add new pathways.
8. Save the model and record the analysis.

## Steps 1–3: Install, load and inspect

### Step 1: Install

Use Python 3.9 or newer, ideally in a virtual environment or Jupyter notebook. COBRApy installs with the open-source GLPK solver; Gurobi or CPLEX are faster for large models if you have a licence.

```bash
pip install cobra
python -c "import cobra; print(cobra.__version__)"
```

### Step 2: Load a model

```python
import cobra
from cobra.io import load_model, read_sbml_model

# Bundled E. coli core model (no download needed)
model = load_model("textbook")

# Or your own genome-scale model, e.g. from BiGG (bigg.ucsd.edu)
# model = read_sbml_model("iML1515.xml")

print(len(model.reactions), len(model.metabolites), len(model.genes))
# 95 72 137
```

### Step 3: Inspect the model

```python
rxn = model.reactions.get_by_id("PGI")
print(rxn.reaction)            # g6p_c <=> f6p_c
print(rxn.bounds)              # (-1000.0, 1000.0)
print(rxn.gene_reaction_rule)  # b4025

print(model.objective.expression)   # biomass reaction
print(model.medium)                 # exchange reactions open for uptake
```

Check before any analysis:

- The objective is the biomass (or other) reaction you intend.
- The medium matches your experiment. In the core model glucose uptake is capped at 10 mmol/gDW/h and oxygen is open.
- Exchange reactions start with `EX_`. A negative flux means uptake, a positive flux means secretion.
- Fluxes are in mmol/gDW/h; growth rate is in 1/h.

## Steps 4–5: FBA and growth conditions

### Step 4: Run flux balance analysis

```python
solution = model.optimize()
print(solution.status)            # optimal
print(solution.objective_value)   # 0.8739  (growth rate, 1/h)
print(solution.fluxes["PGI"])     # flux through one reaction

model.summary()                   # uptake, secretion, objective
```

On glucose with oxygen, the core model grows at 0.874 1/h. It takes up 10 glucose, 21.8 oxygen and 4.77 ammonium, and secretes 22.8 CO2 (all mmol/gDW/h).

Always check `solution.status`. Anything other than `optimal` (for example `infeasible`) means the constraints cannot all be met.

### Step 5: Change the medium and bounds

Use a `with model:` block so changes are undone automatically when the block ends.

```python
# Anaerobic growth
with model:
    medium = model.medium
    medium["EX_o2_e"] = 0.0
    model.medium = medium
    print(model.slim_optimize())   # 0.2117

# Scan glucose uptake
for glc in [2, 5, 10]:
    with model:
        model.reactions.EX_glc__D_e.lower_bound = -glc
        print(glc, model.slim_optimize())
```

| Condition | Growth rate (1/h) |
| --- | --- |
| Aerobic, glucose 10 | 0.874 |
| Aerobic, glucose 5 | 0.416 |
| Aerobic, glucose 2 | 0.141 |
| Anaerobic, glucose 10 | 0.212 |

`slim_optimize()` returns only the objective value and is faster for loops. To set uptake directly, give the `model.medium` dictionary positive uptake limits; to set a bound, give the exchange reaction a negative lower bound.

## Steps 6–7: pFBA and flux variability analysis

FBA usually has many equally optimal flux distributions. These two methods deal with that.

### Step 6: Parsimonious FBA (pFBA)

pFBA keeps optimal growth but picks the solution with the smallest total flux. It is closer to how cells spend enzyme and is a better choice when comparing with measured fluxes.

```python
from cobra.flux_analysis import pfba

pfba_solution = pfba(model)
print(pfba_solution.fluxes.abs().sum())   # 518.4 mmol/gDW/h total
```

### Step 7: Flux variability analysis (FVA)

FVA finds the minimum and maximum flux of each reaction while growth stays at or above a set fraction of the optimum.

```python
from cobra.flux_analysis import flux_variability_analysis

fva = flux_variability_analysis(
    model,
    reaction_list=model.reactions,   # or a subset
    fraction_of_optimum=0.9,         # growth >= 90 % of max
)
print(fva.loc[["ACONTa", "ACKr"]])
#          minimum   maximum
# ACONTa   0.849     8.895
# ACKr    -3.814     0.000
```

How to read it:

- Minimum and maximum both far from 0 with the same sign: the reaction is required.
- Range includes 0: the reaction is optional at that growth level.
- Minimum equals maximum: the flux is fixed.
- A very wide range on reactions in a cycle can mean a thermodynamically impossible loop. Add `loopless=True` to remove those.

## Step 8: Gene and reaction knockouts

Knockouts predict which genes are essential and which deletions redirect flux. COBRApy uses each reaction's gene–protein–reaction (GPR) rule, so deleting one isozyme does not block the reaction.

```python
from cobra.flux_analysis import single_gene_deletion, single_reaction_deletion

wt_growth = model.slim_optimize()

# All single-gene deletions
gene_ko = single_gene_deletion(model)
essential = gene_ko[gene_ko["growth"] < 0.01 * wt_growth]
print(len(essential))          # 5 essential genes on aerobic glucose

# Single-reaction deletions
rxn_ko = single_reaction_deletion(model)

# One targeted knockout: zwf (b1852, glucose-6-phosphate dehydrogenase)
with model:
    model.genes.b1852.knock_out()
    print(model.slim_optimize())   # 0.8638, a small drop from 0.8739
```

Notes:

- Growth below 1 % of wild type is a common cut-off for "essential". State the cut-off you use.
- For double knockouts use `double_gene_deletion(model)`. It scales with the square of the gene number, so run it on a subset for genome-scale models.
- `processes=4` (or more) runs deletions in parallel.
- Compare predicted essential genes with a knockout library such as Keio for *E. coli* to check model quality.

## Steps 9–10: Production envelopes and new pathways

### Step 9: Production envelope (growth vs product)

A production envelope shows the trade-off between growth and product secretion. It tells you the theoretical maximum yield and how much growth you give up to reach it.

```python
from cobra.flux_analysis import production_envelope

env = production_envelope(
    model,
    reactions=["Biomass_Ecoli_core"],   # x-axis
    objective="EX_succ_e",              # product to maximize
    points=5,
)
print(env[["Biomass_Ecoli_core", "flux_minimum", "flux_maximum"]])
```

| Growth rate (1/h) | Max succinate secretion (mmol/gDW/h) |
| --- | --- |
| 0.000 | 16.38 |
| 0.218 | 12.38 |
| 0.437 | 8.31 |
| 0.655 | 4.18 |
| 0.874 | 0.00 |

At maximum growth the wild type secretes no succinate; it can only make it by growing slower. Strain-design tools such as OptKnock look for knockouts that force product out at high growth.

### Step 10: Add a reaction or pathway

Use this to add a heterologous pathway, a transport step or a demand reaction for a product.

```python
from cobra import Reaction, Metabolite

with model:
    rxn = Reaction("DM_succ_c")
    rxn.name = "Succinate demand"
    rxn.bounds = (0, 1000)
    rxn.add_metabolites({model.metabolites.succ_c: -1})
    model.add_reactions([rxn])

    model.objective = "DM_succ_c"
    print(model.slim_optimize())   # 16.73 mmol/gDW/h maximum
```

For a new metabolite, create it with `Metabolite("id_c", formula="...", compartment="c")` before adding it to a reaction. Set `rxn.gene_reaction_rule = "geneA and geneB"` to link genes. Run `rxn.check_mass_balance()` on every new reaction; an empty result means it is balanced.

## Step 11: Save, export and reproducibility

```python
from cobra.io import write_sbml_model, read_sbml_model, save_json_model

write_sbml_model(model, "my_model.xml")    # SBML, the standard exchange format
save_json_model(model, "my_model.json")    # JSON, used by Escher maps

reloaded = read_sbml_model("my_model.xml")
print(reloaded.slim_optimize())             # 0.8739, same as before

# Export results for Excel or plotting
solution.fluxes.to_csv("fba_fluxes.csv")
fva.to_csv("fva_results.csv")
```

For reproducible work, record:

- COBRApy version (`cobra.__version__`) and solver (`model.solver.interface`).
- Model name, version and source (for example BiGG iML1515).
- The exact medium, objective and any bounds you changed.
- Cut-offs used (essential-gene threshold, FVA fraction of optimum).

Re-loading the saved SBML and getting the same growth rate (0.8739 here) is a quick check that the file is complete.

To see fluxes on a pathway map, load the JSON model into [Escher](https://escher.github.io) and colour reactions by the exported flux values.

## Troubleshooting and good practice

| Problem | Likely cause | Fix |
| --- | --- | --- |
| Status `infeasible` | Bounds conflict, e.g. a required nutrient closed or a forced flux too high | Reopen exchanges one at a time; check bounds you set |
| Growth is 0 but status optimal | Missing nutrient in the medium, or a blocked biomass precursor | Compare `model.medium` with your real medium; run FVA on biomass precursors |
| Growth far higher than measured | Uptake bounds left at 1000 | Set measured uptake rates (glucose, oxygen) as bounds |
| Changes carry over between runs | Edits made outside a `with model:` block | Use `with model:` or work on `model.copy()` |
| Huge flux values (near 1000) in loops | Thermodynamically infeasible cycles | Use `pfba`, or `loopless=True` in FVA |
| Slow genome-scale deletions | Single process and GLPK solver | Use `processes=` and a faster solver such as Gurobi or CPLEX |
| SBML fails to load | Non-standard or old file | Check with `cobra.io.validate_sbml_model("file.xml")` |

**Good practice:**

- Start from a curated model (BiGG, BioModels or the organism's published reconstruction), not a draft.
- Set realistic uptake rates from your own experiments before trusting predictions.
- Validate the model: compare predicted growth rates, by-products and essential genes with lab data.
- Keep a script or notebook for every analysis so it can be re-run with a new model version.
- Report fluxes in mmol/gDW/h and growth in 1/h, and say which method (FBA, pFBA, FVA) produced each number.

Sources: [COBRApy documentation](https://cobrapy.readthedocs.io), [BiGG Models](http://bigg.ucsd.edu).
