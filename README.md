# Neuromorphic Control — Formalizing Neuromorphic Control Systems

Code and figures for the paper **"Formalizing Neuromorphic Control Systems: A General Proposal and A Rhythmic Case Study"**.

## Contents
- `src/` — Julia modules: `model.jl`, `analysis.jl`, `utils.jl`, top-level `run_simulation.jl`
- `scripts/` — interactive `run_fig2.ipynb`, `run_fig3_1.ipynb`, `run_fig3_2.ipynb`, `run_fig4.ipynb`
- `figures/` — generated figures for the paper
- `Project.toml` — Julia environment

## Quick start

```bash
git clone https://github.com/TayaMedvedeva/fig-rhythmic-neuromorphic-case-study.git
cd fig-rhythmic-neuromorphic-case-study
julia --project=.
julia> using Pkg; Pkg.instantiate()    # install dependencies from Project.toml
