# Neuromorphic Control — Formalizing Neuromorphic Control Systems

Code and figures for the paper **"Formalizing Neuromorphic Control Systems: A General Proposal and A Rhythmic Case Study"**.

## Contents
- `src/` — Julia modules: `model.jl`, `analysis.jl`, `utils.jl`, top-level `run_simulation.jl`
- `notebooks/` — interactive `main.ipynb`
- `figures/` — generated figures for the paper
- `Project.toml` / `Manifest.toml` — Julia environment

## Quick start

```bash
git clone https://github.com/<your-username>/neuromorphic-control.git
cd neuromorphic-control
julia --project=.
julia> using Pkg; Pkg.instantiate()    # install dependencies from Project.toml
