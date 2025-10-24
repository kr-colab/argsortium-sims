# argsortium-sims

Population genetic simulation pipeline using [stdpopsim](https://stdpopsim.readthedocs.io/), [msprime](https://tskit.dev/msprime/), and [Snakemake](https://snakemake.readthedocs.io/).

## Overview

This pipeline simulates tree sequences across multiple chromosomes using demographic models from the stdpopsim catalog. Each chromosome is simulated independently, allowing for parallel execution across multiple cores.

**Pipeline steps:**
1. **Simulate ancestry**: Generate coalescent tree sequences without mutations
2. **Add mutations**: Overlay mutations on the tree sequences using specified mutation models
3. **Write output**: Write out VCF, applying perturbations like mispolarisation or masking

## Requirements

- Python 3.13+
- [uv](https://docs.astral.sh/uv/) package manager

## Installation

Install dependencies using uv:

```bash
uv sync
```

## Configuration

Edit `config/config.yaml` to configure your simulation:

```yaml
# Simulation parameters
species_name: "HomSap"
model_name: "OutOfAfrica_3G09"

# Sample configuration
samples:
  YRI: 1  # Yoruba
  CEU: 1  # Utah Residents (CEPH)
  CHB: 1  # Han Chinese

# Chromosome list for parallelization
contigs:
  - chr22

# Genetic map
genetic_map: "HapMapII_GRCh38"

# Outgroup configuration
add_outgroup: true
outgroup_time: 215000  # generations ago

# Mutation parameters
mutation_model: "HKY"  # Options: HKY, JC69, GTR, BinaryMutationModel
discrete_genome: true  # Use discrete genome (finite sites)

# Reproducibility
base_seed: 42  # Each contig gets: base_seed + contig_index

# Output
output_dir: "results/simulations"
```

### Example Configurations

- `config/config.yaml` - Default config (chr22 only, for testing)
- `config/examples/ooa_3g09.yaml` - Full genome (all 22 autosomes)

## Running Simulations

### Quick Start

Run a single chromosome simulation:

```bash
uv run snakemake --snakefile workflow/Snakefile --cores 1
```

### Full Genome Simulation

Run all chromosomes in parallel (using 8 cores):

```bash
uv run snakemake --snakefile workflow/Snakefile --configfile config/examples/ooa_3g09.yaml --cores 8
```

### Common Commands

Dry run (see what will execute):
```bash
uv run snakemake --snakefile workflow/Snakefile --cores 1 --dry-run
```

Run with specific number of cores:
```bash
uv run snakemake --snakefile workflow/Snakefile --cores 8
```

Run specific chromosome:
```bash
uv run snakemake --snakefile workflow/Snakefile results/simulations/HomSap_OutOfAfrica_3G09/chr22/sim.mutated.trees --cores 1
```

Run only ancestry simulation (skip mutations):
```bash
uv run snakemake --snakefile workflow/Snakefile results/simulations/HomSap_OutOfAfrica_3G09/chr22/sim.trees --cores 1
```

Visualize workflow DAG:
```bash
uv run snakemake --snakefile workflow/Snakefile --dag | dot -Tpdf > dag.pdf
```

## Output

Each chromosome simulation produces:

- `{output_dir}/{species}_{model}/{contig}/sim.tsz` - Tree sequence (ancestry only, no mutations, compressed with `tszip`)
- `{output_dir}/{species}_{model}/{contig}/sim.log` - Ancestry simulation timing and statistics
- `{output_dir}/{species}_{model}/{contig}/sim.mutated.tsz` - Tree sequence with mutations, compressed with `tszip`
- `{output_dir}/{species}_{model}/{contig}/sim.mutated.log` - Mutation timing and statistics

Example output location:
```
results/simulations/HomSap_OutOfAfrica_3G09/chr22/
├── sim.trees           # Ancestry only
├── sim.log             # Ancestry log
├── sim.mutated.trees   # With mutations (final output)
└── sim.mutated.log     # Mutation log
```

## Mutation Models

The pipeline supports multiple mutation models:

- **HKY** (default): Hasegawa-Kishino-Yano model with transition/transversion rate variation (kappa=2.0)
- **JC69**: Jukes-Cantor model (all mutations equally likely)
- **GTR**: General Time Reversible model with equal rates
- **BinaryMutationModel**: Simple binary (0/1) mutations

Mutation rates are automatically pulled from the stdpopsim species catalog for each chromosome.

## Random Seed Strategy

Each chromosome gets deterministic seeds for reproducibility:

**Ancestry simulation:**
- chr1: `base_seed + 0`
- chr2: `base_seed + 1`
- chr22: `base_seed + 21`

**Mutation simulation:**
- chr1: `base_seed + 1000`
- chr2: `base_seed + 1001`
- chr22: `base_seed + 1021`

This ensures reproducible results while keeping ancestry and mutation seeds distinct.

## Project Structure

```
argsortium-sims/
├── workflow/
│   ├── Snakefile              # Main workflow entry point
│   ├── rules/
│   │   ├── simulate.smk       # Ancestry simulation rules
│   │   └── mutate.smk         # Mutation rules
│   └── scripts/
│       ├── run_simulation.py  # Ancestry simulation script
│       └── add_mutations.py   # Mutation script
├── config/
│   ├── config.yaml            # Main configuration file
│   └── examples/              # Example configurations
├── results/
│   └── simulations/           # Output tree sequences (.trees files)
└── logs/                      # Snakemake execution logs
```

## Available Models

This pipeline uses the stdpopsim catalog. To see available species, models, and genetic maps:

```bash
uv run python -c "import stdpopsim; help(stdpopsim)"
```

Or visit the [stdpopsim catalog documentation](https://stdpopsim.readthedocs.io/en/stable/catalog.html).


