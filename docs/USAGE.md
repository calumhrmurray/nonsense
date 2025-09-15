# Usage Guide

## Quick Start

### 1. Basic LRG Catalog Generation
```bash
cd /home/jovyan/abacus_ia_analysis
python src/test_lrg_hod.py
```

### 2. Add Tidal Shear Measurements
```bash
python src/efficient_tidal_shear.py
```

### 3. Generate Validation Report
```bash
python src/validate_lrg_catalog.py
```

### 4. Full Pipeline (Advanced)
```bash
python src/lrg_hod_pipeline.py
```

## Configuration

### HOD Parameters
Edit the HOD parameters in `src/lrg_hod_pipeline.py`:

```python
self.lrg_hod_params = {
    'logMmin': 13.02,      # log10(Minimum halo mass for central LRGs)
    'sigma_logM': 0.28,    # Scatter in log halo mass
    'logM0': 13.27,        # log10(Cutoff mass for satellites)
    'logM1': 14.08,        # log10(Normalization mass for satellites)
    'alpha': 0.76,         # Power law slope for satellites
}
```

### Tidal Shear Settings
Modify tidal shear computation in `src/efficient_tidal_shear.py`:

```python
tidal_data = compute_efficient_tidal_shear(
    galaxy_pos, halo_pos, halo_mass,
    smoothing_scale=2.0,    # Smoothing scale in Mpc/h
    n_neighbors=30          # Number of nearest neighbors
)
```

## Output Files

- `output/test_lrg_catalog.h5` - Basic LRG catalog
- `output/lrg_catalog_with_tidal_shear.h5` - Complete catalog with IA
- `output/lrg_validation_report.png` - Validation plots
- `output/tidal_shear_analysis.png` - Tidal field analysis

## Data Requirements

The pipeline expects AbacusSummit data in the standard format:
```
/home/jovyan/AbacusSummit/
└── halo_light_cones/
    └── AbacusSummit_base_c000_ph000/
        ├── z0.100/
        ├── z0.500/
        └── ...
```

Each redshift directory should contain:
- `lc_halo_info.asdf` - Halo properties
- `lc_pid_rv.asdf` - Particle data (optional)