#!/usr/bin/env python3
"""
Quick test of LRG HOD Pipeline
"""

import numpy as np
import matplotlib.pyplot as plt
from lrg_hod_pipeline import LRGHODPipeline

def main():
    """Quick test run"""
    print("Testing LRG HOD Pipeline...")

    # Initialize pipeline
    sim_path = "/home/jovyan/AbacusSummit"
    pipeline = LRGHODPipeline(sim_path)

    # Load lightcone data
    halo_data = pipeline.load_lightcone_data("z0.500")

    # Prepare HOD configuration
    config = pipeline.prepare_hod_config(halo_data)

    # Run HOD to generate LRG catalog
    lrg_catalog = pipeline.run_hod(halo_data, config)

    # Apply redshift-space distortions
    lrg_catalog = pipeline.apply_rsd(lrg_catalog)

    # Print summary statistics
    print(f"\nLRG Catalog Summary:")
    print(f"Number of LRGs: {len(lrg_catalog['pos'])}")
    print(f"Redshift: {lrg_catalog['redshift'][0]}")
    print(f"Mean halo mass: {np.mean(lrg_catalog['mass_halo']):.2e} Msun/h")
    print(f"Median halo mass: {np.median(lrg_catalog['mass_halo']):.2e} Msun/h")

    # Calculate some basic clustering statistics
    pos = lrg_catalog['pos']
    print(f"Position range:")
    print(f"  X: [{np.min(pos[:,0]):.1f}, {np.max(pos[:,0]):.1f}] Mpc/h")
    print(f"  Y: [{np.min(pos[:,1]):.1f}, {np.max(pos[:,1]):.1f}] Mpc/h")
    print(f"  Z: [{np.min(pos[:,2]):.1f}, {np.max(pos[:,2]):.1f}] Mpc/h")

    # RSD effect
    if 'pos_rsd' in lrg_catalog:
        rsd_shift = np.mean(np.abs(lrg_catalog['pos_rsd'][:, 2] - lrg_catalog['pos'][:, 2]))
        print(f"Mean RSD shift: {rsd_shift:.2f} Mpc/h")

    # Save catalog
    output_file = pipeline.save_catalog(lrg_catalog, "test_lrg_catalog.h5")

    # Quick visualization
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Real space positions
    sample_idx = np.random.choice(len(pos), size=min(10000, len(pos)), replace=False)
    axes[0].scatter(pos[sample_idx, 0], pos[sample_idx, 1], s=0.5, alpha=0.6)
    axes[0].set_xlabel('x [Mpc/h]')
    axes[0].set_ylabel('y [Mpc/h]')
    axes[0].set_title(f'LRG Positions (Real Space)\n{len(pos)} galaxies')
    axes[0].set_aspect('equal')

    # Mass distribution
    axes[1].hist(np.log10(lrg_catalog['mass_halo']), bins=50, alpha=0.7, edgecolor='black')
    axes[1].set_xlabel('log10(M_halo) [Msun/h]')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Halo Mass Distribution')
    axes[1].axvline(pipeline.lrg_hod_params['logMmin'], color='r', linestyle='--',
                   label=f"logMmin = {pipeline.lrg_hod_params['logMmin']}")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig('output/test_lrg_results.png', dpi=150, bbox_inches='tight')
    plt.show()

    print(f"\nTest completed successfully!")
    print(f"Catalog saved to: {output_file}")
    print(f"Plot saved to: output/test_lrg_results.png")

if __name__ == "__main__":
    main()