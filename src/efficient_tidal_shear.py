#!/usr/bin/env python3
"""
Efficient Tidal Shear Calculation for LRG Catalog
"""

import numpy as np
import h5py
from scipy.spatial import cKDTree
import time


def compute_efficient_tidal_shear(galaxy_pos, halo_pos, halo_mass,
                                smoothing_scale=4.0, n_neighbors=50):
    """
    Efficiently compute tidal shear using k-nearest neighbors

    Parameters:
    -----------
    galaxy_pos : array
        Galaxy positions [N_gal, 3]
    halo_pos : array
        Halo positions [N_halo, 3]
    halo_mass : array
        Halo masses [N_halo]
    smoothing_scale : float
        Smoothing scale in Mpc/h
    n_neighbors : int
        Number of nearest neighbors to use

    Returns:
    --------
    dict : Tidal shear components
    """
    print(f"Computing efficient tidal shear for {len(galaxy_pos)} galaxies...")
    start_time = time.time()

    n_gal = len(galaxy_pos)

    # Use k-d tree for efficient neighbor finding
    print("Building k-d tree...")
    tree = cKDTree(halo_pos)

    # Initialize shear arrays
    gamma1 = np.zeros(n_gal)
    gamma2 = np.zeros(n_gal)
    kappa = np.zeros(n_gal)

    # Process in batches for memory efficiency
    batch_size = 1000
    n_batches = (n_gal + batch_size - 1) // batch_size

    for batch_idx in range(n_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, n_gal)
        batch_gal_pos = galaxy_pos[start_idx:end_idx]

        if batch_idx % 10 == 0:
            print(f"  Processing batch {batch_idx+1}/{n_batches}")

        # Find nearest neighbors within smoothing scale
        distances, indices = tree.query(batch_gal_pos,
                                      k=min(n_neighbors, len(halo_pos)),
                                      distance_upper_bound=smoothing_scale)

        for i, gal_pos in enumerate(batch_gal_pos):
            global_idx = start_idx + i

            # Get valid neighbors (not at infinite distance)
            valid_mask = distances[i] < np.inf
            if np.sum(valid_mask) < 5:  # Need minimum neighbors
                continue

            neighbor_indices = indices[i][valid_mask]
            neighbor_distances = distances[i][valid_mask]

            # Get neighbor properties
            neighbor_pos = halo_pos[neighbor_indices]
            neighbor_mass = halo_mass[neighbor_indices]

            # Relative positions
            r_vec = neighbor_pos - gal_pos
            r_mag = neighbor_distances

            # Gaussian weighting
            weights = neighbor_mass * np.exp(-0.5 * (r_mag / smoothing_scale)**2)

            # Avoid division by zero
            r_mag_safe = np.maximum(r_mag, 1e-10)

            # Compute tidal tensor components (simplified)
            try:
                xx = np.sum(weights * (r_vec[:, 0] / r_mag_safe)**2)
                yy = np.sum(weights * (r_vec[:, 1] / r_mag_safe)**2)
                xy = np.sum(weights * (r_vec[:, 0] * r_vec[:, 1]) / r_mag_safe**2)

                total_weight = xx + yy
                if total_weight > 0:
                    gamma1[global_idx] = (xx - yy) / total_weight
                    gamma2[global_idx] = 2 * xy / total_weight
                    kappa[global_idx] = total_weight / (smoothing_scale**2)
            except:
                continue

    elapsed_time = time.time() - start_time
    print(f"Tidal shear computation completed in {elapsed_time:.2f} seconds")

    # Print statistics
    gamma_mag = np.sqrt(gamma1**2 + gamma2**2)
    print(f"Tidal shear statistics:")
    print(f"  Mean |γ|: {np.mean(gamma_mag):.4f}")
    print(f"  Std |γ|: {np.std(gamma_mag):.4f}")
    print(f"  Mean κ: {np.mean(kappa):.4f}")
    print(f"  Std κ: {np.std(kappa):.4f}")

    return {
        'gamma1': gamma1,
        'gamma2': gamma2,
        'kappa': kappa,
        'gamma_mag': gamma_mag
    }


def test_tidal_shear():
    """Test the efficient tidal shear calculation with LRG catalog"""

    # Load the LRG catalog
    catalog_file = "output/test_lrg_catalog.h5"

    with h5py.File(catalog_file, 'r') as f:
        galaxy_pos = f['pos'][:]
        print(f"Loaded {len(galaxy_pos)} LRGs from catalog")

    # Load halo data (sample for efficiency)
    print("Loading halo data...")
    import asdf
    with asdf.open('/home/jovyan/AbacusSummit/halo_light_cones/AbacusSummit_base_c000_ph000/z0.500/lc_halo_info.asdf') as af:
        # Sample halos for efficiency (use every 10th halo)
        step = 10
        halo_pos = np.array(af['data']['x_L2com'])[::step]
        halo_N = np.array(af['data']['N'])[::step]
        mass_per_particle = af['header']['ParticleMassHMsun']
        halo_mass = halo_N * mass_per_particle

    print(f"Using {len(halo_pos)} halos for tidal field computation")

    # Compute tidal shear
    tidal_data = compute_efficient_tidal_shear(
        galaxy_pos, halo_pos, halo_mass,
        smoothing_scale=2.0,  # Smaller scale for efficiency
        n_neighbors=30
    )

    # Visualize results
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Galaxy positions colored by tidal shear magnitude
    ax = axes[0, 0]
    scatter = ax.scatter(galaxy_pos[:, 0], galaxy_pos[:, 1],
                        c=tidal_data['gamma_mag'],
                        s=0.5, alpha=0.7, cmap='viridis')
    plt.colorbar(scatter, ax=ax, label='|γ|')
    ax.set_xlabel('x [Mpc/h]')
    ax.set_ylabel('y [Mpc/h]')
    ax.set_title('Galaxy Positions (colored by tidal shear)')
    ax.set_aspect('equal')

    # Tidal shear magnitude histogram
    ax = axes[0, 1]
    ax.hist(tidal_data['gamma_mag'], bins=50, alpha=0.7, edgecolor='black')
    ax.set_xlabel('|γ|')
    ax.set_ylabel('Count')
    ax.set_title('Tidal Shear Magnitude Distribution')

    # Convergence histogram
    ax = axes[1, 0]
    ax.hist(tidal_data['kappa'], bins=50, alpha=0.7, edgecolor='black')
    ax.set_xlabel('κ')
    ax.set_ylabel('Count')
    ax.set_title('Convergence Distribution')

    # Shear components correlation
    ax = axes[1, 1]
    ax.scatter(tidal_data['gamma1'], tidal_data['gamma2'],
              s=0.5, alpha=0.5)
    ax.set_xlabel('γ₁')
    ax.set_ylabel('γ₂')
    ax.set_title('Shear Components')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.axvline(0, color='k', linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig('output/tidal_shear_analysis.png', dpi=150, bbox_inches='tight')
    plt.show()

    # Save tidal shear data
    output_file = "output/lrg_catalog_with_tidal_shear.h5"

    # Load existing catalog and add tidal shear
    with h5py.File(catalog_file, 'r') as f_in:
        with h5py.File(output_file, 'w') as f_out:
            # Copy existing data
            for key in f_in.keys():
                f_out.create_dataset(key, data=f_in[key][:])

            # Add tidal shear data
            f_out.create_dataset('gamma1', data=tidal_data['gamma1'])
            f_out.create_dataset('gamma2', data=tidal_data['gamma2'])
            f_out.create_dataset('kappa', data=tidal_data['kappa'])
            f_out.create_dataset('gamma_mag', data=tidal_data['gamma_mag'])

            # Copy attributes
            for key, val in f_in.attrs.items():
                f_out.attrs[key] = val

            f_out.attrs['tidal_shear_computed'] = True
            f_out.attrs['tidal_smoothing_scale'] = 2.0

    print(f"\nTidal shear analysis complete!")
    print(f"Enhanced catalog saved to: {output_file}")
    print(f"Plots saved to: output/tidal_shear_analysis.png")

    return tidal_data


if __name__ == "__main__":
    tidal_data = test_tidal_shear()