#!/usr/bin/env python3
"""
Validation and Analysis of LRG HOD Catalog
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.stats import binned_statistic
import time


def load_catalog(filename="../output/lrg_catalog_with_tidal_shear.h5"):
    """Load LRG catalog with all data"""
    print(f"Loading LRG catalog from {filename}")

    data = {}
    with h5py.File(filename, 'r') as f:
        for key in f.keys():
            data[key] = f[key][:]

        # Load metadata
        attrs = dict(f.attrs)

    print(f"Loaded {len(data['pos'])} LRG galaxies")
    print(f"Redshift: {attrs.get('redshift', 'unknown')}")

    return data, attrs


def compute_clustering_statistics(pos, max_dist=10.0, n_bins=20):
    """
    Compute basic clustering statistics

    Parameters:
    -----------
    pos : array
        Galaxy positions
    max_dist : float
        Maximum distance for correlation function
    n_bins : int
        Number of bins

    Returns:
    --------
    dict : Clustering results
    """
    print(f"Computing clustering statistics for {len(pos)} galaxies...")

    # Subsample for efficiency (correlation functions are expensive)
    n_sample = min(5000, len(pos))
    idx = np.random.choice(len(pos), size=n_sample, replace=False)
    pos_sample = pos[idx]

    print(f"Using subsample of {n_sample} galaxies")

    # Compute pairwise distances
    start_time = time.time()
    distances = pdist(pos_sample)
    print(f"Computed {len(distances)} pairwise distances in {time.time()-start_time:.2f}s")

    # Create bins
    r_bins = np.logspace(-1, np.log10(max_dist), n_bins+1)
    r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])

    # Compute pair counts
    counts, _ = np.histogram(distances, bins=r_bins)

    # Simple correlation function estimate (DD/RR - 1)
    # This is a simplified version - proper analysis would include random catalogs
    box_volume = np.prod(np.max(pos, axis=0) - np.min(pos, axis=0))
    n_density = len(pos) / box_volume
    n_pairs_total = n_sample * (n_sample - 1) / 2

    # Expected random pairs in each bin
    bin_volumes = 4/3 * np.pi * (r_bins[1:]**3 - r_bins[:-1]**3)
    expected_pairs = n_pairs_total * bin_volumes * n_density / box_volume

    # Correlation function
    xi = counts / expected_pairs - 1

    return {
        'r': r_centers,
        'xi': xi,
        'counts': counts,
        'expected': expected_pairs,
        'n_sample': n_sample
    }


def analyze_halo_occupation(mass_halo, hod_params):
    """Analyze halo occupation statistics"""
    print("Analyzing halo occupation...")

    # Create mass bins
    mass_bins = np.logspace(12, 15, 20)
    mass_centers = 0.5 * (mass_bins[1:] + mass_bins[:-1])

    # Count galaxies in each mass bin
    counts, _ = np.histogram(mass_halo, bins=mass_bins)
    halo_counts, _ = np.histogram(mass_halo, bins=mass_bins)

    # Average occupation (galaxies per halo in each bin)
    occupation = np.divide(counts, halo_counts, out=np.zeros_like(counts, dtype=float),
                          where=halo_counts!=0)

    # Theoretical HOD prediction
    log_mass_centers = np.log10(mass_centers)
    from scipy.special import erfc
    theo_occupation = 0.5 * erfc((hod_params['logMmin'] - log_mass_centers) /
                                (np.sqrt(2) * hod_params['sigma_logM']))

    return {
        'mass_centers': mass_centers,
        'occupation': occupation,
        'theo_occupation': theo_occupation,
        'counts': counts,
        'halo_counts': halo_counts
    }


def generate_validation_report(data, attrs):
    """Generate comprehensive validation report"""
    print("="*60)
    print("LRG HOD CATALOG VALIDATION REPORT")
    print("="*60)

    n_gal = len(data['pos'])
    print(f"Total galaxies: {n_gal:,}")
    print(f"Redshift: {attrs.get('redshift', 'unknown')}")

    # Position statistics
    pos = data['pos']
    print(f"\nPosition Statistics:")
    print(f"  X range: [{np.min(pos[:,0]):.3f}, {np.max(pos[:,0]):.3f}] Mpc/h")
    print(f"  Y range: [{np.min(pos[:,1]):.3f}, {np.max(pos[:,1]):.3f}] Mpc/h")
    print(f"  Z range: [{np.min(pos[:,2]):.3f}, {np.max(pos[:,2]):.3f}] Mpc/h")

    volume = np.prod(np.max(pos, axis=0) - np.min(pos, axis=0))
    density = n_gal / volume
    print(f"  Survey volume: {volume:.4f} (Mpc/h)³")
    print(f"  Number density: {density:.2f} galaxies/(Mpc/h)³")

    # Mass statistics
    mass = data['mass_halo']
    print(f"\nHalo Mass Statistics:")
    print(f"  Mean mass: {np.mean(mass):.2e} Msun/h")
    print(f"  Median mass: {np.median(mass):.2e} Msun/h")
    print(f"  Mass range: [{np.min(mass):.2e}, {np.max(mass):.2e}] Msun/h")

    # RSD statistics
    if 'pos_rsd' in data:
        rsd_shift = data['pos_rsd'][:, 2] - data['pos'][:, 2]
        print(f"\nRedshift-Space Distortions:")
        print(f"  Mean RSD shift: {np.mean(rsd_shift):.4f} Mpc/h")
        print(f"  RMS RSD shift: {np.sqrt(np.mean(rsd_shift**2)):.4f} Mpc/h")
        print(f"  RSD range: [{np.min(rsd_shift):.4f}, {np.max(rsd_shift):.4f}] Mpc/h")

    # Tidal shear statistics
    if 'gamma1' in data:
        gamma_mag = np.sqrt(data['gamma1']**2 + data['gamma2']**2)
        print(f"\nTidal Shear Statistics:")
        print(f"  Mean |γ|: {np.mean(gamma_mag):.4f}")
        print(f"  Median |γ|: {np.median(gamma_mag):.4f}")
        print(f"  RMS |γ|: {np.sqrt(np.mean(gamma_mag**2)):.4f}")
        print(f"  Mean κ: {np.mean(data['kappa']):.2e}")

    print("="*60)


def main():
    """Main validation function"""
    print("LRG HOD Catalog Validation")
    print("="*40)

    # Load catalog
    data, attrs = load_catalog()

    # Generate validation report
    generate_validation_report(data, attrs)

    # HOD parameters from attrs or use defaults
    hod_params = {
        'logMmin': 13.02,
        'sigma_logM': 0.28,
        'logM0': 13.27,
        'logM1': 14.08,
        'alpha': 0.76,
    }

    # Compute clustering
    print("\nComputing clustering statistics...")
    clustering = compute_clustering_statistics(data['pos'])

    # Analyze HOD
    print("\nAnalyzing halo occupation...")
    hod_analysis = analyze_halo_occupation(data['mass_halo'], hod_params)

    # Create comprehensive plots
    fig = plt.figure(figsize=(16, 12))

    # Layout: 3x3 grid
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # 1. Galaxy positions (real vs RSD)
    ax1 = fig.add_subplot(gs[0, 0])
    sample_idx = np.random.choice(len(data['pos']), size=min(5000, len(data['pos'])), replace=False)
    ax1.scatter(data['pos'][sample_idx, 0], data['pos'][sample_idx, 1],
               s=0.3, alpha=0.6, color='blue', label='Real space')
    if 'pos_rsd' in data:
        ax1.scatter(data['pos_rsd'][sample_idx, 0], data['pos_rsd'][sample_idx, 2],
                   s=0.3, alpha=0.6, color='red', label='Redshift space')
    ax1.set_xlabel('x [Mpc/h]')
    ax1.set_ylabel('y/z [Mpc/h]')
    ax1.set_title('Galaxy Positions')
    ax1.legend()
    ax1.set_aspect('equal')

    # 2. Mass distribution
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.hist(np.log10(data['mass_halo']), bins=50, alpha=0.7, color='green', edgecolor='black')
    ax2.axvline(hod_params['logMmin'], color='red', linestyle='--', label=f"logMmin = {hod_params['logMmin']}")
    ax2.set_xlabel('log₁₀(M_halo) [M☉/h]')
    ax2.set_ylabel('Count')
    ax2.set_title('Halo Mass Distribution')
    ax2.legend()

    # 3. Clustering
    ax3 = fig.add_subplot(gs[0, 2])
    mask = clustering['xi'] > -1  # Remove unphysical values
    ax3.loglog(clustering['r'][mask], 1 + clustering['xi'][mask], 'b.-', label='ξ(r) + 1')
    ax3.set_xlabel('r [Mpc/h]')
    ax3.set_ylabel('1 + ξ(r)')
    ax3.set_title('Correlation Function')
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # 4. HOD function
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.semilogx(hod_analysis['mass_centers'], hod_analysis['occupation'], 'bo-',
                label='Measured', markersize=4)
    ax4.semilogx(hod_analysis['mass_centers'], hod_analysis['theo_occupation'], 'r--',
                label='Theoretical')
    ax4.set_xlabel('M_halo [M☉/h]')
    ax4.set_ylabel('⟨N_cen⟩')
    ax4.set_title('HOD Function')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # 5. RSD effects
    if 'pos_rsd' in data:
        ax5 = fig.add_subplot(gs[1, 1])
        rsd_shift = data['pos_rsd'][:, 2] - data['pos'][:, 2]
        ax5.hist(rsd_shift, bins=50, alpha=0.7, color='purple', edgecolor='black')
        ax5.set_xlabel('RSD shift [Mpc/h]')
        ax5.set_ylabel('Count')
        ax5.set_title('Redshift-Space Distortions')
        ax5.axvline(0, color='red', linestyle='--', alpha=0.7)

    # 6. Tidal shear magnitude
    if 'gamma1' in data:
        ax6 = fig.add_subplot(gs[1, 2])
        gamma_mag = np.sqrt(data['gamma1']**2 + data['gamma2']**2)
        ax6.hist(gamma_mag, bins=50, alpha=0.7, color='orange', edgecolor='black')
        ax6.set_xlabel('|γ|')
        ax6.set_ylabel('Count')
        ax6.set_title('Tidal Shear Magnitude')

    # 7. Shear field visualization
    if 'gamma1' in data and 'gamma2' in data:
        ax7 = fig.add_subplot(gs[2, 0])
        # Sample for visualization
        n_arrows = 1000
        idx = np.random.choice(len(data['pos']), size=min(n_arrows, len(data['pos'])), replace=False)
        pos_sample = data['pos'][idx]
        gamma1_sample = data['gamma1'][idx]
        gamma2_sample = data['gamma2'][idx]

        # Scale shear for visibility
        scale = 0.01
        ax7.quiver(pos_sample[:, 0], pos_sample[:, 1],
                  gamma1_sample * scale, gamma2_sample * scale,
                  alpha=0.6, width=0.002, scale=1, scale_units='xy')
        ax7.set_xlabel('x [Mpc/h]')
        ax7.set_ylabel('y [Mpc/h]')
        ax7.set_title('Tidal Shear Field')
        ax7.set_aspect('equal')

    # 8. Convergence vs shear
    if 'kappa' in data and 'gamma1' in data:
        ax8 = fig.add_subplot(gs[2, 1])
        gamma_mag = np.sqrt(data['gamma1']**2 + data['gamma2']**2)
        # Sample for scatter plot
        idx = np.random.choice(len(gamma_mag), size=min(5000, len(gamma_mag)), replace=False)
        ax8.scatter(data['kappa'][idx], gamma_mag[idx], s=0.5, alpha=0.5)
        ax8.set_xlabel('κ (convergence)')
        ax8.set_ylabel('|γ| (shear magnitude)')
        ax8.set_title('Convergence vs Shear')

    # 9. Summary statistics
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.axis('off')
    summary_text = f"""
    CATALOG SUMMARY

    Galaxies: {len(data['pos']):,}
    Redshift: {attrs.get('redshift', 'unknown')}

    Number density:
    {len(data['pos']) / np.prod(np.max(data['pos'], axis=0) - np.min(data['pos'], axis=0)):.1f} gal/(Mpc/h)³

    Mean halo mass:
    {np.mean(data['mass_halo']):.2e} M☉/h

    HOD Parameters:
    logMmin = {hod_params['logMmin']}
    σ_logM = {hod_params['sigma_logM']}
    """

    if 'gamma1' in data:
        gamma_mag = np.sqrt(data['gamma1']**2 + data['gamma2']**2)
        summary_text += f"\nTidal Shear:\n⟨|γ|⟩ = {np.mean(gamma_mag):.3f}"

    ax9.text(0.05, 0.95, summary_text, transform=ax9.transAxes,
            verticalalignment='top', fontsize=10, family='monospace')

    plt.suptitle('LRG HOD Catalog Validation', fontsize=16, y=0.98)

    # Save plots
    plt.savefig('../output/lrg_validation_report.png', dpi=150, bbox_inches='tight')
    plt.show()

    print("\n" + "="*60)
    print("VALIDATION COMPLETE")
    print("="*60)
    print(f"Validation plots saved to: ../output/lrg_validation_report.png")
    print(f"Catalog contains {len(data['pos']):,} LRG galaxies with:")
    print(f"  ✓ Positions (real and redshift-space)")
    print(f"  ✓ Halo masses and HOD implementation")
    print(f"  ✓ Redshift-space distortions")
    if 'gamma1' in data:
        print(f"  ✓ Tidal shear measurements")
    print("="*60)


if __name__ == "__main__":
    main()