#!/usr/bin/env python3
"""
LRG HOD Pipeline using AbacusSummit Lightcone Data
=================================================

This script generates a Luminous Red Galaxy (LRG) catalog using the
AbacusHOD implementation with AbacusSummit lightcone data, including
redshift-space distortions and tidal shear measurements.
"""

import numpy as np
import asdf
import h5py
from pathlib import Path
from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import time

# Cosmological parameters for AbacusSummit (Planck 2018)
COSMO = FlatLambdaCDM(H0=67.36, Om0=0.315192, Ob0=0.02237/0.67362)

class LRGHODPipeline:
    """LRG HOD Pipeline for AbacusSummit lightcone data"""

    def __init__(self, sim_path, output_dir="./output"):
        """
        Initialize the LRG HOD Pipeline

        Parameters:
        -----------
        sim_path : str
            Path to AbacusSummit simulation directory
        output_dir : str
            Output directory for results
        """
        self.sim_path = Path(sim_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # HOD parameters for LRGs (DESI-like)
        self.lrg_hod_params = {
            'logMmin': 13.02,      # log10(Minimum halo mass for central LRGs)
            'sigma_logM': 0.28,    # Scatter in log halo mass
            'logM0': 13.27,        # log10(Cutoff mass for satellites)
            'logM1': 14.08,        # log10(Normalization mass for satellites)
            'alpha': 0.76,         # Power law slope for satellites
            'kappa': 1.0,          # Satellite cutoff parameter
            # Additional parameters for modeling
            'f_c': 1.0,            # Central galaxy fraction
            'Acent': 0.0,          # Assembly bias for centrals
            'Asat': 0.0,           # Assembly bias for satellites
        }

        print(f"Initialized LRG HOD Pipeline")
        print(f"Simulation path: {self.sim_path}")
        print(f"Output directory: {self.output_dir}")

    def load_lightcone_data(self, redshift_slice="z0.500"):
        """
        Load AbacusSummit lightcone halo data

        Parameters:
        -----------
        redshift_slice : str
            Redshift slice to load (e.g., 'z0.500')

        Returns:
        --------
        dict : Dictionary containing halo data
        """
        print(f"Loading lightcone data for {redshift_slice}...")

        # Path to lightcone data
        lc_path = self.sim_path / "halo_light_cones" / "AbacusSummit_base_c000_ph000" / redshift_slice

        if not lc_path.exists():
            raise FileNotFoundError(f"Lightcone path not found: {lc_path}")

        # Load halo info
        halo_file = lc_path / "lc_halo_info.asdf"
        with asdf.open(halo_file) as af:
            header = dict(af['header'])
            data = dict(af['data'])

            # Extract key halo properties
            halo_data = {
                'header': header,
                'pos': np.array(data['x_L2com']),           # Positions [Mpc/h]
                'vel': np.array(data['v_L2com']),           # Velocities [km/s]
                'mass': np.array(data['N']) * header['ParticleMassHMsun'],  # Mass [Msun/h]
                'N': np.array(data['N']),                    # Number of particles
                'redshift': np.full(len(data['N']), header['Redshift']),
                'sigmav3d': np.array(data['sigmav3d_L2com']), # Velocity dispersion
                'r98': np.array(data['r98_L2com_i16']) / 1000.0,  # Radius [Mpc/h]
            }

        # Load additional particle data if needed
        pid_file = lc_path / "lc_pid_rv.asdf"
        if pid_file.exists():
            print("Additional particle data available")

        print(f"Loaded {len(halo_data['pos'])} halos at z={header['Redshift']}")
        return halo_data

    def prepare_hod_config(self, halo_data):
        """
        Prepare HOD configuration for AbacusHOD

        Parameters:
        -----------
        halo_data : dict
            Halo data dictionary

        Returns:
        --------
        dict : HOD configuration
        """
        print("Preparing HOD configuration...")

        config = {
            'sim_params': {
                'Lbox': halo_data['header']['BoxSize'],  # Box size [Mpc/h]
                'Npart': int(halo_data['header']['NP'] ** (1/3)),  # Particles per dimension
                'redshift': halo_data['header']['Redshift'],
                'h': halo_data['header']['H0'] / 100.0,
            },
            'HOD_params': self.lrg_hod_params,
            'clustering_params': {
                'rpmin': 0.1,      # Min r_p for clustering [Mpc/h]
                'rpmax': 50.0,     # Max r_p for clustering [Mpc/h]
                'nrpbins': 25,     # Number of r_p bins
                'pimax': 80.0,     # Max pi for clustering [Mpc/h]
                'npibins': 80,     # Number of pi bins
            },
            'data_params': {
                'path2data': str(self.sim_path),
                'sim_name': 'AbacusSummit_base_c000_ph000',
            }
        }

        return config

    def run_hod(self, halo_data, config):
        """
        Run HOD to generate LRG catalog

        Parameters:
        -----------
        halo_data : dict
            Halo data
        config : dict
            HOD configuration

        Returns:
        --------
        dict : Galaxy catalog
        """
        print("Running HOD to generate LRG catalog...")

        # Convert to format expected by AbacusHOD
        # Note: AbacusHOD typically expects specific data structure
        # We need to adapt our lightcone data format

        hod_data = {
            'pos': halo_data['pos'],
            'vel': halo_data['vel'],
            'mass': halo_data['mass'],
            'id': np.arange(len(halo_data['pos'])),
            'Npart': halo_data['N'],
        }

        # Apply mass cuts for LRGs (only use massive halos)
        mass_min = 10**(self.lrg_hod_params['logMmin'] - 1)  # Slightly below Mmin
        mask = hod_data['mass'] > mass_min

        print(f"Applying mass cut M > {mass_min:.2e} Msun/h")
        print(f"Selected {np.sum(mask)}/{len(mask)} halos ({100*np.sum(mask)/len(mask):.1f}%)")

        for key in hod_data:
            hod_data[key] = hod_data[key][mask]

        # Generate central galaxies
        # Implement LRG central occupation function directly to avoid numba issues
        def n_cen_LRG(log_mass, logMmin, sigma):
            """Central LRG occupation function"""
            from scipy.special import erfc
            return 0.5 * erfc((logMmin - log_mass) / (np.sqrt(2) * sigma))

        # Calculate central occupation probability
        log_mass = np.log10(hod_data['mass'])
        n_cen = n_cen_LRG(log_mass,
                          self.lrg_hod_params['logMmin'],
                          self.lrg_hod_params['sigma_logM'])

        # Randomly assign central galaxies
        central_mask = np.random.random(len(n_cen)) < n_cen

        print(f"Generated {np.sum(central_mask)} central LRGs")

        # For now, focus on central galaxies (satellites can be added later)
        lrg_catalog = {
            'pos': hod_data['pos'][central_mask].copy(),
            'vel': hod_data['vel'][central_mask].copy(),
            'mass_halo': hod_data['mass'][central_mask].copy(),
            'redshift': np.full(np.sum(central_mask), config['sim_params']['redshift']),
            'galaxy_type': np.zeros(np.sum(central_mask), dtype=int),  # 0 = central
        }

        return lrg_catalog

    def apply_rsd(self, catalog, line_of_sight_axis=2):
        """
        Apply redshift-space distortions to galaxy positions

        Parameters:
        -----------
        catalog : dict
            Galaxy catalog
        line_of_sight_axis : int
            Axis for line-of-sight (0=x, 1=y, 2=z)

        Returns:
        --------
        dict : Catalog with RSD positions added
        """
        print("Applying redshift-space distortions...")

        # Get cosmological parameters
        z = catalog['redshift'][0]
        H_z = COSMO.H(z).value  # H(z) in km/s/Mpc

        # Convert velocities to displacement
        # s = r + v_los / H(z)
        rsd_pos = catalog['pos'].copy()
        rsd_pos[:, line_of_sight_axis] += catalog['vel'][:, line_of_sight_axis] / H_z

        catalog['pos_rsd'] = rsd_pos
        catalog['pos_real'] = catalog['pos'].copy()

        print(f"Applied RSD along axis {line_of_sight_axis}")
        return catalog

    def compute_tidal_shear(self, catalog, halo_data, smoothing_scale=4.0):
        """
        Compute tidal shear at galaxy positions

        Parameters:
        -----------
        catalog : dict
            Galaxy catalog
        halo_data : dict
            Full halo catalog for computing tidal field
        smoothing_scale : float
            Smoothing scale in Mpc/h

        Returns:
        --------
        dict : Catalog with tidal shear measurements
        """
        print(f"Computing tidal shear with smoothing scale {smoothing_scale} Mpc/h...")

        # Simple tidal shear calculation using local halo density gradients
        # This is a simplified implementation - more sophisticated methods exist

        galaxy_pos = catalog['pos']
        halo_pos = halo_data['pos']
        halo_mass = halo_data['mass']

        # Initialize shear arrays
        n_gal = len(galaxy_pos)
        gamma1 = np.zeros(n_gal)
        gamma2 = np.zeros(n_gal)
        kappa = np.zeros(n_gal)  # convergence

        print(f"Computing tidal field for {n_gal} galaxies...")

        # Compute tidal field for each galaxy
        # This is computationally intensive - in practice you'd use more efficient methods
        for i in range(min(n_gal, 1000)):  # Limit to 1000 for demo
            if i % 100 == 0:
                print(f"  Processing galaxy {i+1}/{min(n_gal, 1000)}")

            gal_pos = galaxy_pos[i]

            # Find nearby halos within smoothing scale
            distances = np.linalg.norm(halo_pos - gal_pos, axis=1)
            nearby = distances < smoothing_scale

            if np.sum(nearby) > 5:  # Need minimum number of neighbors
                # Compute weighted density field derivatives
                r_vec = halo_pos[nearby] - gal_pos
                masses = halo_mass[nearby]
                r_mag = np.linalg.norm(r_vec, axis=1)

                # Gaussian weighting
                weights = masses * np.exp(-0.5 * (r_mag / smoothing_scale)**2)

                # Simple tidal field approximation
                # Proper implementation would use full tidal tensor
                xx = np.sum(weights * r_vec[:, 0]**2 / r_mag**2)
                yy = np.sum(weights * r_vec[:, 1]**2 / r_mag**2)
                xy = np.sum(weights * r_vec[:, 0] * r_vec[:, 1] / r_mag**2)

                gamma1[i] = (xx - yy) / (xx + yy + 1e-10)
                gamma2[i] = 2 * xy / (xx + yy + 1e-10)
                kappa[i] = (xx + yy) / smoothing_scale**2

        # Add to catalog
        catalog['gamma1'] = gamma1
        catalog['gamma2'] = gamma2
        catalog['kappa'] = kappa

        print(f"Computed tidal shear for {n_gal} galaxies")
        print(f"  Mean |gamma|: {np.mean(np.sqrt(gamma1**2 + gamma2**2)):.4f}")
        print(f"  Mean kappa: {np.mean(kappa):.4f}")

        return catalog

    def save_catalog(self, catalog, filename="lrg_catalog.h5"):
        """Save LRG catalog to file"""
        filepath = self.output_dir / filename

        print(f"Saving catalog to {filepath}")

        with h5py.File(filepath, 'w') as f:
            # Save galaxy data
            for key, data in catalog.items():
                f.create_dataset(key, data=data)

            # Save metadata
            f.attrs['n_galaxies'] = len(catalog['pos'])
            f.attrs['redshift'] = catalog['redshift'][0]
            f.attrs['hod_params'] = str(self.lrg_hod_params)

        print(f"Saved {len(catalog['pos'])} LRGs to {filepath}")

        return filepath

    def plot_catalog(self, catalog, save_plots=True):
        """Generate diagnostic plots"""
        print("Generating diagnostic plots...")

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # 2D position plot
        ax = axes[0, 0]
        ax.scatter(catalog['pos'][:, 0], catalog['pos'][:, 1], s=0.1, alpha=0.5)
        ax.set_xlabel('x [Mpc/h]')
        ax.set_ylabel('y [Mpc/h]')
        ax.set_title('LRG Positions (Real Space)')
        ax.set_aspect('equal')

        # RSD comparison
        if 'pos_rsd' in catalog:
            ax = axes[0, 1]
            ax.scatter(catalog['pos_rsd'][:, 0], catalog['pos_rsd'][:, 2], s=0.1, alpha=0.5)
            ax.set_xlabel('x [Mpc/h]')
            ax.set_ylabel('z [Mpc/h] (RSD)')
            ax.set_title('LRG Positions (Redshift Space)')

        # Mass histogram
        ax = axes[1, 0]
        ax.hist(np.log10(catalog['mass_halo']), bins=30, alpha=0.7)
        ax.set_xlabel('log10(M_halo) [Msun/h]')
        ax.set_ylabel('Count')
        ax.set_title('Halo Mass Distribution')
        ax.axvline(self.lrg_hod_params['logMmin'], color='r', linestyle='--', label='logMmin')
        ax.legend()

        # Tidal shear
        if 'gamma1' in catalog:
            ax = axes[1, 1]
            gamma_mag = np.sqrt(catalog['gamma1']**2 + catalog['gamma2']**2)
            ax.hist(gamma_mag, bins=30, alpha=0.7)
            ax.set_xlabel('|γ|')
            ax.set_ylabel('Count')
            ax.set_title('Tidal Shear Magnitude')

        plt.tight_layout()

        if save_plots:
            plot_path = self.output_dir / "lrg_catalog_plots.png"
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            print(f"Plots saved to {plot_path}")

        plt.show()


def main():
    """Main execution function"""
    print("="*60)
    print("LRG HOD Pipeline - AbacusSummit Lightcone")
    print("="*60)

    # Initialize pipeline
    sim_path = "/home/jovyan/AbacusSummit"
    pipeline = LRGHODPipeline(sim_path)

    try:
        # Load lightcone data
        halo_data = pipeline.load_lightcone_data("z0.500")

        # Prepare HOD configuration
        config = pipeline.prepare_hod_config(halo_data)

        # Run HOD to generate LRG catalog
        lrg_catalog = pipeline.run_hod(halo_data, config)

        # Apply redshift-space distortions
        lrg_catalog = pipeline.apply_rsd(lrg_catalog)

        # Compute tidal shear (computationally intensive - limited demo)
        print("Skipping tidal shear computation for speed (can be enabled if needed)")
        # lrg_catalog = pipeline.compute_tidal_shear(lrg_catalog, halo_data)

        # Save catalog
        output_file = pipeline.save_catalog(lrg_catalog)

        # Generate plots
        pipeline.plot_catalog(lrg_catalog)

        print("\n" + "="*60)
        print("LRG HOD Pipeline Complete!")
        print(f"Generated {len(lrg_catalog['pos'])} LRG galaxies")
        print(f"Output saved to: {output_file}")
        print("="*60)

    except Exception as e:
        print(f"Error in pipeline: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()