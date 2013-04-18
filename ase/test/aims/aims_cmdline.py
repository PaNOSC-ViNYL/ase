from ase.cli import run

# warning! parameters are not converged - only an illustration!
atoms = run('-x bcc -a 3.6 Li optimize -c aims -s 0.3 -p ' +
            'kpts=1.5,xc=LDA,sc_accuracy_eev=5.e-2,relativistic=none,' +
            'compute_analytical_stress=True,sc_accuracy_forces=5.e-2')
