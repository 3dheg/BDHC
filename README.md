BDHC
====

Pade fit to the exchange-correlation energy of the 3D homogeneous electron gas

## Use

#### Compiling

#### Functions

DOUBLE PRECISION FUNCTION exc_bdhc(rs,th,xi,exc0,dexc0)

INPUT: rs = Wigner-Seitz radius, scaled by the Bohr radius
       th = Ratio of temperature to Fermi temperature
       xi = Spin-polarization (1 = polarized, 0 = unpolarized)
       exc0 = Exchange-correlation energy for the same system at 0T

RETURNS: exc_bdhc = Exchange-correlation energy (Rydbergs)


DOUBLE PRECISION FUNCTION dexc_bdhc(rs,th,xi,exc0,dexc0)

INPUT: rs = Wigner-Seitz radius, scaled by the Bohr radius
       th = Ratio of temperature to Fermi temperature
       xi = Spin-polarization (1 = polarized, 0 = unpolarized)
       exc0 = Exchange-correlation energy for the same system at 0T
       dexc0 = rs derivative of exc0

RETURNS: dexc_bdhc = rs derivative of exc_bdhc (Rydbergs/rs)

EXAMPLE USAGE:

    PROGRAM main

      USE xc_bdhc_mod
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: rs = 4.0, th = 1.0, xi = 1
      DOUBLE PRECISION, PARAMETER :: exc0 = -0.323254
      DOUBLE PRECISION, PARAMETER :: dexc0 = 0.07655
      DOUBLE PRECISION :: exc, dexc

      exc = exc_bdhc(rs,t,xi,exc0)
      dexc = dexc_bdhc(rs,t,xi,exc0,dexc0)
      write(*,*) exc, dexc

    END PROGRAM main

### Citation

If you use this module in your calculations, please cite both the fit and original Monte Carlo simulation:

  E. W. Brown, J. L. DuBois, M. Holzmann and D. M. Ceperley
  [Exchange-correlation energy for the 3D homogeneous electron gas at arbitrary temperature](http://arxiv.org/abs/1306.1863)
  ArXiv eprints 1306.1863 (cond-mat.str-el)
  Submitted to Phys. Rev. B. Rapid Communications (2013)

  E. W. Brown, B. K. Clark, J. L. DuBois, and D. M. Ceperley
  [Path Integral Monte Carlo simulation of the warm-dense homogeneous electron gas](http://prl.aps.org/abstract/PRL/v110/i14/e146405)  
  Phys. Rev. Lett. 110, 146405 (2013)

Bibtex:

    @ARTICLE{2013arXiv1306.1863B,
       author = {{Brown}, E.~W. and {DuBois}, J.~L. and {Holzmann}, M. and {Ceperley}, D.~M.
            },
        title = "{Exchange-correlation energy for the 3D homogeneous electron gas at arbitrary temperature}",
      journal = {ArXiv e-prints},
    archivePrefix = "arXiv",
       eprint = {1306.1863},
     primaryClass = "cond-mat.str-el",
     keywords = {Condensed Matter - Strongly Correlated Electrons},
         year = 2013,
        month = jun,
       adsurl = {http://adsabs.harvard.edu/abs/2013arXiv1306.1863B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

    @article{PhysRevLett.110.146405,
      title = {Path-Integral Monte~Carlo Simulation of the Warm Dense Homogeneous Electron Gas},
      author = {Brown, Ethan W. and Clark, Bryan K. and DuBois, Jonathan L. and Ceperley, David M.},
      journal = {Phys. Rev. Lett.},
      volume = {110},
      issue = {14},
      pages = {146405},
      numpages = {5},
      year = {2013},
      month = {Apr},
      doi = {10.1103/PhysRevLett.110.146405},
      url = {http://link.aps.org/doi/10.1103/PhysRevLett.110.146405},
      eprint = {1211.6130},
      publisher = {American Physical Society}
    }

