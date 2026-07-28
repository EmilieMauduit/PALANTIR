# About PALANTIR : Prediction Algorithm for star-pLANeT Interactions in Radio

This code has been developed for the construction of an up-to-date and evolutive target catalog, based on observed exoplanet physical parameters, radio emission theory, and magnetospheric physics embedded in scaling laws. It is based on, and extends, previous work by [Griessmeier et al, A&A, 2007.](https://doi.org/10.1051/0004-6361:20077397)

We take advantage of the [exoplanet.eu](https://exoplanet.eu/home/) database to retrieve the maximum parameters possible on every known exoplanet. Using this database allowed us to compute more accurate empirical models, thanks to the increased number of available planets. With this we can compute missing parameters using either empirical models or physical models, to predict the frequency and the radio flux of potential radio emissions.

We consider two types of flow-obstacle interactions :

- Stellar Wind - Planetary Magnetosphere (MS) : these emissions occur at the planet and are therefore limited by the planetary cyclotron frequency, $f_{c,p}^{max} = \frac{eB_p}{2\pi m_e}$. The flow is the stellar wind so we consider the magnetic field it carries and the obstacle is the planetary magnetosphere so we consider its size. We then have the following radio power $P_{mag} = \frac{\beta\pi}{\mu_0}\times v_{eff} B_{SW}^2 R_s^2$,
- Interaction between the host magnetic field and the planet (SPI) : these emissions occur at the star and are therefore limited by the stellar cyclotron frequency, $f_{c,*}^{max} = \frac{eB_*}{2\pi m_e}$. the flow is still the stellar wind of the star, but this time the obstacle is an unmagnetized planet, so we consider the size of its ionosphere. We then have the following radio power $P_{SPI} = \frac{\beta\pi}{\mu_0}\times v_{eff} B_{SW}^2 R_{iono}^2$.

Using PALANTIR, we prepared an updated list of targets of interest for radio emissions. Additionally, we compare our results with previous studies conducted with similar models [Griessmeier, Planetary Radio Emissions VIII, 2017](https://doi.org/10.1553/PRE8s285). 
For the next steps, we aim at improving this code by adding new models and updating those already used. 
There are three papers related to this work, two published and one in writting, along with a PhD manuscript (in French) :

- Mauduit et al, 'PALANTIR: An updated prediction tool for exoplanetary radio emissions', 2023, PRE IX, https://doi.org/10.25546/103092
- Mauduit Emilie, 'Méthodes pour la détection d'exoplanètes en ondes radio basses fréquences : sélection de cibles, identification de contaminations, méthodes de détection et applications à Jupiter', 2024, https://theses.hal.science/tel-04821784v1 
- Duchêne et al, 'Stellar Magnetic Fields and Radio Emissions from Star–Planet Systems', 2026, PRE X, https://doi.org/10.25935/prex-rrvx 
- Mauduit et al, 'PALANTIR: Predicting Star-Planet Interactions in Radio', 2026, to be submitted to Astronomy & Astrophysics

# Installation

Open your terminal. Go to the folder in which you want to use this algorithm. Copy-paste the https link in the 'Clone' button. Clone the repository and initialize the submodules present ([magfieldprediction](https://gitlab.obspm.fr/qduchene/mag_field_prediction)) :

```bash
cd existing_repo
git clone https://github.com/EmilieMauduit/PALANTIR.git
cd PALANTIR
git submodule update --init
```

Then you should load your desired python environment (it should be python3.9 or 3.10, the version of tensorflow used is not compatible with more recent version of python) and install this package and the submodule with :
```bash
pip install .
cd src/palantir/prediction_tools/mag_field_prediction 
pip install .
```

*Warning* If you use an older MacOs system based on Intel processor you might need to manually install some of the dependencies before running the previous installation steps :
```bash
pip install --only-binary=:all: cryptography -v
```

In order to use it from any folder, add the path to the module to your `PYTHONPATH` in your `/.bashrc`.

```bash
export PYTHONPATH=$PYTHONPATH:/mypath/palantir/src
```
Now you can import this package from anywhere.

If updates are pushed on the main branch, you can add them by running the following lines in the repository of the code, inside your python environment:

```bash
git pull
pip install .
```

# Support

Any issue encountered while using this software can be raised on the Issue page of this project. The author will gladly look at them and try to find a fix.

# Contributing

Any suggestion on a missing feature, a way to improve the code will be welcomed. Please do not hesitate to reach out.

# Cite this work

[![DOI](https://zenodo.org/badge/480817410.svg)](https://doi.org/10.5281/zenodo.15599637)

# Authors and acknowledgment

This work has made use of the Extrasolar Planet Encyclopaedia (exoplanet.eu) maintained by [J. Schneider et al](https://doi.org/10.1051/0004-6361/201116713). Philippe Zarka acknowledges funding from the ERC under the European Union's Horizon 2020 research and innovation programme (grant agreement no. 101020459 - Exoradio). Emilie Mauduit acknowledges the precious help of A. Loh for Python usage. This work was supported by the Programme National de Planétologie (PNP) of CNRS/INSU co-funded by CNES and by the Programme National de Physique Stellaire (PNPS) of CNRS/INSU co-funded by CEA and CNES.

# License

[MIT](LICENSE)

# Project status

The project is active. The main branch is maintained up-to-date with the latest working and validated features. Other branches are reserved for developement and are not guaranteed to be working.
