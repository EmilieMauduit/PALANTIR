# About PALANTIR : Prediction Algorithm for star-pLANeT Interactions in Radio

This code has been developed for the construction of an up-to-date and evolutive target catalog, based on observed exoplanet physical parameters, radio emission theory, and magnetospheric physics embedded in scaling laws. It is based on, and extends, previous work by [Griessmeier et al, A&A, 2007.](https://doi.org/10.1051/0004-6361:20077397)

We take advantage of the [exoplanet.eu](https://exoplanet.eu/home/) database to retrieve the maximum parameters possible on every known exoplanet. Using this database allowed us to compute more accurate empirical models, thanks to the increased number of available planets. With this we can compute missing parameters using either empirical models or physical models, to predict the frequency and the radio flux of potential radio emissions.

We consider two types of flow-obstacle interactions :

- Stellar Wind - Planetary Magnetosphere (MS) : these emissions occur at the planet and are therefore limited by the planetary cyclotron frequency, $f_{c,p}^{max} = \frac{eB_p}{2\pi m_e}$. The flow is the stellar wind so we consider the magnetic field it carries and the obstacle is the planetary magnetosphere so we consider its size. We then have the following radio power $P_{mag} = \frac{\beta\pi}{\mu_0}\times v_{eff} B_{SW}^2 R_s^2$,
- Interaction between the host magnetic field and the planet (SPI) : these emissions occur at the star and are therefore limited by the stellar cyclotron frequency, $f_{c,*}^{max} = \frac{eB_*}{2\pi m_e}$. the flow is still the stellar wind of the star, but this time the obstacle is an unmagnetized planet, so we consider the size of its ionosphere. We then have the following radio power $P_{SPI} = \frac{\beta\pi}{\mu_0}\times v_{eff} B_{SW}^2 R_{iono}^2$.

Using PALANTIR, we prepared an updated list of targets of interest for radio emissions. Additionally, we compare our results with previous studies conducted with similar models [Griessmeier, Planetary Radio Emissions VIII, 2017](https://doi.org/10.1553/PRE8s285). 
For the next steps, we aim at improving this code by adding new models and updating those already used. 
There are two papers related to this work, one published and one in writting, along with a PhD manuscript (in French) :

- Mauduit et al, 'PALANTIR: An updated prediction tool for exoplanetary radio emissions', 2023, PRE IX, https://doi.org/10.25546/103092
- Mauduit Emilie, 'Méthodes pour la détection d'exoplanètes en ondes radio basses fréquences : sélection de cibles, identification de contaminations, méthodes de détection et applications à Jupiter', 2024, https://theses.hal.science/tel-04821784v1 
- Mauduit et al, 'PALANTIR: Predicting Star-Planet Interactions in Radio', 2025, to be submitted to Astronomy & Astrophysics

# Installation

Open your terminal. Go to the folder in which you want to use this algorithm. Copy-paste the https link in the 'Clone' button.

```bash
cd existing_repo
git clone https://github.com/EmilieMauduit/PALANTIR.git
```

Then you should load your desired python environment and install this package with :
```bash
pip install .
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

# Usage

## 1 - Defining your parameters for the predictions

The code relies on a variety of models for prediction that can be chosen as wished by the user. These parameters have to be defined in the `run_parameters.csv` file that can be found in the following folder : `~/src/palantir/scripts/input_files/`. The user should only change the last column of this file. In this column, a $0$ means this setting won't be used and $1$ that it will be used, for some parameters real values or a string are expected. We detail this in this sub-section.

- ### Choosing the database :

PALANTIR can use three different databases, it is not yet capable of combining multiple databases. The user can choose between the database from [exoplanet.eu](https://exoplanet.eu/home/), which is the one used in every study published by the author, the NASA exoplanet database, or a custom one.

In each case, a version of the database is present in `~/src/palantir/scripts/input_files/`. The author tries to keep these version updated, but they strongly encourage the user to go onte the corresponding website, dowload and save the file in the above folder to make sure they run with the latest version. As for the custom catalog, an example is also provided in the same folder, the user have to make sure it has the given structure, else the code will not be able to read it.

In the parameter file, the user can choose the database by putting 0 or 1, in the value column of the following setting :

```bash
setting;value
nasa_data;0
exoplanet_data;1
custom_data;0
```

With this, the chosen data base is the one from [exoplanet.eu](https://exoplanet.eu/home/). Note that if more than one is set to `1`, an error will be encountered.

- ### Choose the planetary magnetic moment models :

Different models are available to predict the planetary magnetic moment, $\mathcal{M}_p$ with respect to Jupiter's. They are detailed in [Mauduit, 2024](https://theses.hal.science/tel-04821784v1).
Several models can be selected at the same time, the geometric mean value of every result will be used. As presented in the above citation, there a two types of models available : some that only use the characteristics of the planetary dynamo region and one (`rein-chris`) that additionally uses the apparent luminosity of the planet. When using the latest, we strongly recommend to only use this one because the estimations made with this are one order of magnitude bigger that with the other models, therefore it is better to consider it seprately, as presented in [Mauduit, 2024](https://theses.hal.science/tel-04821784v1).

Here are the two recommended settings :

```bash
setting;value
.
.
blackett;0
busse;1
rein-chris;0
rein-chris-dyn;0
mizu_slow;1
mizu_mod;1
sano;1
```
or 
```bash
setting;value
.
.
blackett;0
busse;0
rein-chris;0
rein-chris-dyn;1
mizu_slow;0
mizu_mod;0
sano;0
```

There are two possibilities for the [Reiners-Christensen]() model, we recommend using the "dynamo" one. In their paper they use an approximated factor to estimate the size of the dynamo region, that is defined so that it is $0.83$ in the case of Jupiter. Since PALANTIR also allows to predict the size of this region, we directely use our prediction in their formula.

- ### Dynamo region density models :

Two models are available based on the litterature, however we recommend the use of the model developed by [Lane-Emden](), since it is more recent and more precise.

```bash
setting;value
.
.
dens_lin;0
dens_laneEmden;1
```

- ### Critical density of the dynamo region, $\rho_{crit}$ :

This is the density below which we consider the core of the planet is not dens enough to produce a magnetic field. It can be set to either $700$ or $1000$ $g.cm^{-3}$.

```bash
setting;value
.
.
rho_crit;700
```

Default value is $700\ g.cm^{-3}$ and this the value we recommend using.

- ### Planetary radius models :

When the planetary radius is not known, we have several models that can be used to predict it. Either the original model present in [Griessmeier et al, A&A, 2007.](https://doi.org/10.1051/0004-6361:20077397) or a new empirical polynomial model can be used. A new model for the computation of the radius expansion, due to the closeness of the planet to its host, is also available. It was also computed empirically with data from the [exoplanet.eu](https://exoplanet.eu/home/) database.

```bash
setting;value
.
.
radius_original;1
radius_polyfit;0
radius_expansion;1
```

- ### Stellar radius models :

For the stellar radius, only one model is available at the moment, it is the mass dependent model presented in [Tout et al, 1996](https://doi.org/10.1093/mnras/281.1.257).

```bash
setting;value
.
.
Tout;1
```

- ### Planetary apperent luminosity, $L_p$ :

As explained above, one of the planetary magnetic moment model requires this parameter. It can be computed from tables. The original one used in [Reiners & Christensen, 2010](https://doi.org/10.1051/0004-6361/201014251) is the one from [Burrows et al, 1997](https://doi.org/10.1086/305002). We added the possibility to use more recent tables from [Baraffe et al, 2008](https://doi.org/10.1051/0004-6361:20079321), that takes, or not, the irradiation from the host star. We recommend using the following parameters :

```bash
setting;value
.
.
Burrows;0
Baraffe_irrad;0
Baraffe_noirrad;1
```

- ### Stellar type code :

This is to set the type of stars to be kept for the prediction. Since most of our models are based on main sequence stars, we recommend using $1$, to keep in the domain of validity of our models.

```bash
setting;value
.
.
sp_type_code;1
```

- ### Stellar magnetic field models :

The prediction of the stellar magnetic field can be done with several models. Y

```bash
setting;value
.
.
Bstar_original;0
Bstar_polyfit;1
Bstar_catalog_only;0
```

- ### Verbose :

This is to decide if you want the code to print information while running. In any case, every logging information is saved into the `palantir_info.log` file, with the rest of the output files.

```bash
setting;value
.
.
talk;0
```

- ### Output path :
The user can specify the path to the folder where the outputs should be saved. 
The pipeline will create a new folder within the given one, with a name of the format `2025-06-25T18h59`, that will correspond to the date and time at which the pipeline started to run. 
```bash
setting;value
.
.
path_outputs;/mypath/myfolder
```

## 2 - Starting the prediction run

```bash
python -m palantir.scripts.main.py 
```

## 3 - Output files :

The following files will be saved in the created folder at the given `output_path` :

  - `catalog_input.csv` : the catalog used will be copied in this file, in order to be able to reproduce results afterwards (values for a given target can be updated with time),
  - `skipped_targets.txt`: a text file containing the list of skipped targets, and the reason why they were skipped,
  - `main_output.csv` : a table containing every predicted parameters, along with the main ones coming from the database,
  - `palantir_info.log`: a logging file, containing every logs of the run. At the begining of the file, there is the version of PALANTIR used for the run, along with a summary of the configuration parameters given for this run (models used, etc)

## 4 - Reading the results 

Once PALANTIR is applied to the given database, the output is saved in a `.csv` file. It can be read and navigated in with the following :

```python
from palantir.output_analysis.data_manipulation import DataManipulation

data = DataManipulation.from_file(
        filename="/mypath/main_output.csv",
        instrument_name='NenuFAR'
        )
```

It can also be visualized with basic `.csv` file readers such as Microsoft Excel, Apple Numbers or LibreOffice.

The `instrument_name` parameters allows to define which instrument will be used for the observation and therefore select only the targets with `ra/dec` coordinates observable with this telescope.

- ### Make frequency-flux plots :

The main goal of these prediction is to obtain the frequency and the flux of the possible radio emissions. Plots of $\Phi_{radio} = f(f_c^{max})$ can be done automatically with the `plot_frequency_flux` function. 
The user can specify which type of interaction to consider ('MS' or 'SPI') and which instruments' sensitivity to overplot.

If `figname`is not given, the plot will not be saved and only displayed. The axis limits can be specified, but if not given they are computed as $0.9 \times min$ and $1.1 \times max$.

```python
data.plot_frequency_flux(
    interaction='MS', 
    instruments = ['NenuFAR', 'LOFAR low', 'LOFAR high','SKA1 low','SKA2 low','GMRT','VLA','UTR-2'],
    figname = '/mypath/Plot_Pmag_fp.png', 
    xmin = 1e-2, xmax = 3e4, 
    ymin = 1e-9, ymax = 1e6)
```

- ### Investigate other parameters dependancy :

It can be interesting to see how one or two parameters evolve with respect to another one. The function allows to display such dependencies in a scatter plot.

```python
data.plot_quantities(x='star_age', y='star_mass', z='star_magfield', zmin = 0, zmax = 1e-3)
```
If only `x` and `y`are given it will be a regular scatter plot, if `z`is also given it will be represented by a colorbar. The ranges and scales for `x,y,z` can be specified as `kwargs`.

- ### Select targets based on various criteria :

For now, only two criteria for selection are available with the function `target_selection`. You can select on instrument sensitivity (default), or specify a custom minimum frequency and flux for your selection. If you specify only one of the parameter, the other one will be set to $0$. This will give you a new `DataManipulation` object.

```python
data_filtered = data.target_selection(fc_min_MHz = 5., flux_min_mJy = 10.)
```

We plan on adding the possibility to select on observability with the chosen telescope.

# Tutorials

A variety of tutorials, in the form of JupyterNotebook files are available for the user to get familiar with the various modules present in this package.

# Support

Any issue encountered while using this software can be raised on the Issue page of this project. The author will gladly look at them and try to find a fix.

# Contributing

Any suggestion on a missing feature, a way to improve the code will be welcomed. Please do not hesitate to reach out.

# Authors and acknowledgment

This work has made use of the Extrasolar Planet Encyclopaedia (exoplanet.eu) maintained by [J. Schneider et al](https://doi.org/10.1051/0004-6361/201116713). Philippe Zarka acknowledges funding from the ERC under the European Union's Horizon 2020 research and innovation programme (grant agreement no. 101020459 - Exoradio). Emilie Mauduit acknowledges the precious help of A. Loh for Python usage. This work was supported by the Programme National de Planétologie (PNP) of CNRS/INSU co-funded by CNES and by the Programme National de Physique Stellaire (PNPS) of CNRS/INSU co-funded by CEA and CNES.

# License

[MIT](LICENSE)

# Project status

The project is active. The main branch is maintained up-to-date with the latest working and validated features. Other branches are reserved for developement and are not guaranteed to be working.
