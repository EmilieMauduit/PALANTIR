.. palantir documentation master file, created by
   sphinx-quickstart on Mon Jun 30 16:26:15 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PALANTIR documentation !
=================================================

This code has been developed for the construction of an up-to-date and evolutive target catalog, based on observed exoplanet physical parameters, radio emission theory, and magnetospheric physics embedded in scaling laws. It is based on, and extends, previous work by [Griessmeier et al, A&A, 2007.](https://doi.org/10.1051/0004-6361:20077397)

We take advantage of the [exoplanet.eu](https://exoplanet.eu/home/) database to retrieve the maximum parameters possible on every known exoplanet. Using this database allowed us to compute more accurate empirical models, thanks to the increased number of available planets. With this we can compute missing parameters using either empirical models or physical models, to predict the frequency and the radio flux of potential radio emissions.

We consider two types of flow-obstacle interactions :

- Stellar Wind - Planetary Magnetosphere (MS) : these emissions occur at the planet and are therefore limited by the planetary cyclotron frequency, $f_{c,p}^{max} = \frac{eB_p}{2\pi m_e}$. The flow is the stellar wind so we consider the magnetic field it carries and the obstacle is the planetary magnetosphere so we consider its size. We then have the following radio power $P_{mag} = \frac{\beta\pi}{\mu_0}\times v_{eff} B_{SW}^2 R_s^2$,
- Interaction between the host magnetic field and the planet (SPI) : these emissions occur at the star and are therefore limited by the stellar cyclotron frequency, $f_{c,*}^{max} = \frac{eB_*}{2\pi m_e}$. the flow is still the stellar wind of the star, but this time the obstacle is an unmagnetized planet, so we consider the size of its ionosphere. We then have the following radio power $P_{SPI} = \frac{\beta\pi}{\mu_0}\times v_{eff} B_{SW}^2 R_{iono}^2$.

Using PALANTIR, we prepared an updated list of targets of interest for radio emissions. Additionally, we compare our results with previous studies conducted with similar models [Griessmeier, Planetary Radio Emissions VIII, 2017](https://doi.org/10.1553/PRE8s285). 
For the next steps, we aim at improving this code by adding new models and updating those already used. 
There are three papers related to this work, one published and one in writting, along with a PhD manuscript (in French) :

- Mauduit et al, 'PALANTIR: An updated prediction tool for exoplanetary radio emissions', 2023, PRE IX, https://doi.org/10.25546/103092
- Mauduit Emilie, 'Méthodes pour la détection d'exoplanètes en ondes radio basses fréquences : sélection de cibles, identification de contaminations, méthodes de détection et applications à Jupiter', 2024, https://theses.hal.science/tel-04821784v1 
- Duchêne et al, 'Stellar Magnetic Fields and Radio Emissions from Star–Planet Systems', 2026, PRE X, https://doi.org/10.25935/prex-rrvx
- Mauduit et al, 'PALANTIR: Predicting Star-Planet Interactions in Radio', 2026, to be submitted to Astronomy & Astrophysics


.. note::
   By default, `logging <https://docs.python.org/3/library/logging.html>`_ is set to ``DEBUG`` level. However, this can be changed dynamically by the user, for e.g.:

   .. code-block:: python

      >>> import palantir
      >>> import logging
      >>> logging.getLogger('palantir').setLevel(logging.INFO)

.. note::
   Users are most welcome to signal bugs or ask for complementary functionalities via sending a 
   `GitLab issue <https://github.com/EmilieMauduit/PALANTIR/issuess>`_.

.. note::
   [![DOI](https://zenodo.org/badge/480817410.svg)](https://doi.org/10.5281/zenodo.15599637)

   BibTeX citation:

   .. code-block:: bash

      @software{emiliemauduit_2025_15600127,
      author       = {EmilieMauduit},
      title        = {EmilieMauduit/PALANTIR: Release v0.1.1},
      month        = jun,
      year         = 2025,
      publisher    = {Zenodo},
      version      = {v0.1.1},
      doi          = {10.5281/zenodo.15600127},
      url          = {https://doi.org/10.5281/zenodo.15600127},
      swhid        = {swh:1:dir:97ea3dbf4f5e732292dcf1f9ab592503e9e00635
                        ;origin=https://doi.org/10.5281/zenodo.15599637;vi
                        sit=swh:1:snp:0a9f61d3b5b345f321a756000d9143f1cace
                        33fe;anchor=swh:1:rel:1161f29eb661a467778dc05b2ef6
                        23fb19aaa947;path=EmilieMauduit-PALANTIR-aa71d2d
                        },
      }

.. _getting-started:

.. toctree::
   :caption: Getting Started
   :maxdepth: 1

   install

.. _prediction-tools:

.. toctree::
   :caption: Prediction tools
   :maxdepth: 1

   prediction_tools/planet
   prediction_tools/star
   prediction_tools/stellar_wind
   prediction_tools/dynamo_region
   prediction_tools/magnetic_moment
   prediction_tools/emission

.. _output-analysis:

.. toctree::
   :caption: Output analysis
   :maxdepth: 1

   output_analysis/data_manipulation

.. _prediction-pipeline:

.. toctree::
   :caption: Prediction pipeline
   :maxdepth: 1

   scripts/main_prediction

.. _contents:

.. toctree::
   :maxdepth: 3
   :caption: Contents

   modules

.. _changelog:

.. toctree::
   :caption: Changelog
   :maxdepth: 1

   changelog


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
