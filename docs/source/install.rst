Installation
============

Installing with pip
-------------------

``Python3.8`` or higher version is required and no support will be provided for ``Python2``.
In your desired python environment :

.. code-block:: bash

   $ cd existing_repo
   $ git clone https://github.com/EmilieMauduit/PALANTIR

To get the latest working features, we recommend only using the main branch as the other branches are reserved for development.
Then you should load your desired python environment and install this package with :

.. code-block:: bash

   $ pip install .

Updates need to be regularly checked for as this package is still in developpment. To update :

.. code-block:: bash

   $ cd existing_repo/PALANTIR
   $ git pull
   $ pip install .

.. note:: 

   There is no support for a ``conda`` insatallation yet.

Load ``palantir`` as a package from anywhere
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the package was installed by cloning this repository, you can modify your ``PYTHONPATH`` to be able to import this package from anywhere on your machine. To do so, add this line in your ``.bashrc`` :

.. code-block:: bash

   $ export PYTHONPATH=$PYTHONPATH:/mypath/PALANTIR/src

Dependencies
------------

If not already installed, these dependencies will be installed along with the package when you first install it.

* `astropy <'https://docs.astropy.org/en/stable/'>`_
* `numpy <'https://numpy.org/doc/stable/'>`_
* `python <'http://docs.python.org/3'>`_
* `matplotlib <'https://matplotlib.org/stable/'>`_
* `h5py <'https://docs.h5py.org/en/stable/'>`_
* `pytest <'https://docs.pytest.org/en/stable/'>`_
* `scipy <'https://docs.scipy.org/doc/scipy/'>`_
* `pandas <'https://pandas.pydata.org/docs/'>`_
* `pyvo <'https://pyvo.readthedocs.io/en/latest/'>`_
* `atroquery <''>`_
