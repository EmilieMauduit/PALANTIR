Changelog
=========

1.0.0 (2026-09-01)
^^^^^^^^^^^^^^^^^^

Second official release for the code paper.

* New models for stellar mass-loss rate prediction.
* New flag system encoded in bits to identify when strong assumptions are made.
* Implemented unit tests.
* Added the use of astropy units to avoid unhomogeneous equations.
* Added complete documentation on readthedocs.


0.2.0 (2026-03-20)
^^^^^^^^^^^^^^^^^^
Major refactoring in the frame of SKA SPI interaction chapter and PRE X proceeding paper.

* New models for coronal temperature, stellar rotation period, 
* 


0.1.1 (2025-05-06)
^^^^^^^^^^^^^^^^^^

First release on zenodo.

* Creation of abstract classes to use the processing pipeline with any instrument. :class:`~spyburst.instruments.instrument_utils.ReadInstrumentData`, :class:`~spyburst.instruments.instrument_utils.ProcessInstrumentData` and :class:`~spyburst.instruments.instrument_utils.InstrumentDataPipeline`.
* Modification of classes concerning the NDA :module:`~spyburst.instruments.NDA` and NenuFAR :module:`~spyburst.instruments.NenuFAR`, they can be used as example for other instruments.
* Fixed bugs in the computation of `slicing_parameters` for the NDA.

0.1.0 (2025-01-05)
^^^^^^^^^^^^^^^^^^

First version of the palantir package used to produce the results in Mauduit Emilie, 'Méthodes pour la détection d'exoplanètes en ondes radio basses fréquences : sélection de cibles, identification de contaminations, méthodes de détection et applications à Jupiter', 2024, https://theses.hal.science/tel-04821784v1 
