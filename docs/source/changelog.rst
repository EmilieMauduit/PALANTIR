Changelog
=========

0.3.0 (2026-03-20)
^^^^^^^^^^^^^^^^^^

* Renamed the package to `spyburst`.
* Implemented various tests.
* Edited the documentation and tutorials.
* The NenuFAR :module:`~spyburst.instruments.NenuFAR` module can work with both GPU or CPU depending on what is available on the machine. It is fully compatible with the MurMuRe pipeline.

0.2.0 (2025-01-11)
^^^^^^^^^^^^^^^^^^

Refactoring of the instrument related functions.

* Creation of abstract classes to use the processing pipeline with any instrument. :class:`~spyburst.instruments.instrument_utils.ReadInstrumentData`, :class:`~spyburst.instruments.instrument_utils.ProcessInstrumentData` and :class:`~spyburst.instruments.instrument_utils.InstrumentDataPipeline`.
* Modification of classes concerning the NDA :module:`~spyburst.instruments.NDA` and NenuFAR :module:`~spyburst.instruments.NenuFAR`, they can be used as example for other instruments.
* Fixed bugs in the computation of `slicing_parameters` for the NDA.

0.1.0 (2025-01-05)
^^^^^^^^^^^^^^^^^^

First version of the taska_a2 package.

* Complete migration from IDL to Python, with succesfull tests.
* Implementation of module :class:`~spyburst.instruments.NDA.read_JunoN_data.ReadJunonData` for the NDA, to read and process its high-resolution data.
* Added a module for output reading and analysis.
* Successfully recomputed probability maps for Jovian radio emissions that were first obtained from the catalogs associated to [`Marques et al, 2017, A&A <10.1051/0004-6361/201630025>`_] that are availble on [`Vizier <https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/604/A17>`_].