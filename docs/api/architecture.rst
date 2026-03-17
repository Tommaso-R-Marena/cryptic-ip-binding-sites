Architecture diagram
====================

Module relationships used by the pipeline:

.. code-block:: text

   cryptic_ip.cli
      ├── cryptic_ip.analysis
      │    ├── analyzer.py
      │    ├── scorer.py
      │    ├── ml_classifier.py
      │    └── comparative_analysis.py
      ├── cryptic_ip.database
      │    ├── downloader.py
      │    ├── manager.py
      │    └── batch_processing.py
      └── cryptic_ip.validation
           ├── adar2.py
           ├── structure_validator.py
           ├── results_validator.py
           └── md_validation.py
