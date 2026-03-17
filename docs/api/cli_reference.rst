CLI command reference
=====================

``cryptic-ip`` root command
---------------------------

.. autofunction:: cryptic_ip.cli.main

``cryptic-ip analyze``
----------------------

.. autofunction:: cryptic_ip.cli.analyze

Example::

   cryptic-ip analyze data/1ZY7.pdb --output results/adar2.csv --score-threshold 0.65

``cryptic-ip validate``
-----------------------

.. autofunction:: cryptic_ip.cli.validate

Example::

   cryptic-ip validate --structure data/1ZY7.pdb --results results/adar2.csv

``cryptic-ip validate-suite``
-----------------------------

.. autofunction:: cryptic_ip.cli.validate_suite

Example::

   cryptic-ip validate-suite

``cryptic-ip download``
-----------------------

.. autofunction:: cryptic_ip.cli.download

Example::

   cryptic-ip download yeast --data-dir data/structures

``cryptic-ip screen``
---------------------

.. autofunction:: cryptic_ip.cli.screen

Example::

   cryptic-ip screen data/structures/yeast --output results/yeast_hits.csv --max-structures 200 --use-ml-model

``cryptic-ip md-validate``
--------------------------

.. autofunction:: cryptic_ip.cli.md_validate

Example::

   cryptic-ip md-validate results/combined_hits.csv --output-dir results/md_validation --top-n 20
