Public API examples
===================

This page provides an **Examples** section for every public symbol exported by the package's ``__all__`` interfaces.

Analysis
--------

``ProteinAnalyzer``::

   from cryptic_ip.analysis import ProteinAnalyzer
   analyzer = ProteinAnalyzer("data/1ZY7.pdb")
   analyzer.score_all_pockets().head()

``PocketScorer``::

   from cryptic_ip.analysis import PocketScorer
   scorer = PocketScorer()

``CrypticSiteMLClassifier`` / ``MLPocketScorer``::

   from cryptic_ip.analysis import CrypticSiteMLClassifier
   model = CrypticSiteMLClassifier()

``StatisticalValidation``::

   from cryptic_ip.analysis import StatisticalValidation
   qvals = StatisticalValidation.benjamini_hochberg([0.01, 0.2, 0.04])

``ComparativeIPAnalysis``::

   from cryptic_ip.analysis import ComparativeIPAnalysis
   analysis = ComparativeIPAnalysis(confidence_level=0.95)

``ElectrostaticsCalculator``::

   from cryptic_ip.analysis import ElectrostaticsCalculator
   calc = ElectrostaticsCalculator()

Database
--------

``ProteomeDownloader``::

   from cryptic_ip.database import ProteomeDownloader
   downloader = ProteomeDownloader("data/structures")

``ProteomeManager``::

   from cryptic_ip.database import ProteomeManager
   manager = ProteomeManager("data/structures/yeast")
   manager.build_catalog().head()

``AlphaFoldBatchDownloader`` / ``ParallelProcessor`` / ``AnalysisCache``::

   from cryptic_ip.database import AlphaFoldBatchDownloader, ParallelProcessor, AnalysisCache

Validation
----------

``validate_adar2``::

   from cryptic_ip.validation import validate_adar2
   report = validate_adar2(use_alphafold=True)

``ValidationSuite``::

   from cryptic_ip.validation import ValidationSuite
   ValidationSuite().run_full_validation()

``OpenMMMDValidationPipeline``::

   from cryptic_ip.validation import OpenMMMDValidationPipeline
   pipeline = OpenMMMDValidationPipeline(output_dir="results/md")

``StructureValidator`` / ``ResultsValidator``::

   from cryptic_ip.validation import StructureValidator, ResultsValidator
   StructureValidator().validate("data/1ZY7.pdb")
