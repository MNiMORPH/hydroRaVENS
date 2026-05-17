Data-Driven Priors
==================

:func:`~hydroravens.suggest_priors` combines
:class:`~hydroravens.BrutsaertNieber` recession analysis with
:class:`~hydroravens.HydrographSeparation` to produce a coherent set of
parameter starting points before any model run or calibration.

See the :doc:`tutorial` for a full worked example.

.. autofunction:: hydroravens.suggest_priors

.. autoclass:: hydroravens.Priors
   :members: summary, to_yaml_snippet
   :member-order: bysource
