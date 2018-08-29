.. _applications:

***************
Applications
***************

Past Applications
====================

The previous version, GLOFRIM 1.0, was successfully applied for replacing the routing scheme of PCR [Hoch2017a]_ and for benchmarking the
hydrodynamic models DFM and LFP [Hoch2017b]_.

Results of a test case in the Amazon River basin show that by applying the full shallow water equations instead of the kinematic wave approximation 
increases the accuracy of simulated discharge:

.. image:: _images/fig_06.png
   :width: 1500px
   :height: 750px
   :scale: 50 %
   :alt: alternate text
   :align: center

Benchmarking the hydrodynamic models revealed similar performace with respect to simulated discharge:

.. image:: _images/Fig5_qsim_DFMvsLFP.png
   :width: 1500px
   :height: 750px
   :scale: 50 %
   :alt: alternate text
   :align: center

.. todo::
    align layout of figures

However, simulated flood extent differs locally due to the different gridding schemes applied, yielding a low critical success index *C*:

+------------------------+------------------------+------------------------+
| Hit rate               | False alarm ration     | Critical Succes Index  |
+========================+========================+========================+
| 0.85                   | 0.50                   | 0.46                   |
+------------------------+------------------------+------------------------+

Current developments and future applications
============================================

We are planning to some more crazy stuff with GLOFRIM.

These include amongst others:

* Coupling PCR with CMF at a global scale
* Further becnhmarking of models
* Using nested models for compound flood modelling.

Stay tuned!

