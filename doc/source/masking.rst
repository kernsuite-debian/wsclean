Masks and auto-masking
======================

Masks are often used to be able to clean up to/below the noise level. Deeper cleaning leaves fewer undeconvolved residuals behind, hence the image quality is higher, and it decreases the self-cal bias when used inside a self-cal loop.

Providing a mask
----------------

WSClean accepts masks in CASA format and in ``fits`` file format. A mask is a normal, single polarization image file, where all zero values are interpreted as being not masked, and all non-zero values are interpreted as masked. In the case of a ``fits`` file, the file may either contain a single frequency or it may contain a cube of images. In the first case, the same mask is used for all images, and in the second case, each imaging output channel will be performed with the matching channel in the ``fits`` file. In the second case, the frequency dimension in the ``fits`` file should therefore match the number of imaged output channels.

For spectral imaging (in particular HI obervations), a typical scenario is to use a source detector such as `SoFiA <https://arxiv.org/abs/1501.03906>`_. These kind of spectral source detectors can produce a three dimensional mask where sources have been marked in the masks. This mask can then be used to constrain the deconvolution and perform very deep cleaning.

A basic example:

.. code-block:: bash

    wsclean -fits-mask mymask.fits -niter 1000 -mgain 0.8 \
      -size 1024 1024 -scale 10asec myobservation.ms

A model image can be used as a mask file. After a standard Högbom clean, this can be used to limit the clean components to significant features. For example, one can first clean to 3 sigma (with the :doc:`-auto-threshold option <basic_cleaning>`), then use the model as fitsmask and clean to e.g. 0.3 sigma. One issue with this is that it requires two runs. Also, while this same method can be used for multi-scale cleaning, in general a multi-scale clean will make large parts of the model image non-zero, and cleaning will therefore not be very well constrained. To overcome both issues, WSClean implements a technique called auto-masking.

Auto-masking
------------

Auto-masking allows automated deep cleaning and solves the two problems mentioned above:

 * Only one run of wsclean is required;
 * It maintains scale-dependent masks, which improves multi-scale cleaning.

Auto-masking works in all modes since :doc:`WSClean 2.2 <changelogs/v2.2>`. The general syntax is as follows:

.. code-block:: bash

    wsclean -multiscale -auto-mask 3 -auto-threshold 0.3 \
      -niter 1000000 -mgain 0.8 -size 8000 8000 -scale 2asec \
      obs.ms
    
This will start multi-scale cleaning first towards a threshold of 3 sigma. During multi-scale cleaning, the positions and scale of each component is recorded and put in a scale-dependent mask. Only the center pixel of the kernel is placed in the mask. Once the 3 sigma level is reached, cleaning will continue towards the final threshold of 0.3 sigma. During that stage, the scale-dependent masks are used to constrain the cleaning.

A scale-dependent mask makes sure that when a certain scale-kernel size was cleaned, only that same scale is allowed at that position when the mask is used.

The combination ``-auto-mask 3 -auto-threshold 0.3`` seems like a good general setting, which normally leaves almost no residuals behind.

In cases where the RMS varies strongly over the field of view, for example because of calibration artefacts, it might be useful to use [local RMS thresholding](LocalRMSThresholding), instead of thresholding relative to the global RMS.

Combining auto-masking with a regular mask
------------------------------------------

Auto-masking works in combination with normal masking, so that WSClean can be forced to to find components only in a certain area, and after reaching the auto-masking threshold it would continue with the automask. Since the automask will only contain components within the fitsmask, no components are found outside the fitsmask. Auto-masking and regular masking can be combined by specifying both ``-fits-mask <file>`` and ``auto-mask <threshold>`` on the command line.

Specifying an absolute auto-masking threshold
---------------------------------------------

An absolute auto-masking threshold can be specified with ``-abs-auto-mask <value>``. Like the ``-abs-threshold`` parameter (see :doc:`basic cleaning <basic_cleaning>`), this value is by default in Janskys, but also allows specifying its unit, e.g. ``-abs-auto-mask 12mJy``.

When both ``-auto-mask`` and ``-abs-auto-mask`` are specified, the (highest) threshold that is reached first will stop the stage in which the auto-mask is generated.

Availability
------------

Auto-masking is available in multi-scale since :doc:`WSClean version 2.1 <changelogs/v2.1>`, and is available in all modes since :doc:`WSClean 2.2 <changelogs/v2.2>`. The absolute auto-mask threshold is available since :doc:`WSClean version 3.4 <changelogs/v3.4>`.

References
----------
The auto-masking algorithm is described and demonstrated in `Offringa and Smirnov (2017)  <http://arxiv.org/abs/1706.06786>`_.
