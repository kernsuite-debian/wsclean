#! /usr/bin/python

# This file demonstrates a simple deconvolution approach, which removes the strongest
# components one by one (similar to Clean with a gain of 1).

# run this e.g. with:
# wsclean -niter 10000 -save-first-residual -auto-threshold 3 -mgain 0.8 -interval 10 11 -size 1024 1024 -scale 1amin -python-deconvolution ~/projects/wsclean/scripts/deconvolution-example.py 1052736496-averaged.ms/

import numpy


def deconvolve(residual, model, psf, meta):
    nchan, npol, height, width = residual.shape
    print(
        "Python deconvolve() function was called for "
        + f"{width} x {height} x {npol} (npol) x {nchan} (chan) dataset"
    )

    # residual and model are numpy arrays with dimensions nchan x npol x height x width
    # psf is a numpy array with dimensions nchan x height x width

    # This file doesn't support multiple channels or polarizations:
    if nchan != 1 or npol != 1:
        raise NotImplementedError("nchan and npol must be one")

    # meta contains several useful meta data:
    # meta.channels is an array, each element has properties 'frequency' and 'weight'
    # meta.square_joined_channels: boolean, true if joined channels are to be squared.
    # meta.spectral_fitter is an object that knows how the user requested spectral fitting,
    # and is also aware of the channel frequencies / weights. It
    # has methods 'fit' and 'fit_and_evaluate'. Example:
    #
    # values = numpy.zeros(nchan, dtype=numpy.float64)
    # coefficients = meta.spectral_fitter.fit(values, x, y)
    # or:
    # values = meta.spectral_fitter.fit_and_evaluate(values, x, y)
    #
    # Furthermore, the meta class has properties major_iter_threshold, final_threshold,
    # mgain, iteration_number and max_iterations. These are demonstrated below.
    # iteration_number can be modified, and will keep its value when deconvolve()
    # is called again.

    # find the largest peak
    peak_index = numpy.unravel_index(numpy.argmax(residual), residual.shape)
    peak_value = residual[peak_index]

    mgain_threshold = peak_value * (1.0 - meta.mgain)
    first_threshold = max(
        meta.major_iter_threshold, meta.final_threshold, mgain_threshold
    )

    print(
        f"Starting iteration {meta.iteration_number}, peak={peak_value}, first threshold={first_threshold}"
    )
    while (
        peak_value > first_threshold
        and meta.iteration_number < meta.max_iterations
    ):
        model[peak_index] += peak_value

        psf_shift = (peak_index[2] + height // 2, peak_index[3] + width // 2)
        residual = residual - peak_value * numpy.roll(
            psf, psf_shift, axis=(1, 2)
        )

        peak_index = numpy.unravel_index(
            numpy.argmax(residual), residual.shape
        )
        peak_value = residual[peak_index]

        meta.iteration_number = meta.iteration_number + 1

    print(
        f"Stopped after iteration {meta.iteration_number}, peak={peak_value}"
    )

    # Fill a dictionary with values that wsclean expects:
    result = dict()
    result["residual"] = residual
    result["model"] = model
    result["level"] = peak_value
    result["continue"] = (
        peak_value > meta.final_threshold
        and meta.iteration_number < meta.max_iterations
    )

    print("Finished deconvolve()")
    return result
