from wsgiref import validate
import pytest
import os, glob
import sys
from astropy.io import fits
import numpy as np
from utils import validate_call, compute_rms


# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf

"""
Test script containing a collection of wsclean commands, tested on big MWA/SKA
measurement sets. Tests contained in this file can be invoked via various routes:

- execute "make longsystemcheck"  in your build directory
- execute "[python3 -m] pytest [OPTIONS] source/long_system_checks.py::TestLongSystem::<test_name.py>" in your build/tests/python directory
"""


def name(name: str):
    return os.path.join(tcf.RESULTS_DIR, name)


@pytest.mark.usefixtures("prepare_large_ms")
class TestLongSystem:
    """
    Collection of long system tests.
    """

    def test_dirty_image(self):
        # Make dirty image
        s = f"{tcf.WSCLEAN} -name {name('test-dirty')} {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_clean_rectangular_unpadded_image(self):
        # Clean a rectangular unpadded image
        s = f"{tcf.WSCLEAN} -name {name('clean-rectangular')} -padding 1 -local-rms \
              -auto-threshold 5 -mgain 0.8 -niter 100000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_automask_multiscale_clean(self):
        # Auto-masked multi-scale clean
        s = f"{tcf.WSCLEAN} -name {name('multiscale-automasked')} -auto-threshold 0.5 -auto-mask 3 \
              -mgain 0.8 -multiscale -niter 100000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multiple_intervals(self):
        # Multiple intervals
        s = f"{tcf.WSCLEAN} -name {name('intervals')} -intervals-out 3 \
            {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multiple_intervals_and_channels(self):
        # Multiple intervals + multiple channels with some cleaning
        s = f"{tcf.WSCLEAN} -name {name('intervals-and-channels')} -intervals-out 3 \
            -channels-out 2 -niter 1000 -mgain 0.8 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multiple_intervals_and_facets(self):
        # Multiple intervals + multiple facets with some cleaning
        s_base = f"{tcf.WSCLEAN} -name {name('intervals-and-facets')} -intervals-out 3 \
            -facet-regions {tcf.FACETFILE_4FACETS}"
        s = f"{s_base} -niter 1000 -mgain 0.8 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        validate_call(s.split())

        # Run predict, using the model generated above.
        s = f"{s_base} -predict {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multifrequency_hogbom(self):
        # Multi-frequency Högbom clean, no parallel gridding
        s = f"{tcf.WSCLEAN} -name {name('mfhogbom')} -channels-out 4 -join-channels -auto-threshold 3 \
            -mgain 0.8 -niter 1000000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multifrequency_without_joining_pol(self):
        # Multi-frequency clean, no joining of pols (reproduces bug #128)
        s = f"{tcf.WSCLEAN} -name {name('mf-no-join-pol')} -pol iv -channels-out 2 -join-channels -niter 1 -interval 10 13 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_multifrequency_hogbom_spectral_fit(self):
        # Multi-frequency Högbom clean with spectral fitting
        s = f"{tcf.WSCLEAN} -name {name('mfhogbom-fitted')} -channels-out 4 -join-channels -parallel-gridding 4 \
           -fit-spectral-pol 2 -auto-threshold 3 -mgain 0.8 \
               -niter 1000000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_mutifrequency_multiscale_parallel(self):
        # Multi-frequency multi-scale clean with spectral fitting, pallel gridding & cleaning
        s = f"{tcf.WSCLEAN} -name {name('mfms-fitted')} -channels-out 4 -join-channels -parallel-gridding 4 \
             -parallel-deconvolution 1000 -fit-spectral-pol 2 -multiscale -auto-threshold 0.5 \
                  -auto-mask 3 -mgain 0.8 -niter 1000000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_save_components(self):
        test_name = name("mfms-components")

        # Remove old component list if it exists
        component_file = test_name + "-sources.txt"
        if os.path.exists(component_file):
            os.remove(component_file)

        # Save the list of components
        s = f"{tcf.WSCLEAN} -name {test_name} -save-source-list -channels-out 4 \
            -join-channels -parallel-gridding 4 -fit-spectral-pol 2 \
                -auto-threshold 0.5 -auto-mask 3 -mgain 0.8 -niter 1000000 \
                    -multiscale -parallel-deconvolution 1000 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        validate_call(s.split())

        # Check whether source files is generated
        assert os.path.isfile(component_file)

    def test_linear_joined_polarizations(self):
        # Linear joined polarizations with 4 joined channels
        s = f"{tcf.WSCLEAN} -name {name('linearpol')} -niter 1000000 -auto-threshold 3.0 \
             -pol XX,YY,XY,YX -join-polarizations -join-channels -mgain 0.85 \
                 -channels-out 4 -parallel-gridding 16 -gridder wstacking {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_two_timesteps(self):
        # Image two timesteps
        s = f"{tcf.WSCLEAN} -name {name('two-timesteps')} -niter 1000000 -auto-threshold 3.0 \
            -intervals-out 2 -interval 20 22 -mgain 0.85 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_stop_on_negative_components(self):
        # Stop on negative components
        s = f"{tcf.WSCLEAN} -name {name('stop-on-negatives')} -stop-negative -niter 100000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_save_imaging_weights(self):
        s = f"{tcf.WSCLEAN} -name {name('store-imaging-weights')} -no-reorder -store-imaging-weights {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    @pytest.mark.parametrize(
        "gridder, test_name", (["", "shift-ws"], ["-use-wgridder", "shift-wg"])
    )
    def test_shift_image(self, gridder, test_name):
        # Shift the image with w-stacking and w-gridder gridder
        s = f"{tcf.WSCLEAN} {gridder} -name {name(test_name)} -mgain 0.8 -auto-threshold 5 -niter 1000000 -make-psf {tcf.DIMS_RECTANGULAR} -shift 08h09m20s -39d06m54s -no-update-model-required {tcf.MWA_MS}"
        validate_call(s.split())

    def test_shifted_source_list(self):
        # Shift the image and check coordinates in source list
        s = f"{tcf.WSCLEAN} -name {name('shifted-source-list')} -niter 1 {tcf.DIMS_RECTANGULAR} -shift 08h09m20s -39d06m54s -save-source-list {tcf.MWA_MS}"
        validate_call(s.split())
        source_file = f"{name('shifted-source-list')}-sources.txt"
        assert os.path.isfile(source_file)
        with open(source_file) as f:
            lines = f.readlines()
            # There should be a header line and a single source line in the file
            assert len(lines) == 2
            # 3rd and 4th column contain ra and dec
            cols = lines[1].split(",")
            assert len(cols) >= 4
            ra_str = cols[2]
            dec_str = cols[3]
            assert ra_str[0:5] + " " + dec_str[0:6] == "07:49 -44.12"

    def test_missing_channels_in_deconvolution(self):
        # The test set has some missing MWA subbands. One MWA subband is 1/24 of the data (32/768 channels), so
        # by imaging with -channels-out 24, it is tested what happens when an output channel has no data.
        s = f"{tcf.WSCLEAN} -name {name('missing-channels-in-deconvolution')} -use-wgridder {tcf.DIMS_LARGE} -baseline-averaging 2.0 -no-update-model-required -niter 150000 -auto-threshold 2.0 -auto-mask 5.0 -mgain 0.9 -channels-out 24 -join-channels -fit-spectral-pol 4 {tcf.MWA_MS}"
        validate_call(s.split())

    def test_grid_with_beam(self):
        """Requires that WSClean is compiled with IDG and EveryBeam"""
        name = "idg-beam"

        # Remove existing component files if present
        for source_file in ["sources", "sources-pb"]:
            component_file = name + "-" + source_file + ".txt"
            if os.path.exists(component_file):
                os.remove(component_file)

        s = f"{tcf.WSCLEAN} -name {name} -use-idg -grid-with-beam -save-source-list -mgain 0.8 -auto-threshold 5 -niter 1000000 -interval 10 14 {tcf.DIMS_LARGE} -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())
        for image_type in [
            "psf",
            "beam",
            "dirty",
            "image",
            "image-pb",
            "model",
            "model-pb",
            "residual",
            "residual-pb",
        ]:
            image_name = name + "-" + image_type + ".fits"
            assert os.path.isfile(image_name)
        # Check whether source files are generated
        for source_file in ["sources", "sources-pb"]:
            assert os.path.isfile(name + "-" + source_file + ".txt")

    def test_two_facets(self):
        # Apply the facet to the image
        s = f"{tcf.WSCLEAN} -name {name('two-facets')} -facet-regions {tcf.FACETFILE_2FACETS} \
            {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    def test_nfacets_pol_xx_yy(self):
        # Request two polarizations on approximately 25 facets
        s = f"{tcf.WSCLEAN} -name {name('nfacets-XX_YY')} -pol XX,YY \
            -facet-regions {tcf.FACETFILE_NFACETS} {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
        validate_call(s.split())

    @pytest.mark.parametrize("npol", (2, 4))
    def test_facet_h5solution(self, npol):
        # Test facet-based imaging and applying h5 solutions
        # where the polarization axis in the h5 file has size npol
        h5download = (
            f"wget -N -q {tcf.WSCLEAN_DATA_URL}/mock_soltab_{npol}pol.h5"
        )
        validate_call(h5download.split())

        name = f"facet-h5-{npol}pol"
        s = f"{tcf.WSCLEAN} -use-wgridder -name {name} -apply-facet-solutions mock_soltab_{npol}pol.h5 ampl000,phase000 -pol xx,yy -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_LARGE} -join-polarizations -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 {tcf.MWA_MS}"
        validate_call(s.split())

        # Check for output images
        assert os.path.isfile(f"{name}-psf.fits")
        for pol in ["XX", "YY"]:
            trunk = name + "-" + pol
            for image_type in [
                "image",
                "image-pb",
                "dirty",
                "model",
                "model-pb",
                "residual",
                "residual-pb",
            ]:
                image_name = trunk + "-" + image_type + ".fits"
                assert os.path.isfile(image_name)

    def test_facet_beam(self):
        # Test facet beam, using 4 polarizations
        s = f"{tcf.WSCLEAN} -name {name('nfacets-iquv-facet-beam')} -interval 10 14 -apply-facet-beam -pol iquv \
            -facet-regions {tcf.FACETFILE_NFACETS} {tcf.DIMS_RECTANGULAR} \
                -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())

    def test_mpi_join_channels(self):
        # Test wsclean-mp command
        s = f"mpirun {tcf.WSCLEAN_MP} -name {name('mpi-join')} {tcf.DIMS_RECTANGULAR} -channels-out 2 -join-channels -niter 1000000 -mgain 0.8 -auto-threshold 5 -multiscale -no-update-model-required {tcf.MWA_MS}"
        validate_call(s.split())

    def test_mpi_split_channels(self):
        s = f"mpirun {tcf.WSCLEAN_MP} -name {name('mpi-split')} {tcf.DIMS_RECTANGULAR} -channels-out 2 -niter 1000000 -mgain 0.8 -auto-threshold 5 -multiscale -no-update-model-required {tcf.MWA_MS}"
        validate_call(s.split())

    def test_idg_with_reuse_psf(self):
        # Test for issue #81: -reuse-psf gives segmentation fault in IDG
        # First make sure input files exist:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-psf-A')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())
        # Model image A is copied to B-model-pb corrected image, to avoid
        # issues due to NaN values in the A-model-pb.fits file.
        # As such, this test is purely illlustrative.
        os.rename(
            name("idg-reuse-psf-A") + "-model.fits",
            name("idg-reuse-psf-B") + "-model-pb.fits",
        )
        os.rename(
            name("idg-reuse-psf-A") + "-beam.fits",
            name("idg-reuse-psf-B") + "-beam.fits",
        )
        # Now continue:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-psf-B')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -continue -reuse-psf {name('idg-reuse-psf-A')} -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())

    @pytest.mark.skip(
        reason="-reuse-dirty and -grid-with-beam options conflict due to average beam computation (AST-995)"
    )
    def test_idg_with_reuse_dirty(self):
        # Test for issue #80: -reuse-dirty option fails (#80)
        # First make sure input files exist:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-dirty-A')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())
        # Model image A is copied to B-model-pb corrected image, to avoid
        # issues due to NaN values in the A-model-pb.fits file.
        # As such, this test is purely illlustrative.
        os.rename(
            name("idg-reuse-dirty-A") + "-model.fits",
            name("idg-reuse-dirty-B") + "-model-pb.fits",
        )
        # Now continue:
        s = f"{tcf.WSCLEAN} -name {name('idg-reuse-dirty-B')} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu -grid-with-beam -interval 10 14 -mgain 0.8 -niter 1 -continue -reuse-dirty {name('idg-reuse-dirty-A')} -mwa-path . {tcf.MWA_MS}"
        validate_call(s.split())

    def test_masked_parallel_deconvolution(self):
        # Test for two issues:
        # - issue #96: Source edges in restored image after parallel deconvolution
        # - issue #31: Model images are masked in parallel cleaning
        # The result of this test should be a model image with an unmasked Gaussian and a
        # properly residual. Before the fix, the Gaussian was masked in the model, and
        # therefore only a single pixel was visible, and residual would only be updated
        # on the place of the pixel.

        # First create a mask image with one pixel set:
        s = f"{tcf.WSCLEAN} -name {name('masked-parallel-deconvolution-prepare')} -size 256 256 -scale 1amin -interval 10 14 -niter 1 {tcf.MWA_MS}"
        validate_call(s.split())
        # Now use this as a mask, and force a Gaussian on the position
        s = f"{tcf.WSCLEAN} -name {name('masked-parallel-deconvolution')} -size 256 256 -scale 1amin -fits-mask {name('masked-parallel-deconvolution-prepare')}-model.fits -interval 10 14 -niter 10 -parallel-deconvolution 128 -multiscale -multiscale-scales 10 {tcf.MWA_MS}"
        validate_call(s.split())
        for f in glob.glob(
            f"{name('masked-parallel-deconvolution-prepare')}*.fits"
        ):
            os.remove(f)

    @pytest.mark.parametrize("use_beam", (False, True))
    def test_idg_predict(self, use_beam):
        # Check whether primary beam corrected image is used in -predict
        # First make sure model images exist
        run_name = name("idg-predict")
        grid_with_beam = "-grid-with-beam" if use_beam else ""
        s0 = f"{tcf.WSCLEAN} -name {run_name} {tcf.DIMS_LARGE} -use-idg -idg-mode cpu {grid_with_beam} -interval 10 12 -mgain 0.8 -niter 1 -mwa-path . {tcf.MWA_MS}"
        validate_call(s0.split())

        # Remove the model image that shouldn't be needed for the predict
        if use_beam:
            # Move model.fits to model-pb.fits file. Formally, the model-pb.fits file
            # should be used directly, but as this file can contain NaN values, a predict run
            # can bail out on these NaN values.
            os.rename(f"{run_name}-model.fits", f"{run_name}-model-pb.fits")

        s1 = f"{tcf.WSCLEAN} -name {run_name} {tcf.DIMS_LARGE} -predict -use-idg -idg-mode cpu {grid_with_beam} -interval 10 12 -mwa-path . {tcf.MWA_MS}"
        validate_call(s1.split())

    def test_catch_invalid_channel_selection(self):
        # Invalid selection: people often forget the second value of -channel-range is an open interval end (i.e. excluded the value itself).
        s = f"{tcf.WSCLEAN} -name {name('test-caught-bad-selection')} -channels-out 256 -channel-range 0 255 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        with pytest.raises(Exception):
            validate_call(s.split())

    def test_catch_invalid_channel_selection_with_gaps(self):
        s = f"{tcf.WSCLEAN} -name {name('test-caught-bad-selection')} -gap-channel-division -channels-out 256 -channel-range 0 255 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        with pytest.raises(Exception):
            validate_call(s.split())

    def test_catch_invalid_channel_selection_with_division(self):
        s = f"{tcf.WSCLEAN} -name {name('test-caught-bad-selection')} -channel-division-frequencies 145e6 -channels-out 256 -channel-range 0 255 {tcf.DIMS_LARGE} {tcf.MWA_MS}"
        with pytest.raises(Exception):
            validate_call(s.split())

    def test_multiband_no_mf_weighting(self):
        # Tests issue #105: Segmentation fault (core dumped), when grouping spectral windows + no-mf-weighting Master Branch
        # The issue was caused by invalid indexing into the BandData object.
        s = f"{tcf.WSCLEAN} -name {name('vla-multiband-no-mf')} -size 768 768 -scale 0.05arcsec -pol QU -mgain 0.85 -niter 1000 -auto-threshold 3 -join-polarizations -squared-channel-joining -no-update-model-required -no-mf-weighting {tcf.JVLA_MS}"
        validate_call(s.split())
        for f in glob.glob(f"{name('vla-multiband-no-mf')}*.fits"):
            os.remove(f)

    def test_spectrally_fitted_with_joined_polarizations(self):
        s = f"{tcf.WSCLEAN} -name {name('iv-jointly-fitted')} {tcf.DIMS_LARGE} -parallel-gridding 4 -channels-out 4 -join-channels -fit-spectral-pol 2 -pol i,v -join-polarizations -niter 1000 -auto-threshold 5 -multiscale -mgain 0.8 {tcf.MWA_MS}"
        validate_call(s.split())

    def test_direction_dependent_psfs(self):
        """Tests direction-dependent PSFs.
        Checks that the PSF generated which lies close to the source point is more similar to the dirty image than the one lying further away.
        """

        def get_peak_centered_normalized_subimage(img, subimage_size):
            """Get a subimage centered at the pixel with the heighest value
            Image is normalized by the heighest pixel value
            """
            # Get coordinates of the peak
            center_point_x, center_point_y = np.unravel_index(
                np.argmax(img, axis=None), img.shape
            )

            return (
                img[
                    center_point_x
                    - subimage_size // 2 : center_point_x
                    - subimage_size // 2
                    + subimage_size,
                    center_point_y
                    - subimage_size // 2 : center_point_y
                    - subimage_size // 2
                    + subimage_size,
                ]
                / img[center_point_x, center_point_y]
            )

        # Make template model image
        s = f"{tcf.WSCLEAN} -name {name('DD-PSFs')} -no-reorder -size 4800 4800 -scale 5asec -weight briggs -1 -padding 1.2 -gridder idg -grid-with-beam -beam-mode array_factor -aterm-kernel-size 15 -beam-aterm-update 120 {tcf.SKA_MS}"
        validate_call(s.split())

        # Fill model images with grid of point sources
        f_image = fits.open(name("DD-PSFs-image.fits"))
        f_beam = fits.open(name("DD-PSFs-beam.fits"))
        image_size = f_image[0].data.shape[-1]
        PSF_GRID_SIZE_1D = 3
        point_source_spacing = image_size // PSF_GRID_SIZE_1D
        position_range_1d = (
            point_source_spacing // 2
            + point_source_spacing * np.arange(PSF_GRID_SIZE_1D)
        )
        f_image[0].data[:] = 0.0
        for i in position_range_1d:
            for j in position_range_1d:
                f_image[0].data[0, 0, i, j] = 1.0 / f_beam[0].data[0, 0, i, j]
        f_image.writeto(name("DD-PSFs-model-pb.fits"), overwrite=True)

        # Predict visibilites using wsclean
        # predicted visibiliies are written to the MODEL_DATA column
        s = f"{tcf.WSCLEAN} -name {name('DD-PSFs')} -no-reorder -predict -padding 1.2 -gridder idg -grid-with-beam  -beam-mode array_factor  -beam-aterm-update 120 {tcf.SKA_MS}"
        validate_call(s.split())

        # Location of a python implementation of a "deconvolution algorithm" that does
        # nothing except storing its input images to disk
        deconvolution_script = os.path.join(
            os.path.dirname(__file__), "test_deconvolution_write_input.py"
        )

        # Generate the regular (direction independent) psf
        s = f"{tcf.WSCLEAN} -name {name('NO-DD-PSFs')} -make-psf-only -data-column MODEL_DATA -no-reorder -size 4800 4800 -scale 5asec -weight briggs -1 -padding 1.2 -gridder idg -grid-with-beam -beam-mode array_factor -aterm-kernel-size 15 -beam-aterm-update 120 {tcf.SKA_MS}"
        validate_call(s.split())

        # Generate dirty image and PSF_GRID_SIZE_1D x PSF_GRID_SIZE_1D direction-dependent PSFs
        s = f"{tcf.WSCLEAN} -name {name('DD-PSFs')} -data-column MODEL_DATA -parallel-deconvolution 1600 -no-reorder -size 4800 4800 -scale 5asec -mgain 0.8 -niter 10000000 -threshold 10.0mJy -auto-mask 5.0 -weight briggs -1 -padding 1.2 -gridder idg -grid-with-beam -beam-mode array_factor -aterm-kernel-size 15 -beam-aterm-update 120 -dd-psf-grid 3 3 -nmiter 1 -python-deconvolution {deconvolution_script} {tcf.SKA_MS}"
        validate_call(s.split())

        # Check whether the restoring beam for dd-psf and regular imaging is the same
        reference_header = fits.getheader(name("NO-DD-PSFs" + "-psf.fits"))
        ddpsf_header = fits.getheader(name("DD-PSFs" + "-image.fits"))
        assert np.isclose(
            reference_header["BMAJ"], ddpsf_header["BMAJ"], rtol=1e-3, atol=0.0
        )
        assert np.isclose(
            reference_header["BMIN"], ddpsf_header["BMIN"], rtol=1e-3, atol=0.0
        )
        assert np.isclose(
            reference_header["BPA"], ddpsf_header["BPA"], rtol=0.0, atol=0.1
        )

        SUBIMAGE_SIZE = 40
        dirty = get_peak_centered_normalized_subimage(
            fits.open(f"{name('DD-PSFs')}-dirty.fits")[0].data.squeeze()[
                :1600, :1600
            ],
            SUBIMAGE_SIZE,
        )
        psf_on_source = get_peak_centered_normalized_subimage(
            fits.open(f"{name('DD-PSFs')}-d0000-psf.fits")[0].data.squeeze(),
            SUBIMAGE_SIZE,
        )
        psf_off_source = get_peak_centered_normalized_subimage(
            fits.open(f"{name('DD-PSFs')}-d0004-psf.fits")[0].data.squeeze(),
            SUBIMAGE_SIZE,
        )

        # Verify that the psf generated at the location of a point source
        # is indeed a better match then a psf for a diffetent location
        expected_improvement_factor = 0.3
        assert np.sqrt(
            np.mean(np.square(dirty - psf_on_source))
        ) < expected_improvement_factor * np.sqrt(
            np.mean(np.square(dirty - psf_off_source))
        )

        num_psfs = PSF_GRID_SIZE_1D * PSF_GRID_SIZE_1D

        dirty = []
        psf = []

        # Load the psfs and dirty images that were stored by the dummy deconvolution algorithm
        for i in range(num_psfs):
            dirty.append(
                get_peak_centered_normalized_subimage(
                    np.load(
                        f"{name(f'test-deconvolution-write-input-dirty-{i}.npy')}"
                    )[0][0],
                    SUBIMAGE_SIZE,
                )
            )
            psf.append(
                get_peak_centered_normalized_subimage(
                    np.load(
                        f"{name(f'test-deconvolution-write-input-psf-{i}.npy')}"
                    )[0],
                    SUBIMAGE_SIZE,
                )
            )

        # Create a difference matrix of rms differences between all pairs of
        # dirty images and psfs.
        # The best match should occur when psf and dirty image index match,
        # i.e. on the diagonal of the difference matrix
        diff = np.zeros((num_psfs, num_psfs))
        for i in range(num_psfs):
            for j in range(num_psfs):
                diff[i, j] = np.sqrt(np.mean(np.square(dirty[i] - psf[j])))

        # Assert that for each dirty image, the best match is indeed the psf
        # with the same index
        assert all(np.argmin(diff, axis=0) == np.arange(num_psfs))

    def test_read_only_ms(self):
        chmod = f"chmod a-w -R {tcf.MWA_MS}"
        validate_call(chmod.split())
        try:
            # When "-no-update-model-required" is specified, processing a read-only measurement set should be possible.
            s = f"{tcf.WSCLEAN} -interval 10 20 -no-update-model-required -name {name('readonly-ms')} -auto-threshold 0.5 -auto-mask 3 \
                -mgain 0.95 -nmiter 2 -multiscale -niter 100000 {tcf.DIMS_RECTANGULAR} {tcf.MWA_MS}"
            validate_call(s.split())
        finally:
            chmod = f"chmod u+w -R {tcf.MWA_MS}"
            validate_call(chmod.split())

    def test_rr_polarization(self):
        s = f"{tcf.WSCLEAN} -pol rr -name {name('gmrt-rr')} -mgain 0.8 -niter 1 -size 512 512 -scale 10asec -gridder wstacking {tcf.GMRT_MS}"
        validate_call(s.split())
        rms_dirty = compute_rms(f"{name('gmrt-rr')}-dirty.fits")
        rms_image = compute_rms(f"{name('gmrt-rr')}-image.fits")
        # This was 0.215 when measured
        assert rms_dirty > 0.2 and rms_dirty < 0.22
        assert rms_dirty > rms_image

    def test_gmrt_beam(self):
        s = f"{tcf.WSCLEAN} -pol rr -name {name('gmrt-beam')} -apply-primary-beam -mgain 0.8 -size 512 512 -scale 10asec -gridder wstacking {tcf.GMRT_MS}"
        validate_call(s.split())
        rms_beam = compute_rms(f"{name('gmrt-beam')}-beam-0.fits")
        # This was measured at 0.6306
        assert rms_beam > 0.61 and rms_beam < 0.65
        rms_corrected = compute_rms(f"{name('gmrt-beam')}-image-pb.fits")
        # Measured at 0.45849
        assert rms_corrected > 0.42 and rms_corrected < 0.49

    def test_mf_full_polarization_beam_correction(self):
        prefix = name("mf-full-pol-beam")
        s = f"{tcf.WSCLEAN} -name {prefix} {tcf.DIMS_LARGE} -interval 10 12 -mwa-path . -channels-out 2 -apply-primary-beam -pol iquv -link-polarizations i -mgain 0.8 -niter 1000 -auto-threshold 6 -size 512 512 -scale 2amin {tcf.MWA_MS}"
        validate_call(s.split())

        assert os.path.isfile(prefix + "-0000-psf.fits")
        assert os.path.isfile(prefix + "-0001-psf.fits")
        assert os.path.isfile(prefix + "-MFS-psf.fits")
        for image_type in [
            "dirty",
            "image",
            "image-pb",
            "model",
            "model-pb",
            "residual",
            "residual-pb",
        ]:
            for pol_type in ["I", "Q", "U", "V"]:
                postfix = pol_type + "-" + image_type + ".fits"
                image_name = prefix + "-0000-" + postfix
                assert os.path.isfile(image_name)
                image_name = prefix + "-0001-" + postfix
                assert os.path.isfile(image_name)
                image_name = prefix + "-MFS-" + postfix
                assert os.path.isfile(image_name)
