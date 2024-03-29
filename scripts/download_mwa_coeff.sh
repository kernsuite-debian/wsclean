#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading the MWA coefficients. Coefficients
# file will only be downloaded if not yet present

set -e

[ "$1" = false ] && exit 0;

mkdir -p test_data
cd test_data/

MWA_COEFF_ARCHIVE=MWA_COEFF.tar.bz2
MWA_COEFF_H5=mwa_full_embedded_element_pattern.h5

if [ ! -f ${MWA_COEFF_H5} ] ; then

    wget -q https://support.astron.nl/software/ci_data/EveryBeam/mwa_full_embedded_element_pattern.tar.bz2 -O ${MWA_COEFF_ARCHIVE}

    tar -xjf $MWA_COEFF_ARCHIVE
    rm $MWA_COEFF_ARCHIVE
fi
