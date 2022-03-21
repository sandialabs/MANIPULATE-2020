#!/bin/bash
#
# manipulate-2020_regression Unix shell
#
# Note: Run regression tests for manipulate-2020
#       Later, add a comparison of output files
#
# check number of input arguments
if [[ $# > 0 ]]; then
  echo "manipulate-2020_regression called with Invalid number of arguments: $#"
  echo "Usage: manipulate-2020_regression"
  (exit 1) && true
else
  perl manipulate-2020.prl example_NJOY_groupr_convert
  perl manipulate-2020.prl example_cross_section_covariance_verification
  perl manipulate-2020.prl example_NJOY_groupr_response_convert
  perl manipulate-2020.prl example_resp_unc_spectrum_averaged_response
  perl manipulate-2020.prl example_NJOY_groupr_spectrum_convert
  perl manipulate-2020.prl example_response_fold
  perl manipulate-2020.prl example_NJOY_groupr_xsec_convert
  perl manipulate-2020.prl example_spct_unc_spectrum_averaged_response
  perl manipulate-2020.prl example_NJOY_response_combination
  perl manipulate-2020.prl example_composite_uncertainty
fi
