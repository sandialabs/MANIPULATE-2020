@echo off
rem manipulate-2020_regression windows shell
rem
rem Note: Run regression tests for manipulate-2020
rem       Later, add a comparison of output files
rem
set argC=0
for %%x in (%*) do Set /A argC+=1
if %argC% gtr 0 (
  echo "ERROR: manipulate-2020_regression called with Invalid number of arguments: %argC%"
  echo "Usage: manipulate-2020_regression"
  goto END
) else (
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
)
:END