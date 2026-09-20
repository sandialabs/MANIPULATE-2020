Note:

To do a clean install on a Ubuntu Linux system, go into the source subdirectory and type "./build_manipulate.sh --Clean-install"

To validate install, run regression tests by going into the input subdirectory and type "./manipulate-2020_regression.sh"

Note: If using other self-shield funcions, these must have been placed in the self-shield subdirectory.
      Tf this is not done, you will get a "STOP Open" error.

To run a calculation, go into the input directory and type: "perl manipulate-2020.prl xxx", where "xxx" is the filename.

Needed/Future changes:
1) add a diff step to compare new installed results to testcase outputs
2) add more regression tests for new MANIPULATE-2020 options
3) add regression tests for legacy ASTM standard support calculations
