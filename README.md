## How to run
- environment.yml : file for setup and installation of correct versions of packages to run pipeline.
- control_df.csv : control file, input correct variables for calculation of trajectory; see examples in this csv.
- 2508_tspacepipeline_v3.R : Main pipeline, run this from jobscript.
- tSpace.R - tspace function rewritten slightly to work with the version of different packages from the yml file.
- tspace_helperfunctions.R - different helper functions from original package, edited to fit the package versioned in yml file.
