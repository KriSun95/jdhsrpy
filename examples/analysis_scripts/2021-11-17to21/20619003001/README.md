# 2021 November 17: 3rd Orbit

Data downloaded with the contents of [./01-browse_download_script.txt](./01-browse_download_script.txt) in the terminal.

Created the science ready files by running [./02-nu_pipe_run_all.sh](./02-nu_pipe_run_all.sh) in the terminal.

Sun-position ("sunpos") files were created using [./03-to_solar.py](./03-to_solar.py).

To produce a time filtered EVT file and GTI file, [./04-inspect_evt.py](./04-inspect_evt.py) was used to plot the time profile and select times then create the files.

- Times of interest are then manually stored in the [./times_of_interest.yaml](./times_of_interest.yaml) file.

To obtain a region file, [SAOImageDS9](https://sites.google.com/cfa.harvard.edu/saoimageds9?pli=1&authuser=0) was used.

Using the contents of [./05-nup_script.txt](./05-nup_script.txt), the event files (_not_ the "sunpos" ones) were screened for grade 0 then time + spatially filterd.

The file _./06\*-spectral_fitting\*.py_ was used to perform spectral fitting and the associated plots.
