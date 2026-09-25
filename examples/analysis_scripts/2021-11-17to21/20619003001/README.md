# 2021 November 17: 3rd Orbit

Data downloaded with the contents of [./browse_download_script.txt](./browse_download_script.txt) in the terminal.

Created the science ready files by running [./nu_pipe_run_all.sh](nu_pipe_run_all.sh) in the terminal.

Sun-position ("sunpos") files were created using [./toSolar.py](./toSolar.py).

To produce a time filtered EVT file and GTI file, [./inspect_evt.py](./inspect_evt.py) was used to plot the time profile and select times then create the files.

To obtain a region file, [SAOImageDS9](https://sites.google.com/cfa.harvard.edu/saoimageds9?pli=1&authuser=0) was used.

Using the contents of [./nup_script.txt](./nup_script.txt), the event files (_not_ the "sunpos" ones) were screened for grade 0 then time + spatially filterd.
