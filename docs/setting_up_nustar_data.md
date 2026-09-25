# Downloading and Processing NuSTAR Data

A repository that might help is another I own called [heasoft-test-base](https://github.com/KriSun95/heasoft-test-base). I will heavily pull from information there as it's purpose is to describe how to get NuSTAR data then process it for spectral fitting using terminal commands and Python.

## Downloading NuSTAR data

Visit the [NuSTAR Catalog webite](https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dnumaster&Action=More+Options). You can search for the observation you want with any field but the solar ones will begin with ``"Sol"`` in the ``"name"`` field.

The easiest might be to search using the ``"obsid"`` field and you can find a nice summary of all NuSTAR solar observations and corresponding IDs from [Iain Hannah's work](https://ianan.github.io/nsigh_all/).

Once you have found the data you want, hit ``Retrieve`` and download the file from the link.

When you get the data, it will be compressed so you'll have to do something similar to ``tar -xf`` then ``gunzip -r`` to that output to get workable files.

### Directly download with the OBSID

If you know the exact NuSTAR observation ID you want to download then you can use something like

```bash
wget -q -nH --no-check-certificate --cut-dirs=6 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks https://heasarc.gsfc.nasa.gov/FTP/nustar/data/obs/06/2//20619003001/
```

where I've used the observation ID ``20619003001`` in the URL. The ``06/2//`` section in the URL also comes from the ID too. The ``2//`` is the first number and the ``06/`` come from the second and third digit.

The a user only has to run ``gunzip -r 20619003001`` to get the files ready to work with.

## Processing NuSTAR Data

Processing NuSTAR data can have a lot of steps to it so let's get started.

### HEASoft install

In order to process NuSTAR data, you will require a working ``HEASoft`` install. I have a walkthrough in [heasoft-test-base/heasoft_install_instructions](https://github.com/KriSun95/heasoft-test-base/blob/main/heasoft_install_instructions/README.md) but you can also use the [official instructions](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/download.html) which look _a lot_ nicer than when I went through them a while ago (there are videos and Docker containers now).

Make sure to source the executable for your HEASoft install anytime you want to use it if necessary.

#### Troubles

I have had many errors installing and using HEASoft over the years. I'll try and leave some useful tips:

- The absolute paths to your files are too long.
  - I have this account for a variety of cryptic errors and the solution is moving whatever folder your data is in to a higher level on your machine.

### Get useable science files

Now we can use HEASoft to generate some event files in the ``event_cl`` folder. The following is what I have used:

```bash
nupipeline obsmode=SCIENCE_SC indir=./$OBSID steminputs=nu$OBSID outdir=event_cl entrystage=1 exitstage=2 pntra=OBJECT pntdec=OBJECT statusexpr=STATUS==b0000xx00xx0xx000 cleanflick=no hkevtexpr=NONE clobber=yes runsplitsc=yes splitmode=STRICT
```

where I'm assuming I'm running this command in the directory above the ``$OBSID`` folder, this being the same number as the ``"obsid"`` field mentioned in the downloading data instructions.

_We should now have a bunch of files we can do a lot with so you might want to pick and choose what bits you need from this from now on._

### Get grade 0 filtered files

To get files that have been filtered to only include grade 0 files, try:

```bash
nuscreen infile=nu"$OBSID"A06_cl.evt gtiscreen=no evtscreen=yes gtiexpr=NONE gradeexpr=0 statusexpr=NONE outdir=./ hkfile=./nu"$OBSID"A_fpm.hk outfile=nu"$OBSID"A06_cl_grade0.evt
```

where you would want to run this for FPMB too as the ``"A"`` in all the file names above represent this is using the FPMA files.

### Get the PHA, ARF, and RMF spectral files

The finale step is to get the spectral files we need for spectral fitting. For FPMA again, something like the following will do

```bash
nuproducts indir=./ instrument=FPMA steminputs=nu"$OBSID" outdir=./ extended=no runmkarf=yes runmkrmf=yes infile=nu"$OBSID"A06_cl_grade0.evt bkgextract=no srcregionfile="$REGION_FILEA" attfile=./nu"$OBSID"_att.fits hkfile=./nu"$OBSID"A_fpm.hk usrgtifile="$TIME_INTERVAL_FILE"
```

where ``REGION_FILEA`` and ``$TIME_INTERVAL_FILE`` represent a ``.reg`` file and a ``gti.fits`` file for region and time selection, respectively.

Once the above has run, you should now have PHA, ARF, and RMF files created from the ``nu"$OBSID"A06_cl_grade0.evt`` file you provided for the ``infile`` argument.

#### Getting a region file

You can obtain a region file using FITS viewing software like [SAOImageDS9](https://sites.google.com/cfa.harvard.edu/saoimageds9?pli=1&authuser=0).

#### Getting a good time interval file

You can use the function ``~jdhsrpy.screening.make_gti_file``.
