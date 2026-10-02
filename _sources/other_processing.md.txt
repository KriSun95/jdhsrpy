# Additional Data Processing

A lot of case will be handled well in the [Downloading and Processing NuSTAR Data](https://krisun95.github.io/jdhsrpy/setting_up_nustar_data.html) section but certain scenarios might arise that need extra steps and care.

## Pileup

Many solar observations are subject to pulse pileup and require an additional processing step before their data should be used for science.
**Not applying this correction will result in erroneous removal of piled-up counts for any image-based product, *including spatially-filtered spectra*.**

Applying the correction is simple.
You need the following two `pixpos` files, obtainable from [here]():

```bash
nuApixpos20100101v007_pileup.fits
nuBpixpos20100101v007_pileup.fits
```

These files tell `nupipeline` how to properly assign position coordinates to events with grades 21-24.
They must be passed to `nupipeline` as keywords, e.g.

```bash
nupipeline \
<<other parameters>> \
fpma_pixposfile=/path/to/nuApixpos20100101v007_pileup.fits \
fpmb_pixposfile=/path/to/nuBpixpos20100101v007_pileup.fits \
<<other parameters>>
```

(see ["Get useable science files" section from here](setting_up_nustar_data.md#get-useable-science-files) for an example of a full `nupipeline` call)

## Fe 6.7 keV complex gain shift
