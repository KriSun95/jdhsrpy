"""Essentially a wrapper for sunkit-spex things."""

from sunkit_spex.legacy.fitting.fitter import Fitter

__all__ = ["get_fitter_object", 
           "multi_or_single_fit_model", 
           "set_single_thermal_model", 
           "set_double_thermal_model", 
           "set_triple_thermal_model",
           "set_single_thermal_single_nonthermal_model",
           "set_double_thermal_single_nonthermal_model",
           ]

def get_fitter_object(*args, arf_files=None, rmf_files=None):
    """ Returns the ``sunkit-spex.legacy.fitting.fitter.Fitter`` object.
    
    The fitter object has been initialised with the given strings for 
    NuSTAR PHA files.

    If the ARF and RMF files are in the same directory as the PHA with 
    the same naming structure then they are not necessary to be given, 
    else give a corresponding ARF and/or RMF file for each PHA given.
    """

    return Fitter(pha_file=[*args, ], arf_file=arf_files, rmf_file=rmf_files)

def multi_or_single_fit_model(base_model, fitter):
    """Adds a constant to the model if multiple data-sets exist.
    
    So return `C*base_model` if `len(fitter.data)>1` else `base_model`.
    """
    if len(fitter.data)>1:
        return f"C*({base_model})"
    return f"{base_model}"

def set_single_thermal_model(fitter):
    """Set `fitter.model` with the `f_vth`model."""
    fitter.model = multi_or_single_fit_model("f_vth", fitter)

def set_double_thermal_model(fitter):
    """Set `fitter.model` with the `f_vth+f_vth`model."""
    fitter.model = multi_or_single_fit_model("f_vth+f_vth", fitter)

def set_triple_thermal_model(fitter):
    """Set `fitter.model` with the `f_vth+f_vth+f_vth`model."""
    fitter.model = multi_or_single_fit_model("f_vth+f_vth+f_vth", fitter)

def set_single_thermal_single_nonthermal_model(fitter):
    """Set `fitter.model` with the `f_vth+thick_fn`model."""
    fitter.model = multi_or_single_fit_model("f_vth+thick_fn", fitter)

def set_double_thermal_single_nonthermal_model(fitter):
    """Set `fitter.model` with the `f_vth+f_vth+thick_fn`model."""
    fitter.model = multi_or_single_fit_model("f_vth+f_vth+thick_fn", fitter)