import nustar_pysolar

__all__ = ["sunpos_evt"]

def sunpos_evt(file, load_path=None):
	"""Convert a `.evt` NuSTAR file to a `_sunpos.evt` file
	
	The new file will contain Solar-X/-Y coordinates for each event.
	"""
	load_path = "./" if load_path is None else load_path
	nustar_pysolar.convert.convert_file(file, load_path=load_path)
