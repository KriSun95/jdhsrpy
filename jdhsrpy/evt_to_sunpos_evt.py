import nustar_pysolar

def evt_to_sunpos_evt(file, load_path=None):
	load_path = "./" if load_path is None else load_path
	nustar_pysolar.convert.convert_file(file, load_path=load_path)
