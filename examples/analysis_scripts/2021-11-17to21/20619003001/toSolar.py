import os

import nustar_pysolar

obs_id = "20619003001"

file_dir = "/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/event_cl/"

orig_files = [os.path.join(file_dir, f"nu{obs_id}A06_cl.evt"), 
              os.path.join(file_dir, f"nu{obs_id}B06_cl.evt")]

for c,f in enumerate(orig_files):
	nustar_pysolar.convert.convert_file(f, load_path="./")
	print(f"\rDone file {c+1} of {len(orig_files)}", end="")