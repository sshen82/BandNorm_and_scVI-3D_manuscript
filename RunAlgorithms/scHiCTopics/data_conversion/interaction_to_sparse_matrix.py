from utils import *
import os,sys
import numpy as np
import pandas as pd

def main(argv):
	fileDir = sys.argv[1]
	outDir = sys.argv[2]
	libname = sys.argv[3]
	resolution = int(sys.argv[4])
	dist = int(sys.argv[5]) # number of bins
	
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	
	assembly = 'hg19'
	chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']
	chr_list = ["chr" + i for i in chr_list]
	valid_LPs = []

	for this_chr in chr_list:
		_, midPoints = generateChrMidpoints(chrs = [this_chr], assembly = assembly, resolution = resolution)
		offset = 0
		for mid_point in midPoints:
			for k in range(dist):
				if (offset+k) < len(midPoints):
					valid_LPs.append([this_chr, mid_point, this_chr, midPoints[k+offset]])
			offset += 1

	LP_list_saving = [LP[0] + ':' + str(LP[1]) + '-' + str(LP[3]) for LP in valid_LPs]
	set_lp_list = set(LP_list_saving)
	
	out_LP_name = os.path.join(outDir, libname + "_" + str(dist) + "_" + str(resolution) + "_LPnames.txt")
	np.savetxt(out_LP_name, LP_list_saving, fmt='%s', delimiter="\t", newline="\n")		

	cell_folders = [fileDir + "/" + i for i in os.listdir(fileDir)]
	filenames = list()
	names = list()
	for folder in cell_folders:
		new = [folder + "/" + i for i in os.listdir(folder)]
		filenames = filenames + new
		names = names + os.listdir(folder)
	cell_number = 0
	out = np.empty([0, 3])
	for filename in filenames:
		mat_file = open(filename)
		cell_save = []
		for line in mat_file:
			target = line.split()
			if target[0] == target[2]:
				search_lp = target[0] + ":" + str(int(target[1]) + 500000) + "-" + str(int(target[3]) + 500000)
			if search_lp in set_lp_list:
				cell_save.append([names[cell_number], LP_list_saving.index(search_lp), int(target[4])])

		mat_file.close()
		cell_save = np.asarray(cell_save)
		cell_save = cell_save[cell_save[:,1].argsort()]
		cell_number += 1
		out = np.concatenate((out, cell_save), axis=0)
	outs = pd.DataFrame(out)
	outs.to_csv(outDir + "/" + "cs_mat.csv", index=False)

if __name__ == "__main__":
        main(sys.argv)
