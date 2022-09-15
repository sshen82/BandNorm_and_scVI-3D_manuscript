import sys,os,subprocess

inputfolder=sys.argv[1]
mapbed=sys.argv[2]
resolution=int(sys.argv[3])
outputdir = sys.argv[4]

# inputfolder="C:/Users/solei/Downloads/cstopictesting/cistopic_1mb"
# mapbed="C:/Users/solei/Downloads/cstopictesting/bed_file.bed"
# resolution=1000000
# outputdir="C:/Users/solei/Downloads/cstopictesting/cistopic_step2"


mapbedfh = open(mapbed)
bintocoord = {}
for line in mapbedfh :
	tokens = line.split()
	bintocoord[tokens[2]] = (tokens[0],str(int(tokens[1])+resolution/2))

inputfolder = inputfolder + "/"
all_files = [str(i) for i in range(1,401)]
filenames = [inputfolder + s for s in all_files]
filenames = [i + ".txt" for i in filenames]

outputdir = outputdir + "/"
all_files = [str(i) for i in range(1,401)]
filenames_out = [outputdir + s for s in all_files]
filenames_out = [i + ".txt" for i in filenames_out]

j=0
for matrixfilename in filenames :
	if matrixfilename != '' :
		matrixfilefh = open(matrixfilename)
		outfilefh = filenames_out[j]
		for line in matrixfilefh :
			tokens = line.split()
			firstcoor = bintocoord[tokens[0]]; secondcoor = bintocoord[tokens[1]] 
			outtokens = (firstcoor[0], firstcoor[1], secondcoor[0], secondcoor[1], tokens[2], tokens[3] )
			outline = '\t'.join(outtokens)
			print >> outfilefh, outline
		outfilefh.close()
	j += 1


