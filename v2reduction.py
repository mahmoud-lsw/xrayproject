import os
import sys
import numpy as np

__author__ = "Sulistiyowati, Febrie A. Azizi"
__version__ = "0.1.0@2016.4.4"


def read_data(ifile):
	"""Read input file (for this particular code) manually."""

	listdata = []
	num = 0
	with open(ifile, 'r') as data:
		for line in data:
			line = line.strip()
			if len(line) != 0 and line[0] != '#':
				iline = line.split()
				listdata.append(iline)
				num += 1
	listdata = np.array(listdata)
	return listdata, num


if __name__== "__main__":

	old_dir = os.getcwd() #get and store current working directory

	#Reading input file
	if len(sys.argv) == 1:
		input_file = str(raw_input("Write your reduction input filename: "))
	else:
		option = sys.argv[1]; del sys.argv[1]
		if option == '-f':
			input_file = sys.argv[1]; del sys.argv[1]
		else:
			print 'Error[1] no option available! -f is expected'

	input_object, nobject = read_data(input_file)

	for i in range(0,nobject):
		ra = input_object[:,0][i] # RA source in hhmmss.ss
		ra = ra[0:2]+" "+ra[2:4]+" "+ra[4:]
		dec = input_object[:,1][i] # Dec source in ddmmss.ss
		dec = dec[0:3]+" "+dec[3:5]+" "+dec[5:]
		ulx = input_object[:,2][i] # ULX name or ID
		host_name = input_object[:,3][i] # host name as written in the directory
		obj_path = input_object[:,4][i] # directory listing all objects
		comp_path = obj_path+"/"+host_name+"/"


		clean_dir = obj_path+"/"+host_name+"_cl/"+ulx+"_cl"
		cleans_dir = obj_path+"/"+host_name+"_cls/"+ulx+"_cls"

		os.makedirs(clean_dir)
		os.makedirs(cleans_dir)

		dirlist = os.listdir(comp_path)

		for datadir in dirlist:
			if len(datadir) != 11:
				continue
			else:
				new_dir = comp_path+datadir
				os.chdir(new_dir)		
				odir=clean_dir+"/"+datadir+"_"+ulx+"_cl"
				adir=cleans_dir+"/"+datadir+"_"+ulx+"_cls"
				if not os.path.exists(odir):
					os.makedirs(odir)
				if not os.path.exists(adir):
					os.makedirs(adir)
				os.system("xrtpipeline indir=./ outdir="+odir+" steminputs=sw"+datadir+" srcra='"+ra+"' srcdec='"+dec+\
					"' datamode=PC gtiexpr=\"(ELV>=45||BR_EARTH>=120)&&(SUN_ANGLE>=45&&ANG_DIST<=0.08&&MOON_ANGLE>=14)\" cleanup=yes chatter=5 clobber=yes > "+odir+"/"+datadir+"_"+ulx+"_pipe.log")
				os.chdir(odir)
				os.system("cp *cl.evt "+adir)
				os.system("cp *ex.img "+adir)
				os.system("cp *arf "+adir)
				print datadir+" done"
				if len(os.listdir(adir)) != 0:
					continue
				else:
					os.rmdir(adir)
				os.chdir("..")

		dirlist = os.listdir(cleans_dir)

		print 'Processes completed!'
		print ulx+" has "+str(len(dirlist))+" OBSIDs ready for the next steps :)"	
		os.chdir(old_dir)
