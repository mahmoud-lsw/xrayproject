#!/usr/bin/env python -tt
import re,os,sys,math,scipy,glob, datetime
import numpy as np
from xspec import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

__author__ = "Sulistiyowati, Febrie A. Azizi"
__version__ = "0.1.0@2015.11.25"


def regfilt(name,ra='',dec='',sreg='20',breg='50'):
	"""Region maker."""

	sreg = sreg*2.36
	breg = breg*2.36
	bkgra  = str(rab)
	bkgdec = str(decb)
	name = name[:-7]
	with open(name+"_src.reg","w") as fout:
		fout.write("# Region file format: DS9 version 7.2\n")
		fout.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' ")
		fout.write("select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
		fout.write("fk5\n")
		fout.write("circle(%s:%s:%s,%s:%s:%s,%f\")" %(ra[0:2],ra[2:4],ra[4:],dec[0:3],dec[3:5],dec[5:],sreg))
	with open(name+"_bkg.reg","w") as fout:
		fout.write("# Region file format: DS9 version 7.2\n")
		fout.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" ")
		fout.write("select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
		fout.write("fk5\n")
		fout.write("circle(%s:%s:%s,%s:%s:%s,%f\") # background" %(bkgra[0:2],bkgra[2:4],bkgra[4:],bkgdec[0:3],\
		bkgdec[3:5],bkgdec[5:],breg))
		

def sextodeg(ra,dec):
	"""Coordinate conversion to degree."""

	ra_d = (float(ra[0:2])+float(ra[2:4])/60+float(ra[4:])/3600)*15
	dec_d = float(dec[1:3])+float(dec[3:5])/60+float(dec[4:])/3600
	dec = dec[0]+str(dec_d)
	return str(ra), dec

	
def Xselspectra(name):
	"""Spectra filtering and extraction."""

	with open("xselectinputspec","w") as fout:
		fout.write("xselect\n")
		fout.write("read events "+str(name)+"_cl.evt\n")
		fout.write("./\n")
		fout.write("yes\n")
		fout.write("extract image\n")
		name = name.replace(".evt","",1)
		fout.write("save image"+name+"_all.img\n")
		fout.write("filter region "+name+"_src.reg\n")
		fout.write("extract all\n")
		fout.write("save spectrum "+str(name)+"_src.pha\n")
		fout.write("save curve "+str(name)+"_src.lc\n")
		fout.write("clear region\n")
		fout.write("filter region "+name+"_bkg.reg\n")
		fout.write("extract all\n")
		fout.write("save spectrum "+str(name)+"_bkg.pha\n")
		fout.write("save curve "+str(name)+"_bkg.lc\n")
		fout.write("clear region\n")
		fout.write("exit\n")
		fout.write("no\n")
	spec = name+"_xsel.log"
	os.system("xselect < xselectinputspec > "+str(spec))
		

def prefit(count,olddir,name=''):
	"""Generating ancillary file, grouping spectra min 20.
	Appropriate rmf file from http://www.swift.ac.uk/analysis/xrt/rmfarf.php.
	Updated 2014 July 2."""

	grup = 0
	os.system("ftlist "+str(name)+"_src.pha k include='MJD-OBS' outfile=keywords.txt clobber=yes")
	os.system("sed -i '1!d' keywords.txt")
	with open('keywords.txt') as f:
		mjd_obs = str(f.readlines())
		mjd_obs = float(re.search("= (.+?) /",mjd_obs).group(1))
	if(mjd_obs < 54101.): #2007Jan01 at 00:00:00
		resp = olddir+"/rmf/swxpc0to12s0_20010101v012.rmf"
		grup = 1
	elif(mjd_obs < 54343.): #2007Aug31 at 00:00:00
		resp = olddir+"/rmf/swxpc0to12s0_20070101v012.rmf"
		grup = 2
	elif(mjd_obs < 54832.): #2009Jan01 at 00:00:00
		resp = olddir+"/rmf/swxpc0to12s6_20010101v014.rmf"
		grup = 3
	elif(mjd_obs < 55562.): #2011Jan01 at 00:00:00
		resp = olddir+"/rmf/swxpc0to12s6_20090101v014.rmf"
		grup = 4
	elif(mjd_obs < 56293.): #2013Jan01 at 00:00:00
		resp = olddir+"/rmf/swxpc0to12s6_20110101v014.rmf"
		grup = 5
	else:
		resp = olddir+"/rmf/swxpc0to12s6_20130101v014.rmf"
		grup = 6
	os.system("xrtmkarf outfile=./"+name+"_src.arf phafile=./"+name+"_src.pha srcx=0 srcy=0 expofile="+name+\
	"_ex.img rmffile="+resp+" psfflag=yes cleanup=yes clobber=yes history=yes")
	arf = name+"_src.arf"
	os.system("grppha "+name+"_src.pha "+name+"_srcgrp.pha comm='sys_err 0.01 & group min 20 & chkey backfile "\
	+name+"_bkg.pha & chkey ancrfile "+arf+" & chkey respfile "+resp+" & exit' clobber=yes chatter=0")
	return grup


def pha_stats(name):
	"""Extracting exposure time, MJD, and total count."""

	os.system("ftlist "+str(name)+"_src.pha k include='TOTCTS,OBJECT,MJD-OBS,EXPOSURE' outfile=keywords.txt clobber=yes")
	with open('keywords.txt','r') as f:
		lines = f.readlines()
		OBJECT = (re.search("'(.+?)'",lines[0]).group(1))
		EXP = float(re.search("= (.+?) /",lines[1]).group(1))
		MJD_OBS = (float(re.search("= (.+?) /",lines[2]).group(1)))
		TOTCTS = (float(re.search("= (.+?) /",lines[3]).group(1)))
	return OBJECT,EXP,MJD_OBS,TOTCTS


def division(valx,valy):
	"""Division, used to calculate HR."""

	try:
		value = valy/valx
		return value
	except:
		value = float("Inf")
		return value


def error_calc(valx,evalx,valy,evaly,valz):
	"""Error propagation, used to calculate error_HR."""

	try:
		evalz = math.sqrt((evalx/valx)**2+(evaly/valy)**2)*valz
		return evalz
	except:
		return 0


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


def prepare_add(name,formerdir,laterdir,j):
	"""Preparing files to be added with addspec.
	Copy and rename source and background spectra, arf file to the addspec directory."""

	os.system('cp '+formerdir+'/'+name+'_src.pha '+laterdir+'/'+str(j)+'.pha')
	os.system('cp '+formerdir+'/'+name+'_src.arf '+laterdir+'/'+str(j)+'.arf')
	os.system('cp '+formerdir+'/'+name+'_bkg.pha '+laterdir+'/'+str(j)+'b.pha')


def filestoadd(addloc,groupname):
	"""Classifying spectra based on rmf and HR.
	Creating input files for addspec and addarf command."""

	var = np.array(eval(groupname))
	if len(var) != 0:
		np.savetxt(addloc+'/'+groupname+'_src',var[:,0],fmt='%s')
		np.savetxt(addloc+'/'+groupname+'_arf',var[:,1],fmt='%s')
		np.savetxt(addloc+'/'+groupname+'_bkg',var[:,2],fmt='%s')

if __name__ == "__main__":

	# Reading input file
	if len(sys.argv) == 1:
		input_file = str(raw_input("Write your extraction input filename: "))
	else:
		option = sys.argv[1]; del sys.argv[1]
		if option == '-f':
			input_file = sys.argv[1]; del sys.argv[1]
		else:
			print 'Error[1] no option available! -f is expected'
	
	input_object, nobject = read_data(input_file)
	
	for i in range(0,nobject):
		ra = input_object[:,0][i] # RA source in hhmmss.ss
		dec = input_object[:,1][i] # Dec source in ddmmss.ss
		rab = input_object[:,2][i] # RA background in hhmmss.ss
		decb = input_object[:,3][i] # Dec background in hhmmss.ss
		dist = input_object[:,4][i] # distance of the source/host galaxy  in kpc
		ulx = input_object[:,5][i] # ULX name or ID
		comp_path = input_object[:,6][i] # ULX directory, containing cls directory produced from reduction process
	
		# All possible situation, number indicating which rmf file is used:
		# 1  swxpc0to12s0_20010101v012.rmf
		# 2  swxpc0to12s0_20070101v012.rmf
		# 3  swxpc0to12s6_20010101v014.rmf
		# 4  swxpc0to12s6_20090101v014.rmf
		# 5  swxpc0to12s6_20110101v014.rmf
		# 6  swxpc0to12s6_20130101v014.rmf
		# first character indicating: s for softstate (HR < 1)
		#                             h for hardstate (HR >= 1)

		#groups = ['sgrup1','sgrup2','sgrup3','sgrup4','sgrup5','sgrup6',\
		#		'hgrup1','hgrup2','hgrup3','hgrup4','hgrup5','hgrup6']

		#sgrup1 = []
		#sgrup2 = []
		#sgrup3 = []
		#sgrup4 = []
		#sgrup5 = []
		#sgrup6 = []
		#hgrup1 = []
		#hgrup2 = []
		#hgrup3 = []
		#hgrup4 = []
		#hgrup5 = []
		#hgrup6 = []
		
	    # Initial parameters		

		low_enrg = 0.3 # definition of low energy limit to calculate HR
		med_enrg = 1.5 # definition of medium energy limit to calculate HR
		hig_enrg = 10. # definition of high energy limit to calculate HR
		#theta   = 0.
		#nh	= 1.0
		#fre = 'n'
		
		sreg= 20 # radius (in pixel) of source
		breg= 50 # radius (in pixel) of background

		old_dir = os.getcwd() # current working directory (finale)
		#dir_new = 'addspec_'+ulx # name of directory for everything regarding spectral addition for each ULX
		#os.makedirs(comp_path+dir_new) # make directory for spectral addition for each ULX
		#os.system('cp '+old_dir+'/rmf/swxpc0to12s0_20010101v012.rmf '+comp_path+'/'+dir_new+'/1.rmf') # copy and change name of 
		#os.system('cp '+old_dir+'/rmf/swxpc0to12s0_20070101v012.rmf '+comp_path+'/'+dir_new+'/2.rmf') # corresponding rmf file
		#os.system('cp '+old_dir+'/rmf/swxpc0to12s6_20010101v014.rmf '+comp_path+'/'+dir_new+'/3.rmf')
		#os.system('cp '+old_dir+'/rmf/swxpc0to12s6_20090101v014.rmf '+comp_path+'/'+dir_new+'/4.rmf')
		#os.system('cp '+old_dir+'/rmf/swxpc0to12s6_20110101v014.rmf '+comp_path+'/'+dir_new+'/5.rmf')
		#os.system('cp '+old_dir+'/rmf/swxpc0to12s6_20130101v014.rmf '+comp_path+'/'+dir_new+'/6.rmf')
		#document = open(comp_path+'/'+dir_new+'/'+dir_new+'_doc','w') # creating documentation of which OBSID corresponds to which shortened filename
		os.chdir(comp_path) # entering ULX directory containing cls directories
		dirlist = glob.glob('*cls') # list all cls directories from all OBSIDs
		now = datetime.datetime.now()
		summary_extract = open('extract_product_'+ulx,'w') # create file to give brief summary on below parameters
		summary_extract.write('%s\n' %('#'+str(now.strftime("Y%-%m-%d %H:%M"))))
		summary_extract.write('%11s%15s%15s%15s%17s%15s%15s%12s\n' %('#OBSID','MJD','EXPTIME','Count rate','err_Count rate','HR','err_HR','Group rmf'))
		n = 0 # counter extracted OBSIDs
		for datadir in dirlist:
			ncleanfile = len(glob.glob(datadir+'/*cl*')) # check if clean, exposure image, and ancillary file available
			nexpimage = len(glob.glob(datadir+'/*ex.img*'))
			narffile = len(glob.glob(datadir+'/*arf*'))
			obsid = datadir[:11]
			
			if ncleanfile == 1 and nexpimage == 1 and narffile == 1:
				new_dir = comp_path+datadir 
				os.chdir(new_dir) # enter each cls directory
				cl_file = glob.glob('*cl*')[0] 
				ex_imag = glob.glob('*ex.img*')[0]
				ar_file = glob.glob('*arf*')[0]
				
				regfilt(cl_file,ra,dec,sreg,breg)
				basename = cl_file[:-7] # clean filename, cut at "_cl.evt"
				
				Xselspectra(basename) # extract and filter clean file
				ob_na,exp,mjd,tot = pha_stats(basename) 
				grup = prefit(tot,old_dir,basename)
				Xset.chatter = 10
				logFile = Xset.openLog(cl_file+"_"+ulx+"_xspec_"+".log") 
				logFile = Xset.log
				s       = Spectrum(basename+"_srcgrp.pha") # pha file to be extracted using PyXspec
				color   = [] # array to store low and high energy intesity (represented by count rate), including the error
				s.ignore("**-"+str(low_enrg)+","+str(med_enrg)+"-**") # low energy range (0.3-1.5 keV)
				color.append([s.rate[0],s.rate[1]]) # low energy count rate and its error
				s.notice("all")
				s.ignore("**-"+str(med_enrg)+","+str(hig_enrg)+"-**") # high energy range (1.5-10. keV)
				color.append([s.rate[0],s.rate[1]]) # high energy count rate and its error
				s.notice("all") 
				s.ignore("**-"+str(low_enrg)+","+str(hig_enrg)+"-**") 
				color.append([s.rate[0],s.rate[1]]) # all energy range (0.3-10. keV) and its error
				s.notice("all")
	#Open file to write:
				softctr		= color[0][0]
				err_softctr = color[0][1]
				hardctr 	= color[1][0]
				err_hardctr = color[1][1]
				allctr 		= color[2][0]
				err_allctr 	= color[2][1]
				HR   		= division(softctr,hardctr) # calculate HR by dividing hard over soft count rate
				err_HR		= error_calc(softctr,color[0][1],hardctr,color[1][1],HR) # calculate error_HR
				summary_extract.write('%11s%15.5E%15.5E%15.5E%17.5E%15.5E%15.5E%12i\n' %(obsid,mjd,exp,allctr,err_allctr,HR,err_HR,grup))
				n += 1
		#		document.write('%11s%7s\n' %(obsid,str(n)))
		#		previous = new_dir # present cls directory
		#		next = comp_path+dir_new # addspec directory inside present cls directory
		#		prepare_add(basename,previous,next,n) # copy source and background pha files, and arf file to addspec directory
				
				#source, arf, background, not including those with HR = float('INF') or negative value
		#		if 0. < HR < 1.:
		#			if grup == 1:
		#				sgrup1.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])
		#			elif grup == 2:
		#				sgrup2.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])
		#			elif grup == 3:
		#				sgrup3.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])
		#			elif grup == 4:
		#				sgrup4.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])
		#			elif grup == 5:
		#				sgrup5.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])
		#			else:
		#				sgrup6.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])		
		#		elif HR >= 1. and HR != float('INF'):
		#			if grup == 1:
		#				hgrup1.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])
		#			elif grup == 2:
		#				hgrup2.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])
		#			elif grup == 3:
		#				hgrup3.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])
		#			elif grup == 4:
		#				hgrup4.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])
		#			elif grup == 5:
		#				hgrup5.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])
		#			else:
		#				hgrup6.append([str(n)+'.pha',str(n)+'.arf',str(n)+'b.pha'])					
				os.chdir(comp_path)
				print ulx+': '+str(n)+' of '+str(len(dirlist))+' listed OBSIDs are extracted'
			else:
				continue
			
		summary_extract.close()
		#document.close()
	

	#Preparing to add spectra
	#	dir_new2 = comp_path+'/'+dir_new
	#	for shgrup in groups:
	#	 	filestoadd(dir_new2,shgrup)

		os.chdir(old_dir)
		print "Extraction completed!"

	# #os.chdir()
	# #Creating plots
	# 	print "Creating light curve and HR diagram"

	# 	host = input_file[:-4]
	# 	#list and read extract file
	# 	file_extract = []
	# 	for file in os.listdir("."):
	# 		if file.startswith("extract"):
	# 			file_extract.append(file)

	# 	majorLocator = MultipleLocator(1000)
	# 	majorFormatter = FormatStrFormatter("%d")
	# 	minorLocator = MultipleLocator(100)

	# 	#lightcurve -> fig, HR -> fig2
	# 	#=======================================================================================
	# 	fig1 = plt.figure(figsize=(16.,12.))
	# 	ax1 = plt.subplot()

	# 	fig2 = plt.figure(figsize=(16.,12.))
	# 	ax2 = plt.subplot()

	# 	for extract_ulx in file_extract:
	# 		alldata = np.genfromtxt(extract_ulx)
	# 		mjd = alldata[:,1]
	# 		ctr = alldata[:,3]
	# 		err_ctr = alldata[:,4]
	# 		hr = alldata[:,5]
	# 		err_hr = alldata[:,6]
	# 		print extract_ulx
	# 		ax1.errorbar(mjd, ctr, yerr=err_ctr, fmt='o', label=extract_ulx[-2:], lw=2)
	# 		ax2.errorbar(ctr, hr, xerr=err_ctr, yerr=err_hr, fmt='o', label=extract_ulx[-2:], lw=2)


	# 	ax1.legend(loc=0,ncol=len(file_extract), mode="expand", borderaxespad=0., fontsize = 24)
	# 	ax2.legend(loc=0,ncol=len(file_extract), mode="expand", borderaxespad=0., fontsize = 24)

	# 	ax1.set_xlabel("MJD", fontsize=22)
	# 	ax1.set_ylabel("Count rate (ctr/s)", fontsize=22)
	# 	ax1.set_title(host, fontsize=28)

	# 	ax2.set_xlabel("Count rate (ctr/s)", fontsize=22)
	# 	ax2.set_ylabel("Hardness ratio", fontsize=22)
	# 	ax2.set_title(host, fontsize=28)

	# 	ax1.xaxis.set_major_locator(majorLocator)
	# 	ax1.xaxis.set_major_formatter(majorFormatter)
	# 	ax1.tick_params(axis='x', labelsize=22, length = 12)
	# 	ax1.tick_params(axis='y', labelsize=22, length = 12)
	# 	ax1.set_xlim([53000.,58000.])

	# 	ax1.xaxis.set_minor_locator(minorLocator)

	# 	ax2.tick_params(axis='x', labelsize=22, length = 12)
	# 	ax2.tick_params(axis='y', labelsize=22, length = 12)

	# 	fig1.savefig("lc_"+host)
	# 	fig2.savefig("hr_"+host)
	# 	os.chdir(old_dir)
		#=======================================================================================
