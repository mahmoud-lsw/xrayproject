#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  final.py
#  
#  Copyright 2014 febrie <febrimaru13@gmail.com>
#  Copyright 2016 sulis <sulis.astro08@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
from os import *
import sys

__author__ = "Sulistiyowati, Febrie A. Azizi"
__version__ = "0.1.0@2016.4.4"



#manual one by one
#=======================================
if len(sys.argv) == 1:
	print "Welcome to AutoFit"
	print "Please choose mode"
	print "1. Reduction process"
	print "2. Extraction process"
	print "3. Spectral Co-addition"
	#print "4. Create All ObsID light curve"

	mode = int(raw_input("Mode: "))

	if (mode == 1):
		system("python v2reduction.py")
		#system("echo 'Executing!'")
		#chdir(".")

	elif (mode == 2):
		system("python v2extraction.py")

	elif (mode == 3):
		system("python v2coaddspectra.py")

	else:
		print "Error[0] mode is not available, expected 1, 2, or 3"
else:
	toplot = ""
	while (len(sys.argv) > 1):
		option = sys.argv[1]; del sys.argv[1]
		if option == '-m':
			mode = sys.argv[1]; del sys.argv[1]
		elif option == '-f':
			filename = sys.argv[1]; del sys.argv[1]
		elif option == '-p':
			toplot = ' -p'

	if (mode == "1"):
		system("python v2reduction.py -f "+filename)
	elif (mode == "2"):
		system("python v2extraction.py -f "+filename)
	elif (mode == "3"):
		system("python v2coaddspectra.py -f "+filename+toplot)
	else:
		print "Error[0] mode is not available! 1, 2, or 3 is expected"