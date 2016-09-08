# xrayproject
Scripts to analyze Swift-XRT data


Four scripts are available in the package:
1. v2final.py
2. v2reduction.py
3. v2extraction.py
4. v2coaddspectra.py

1. v2final.py contains processing mode, consists of: reduction, extraction, co-addition (with certain condition, read more for detail)

2. v2reduction.py provides you facilities to do data reduction through XRTPIPELINE command in Heasoft package. You are allowed to do data reduction for multiple objects by listing them in the input file. You are required to make an input file using any text editor you are familiar with. The input file should be written as below:

Example of input file:
#ra			dec			rab			decb		dist	ulx	location
031820.00	-662910.97	031854.57	-663228.06	3700.	X1	/home/alien/Document/myOBSIDs/
051822.26	-663603.56	051854.57	-663228.06	5700.	X2	/home/alien/Document/yourOBSIDs/

Explanation for each column:
1. ra: object's RA position in hhmmss.ss
2. dec: object's Dec position in ddmmss.ss
3. rab: background's RA position in hhmmss.ss
4. decb: background's dec position in ddmmss.ss
5. dist: distance in kpc
6. ulx: ULX ID, it is practically up to you, but suggested to use ID from catalog
7. location: directory containing OBSID directories

If everything goes well, by the end of reduction process, there will be <ULX ID>_cl (containing <OBSID>_cl directories) and <ULX ID>_cls (containing <OBSID>_cls directories) directories added in your OBSID directory. <ULX ID>_cls directories and everything inside them are needed for extraction process.

3. v2extraction.py allows you to extract information from previously reduced data using XSELECT and XSPEC. During extraction, you actually perform:

- region filtering
- energy filtering
- HR calculation

For region filtering:
Source region is defined as a circle enclosed by 20 pixels radius centered at the given RA and dec position (in the input file).
Background region is defined as a circle enclosed by 50 pixels radius centered at the given RA and dec position (in the input file).

For energy filtering and HR calculation:
Low energy is defined to be 0.3 keV
Medium energy is defined to be 1.5 keV
High energy is defined to be 10 keV
HR is calculated by comparing detected count rate in high energy band (1.5 - 10 keV) and low energy band (0.3 - 1.5 keV). The error propagation is computed with simple error propagation from division between two count rates in each enegy band.

Grouping with grppha is also performed with minimum grouping 20. Appropriate rmf file is applied during grouping and arf file for each spectra is generated.

After extraction, inside <ULX ID>_cls, there will appear:
- A directory called addspec_<ULX ID>, which contains:
	- each ungrouped ("group" here means using GRPPHA in Heasoft) spectra with its correspondence background, arf, and rmf file.
	- list files of each group of spectrums:
		<x>grup<#>_src for source, <x>grup<#>_bkg for background, <x>grup<#>_arf for ancillary
		with: x -> s for soft, HR < 1
			  x -> h for hard, HR >= 1
			  # -> 1-6, depends on MJD of the observation (currently, Apr 4th 2016, 6 rmf groups are available for this code)
	  This list files will be used for spectral coaddition.
- A file called extract_product_<ULX ID>, which contains tabulation of: OBSID, MJD, EXPTIME, Count rate, Count rate error, HR, HR error, the correspondence rmf group of the spectra. This file can be utilized to build light curve and hardness ratio diagram.

4. v2coaddspectra.py is used for spectral co-addition, including source spectrums, background spectrums, and weighted arf files based on exposure time. Source and background spectrums are added by employing MATHPHA while arf files are combined with ADDARF. Since MATHPHA is used, it is possible to add as many 1000 source and background spectrums. Prior to co-addition, spectrums are grouped based on HR and rmf as processed during extraction.
After co-addition, inside <ULX ID>_cls, there will appear:
- A file called coadd_summary_<ULX ID>
While inside addspec_<ULX ID> directory, there will appear:
- A bunch of ASCII files also under the name of add* which are used as input file during spectral co-addition
- Coadded spectrum files under the name of add*.pha for each spectral group

(Possibly) (u)Useful information for files nomenclature (* shows anything, as in most of linux terminals):
- Files add*arfw.arf -> weighted arf file
- Files add*bkg.pha  -> co-added background spectrum
- Files add*src.pha  -> co-added ungrouped final spectrum
- Files add*grp.pha  -> co-added grouped final spectrum using GRPPHA, with minimum grouping 20. The grouping is carried out using the corresponding add*src.pha, add*bkg.pha, add*arfw.arf, and rmf file.
