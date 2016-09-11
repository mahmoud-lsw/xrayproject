#!/usr/bin/python -tt
import re,os,sys,math,scipy,glob, shutil, datetime
import numpy as np
from xspec import *
#import matplotlib
#backend = 'Agg'
#matplotlib.use(backend)
#import matplotlib.pyplot as plt
#matplotlib.get_backend()
import matplotlib.pyplot as plt
import v2extraction as ex
from scipy.optimize import curve_fit


def addspectra(arflist,srclist,bkglist):
    """Add spectra."""

    with open("addarfinput_"+arflist,"w") as fout:
        fout.write("@"+arflist+"\n")
        fout.write(arflist+".arf\n")
        fout.write("no\n")
        fout.write("no\n")

    with open("addspectrainput_"+srclist,"w") as fout:
        fout.write("@"+srclist+"\n")
        fout.write("C\n")
        fout.write(srclist+".pha\n")
        fout.write("CALC\n")
        fout.write("%\n")
        fout.write("\n")
        fout.write("\n")
        #fout.write("no\n")
        #fout.write("no\n")

    with open("addbkginput_"+bkglist,"w") as fout:
        fout.write("@"+bkglist+"\n")
        fout.write("C\n")
        fout.write(bkglist+".pha\n")
        fout.write("CALC\n")
        fout.write("%\n")
        fout.write("\n")
        fout.write("\n")
        #fout.write("no\n")
        #fout.write("no\n")

    os.system("addarf < addarfinput_"+arflist)
    os.system("mathpha < addspectrainput_"+srclist)
    os.system("mathpha < addbkginput_"+bkglist)


def plot_hidcut(time, countr, err_cr, HR, err_HR, HR_est, sigma_est, popt, host, ulx, show=True, savefig=True, filename='HID.png'):
    """Plot """
    indexmin = np.argmin(countr)
    indexmax = np.argmax(countr)
    xmin = countr[indexmin]-err_cr[indexmin]
    xmax = countr[indexmax]+err_cr[indexmax]
    width = xmax-xmin
    xmin = xmin - 0.05*width
    xmax = xmax + 0.05*width
    xfine = np.linspace(xmin,xmax,100)


    GOLDEN_RATIO = 0.5*(1. + np.sqrt(5))
    xSize = 8
    ySize = xSize/GOLDEN_RATIO

    fig = plt.figure(figsize=(xSize, ySize)) 
    ax = fig.add_subplot(111)
    ax.errorbar(countr, HR, xerr=err_cr, yerr=err_HR, fmt='b.', ecolor='b')
    ax.plot(xfine, line(xfine, popt[0], popt[1]), 'r-', linewidth=1.5)
    ax.hlines(HR_est, xmin=xmin, xmax=xmax, color='g', linewidth=1.5)
    #ax.hlines(HR_est+sigma_est, xmin=xmin, xmax=xmax, color='r', linestyle='--', linewidth=1)
    #ax.hlines(HR_est-sigma_est, xmin=xmin, xmax=xmax, color='r', linestyle='--', linewidth=1)
    xfill = np.linspace(xmin, xmax, 100)
    ax.fill_between(xfill, HR_est-sigma_est, HR_est+sigma_est, color='green', alpha=0.3)
    ax.set_xlim([xmin,xmax])
    ax.set_xlabel('Count rate')
    ax.set_ylabel('HR')
    ax.annotate(host+'-'+ulx, xy = (0., 0.), xycoords = 'axes fraction', xytext=(0., 1.05), textcoords='axes fraction', fontsize=8.)
    ax.annotate(time, xy = (0., 0.5), xycoords = 'axes fraction', xytext=(1., 1.05), textcoords='axes fraction', horizontalalignment='right', fontsize=8.)
    ax.annotate('Co = '+str(float("{0:.3f}".format(HR_ests)))+' +/- '+str(float("{0:.3f}".format(sigma_ests)))+',   Li = ['+str(float("{0:.3f}".format(popt[1])))\
        +' +/- '+str(float("{0:.3f}".format(pcov[1,1]**0.5)))+', '+str(float("{0:.3f}".format(popt[0])))+' +/- '+str(float("{0:.3f}".format(pcov[0,0]**0.5)))+']',\
         xy = (0., 0.6), xycoords = 'axes fraction', xytext=(0.5, 1.05), textcoords='axes fraction', horizontalalignment='center', fontsize=8.)

    if savefig:
        plt.savefig(filename, bbox_inches='tight', dpi=200)
    
    if show:
        plt.show()
    
    plt.close()


def HR_group(HR, err_HR):
    w = 1./(err_HR*err_HR)
    HR_est = (w*HR).sum()/w.sum()
    sigma_est = 1./np.sqrt(w.sum())
    return HR_est, sigma_est

def line(countrs, a, b):
    return a*countrs + b


def convertovariable(extract):
    obsid = extract[:,0]
    mjd = extract[:,1].astype(np.float)
    texp = extract[:,2].astype(np.float)
    countr = extract[:,3].astype(np.float)
    err_countr = extract[:,4].astype(np.float)
    HR = extract[:,5].astype(np.float)
    err_HR = extract[:,6].astype(np.float)
    rmf = extract[:,7].astype(np.int)
    return obsid, mjd, texp, countr, err_countr, HR, err_HR, rmf


    
if __name__ == "__main__":
    """Adding weighted spectra in one rmf group and spectral characteristic (hard or soft state)."""
    
    plot = False
    # Reading input file
    if len(sys.argv) == 1:
        input_file = str(raw_input("Write your coaddition input filename: "))
        input_plot = str(raw_input("Do you want to create HID and lc plot? (y/n): "))
        if input_plot == "y" or input_plot == "yes":
            plot = True
        else:
            print "Plots are not going to be created."

    else:
        while (len(sys.argv) > 1):
            option = sys.argv[1]; del sys.argv[1]
            if option == '-f':
                input_file = sys.argv[1]; del sys.argv[1]
            elif option == '-p':
                plot = True
            else:
                print 'Error[1] no option available! -f is expected'

    input_object, nobject = ex.read_data(input_file)

    now = datetime.datetime.now()
    time = now.strftime("%Y-%m-%d %H:%M")

    for i in range(0,nobject):
        ulx = input_object[:,0][i] # ULX name or ID, the same as extraction
        comp_path = input_object[:,1][i] # ULX directory, containing cls directory produced from reduction process, the same as extraction
        host = input_object[:,2][i] # galaxy host name

        #addspec_dir = comp_path+'addspec_'+ulx
        old_dir = os.getcwd()
        coaddir = comp_path+'coadd'+host+ulx

        #os.makedirs(coaddir)

        if not os.path.exists(coaddir):
            os.makedirs(coaddir)
        else:
            shutil.rmtree(coaddir)
            os.makedirs(coaddir)
    
        os.system('cp '+old_dir+'/rmf/swxpc0to12s0_20010101v012.rmf '+coaddir+'/1.rmf') # copy and change name of 
        os.system('cp '+old_dir+'/rmf/swxpc0to12s0_20070101v012.rmf '+coaddir+'/2.rmf') # corresponding rmf file
        os.system('cp '+old_dir+'/rmf/swxpc0to12s6_20010101v014.rmf '+coaddir+'/3.rmf')
        os.system('cp '+old_dir+'/rmf/swxpc0to12s6_20090101v014.rmf '+coaddir+'/4.rmf')
        os.system('cp '+old_dir+'/rmf/swxpc0to12s6_20110101v014.rmf '+coaddir+'/5.rmf')
        os.system('cp '+old_dir+'/rmf/swxpc0to12s6_20130101v014.rmf '+coaddir+'/6.rmf')

        os.chdir(comp_path)

        extract, nextract = ex.read_data(comp_path+'extract_product_'+ulx)
        countr = extract[:,3].astype(np.float)
        err_countr = extract[:,4].astype(np.float)
        SN = countr/err_countr

        SN_HID = extract[SN >= 3.]
        obsids, mjds, texps, countrs, err_countrs, HRs, err_HRs, rmfs = convertovariable(SN_HID)
        SNs = countrs/err_countrs
        HR_ests, sigma_ests = HR_group(HRs, err_HRs)
        popt, pcov = curve_fit(line, countrs, HRs, sigma=err_HRs)

        forhid = open('hid_'+host+ulx,'w')
        forhid.write('%s\n' %('#HID '+host+' '+ulx))
        forhid.write('%s%6.3f %s%6.3f   %s%6.3f %s%6.3f %s%6.3f %s%6.3f %s\n\n' %('#Co=', HR_ests, '+/-', sigma_ests, 'Li=[', popt[1], '+/-', pcov[1,1]**0.5, ',', popt[0], '+/-', pcov[0,0]**0.5, ']'))
        forhid.write('%12s%15s%15s%15s%17s%15s%15s%12s%15s\n' %('#OBSID','MJD','EXPTIME','Count rate','err_Count rate','HR','err_HR','Group rmf', 'S/N'))
        for h in range(len(SN_HID)):
            forhid.write('%12s%15.5E%15.5E%15.5E%17.5E%15.5E%15.5E%12i%15.5f\n' %(obsids[h], mjds[h], texps[h], countrs[h], err_countrs[h], HRs[h], err_HRs[h], rmfs[h], SNs[h]))
        forhid.close()

        


        if plot:
            plot_hidcut(time, countrs, err_countrs, HRs, err_HRs, HR_ests, sigma_ests, popt, host, ulx, show=False, savefig=True, filename='HID_'+host+ulx+'.png')

        
        extract = extract[SN >= 10.]
        for i in range(6):
            allext = extract[extract[:,7] == str(i+1)]
            

            if len(allext) >= 10:
                obsid, mjd, texp, countr, err_countr, HR, err_HR, rmf = convertovariable(allext)
                
                HR_est, sigma_est = HR_group(HR, err_HR)
                
                addsrc = [] 
                addbkg = [] 
                addarf = [] 

                grouprmf = open('coad_'+host+ulx+'gr'+str(i+1),'w')
                grouprmf.write('%s\n' %('#coadd '+host+' '+ulx))
                grouprmf.write('%7s%15s%11.5E%15.5E\n\n' %('#HR_est', 'sigma_est = ', HR_est, sigma_est))
                grouprmf.write('%12s%15s%15s%15s%17s%15s%15s%12s%15s\n' %('#OBSID','MJD','EXPTIME','Count rate','err_Count rate','HR','err_HR','Group rmf', 'S/N'))
                
                top = HR_est+sigma_est
                bottom = HR_est-sigma_est
                for j in range(len(allext)):
                    SN = countr[j]/err_countr[j]
                    if ((HR[j]-err_HR[j]) > top or ((HR[j]+err_HR[j])) < bottom):
                        grouprmf.write('%12s%15.5E%15.5E%15.5E%17.5E%15.5E%15.5E%12i%15.5f\n' %('#'+obsid[j], mjd[j], texp[j], countr[j], err_countr[j], HR[j], err_HR[j], rmf[j], SN))
                    else:
                        grouprmf.write('%12s%15.5E%15.5E%15.5E%17.5E%15.5E%15.5E%12i%15.5f\n' %(obsid[j], mjd[j], texp[j], countr[j], err_countr[j], HR[j], err_HR[j], rmf[j], SN))
                        formerdir = glob.glob(obsid[j]+'*')[0]
                        os.system('cp '+formerdir+'/'+'*src.pha '+coaddir+'/'+str(j)+'.pha')
                        os.system('cp '+formerdir+'/'+'*bkg.pha '+coaddir+'/'+str(j)+'b.pha')
                        os.system('cp '+formerdir+'/'+'*src.arf '+coaddir+'/'+str(j)+'.arf')
                        addsrc.append(str(j)+'.pha+')
                        addbkg.append(str(j)+'b.pha+')
                        addarf.append([str(j)+'.arf', SN])
                
                os.chdir(coaddir)

                srclist = 'addsrc'+str(i+1)
                addsrc[-1] = addsrc[-1][:-1]
                addsrc = np.array(addsrc)
                np.savetxt(srclist,addsrc,fmt='%s')
                
                bkglist = 'addbkg'+str(i+1)
                addbkg[-1] = addbkg[-1][:-1]
                addbkg = np.array(addbkg)
                np.savetxt(bkglist,addbkg,fmt='%s')
                
                arflist = 'addarf'+str(i+1)
                addarf = np.array(addarf)
                np.savetxt(arflist,addarf,fmt='%s')
                addarf = np.genfromtxt(arflist,dtype='str')
                weight = addarf[:,1].astype(np.float)
                weight = weight/weight.sum()
                farf = open(arflist,'w')
                for k in range(len(weight)):
                    farf.write('%s %g\n' %(addarf[k,0],weight[k]))
                farf.close()    

                grouprmf.close()

                addspectra(arflist,srclist,bkglist)

                source = srclist+'.pha'
                outputgrp = srclist+'grp20.pha'
                backgr = bkglist+'.pha'
                arf = arflist+'.arf'
                rmf = str(i+1)+'.rmf'
                os.system("grppha "+source+" "+outputgrp+" comm='sys_err 0.01 & group min 20 & chkey backfile "+backgr+" & chkey ancrfile "+arf+" & chkey respfile "+rmf+" & exit' clobber=yes chatter=0")

                os.chdir(comp_path)
                

        os.chdir(old_dir)
        print popt
