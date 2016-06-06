import healpy as hp
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as py
from matplotlib import rcParams
import ntpath

#import special colormaps as scinet python version does not have the new good colorbars
import colormaps as cmaps

#set up colormaps
py.register_cmap(name='viridis', cmap=cmaps.viridis)
py.register_cmap(name='inferno', cmap=cmaps.inferno)
py.register_cmap(name='magma', cmap=cmaps.magma)

cmapinf = cmaps.magma#inferno

#USAGE
if len(sys.argv)<3:
    sys.exit("USAGE: python <filein> <flatsky>")

rgbab    = cmapinf(0.51) #border colour
rgbaf    = cmapinf(0.9)  #font colour


#Update label sizes, etc to look good
if int(sys.argv[2])==0:
    py.rcParams.update({'axes.labelsize': 20, 'font.size': 22,
                        'legend.fontsize': 20, 'xtick.labelsize': 32,
                        'ytick.labelsize': 22, 
                        'axes.edgecolor': rgbab,'axes.labelcolor': rgbaf,
                        'axes.linewidth': 6, 'xtick.color': rgbaf,
                        'xtick.labelsize': 28})
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 24}
matplotlib.rc('font', **font)


#change font
py.rc('text', usetex=True)
py.rc('font',**{'family':'serif','serif':['Palatino']})

#filein for flatsky
infile=open(sys.argv[1],"rb")

#=====================HEALPIX FULLSKY=====================================
if int(sys.argv[2]) == 0:
    
    #setup fig
    dpi          = 300
    figsize_inch = 24, 16
    fig          = py.figure(figsize=figsize_inch,dpi=dpi)

    y = hp.read_map(sys.argv[1])
    y = np.abs(y)
    y = np.log10(y+1e-9)

    print "min, max, mean = ", np.min(y), np.max(y), np.mean(y)

    #Plot
    hp.mollview(y,fig=fig.number, xsize=figsize_inch[0]*dpi,cmap=cmapinf,cbar=None,title="", min=-1.8-6,max=2.48-6) 

    #make border and titles of plot
    autoAxis = py.gca().axis()
    rec = py.Rectangle( (autoAxis[0],autoAxis[2]),(autoAxis[1]-autoAxis[0]),
                    (autoAxis[3]-autoAxis[2]),fill=False,lw=18,color=rgbab) 
    rec = py.gca().add_patch(rec)
    rec.set_clip_on(False)
    
    py.annotate(r'0.00 $<$ z $<$ 1.25', xy=(0.01, 0.93),  xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='bottom',
                color = rgbaf, fontsize = 36
                )
    py.annotate(r'8Gpc, 4096$^3$ Box', xy=(0.01, 0.88),  xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='bottom',
                color = rgbaf, fontsize = 36
                )
    py.annotate(r'N = 6.5$\times$10$^6$', xy=(0.01, 0.83),  xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='bottom',
                color = rgbaf, fontsize = 36
                )
    py.annotate(r'tSZ', xy=(0.96, 0.83),  xycoords='axes fraction',
                horizontalalignment='right', verticalalignment='bottom',
                color = rgbaf, fontsize = 98
                )    
    
    #cut off extra tiny bit of white space on edges
    py.gca().set_axis_off()
    py.subplots_adjust(top = autoAxis[3], bottom = autoAxis[2], 
                    right = autoAxis[1], left = autoAxis[0], 
                    hspace = 0, wspace = 0)
    py.margins(0,0)
    py.gca().xaxis.set_major_locator(py.NullLocator())
    py.gca().yaxis.set_major_locator(py.NullLocator())


    #Add colorbar    
    ax = py.gca()
    image = ax.get_images()[0]
    cbaxes=fig.add_axes([0.5,0.17,0.47,0.03])
    cb = py.colorbar(image,cax=cbaxes,orientation="horizontal")#, fraction=0.05)
    cb.set_label(r'log $\Delta T$', labelpad=-100, fontsize=36)
    
    #Output image
    fileout = "Images/Fullsky_tSZ_"+ntpath.basename(sys.argv[1])[:-4]+".pdf"

    py.savefig(fileout,bbox_inches='tight',pad_inches=0)




#=====================FLATSKY=====================================
elif int(sys.argv[2]) == 1:

    zoom = 4096

    #get field of view and pixel size
    npix = np.fromfile(infile,dtype=np.int32,count=1)
    npix = np.fromfile(infile,dtype=np.int32,count=1)
    fov  = np.fromfile(infile,dtype=np.float32,count=1)
    fov = np.fromfile(infile,dtype=np.float32,count=1)
    fov = float(fov/2/np.pi*360 * zoom/4096)

    #get data and reshape to square
    y=np.fromfile(sys.argv[1],dtype=np.float32,count=npix**2)
    y = np.abs(y)
    y=np.reshape(y, (npix,npix))
    y = y[:zoom,-zoom:]
    y = np.log10(y+1e-9)

    print "min max mean y = ", np.min(y), np.max(y), np.mean(y)
    #set up and plot figure
    fig = py.figure(figsize=(12,12))
    ax = py.subplot(111)

    im = ax.imshow(y, vmin=-1.8-6,vmax=2.48-6,cmap=cmapinf,extent=[-fov/2,fov/2,-fov/2,fov/2],origin="upper")

    #set up colorbar
    cb = py.colorbar(im,shrink=0.78)#fraction=0.03)#, cax=cbaxes)
    cb.set_label('log Compton-y')

    #set labels
    ax.set_xlabel("Degrees")
    ax.set_ylabel("Degrees")

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
        ax.tick_params('both', length=10, width=4, which='major')
        ax.tick_params('both', length=6, width=3, which='minor')

    #output
    fileout = "Images/Flatsky_tSZ_"+ntpath.basename(sys.argv[1])[:-4]+"_fov_"+str(fov)+".pdf"
    py.savefig(fileout, bbox_inches="tight")


