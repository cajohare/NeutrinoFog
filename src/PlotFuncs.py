#================================PlotFuncs.py==================================#
# Created by Ciaran O'Hare 2020

# Description:
# This file has many functions which are used throughout the project, but are
# all focused around the bullshit that goes into making the plots

#==============================================================================#

from numpy import *
from numpy.random import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib import colors
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as pe


pltdir = '../plots/'
pltdir_png = pltdir+'plots_png/'
def MySaveFig(fig,pltname,pngsave=True):
    fig.savefig(pltdir+pltname+'.pdf',bbox_inches='tight')
    if pngsave:
        fig.savefig(pltdir_png+pltname+'.png',bbox_inches='tight')

def cbar(mappable,extend='neither',minorticklength=8,majorticklength=10,\
            minortickwidth=2,majortickwidth=2.5,pad=0.2,side="right",orientation="vertical"):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(side, size="5%", pad=pad)
    cbar = fig.colorbar(mappable, cax=cax,extend=extend,orientation=orientation)
    cbar.ax.tick_params(which='minor',length=minorticklength,width=minortickwidth)
    cbar.ax.tick_params(which='major',length=majorticklength,width=majortickwidth)
    cbar.solids.set_edgecolor("face")

    return cbar

def line_background(lw,col):
    return [pe.Stroke(linewidth=lw, foreground=col), pe.Normal()]


# from matplotlib.backends.backend_pgf import FigureCanvasPgf
# mpl.backend_bases.register_backend('pdf', FigureCanvasPgf)


#==============================================================================#
def MySquarePlot(xlab='',ylab='',\
                 lw=2.5,lfs=45,tfs=25,size_x=13,size_y=12,Grid=False):
    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']

    fig = plt.figure(figsize=(size_x,size_y))
    ax = fig.add_subplot(111)

    ax.set_xlabel(xlab,fontsize=lfs)
    ax.set_ylabel(ylab,fontsize=lfs)

    ax.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    if Grid:
        ax.grid()
    return fig,ax

def MyDoublePlot(xlab1='',ylab1='',xlab2='',ylab2='',\
                 wspace=0.25,lw=2.5,lfs=45,tfs=25,size_x=20,size_y=8,Grid=False,tick_pad=7):
    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']
    fig, axarr = plt.subplots(1, 2,figsize=(size_x,size_y))
    gs = gridspec.GridSpec(1, 2)
    gs.update(wspace=wspace)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax1.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=tick_pad)
    ax1.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    ax2.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=tick_pad)
    ax2.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax1.set_xlabel(xlab1,fontsize=lfs)
    ax1.set_ylabel(ylab1,fontsize=lfs)

    ax2.set_xlabel(xlab2,fontsize=lfs)
    ax2.set_ylabel(ylab2,fontsize=lfs)

    if Grid:
        ax1.grid()
        ax2.grid()
    return fig,ax1,ax2


def MyDoublePlot_Vertical(xlab1='',ylab1='',xlab2='',ylab2='',\
                     hspace=0.25,lw=2.5,lfs=45,tfs=35,size_x=20,size_y=10,Grid=False):
        plt.rcParams['axes.linewidth'] = lw
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=tfs)
        mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']


        fig, axarr = plt.subplots(2,1,figsize=(size_x,size_y))
        gs = gridspec.GridSpec(2, 1)
        gs.update(hspace=hspace)
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])

        ax1.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
        ax1.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

        ax2.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
        ax2.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

        ax1.set_xlabel(xlab1,fontsize=lfs)
        ax1.set_ylabel(ylab1,fontsize=lfs)

        ax2.set_xlabel(xlab2,fontsize=lfs)
        ax2.set_ylabel(ylab2,fontsize=lfs)


        if Grid:
            ax1.grid()
            ax2.grid()
        return fig,ax1,ax2


def MyTriplePlot(xlab1='',ylab1='',xlab2='',ylab2='',xlab3='',ylab3='',\
                 wspace=0.25,lw=2.5,lfs=45,tfs=25,size_x=20,size_y=7,Grid=False,width_ratios=[1,1,1]):

    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']

    fig, axarr = plt.subplots(1, 3,figsize=(size_x,size_y))
    gs = gridspec.GridSpec(1, 3,width_ratios=width_ratios)
    gs.update(wspace=wspace)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])

    ax1.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax1.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax2.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax2.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax3.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax3.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax1.set_xlabel(xlab1,fontsize=lfs)
    ax1.set_ylabel(ylab1,fontsize=lfs)

    ax2.set_xlabel(xlab2,fontsize=lfs)
    ax2.set_ylabel(ylab2,fontsize=lfs)

    ax3.set_xlabel(xlab3,fontsize=lfs)
    ax3.set_ylabel(ylab3,fontsize=lfs)

    if Grid:
        ax1.grid()
        ax2.grid()
        ax3.grid()
    return fig,ax1,ax2,ax3
#==============================================================================#


import matplotlib.patches as patches

def CurvedArrow(x0,x1,y0,y1,alpha=0.7,color='orangered',connectionstyle="arc3,rad=-0.3",\
                style = "Simple, tail_width=2, head_width=12, head_length=12",zorder=10,**kw):
    kw = dict(arrowstyle=style, color=color,alpha=alpha)
    a1 = patches.FancyArrowPatch((x0, y0), (x1, y1),connectionstyle=connectionstyle,zorder=zorder,**kw)
    plt.gca().add_patch(a1)
    return


def SILimits(ax,Annotations=True,facecolor=[0.0, 0.62, 0.38],edgecolor='darkgreen',alph=1):
    ymax = ax.get_ylim()[1]
    plt.sca(ax)
    pek = line_background(5,'k')

    # Expt limits
    CRESST = loadtxt("../data/WIMPLimits/SI/CRESST.txt")
    plt.fill_between(CRESST[:,0], CRESST[:,1],edgecolor=None,y2=ymax,facecolor=[0.74, 0.56, 0.56])
    plt.plot(CRESST[:,0], CRESST[:,1],color=[0.8, 0.25, 0.33],linewidth=3,path_effects=pek)

    CDMSLite = loadtxt("../data/WIMPLimits/SI/CDMSLite.txt")
    plt.fill_between(CDMSLite[:,0], CDMSLite[:,1],edgecolor=None,y2=ymax,facecolor=[0.27, 0.51, 0.71])
    plt.plot(CDMSLite[:,0], CDMSLite[:,1],color="blue",linewidth=3,path_effects=pek)

    DarkSide = loadtxt("../data/WIMPLimits/SI/DarkSide.txt")
    plt.fill_between(DarkSide[:,0], DarkSide[:,1],edgecolor=None,y2=ymax,facecolor="forestgreen",alpha=0.5)
    plt.plot(DarkSide[:,0], DarkSide[:,1],color="green",linewidth=3,path_effects=pek)

    PandaX = loadtxt("../data/WIMPLimits/SI/PandaX.txt")
    plt.fill_between(PandaX[:,0], PandaX[:,1],edgecolor=None,y2=ymax,facecolor="teal")
    plt.plot(PandaX[:,0], PandaX[:,1],color="navy",linewidth=3,path_effects=pek)

    XENON1T = loadtxt("../data/WIMPLimits/SI/XENON1T.txt")
    plt.fill_between(XENON1T[:,0], XENON1T[:,1],edgecolor=None,y2=ymax,facecolor=facecolor,alpha=0.9)
    plt.plot(XENON1T[:,0], XENON1T[:,1],color=edgecolor,linewidth=3,path_effects=pek)

    dat = loadtxt('../data/WIMPLimits/SI/LUX.txt')
    plt.plot(dat[:,0], dat[:,1],color='crimson',linewidth=3,path_effects=pek)


    EDELWEISS = loadtxt("../data/WIMPLimits/SI/EDELWEISS.txt")
    plt.plot(EDELWEISS[:,0], EDELWEISS[:,1],color=[0.67, 0.31, 0.32],linewidth=3,path_effects=pek)

    PICO60 = loadtxt("../data/WIMPLimits/SI/PICO60.txt")
    plt.plot(PICO60[:,0], PICO60[:,1],color=[0.5, 0.0, 0.13],linewidth=3,path_effects=pek)

    PICO2L = loadtxt("../data/WIMPLimits/SI/PICO2L.txt")
    plt.plot(PICO2L[:,0], PICO2L[:,1],color=[0.5, 0.0, 0.13],linewidth=3,path_effects=pek)

    DAMA1 = loadtxt("../data/WIMPLimits/SI/DAMA1.txt")
    DAMA2 = loadtxt("../data/WIMPLimits/SI/DAMA2.txt")
    plt.fill_between(DAMA1[:,0], DAMA1[:,1],edgecolor=None,y2=1.0e-45,facecolor='forestgreen')
    plt.fill_between(DAMA2[:,0], DAMA2[:,1],edgecolor=None,y2=1.0e-45,facecolor='forestgreen')
    plt.plot(DAMA1[:,0], DAMA1[:,1],color='darkslategray',linewidth=3,path_effects=pek)
    plt.plot(DAMA2[:,0], DAMA2[:,1],color='darkslategray',linewidth=3,path_effects=pek)

    COSINE = loadtxt("../data/WIMPLimits/SI/COSINE-100.txt")
    plt.plot(COSINE[:,0], COSINE[:,1],color="gold",linewidth=3,path_effects=pek)

    dat = loadtxt('../data/WIMPLimits/SI/DEAP-3600.txt')
    plt.plot(dat[:,0], dat[:,1],color='#4ff09d',linewidth=3,path_effects=pek)

    dat = loadtxt('../data/WIMPLimits/SI/XENON1T-Migdal.txt')
    plt.plot(dat[:,0],dat[:,1],':',lw=3,color='darkgreen',alpha=0.8)

    dat = loadtxt('../data/WIMPLimits/SI/NEWS-G.txt')
    plt.plot(dat[:,0], dat[:,1],color='m',linewidth=3,path_effects=pek)
    return


def MakeLimitPlot_SI(Annotations=True,Collected=False,\
                     xmin=0.1,xmax=1.0e4,ymax=1.0e-36,ymin=1.0e-51,\
                     facecolor=[0.0, 0.62, 0.38],edgecolor='darkgreen',edgecolor_collected='darkgray',\
                     alph=0.5,lfs=35,tfs=25,\
                     xlab=r"DM mass [GeV$/c^2$]",ylab=r"SI DM-nucleon cross section [cm$^2$]"):
    pek = line_background(5,'k')

    fig,ax = MySquarePlot(xlab,ylab,lfs=lfs,tfs=tfs)


    if Collected:
        AllLimits = loadtxt("../data/WIMPLimits/SI/AllLimits-2021.txt")
        plt.fill_between(AllLimits[:,0], AllLimits[:,1],edgecolor=None,y2=ymax,facecolor=facecolor,alpha=alph,zorder=0)
        plt.plot(AllLimits[:,0], AllLimits[:,1],color=edgecolor_collected,linewidth=3,alpha=alph,zorder=0.01)
        for lim in ['CRESST','CDMSLite','DarkSide','PandaX','Xenon1T','LUX','DEAP-3600','EDELWEISS','PICO60','PICO2L','COSINE-100','NEWS-G']:
            dat = loadtxt("../data/WIMPLimits/SI/"+lim+".txt")
            plt.plot(dat[:,0], dat[:,1],color=edgecolor,linewidth=3,alpha=alph,zorder=0)

        DAMA1 = loadtxt("../data/WIMPLimits/SI/DAMA1.txt")
        DAMA2 = loadtxt("../data/WIMPLimits/SI/DAMA2.txt")
        plt.plot(DAMA1[:,0], DAMA1[:,1],color=edgecolor,linewidth=3,alpha=alph,zorder=0)
        plt.plot(DAMA2[:,0], DAMA2[:,1],color=edgecolor,linewidth=3,alpha=alph,zorder=0)

    else:
        if Annotations:
            plt.text(0.12,1.7e-38,r"{\bf CRESST}",color=[0.8, 0.25, 0.33],fontsize=24,rotation=0)
            plt.text(0.4,1.0e-40,r"{\bf CDMSlite}",color="blue",fontsize=22,rotation=0)
            plt.text(1.65,1e-41,r"{\bf DarkSide}",color="green",fontsize=22,rotation=0,ha='right')
            plt.text(850.0,2.e-45,r"{\bf PandaX}",color="navy",fontsize=22,rotation=19)
            plt.text(1800.0,0.57e-45,r"{\bf XENON1T}",color="darkgreen",fontsize=22,rotation=19.3)
            plt.text(6.5,2e-44,r"{\bf EDELWEISS}",color=[0.67, 0.31, 0.32],fontsize=18,rotation=-25)
            plt.text(2000.0,5.4e-43,r"{\bf PICO60}",color=[0.5, 0.0, 0.13],fontsize=22,rotation=18)
            plt.text(2000.0,5.9e-41,r"{\bf PICO2L}",color=[0.5, 0.0, 0.13],fontsize=22,rotation=19)
            plt.text(21.0,1e-39,r"{\bf DAMA}",color='darkslategray',fontsize=22)
            plt.text(1200.0,0.65e-41,r"{\bf COSINE-100}",color="gold",fontsize=22,rotation=19)
            plt.text(1.5e3,3.5e-44,r'{\bf DEAP-3600}',color='#4ff09d',fontsize=22,rotation=19)
            plt.text(4000.0,8.7e-45,r"{\bf LUX}",color="crimson",fontsize=21,rotation=19)
            plt.text(0.86,0.4e-39,'XENON1T \n (Migdal)',alpha=0.8,color='darkgreen',fontsize=18,ha='right')
            plt.text(1.8,0.5e-38,r"{\bf NEWS-G}",color="m",fontsize=22,rotation=-16)

            plt.arrow(0.48, 0.83, -0.04, -0.04, transform=ax.transAxes,
                  length_includes_head=True,
                  head_width=0.01, head_length=0.01, overhang=0.4,
                  edgecolor='darkslategray',facecolor='darkslategray')
            plt.arrow(0.54, 0.83, 0.04, -0.08, transform=ax.transAxes,
                  length_includes_head=True,
                  head_width=0.01, head_length=0.01, overhang=0.4,
                  edgecolor='darkslategray',facecolor='darkslategray')
            CurvedArrow(2.8e3,4e3,0.9e-44,0.5e-44,alpha=1,color='navy',connectionstyle="arc3,rad=-0.3",\
            style = "Simple, tail_width=2, head_width=6, head_length=8")

        SILimits(ax,facecolor=facecolor,edgecolor=edgecolor,alph=alph)

    # Labels
    plt.yscale('log')
    plt.xscale('log')
    plt.yticks(10.0**arange(-51,-30,1),fontsize=25)
    ax.tick_params(which='major',pad=10)
    ax.set_xlim(left=xmin, right=xmax)
    ax.set_ylim(bottom=ymin, top=ymax)
    return fig,ax

def SDpLimits(ax,Annotations=True,facecolor=[0.7111, 0.24352, 0.1755],edgecolor='darkred',alph=1,NeutrinoLimits=False):
    ymax = ax.get_ylim()[1]
    plt.sca(ax)
    pek = line_background(5,'k')

    PICASSO = loadtxt("../data/WIMPLimits/SDp/PICASSO.txt")
    plt.fill_between(PICASSO[:,0], PICASSO[:,1],edgecolor=None,y2=ymax,facecolor='tomato',zorder=0)
    plt.plot(PICASSO[:,0], PICASSO[:,1],color='orangered',linewidth=3,path_effects=pek,zorder=0)

    PICO60 = loadtxt("../data/WIMPLimits/SDp/PICO60.txt")
    plt.fill_between(PICO60[:,0], PICO60[:,1],edgecolor=None,y2=ymax,facecolor=facecolor,zorder=0)
    plt.plot(PICO60[:,0], PICO60[:,1],color=[0.5, 0.0, 0.13],linewidth=3,path_effects=pek,zorder=0)

    PICO2L = loadtxt("../data/WIMPLimits/SDp/PICO2L.txt")
    plt.plot(PICO2L[:,0], PICO2L[:,1],color=[0.5, 0.0, 0.13],linewidth=3,path_effects=pek,zorder=0)

    COUPP = loadtxt("../data/WIMPLimits/SDp/COUPP.txt")
    plt.plot(COUPP[:,0], COUPP[:,1],color=[1.0, 0.55, 0.41],linewidth=3,path_effects=pek,zorder=0)

    KIMS = loadtxt("../data/WIMPLimits/SDp/KIMS.txt")
    plt.plot(KIMS[:,0], KIMS[:,1],color=[0.81, 0.44, 0.69],linewidth=3,path_effects=pek,zorder=0)


    if NeutrinoLimits:
        IceCube = loadtxt("../data/WIMPLimits/SDp/IceCube-tt.txt")
        plt.plot(IceCube[:,0], IceCube[:,1],'--',color=[0.55, 0.71, 0.0],linewidth=3,path_effects=pek)

        SK = loadtxt("../data/WIMPLimits/SDp/SuperK-tt.txt")
        plt.plot(SK[:,0], SK[:,1],'--',color=[0.24, 0.71, 0.54],linewidth=3,path_effects=pek)

    return

def SDnLimits(ax,Annotations=True,facecolor=[0.0, 0.62, 0.38],edgecolor=[0.5, 0.0, 0.13],alph=1):
    ymax = ax.get_ylim()[1]
    plt.sca(ax)
    pek = line_background(5,'k')

    dat = loadtxt('../data/WIMPLimits/SDn/CRESST.txt')
    plt.fill_between(dat[:,0], dat[:,1],edgecolor=None,y2=ymax,facecolor=col_alpha([0.8, 0.25, 0.33],0.5))
    plt.plot(dat[:,0],dat[:,1],color=[0.8, 0.25, 0.33],linewidth=3,path_effects=pek)

    dat = loadtxt('../data/WIMPLimits/SDn/CDMSlite.txt')
    plt.fill_between(dat[:,0], dat[:,1],edgecolor=None,y2=ymax,facecolor=col_alpha('blue',0.5))
    plt.plot(dat[:,0],dat[:,1],color='blue',linewidth=3,path_effects=pek)

    dat = loadtxt('../data/WIMPLimits/SDn/CDMS.txt')
    plt.plot(dat[:,0],dat[:,1],color='blue',linewidth=3,path_effects=pek)

    PICASSO = loadtxt("../data/WIMPLimits/SDn/PICASSO.txt")
    plt.plot(PICASSO[:,0], PICASSO[:,1],color='darkred',linewidth=3,path_effects=pek)

    dat = loadtxt("../data/WIMPLimits/SDn/PandaX.txt")
    plt.fill_between(dat[:,0], dat[:,1],edgecolor=None,y2=ymax,facecolor=col_alpha('navy',0.5))
    plt.plot(dat[:,0], dat[:,1],color='navy',linewidth=3,path_effects=pek)

    dat = loadtxt("../data/WIMPLimits/SDn/LUX.txt")
    plt.plot(dat[:,0], dat[:,1],color='crimson',linewidth=3,path_effects=pek)

    dat = loadtxt("../data/WIMPLimits/SDn/XENON100.txt")
    plt.plot(dat[:,0], dat[:,1],color='darkgreen',linewidth=3,path_effects=pek)

    dat = loadtxt("../data/WIMPLimits/SDn/XENON1T.txt")
    plt.fill_between(dat[:,0], dat[:,1],edgecolor=None,y2=ymax,facecolor=facecolor)
    plt.plot(dat[:,0], dat[:,1],color='darkgreen',linewidth=3,path_effects=pek)

    dat = loadtxt('../data/WIMPLimits/SI/XENON1T-Migdal.txt')
    C = (4/3)*(0.264*(0.5+1)/0.5*(0.293)**2+0.2129*(1.5+1)/1.5*(0.242)**2)
    plt.plot(dat[:,0],dat[:,1]*(131**2/C),':',lw=3,color='darkgreen',alpha=0.8)
    return



def MakeLimitPlot_SDn(Annotations=True,Collected=False,\
                     xmin=0.1,xmax=1.0e4,ymax = 1.0e-32,ymin = 1.0e-46,\
                     facecolor=[0.0, 0.62, 0.38],edgecolor='darkgreen',edgecolor_collected='darkgray',\
                     alph=0.5,lfs=35,tfs=25,\
                     xlab=r"DM mass [GeV$/c^2$]",ylab=r"SD DM-neutron cross section [cm$^2$]",NeutrinoLimits=False):
    pek = line_background(5,'k')

    fig,ax = MySquarePlot(xlab,ylab,lfs=lfs,tfs=tfs)

    if Collected:
        AllLimits = loadtxt("../data/WIMPLimits/SDn/AllLimits-2021.txt")
        plt.fill_between(AllLimits[:,0], AllLimits[:,1],edgecolor=None,y2=ymax,facecolor=facecolor,alpha=alph,zorder=0)
        plt.plot(AllLimits[:,0], AllLimits[:,1],color=edgecolor_collected,linewidth=3,alpha=alph,zorder=0.01)
        for lim in ['PandaX','LUX','XENON100','XENON1T','CRESST','CDMSlite','PICASSO','CDMS']:
            dat = loadtxt("../data/WIMPLimits/SDn/"+lim+".txt")
            plt.plot(dat[:,0], dat[:,1],color=edgecolor,linewidth=3,alpha=alph,zorder=0)
    else:
        if Annotations:
            plt.text(850.0,2.7e-40,r"{\bf PandaX}",color='navy',fontsize=22,rotation=17)
            plt.text(4000.0,1.3e-39,r"{\bf LUX}",color='crimson',fontsize=22,rotation=17)
            plt.text(300.0,1e-39,r"{\bf XENON100}",color='darkgreen',fontsize=22,rotation=16.5)
            plt.text(1300.0,5.5e-41,r"{\bf XENON1T}",color='darkgreen',fontsize=22,rotation=18)
            plt.text(0.9,3.5e-35,'XENON1T \n (Migdal)',alpha=0.8,color='darkgreen',fontsize=18,ha='right')
            plt.text(0.3,1.9e-33,r"{\bf CRESST}",color=[0.8, 0.25, 0.33],fontsize=24,rotation=0)
            plt.text(0.6,3.0e-36,r"{\bf CDMSlite}",color="blue",fontsize=22,rotation=0)
            plt.text(2000,5.5e-37,r"{\bf CDMS}",color="blue",fontsize=22,rotation=17)
            plt.text(1000.0,3.6e-34,r"{\bf PICASSO}",color='darkred',fontsize=22,rotation=17)
            CurvedArrow(2.8e3,4e3,1.3e-39,0.8e-39,alpha=1,color='navy',connectionstyle="arc3,rad=-0.3",\
                    style = "Simple, tail_width=2, head_width=6, head_length=8")

        SDnLimits(ax,facecolor=facecolor,edgecolor=edgecolor,alph=alph)

    # Labels
    plt.yscale('log')
    plt.xscale('log')
    plt.yticks(10.0**arange(-51,-20,1),fontsize=25)
    ax.tick_params(which='major',pad=10)
    ax.set_xlim(left=xmin, right=xmax)
    ax.set_ylim(bottom=ymin, top=ymax)
    return fig,ax


def MakeLimitPlot_SDp(Annotations=True,Collected=False,\
                     xmin=0.1,xmax=1.0e4,ymax = 1.0e-32,ymin = 1.0e-46,\
                     facecolor=[0.7111, 0.24352, 0.1755],edgecolor='darkred',edgecolor_collected='darkgray',\
                     alph=0.5,lfs=35,tfs=25,\
                     xlab=r"DM mass [GeV$/c^2$]",ylab=r"SD DM-proton cross section [cm$^2$]",NeutrinoLimits=False):
    pek = line_background(5,'k')

    fig,ax = MySquarePlot(xlab,ylab,lfs=lfs,tfs=tfs)

    if Collected:
        AllLimits = loadtxt("../data/WIMPLimits/SDp/AllLimits-2021.txt")
        plt.fill_between(AllLimits[:,0], AllLimits[:,1],edgecolor=None,y2=ymax,facecolor=facecolor,alpha=alph,zorder=0)
        plt.plot(AllLimits[:,0], AllLimits[:,1],color=edgecolor_collected,linewidth=3,alpha=alph,zorder=0.01)
        for lim in ['PICASSO','PICO60','PICO2L','COUPP','KIMS']:
            dat = loadtxt("../data/WIMPLimits/SDp/"+lim+".txt")
            plt.plot(dat[:,0], dat[:,1],color=edgecolor,linewidth=3,alpha=alph,zorder=0)
    else:
        if Annotations:
            plt.text(200.0,8.0e-38,r"{\bf PICASSO}",color='tomato',fontsize=21,rotation=19)
            plt.text(2500.0,1.7e-39,r"{\bf PICO60}",color='darkred',fontsize=19,rotation=20)
            plt.text(2500.0,2.5e-38,r"{\bf PICO2L}",color='darkred',fontsize=19,rotation=20)
            plt.text(200.0,3e-39,r"{\bf COUPP}",color=[1.0, 0.55, 0.41],fontsize=19,rotation=18)
            plt.text(3000.0,3.5e-37,r"{\bf KIMS}",color=[0.81, 0.44, 0.69],fontsize=19,rotation=18)

            if NeutrinoLimits:
                plt.text(2000.0,3e-41,r"{\bf IceCube} $\tau\bar{\tau}$",color=[0.55, 0.71, 0.0],fontsize=20,rotation=35)
                plt.text(15.0,5e-41,r"{\bf SK} $\tau\bar{\tau}$",color=[0.24, 0.71, 0.54],fontsize=25,rotation=0)

        SDpLimits(ax,facecolor=facecolor,edgecolor=edgecolor,alph=alph)

    # Labels
    plt.yscale('log')
    plt.xscale('log')
    plt.yticks(10.0**arange(-51,-20,1),fontsize=25)
    ax.tick_params(which='major',pad=10)
    ax.set_xlim(left=xmin, right=xmax)
    ax.set_ylim(bottom=ymin, top=ymax)
    return fig,ax

#==============================================================================#
def reverse_colourmap(cmap, name = 'my_cmap_r'):
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1-t[0],t[2],t[1]))
        reverse.append(sorted(data))

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL)
    return my_cmap_r
#==============================================================================#


#==============================================================================#
def col_alpha(col,alpha=0.1):
    rgb = colors.colorConverter.to_rgb(col)
    bg_rgb = [1,1,1]
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(rgb, bg_rgb)]
#==============================================================================#
