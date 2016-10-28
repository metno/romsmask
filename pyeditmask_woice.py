#!/usr/bin/env python
######################################################
## Edits ROMS masks using a GUI
## Nov 2014
## rsoutelino@gmail.com
## Feb 2016
## editted by Marta Trodahl
## martat@met.no
######################################################
import os
import wx

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as Navbar
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
import	wx.lib.dialogs
from wx.lib.stattext import GenStaticText

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import scipy.io as sp
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap, cm
from scipy import shape

# CHANGES DONE (Marta Trodahl 2016)
#   - changed to plot with basemap
#   - added dialogboxes for user input of plotting parameters ... 
#   - changed the coastline feature : plotting the coastline from GHHSH instead of having a coastfile
#   - created new mask/unmask buttons to fit an editted colormap  
#   - added display of critical grid boxes(a version also for coupling with ice model), boxes are marked in red
#     Coupled version: Jens D. implemented the test for critical points 	
#   - added display of number of critical points left
#   - added a statusbar where additional information is displayed, s.a. new depth when 'unmasking'
#   - added checks on topography, and replacing the depth at new wet points with the surrounding mean!=hmin  	
#   - created a new button which automatically removes critical wetpoints, mask = 1 => 0
#   - implemented toggle/untoggle function for the buttons: released when a new button is chosen
#   - resolved trouble with functions not disconnecting when another is selected
#   - made the main toolbar vertical, and moved the panel tools
#   - mask/unmask: plotting only the selected point and its neighboring points (more efficient than redrawing the entire grid)

# NICE TIP TO DEBUG THIS PROGRAM: ================================
#   - comment out app.MainLoop at the last line of this script
#   - ipython --gui=wx
#   - run pyeditmask.py
#   - trigger the events and check out the objects in the shell
# ================================================================

global currentDirectory
currentDirectory = os.getcwd()

PROJECT_DIR = os.path.abspath(os.path.dirname(__file__))
DEFAULT_VMIN = -5#-1.3 
DEFAULT_VMAX= 1.5#2.1
DEFAULT_CMAP = plt.cm.Spectral#plt.cm.GnBu
DEFAULT_DEPTH_FOR_LAND = -50


# ROMS related objects ---------------------------------------------
class RomsGrid(object):
    """ 
    Stores and manipulates netcdf ROMS grid file information
    """
    def __init__(self,filename):
        self.filename = filename  
	self.modfile = nc.Dataset(filename, mode='r+')  
        self.ncfile = nc.Dataset(filename, mode='r+')
        self.lonr  = self.ncfile.variables['lon_rho'][:]
        self.latr  = self.ncfile.variables['lat_rho'][:]
        self.lonu  = self.ncfile.variables['lon_u'][:]
        self.latu  = self.ncfile.variables['lat_u'][:]
        self.lonv  = self.ncfile.variables['lon_v'][:]
        self.latv  = self.ncfile.variables['lat_v'][:]
        self.h     = self.ncfile.variables['h'][:]
        self.maskr = self.ncfile.variables['mask_rho'][:]
        self.masku = self.ncfile.variables['mask_u'][:]
        self.maskv = self.ncfile.variables['mask_v'][:]     
	self.maskcrit = self.modfile.variables['mask_rho'][:]

def uvp_mask(rfield):
    Mp, Lp = rfield.shape
    M      = Mp - 1
    L      = Lp - 1

    vfield = rfield[0:M,:] * rfield[1:Mp,:]
    ufield = rfield[:,0:L] * rfield[:,1:Lp]
    pfield = ufield[0:M,:] * ufield[1:Mp,:]

    return ufield, vfield, pfield
# -------------------------------------------------------------------


class App(wx.App):
    def OnInit(self):
        self.frame = Interface("PyEditMask 0.1.0", size=(1024,800))
        self.frame.Show()
        return True

class Interface(wx.Frame):
    def __init__(self, title=wx.EmptyString, pos=wx.DefaultPosition, 
                       size=wx.DefaultSize, style=wx.DEFAULT_FRAME_STYLE,
                       *args, **kwargs):
        wx.Frame.__init__(self, None, -1, "PyEditMask 0.1.0", pos=pos, 
                          size=size, style=style, *args, **kwargs)
        
        # Initializing toolbar
        self.toolbar = MainToolBar(self)
	# Initializing statusbar
	self.statusbar = self.CreateStatusBar()
	self.statusbar.SetStatusText('Additional information will be displayed here.')

	# BASIC LAYOUT OF THE NESTED SIZERS ======================
        mplpanel = wx.Panel(self, wx.ID_ANY, style=wx.SUNKEN_BORDER)
        mplpanel.SetBackgroundColour("white")

        # BOX 1 is the main sizer
        box1 = wx.BoxSizer(wx.HORIZONTAL)
        box1.Add(mplpanel, 15, wx.EXPAND)

        # BOX 2 is the inner sizer of the left big control panel
        box2 = wx.BoxSizer(wx.VERTICAL)

        # BOX 3 is the sizer of the right big parent panel, the one that will
        # serve as base for the child panel which will hold the matplotlib canvas
        box3 = wx.BoxSizer(wx.VERTICAL)

        # mplpanel content ========================================
        self.mplpanel = SimpleMPLCanvas(mplpanel)
        box3.Add(self.mplpanel.canvas, 1, flag=wx.CENTER) 

        # FINAL LAYOUT CONFIGURATIONS ============================
        self.SetAutoLayout(True)
        mplpanel.SetSizer(box3)

        self.SetSizer(box1)

        self.InitMenu()
	self.Layout()
        self.Centre()

    def InitMenu(self):
        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        fileMenu.Append(wx.ID_OPEN, u'&Open ROMS grid file')
        fileMenu.Append(wx.ID_OPEN, u'&Open coastline file')
        fileMenu.Append(wx.ID_SAVE, '&Save grid')
        fileMenu.AppendSeparator()

        qmi = wx.MenuItem(fileMenu, wx.ID_EXIT, '&Quit\tCtrl+W')
        opf = wx.MenuItem(fileMenu, wx.ID_OPEN, '&Open\tCtrl+O')
        opc = wx.MenuItem(fileMenu, wx.ID_OPEN, '&Open\tCtrl+O+C')
        svf = wx.MenuItem(fileMenu, wx.ID_SAVE, '&Save\tCtrl+S')
        fileMenu.AppendItem(qmi)
        # fileMenu.AppendItem(svf)

        self.Bind(wx.EVT_MENU, self.OnQuit, qmi)
        self.Bind(wx.EVT_MENU, self.toolbar.OnLoadGrid, opf)
        self.Bind(wx.EVT_MENU, self.toolbar.OnLoadCoastline, opc)
        self.Bind(wx.EVT_MENU, self.toolbar.OnSaveGrid, svf)

        menubar.Append(fileMenu, u'&PyEditMask')
        self.SetMenuBar(menubar)

    def OnQuit(self, e):
        self.Close()
        self.Destroy()

    def OnCloseWindow(self, e):
        self.Destroy()


class SimpleMPLCanvas(object):
    """docstring for SimpleMPLCanvas"""
    def __init__(self, parent):
        super(SimpleMPLCanvas, self).__init__()
        self.parent = parent
        self.plot_properties()
        self.make_navbar()
        
    def make_navbar(self):
        self.navbar = Navbar(self.canvas)   
        self.navbar.SetPosition((180,0)) # this is now working (MT 2016)


    def plot_properties(self):
        # Create matplotlib figure
        self.fig = Figure(facecolor='lightsteelblue', figsize=(12.5,8.5))
        self.canvas = FigureCanvas(self.parent, -1, self.fig)
        # self.fig.patch.set_alpha(0.5)
        self.ax   = self.fig.add_subplot(111)
	self.ax.xaxis.set_ticklabels([])
	self.ax.yaxis.set_ticklabels([])
        # tit = self.ax1.set_title("ROMS mask_rho", fontsize=12, fontweight='bold')
        # tit.set_position([0.9, 1.05])

class MainToolBar(object):
    def __init__(self, parent):
	self.cid = None
        self.currentDirectory = os.getcwd()
        self.parent = parent
        self.toolbar = parent.CreateToolBar(style=wx.TB_VERTICAL, winid=1,
                                            name="Toolbar")
        self.tools_params ={ 
            'load_grid': (load_bitmap('grid.png'), u"Load grid",
                        "Load ocean_grd.nc ROMS grid netcdf file"),
            'load_coastline': (load_bitmap('coast.png'), u"Load coastline",
                        "Load *.mat coastline file [lon / lat polygons]"),
            'save_grid': (load_bitmap('save.png'), u"Apply and save",
                        "Save changes to ocean_grd.nc ROMS grid netcdf file"),
            'set_land': (load_bitmap('L.png'), u"Set land",
                        "Set grid point to land"),
            'set_land_area': (load_bitmap('LA.png'), u"Set land area",
                        "Set poligon area to land"),
            'set_water': (load_bitmap('W.png'), u"Set water",
                        "Set grid point to water"),
            'set_water_area': (load_bitmap('WA.png'), u"Set water area",
                        "Set poligon area to water"),
	    'mask_single_wetpoints': (load_bitmap('point_on_ocean.png'), u"Remove all critical points (new ones may occur)",
                        "When button is toggled, activate by right clicking on the map. This sets all current critical points to land"),
#            'settings': (load_bitmap('settings.png'), u"PyEditMask settings",
#                        "PyEditMask configurations"),
            'quit': (load_bitmap('exit.png'), u"Quit",
                        "Quit PyEditMask"),
        }
        parent.btnDict = dict()

        self.createTool(self.toolbar, self.tools_params['load_grid'], 
                        self.OnLoadGrid)
        self.createTool(self.toolbar, self.tools_params['load_coastline'], 
                        self.OnLoadCoastline)
        self.createTool(self.toolbar, self.tools_params['save_grid'], 
                        self.OnSaveGrid)
        
        self.toolbar.AddSeparator()

	self.mask_tool,ID = self.createTool(self.toolbar, self.tools_params['set_land'], 
                                         self.OnSetLand, isToggle=True)
	# Add a dictionary key-value pair entry for the button.
        btnId = ID	  
	# Create a dictionary entry.
        parent.btnDict[ btnId ] = 'set_land'

	self.mask_area_tool,ID = self.createTool(self.toolbar, 
                                              self.tools_params['set_land_area'], 
                                              self.OnSetLandArea, isToggle=True)
	# Add a dictionary key-value pair entry for the button.
        btnId = ID	  
	# Create a dictionary entry.
        parent.btnDict[ btnId ] = 'set_land_area'

        self.unmask_tool,ID = self.createTool(self.toolbar, self.tools_params['set_water'], 
                                           self.OnSetWater, isToggle=True)
	# Add a dictionary key-value pair entry for the button.
        btnId = ID	  
	# Create a dictionary entry.
        parent.btnDict[ btnId ] = 'set_water'

        self.unmask_area_tool,ID = self.createTool(self.toolbar, 
                                                self.tools_params['set_water_area'], 
                                                self.OnSetWaterArea, isToggle=True)
	# Add a dictionary key-value pair entry for the button.
        btnId = ID	  
	# Create a dictionary entry.
        parent.btnDict[ btnId ] = 'set_water_area'

	self.mask_single_tool,ID = self.createTool(self.toolbar, 
                                                self.tools_params['mask_single_wetpoints'], 
                                                self.OnMaskSingleWetpoint, isToggle=True)
	# Add a dictionary key-value pair entry for the button.
        btnId = ID	  
	# Create a dictionary entry.
        parent.btnDict[ btnId ] = 'mask_single_wetpoints'

        self.toolbar.AddSeparator()

#        self.createTool(self.toolbar, self.tools_params['settings'], 
#                        self.OnSettings)
        self.createTool(self.toolbar, self.tools_params['quit'], 
                        self.parent.OnQuit)

        self.toolbar.Realize()

    def createTool(self, parent, params, evt, isToggle=False):
	ID = wx.NewId()
	tool = parent.AddTool(ID, bitmap=params[0], shortHelpString=params[1],
                    longHelpString=params[2], isToggle=isToggle)
        self.parent.Bind(wx.EVT_TOOL, evt, id=tool.GetId())
        return tool,ID

    def masktest(self,maskr,i,j):
	mask = 1
	if (maskr[j-1,i] + maskr[j+1,i] == 0 or \
	    maskr[j,i-1] + maskr[j,i+1] == 0) :
		mask = 0
	return mask

    def plot_grid(self,grd,lon0):
	mplpanel = app.frame.mplpanel
        ax = mplpanel.ax	
	# Create a Basemap instance.
	lonc=grd.lonr
        latc=grd.latr
	self.m = Basemap(projection=str(self.proj),lon_0=lon0,lat_1=77.5,lat_2=77.5,\
                llcrnrlat=latc[0,0],urcrnrlat=latc[-1,-1],\
                llcrnrlon=lonc[0,0],urcrnrlon=lonc[-1,-1],\
                rsphere=6371200.,resolution='i',ax=mplpanel.ax)
#	self.m = Basemap(projection=str(self.proj),lon_0=lon0,lat_0=90,lat_ts=60,\
#                llcrnrlat=latc[0,0],urcrnrlat=latc[-1,-1],\
#                llcrnrlon=lonc[0,0],urcrnrlon=lonc[-1,-1],\
#                rsphere=6371200.,resolution='i',ax=mplpanel.ax)
        self.lo, self.la = self.m(grd.lonr, grd.latr)
	for i in range(2,len(grd.maskr[0,:])-1):
      		for j in range(2,len(grd.maskr[:,0])-1):
			if (grd.maskr[j,i] == 1 and self.masktest(grd.maskr,i,j)==0):
				grd.maskcrit[j,i]=-10.		
			else:
				grd.maskcrit[j,i]=grd.maskr[j,i]	
	self.pcolor = ax.pcolormesh(self.lo, self.la, grd.maskcrit,vmin=DEFAULT_VMIN, vmax=DEFAULT_VMAX, 
                                   cmap=DEFAULT_CMAP)
        ax.set_aspect('equal')
	a=(grd.maskcrit[:] == -10).sum()
	ax.set_title('Number of critical points: ' + str(a),loc='Left')
        return ax
	
    def plot_maskcrit(self,grd,lon0,line,col):
	# Create a Basemap instance.
    	mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
	for i in range(col-3,col+3):
      		for j in range(line-3,line+3):
			if (grd.maskr[j,i] == 1 and self.masktest(grd.maskr,i,j)==0):
				grd.maskcrit[j,i]=-10.
			else:
				grd.maskcrit[j,i]=grd.maskr[j,i]	
	
	self.pcolor = ax.pcolormesh(self.lo[line-3:line+3,col-3:col+3], self.la[line-3:line+3,col-3:col+3], grd.maskcrit[line-3:line+3,col-3:col+3],vmin=DEFAULT_VMIN,vmax=DEFAULT_VMAX,cmap=DEFAULT_CMAP)
        ax.set_aspect('equal')
	a=(grd.maskcrit[:] == -10).sum()
	ax.set_title('Number of critical points: ' + str(a),loc='Left')
    
    def ask(self,parent, message='', title='',default_value=''):
	dlg = wx.TextEntryDialog(parent,message, title,defaultValue=default_value)
	dlg.ShowModal()	
	result = dlg.GetValue()
	dlg.Destroy()
       	return result

    def OnLoadGrid(self, evt):
	mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
	# Call Dialog, asking user to specify parameters for the basemap
        self.lon0 = self.ask(self.parent,message = 'lon0 (A4km = A20km = 58) \n See "grid_mapping" in grid file',title = 'Please enter plotting parameters')
	self.proj = self.ask(self.parent,message = 'proj (A4km = A20km = stere)',title = 'Please enter plotting parameters')
	openFileDialog = wx.FileDialog(self.parent, "Open grid netcdf file [*.nc]",
                                       "/ops/hindcast/roms", " ",
                                       "netcdf files (*.nc)|*.nc",
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return     # the user changed idea...
        ax = mplpanel.ax
	self.filename = openFileDialog.GetPath()
	grd = RomsGrid(self.filename)        
	self.plot_grid(grd,self.lon0)
	ax.set_xlim([self.lo.min(), self.lo.max()])
        ax.set_ylim([self.la.min(), self.la.max()])
	ax.plot(self.lo,self.la, 'k', alpha=0.2)
        ax.plot(self.lo.transpose(), self.la.transpose(), 'k', alpha=0.2)
	mplpanel.canvas.draw()
	self.grd = grd
        self.grd.hmin = grd.ncfile.variables['h'][:].min()
	
    def OnLoadCoastline(self, evt):
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        lo,la = self.lo, self.la
        self.m.drawcoastlines(linewidth=0.5,color='k')
	mplpanel.canvas.draw()

    def OnSaveGrid(self, evt):
        maskr = self.grd.maskr
        [masku, maskv, maskp] = uvp_mask(maskr)
        self.grd.ncfile.variables['mask_rho'][:] = maskr
        self.grd.ncfile.variables['mask_u'][:]   = masku
        self.grd.ncfile.variables['mask_v'][:]   = maskv
        self.grd.ncfile.variables['mask_psi'][:] = maskp
        self.grd.ncfile.variables['h'][:] = self.grd.h
        self.grd.ncfile.sync()

    def OnSetLand(self, evt):
        mplpanel = app.frame.mplpanel
	if self.cid != None:
		mplpanel.canvas.mpl_disconnect(self.cid)
		self.cid = None
	if self.mask_tool.IsToggled():
		self.cid = mplpanel.canvas.mpl_connect('button_press_event', self.mask)
        	# Deselects whichever button that was previously clicked.
        	for btnId in self.parent.btnDict.iterkeys() :   # Deselect all buttons including those not now selected.
        		self.toolbar.ToggleTool( btnId, 0 )

        	eventId = evt.GetId()     
        	self.toolbar.ToggleTool( eventId, 1 )   	# Show the associated button as being selected.
        	print self.parent.btnDict[ eventId ]            # Extract and print the buttons ID.

        else:
        	mplpanel.canvas.mpl_disconnect(self.cid)
		

    def OnSetLandArea(self, evt):
        mplpanel = app.frame.mplpanel
	if self.cid != None:
		mplpanel.canvas.mpl_disconnect(self.cid)
		self.cid = None 	
	# Deselect whichever button that was previously clicked.
	if self.mask_area_tool.IsToggled():
		self.cid = mplpanel.canvas.mpl_connect('button_press_event', 
        	                                       self.mask_area)        
		for btnId in self.parent.btnDict.iterkeys() :   # Deselect all buttons.
                	self.toolbar.ToggleTool( btnId, 0 )
        	eventId = evt.GetId()    
        	self.toolbar.ToggleTool( eventId, 1 )   	# Show the associated button as being selected.
        	print self.parent.btnDict[ eventId ]            # Extract and print the buttons ID.

        else:
        	mplpanel.canvas.mpl_disconnect(self.cid)
	

    def OnSetWater(self, evt):
	mplpanel = app.frame.mplpanel
	if self.cid != None:
		mplpanel.canvas.mpl_disconnect(self.cid)
		self.cid = None        
	if self.unmask_tool.IsToggled():
        	self.cid = mplpanel.canvas.mpl_connect('button_press_event', self.unmask)
		# Deselect whichever button that was previously clicked.
	        for btnId in self.parent.btnDict.iterkeys() :   # Deselect all buttons.
	        	self.toolbar.ToggleTool( btnId, 0 )
		eventId = evt.GetId()      
        	self.toolbar.ToggleTool( eventId, 1 )   	# Show the associated button as being selected.
        	print self.parent.btnDict[ eventId ]    	# Extract and print the buttons ID.
        else:
        	mplpanel.canvas.mpl_disconnect(self.cid)

    def OnSetWaterArea(self, evt):
        mplpanel = app.frame.mplpanel
	if self.cid != None:
		mplpanel.canvas.mpl_disconnect(self.cid)
		self.cid = None        
	if self.unmask_area_tool.IsToggled():
        	self.cid = mplpanel.canvas.mpl_connect('button_press_event', 
                                                       self.unmask_area)
       		# Deselect whichever button was previously clicked.
        	for btnId in self.parent.btnDict.iterkeys() :   # Deselect all buttons.
            		self.toolbar.ToggleTool( btnId, 0 )
        	eventId = evt.GetId()      
		self.toolbar.ToggleTool( eventId, 1 )   	# Show the associated button as being selected.
	        print self.parent.btnDict[ eventId ]            # Extract and print the buttons ID.

        else:
		mplpanel.canvas.mpl_disconnect(self.cid)


    def OnMaskSingleWetpoint(self, evt):
        mplpanel = app.frame.mplpanel
	if self.cid != None:
		mplpanel.canvas.mpl_disconnect(self.cid)
		self.cid = None   	
	if self.mask_single_tool.IsToggled():
        	self.cid = mplpanel.canvas.mpl_connect('button_press_event', 
                                                       self.mask_single_wetpoints)
        	# Deselect whichever button that was previously clicked.
	        for btnId in self.parent.btnDict.iterkeys() :   # Deselect all buttons.
		        self.toolbar.ToggleTool( btnId, 0 )
		eventId = evt.GetId()      
		self.toolbar.ToggleTool( eventId, 1 )  		# Show the associated button as being selected.
	        print self.parent.btnDict[ eventId ]           	# Extract and print the buttons ID.

        else:
        	mplpanel.canvas.mpl_disconnect(self.cid)

    def mask_single_wetpoints(self,evt):
	if evt.inaxes != app.frame.mplpanel.ax: return
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        x, y = evt.xdata, evt.ydata
	grd = self.grd
        for i in range(2,len(grd.maskr[0,:])-1):
      		for j in range(2,len(grd.maskr[:,0])-1):
			if (grd.maskr[j,i] == 1 and self.masktest(grd.maskr,i,j)==0): 
				grd.maskcrit[j,i] = 0
				grd.maskr[j,i] = 0
				grd.h[j,i] = grd.hmin
	self.plot_grid(grd,self.lon0)
	mplpanel.canvas.draw()
	
    def mask(self, evt):
        if evt.inaxes != app.frame.mplpanel.ax: return
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        x, y = evt.xdata, evt.ydata
	grd = self.grd
        line, col = find_lower_left_node(self.lo, self.la, x, y)
        grd.maskr[line, col] = 0 		# assigning new value
        grd.maskcrit[line,col] = 0
	grd.h[line, col] = self.grd.hmin	
        self.plot_maskcrit(grd,self.lon0,line,col)
	mplpanel.canvas.draw()			

    def mask_area(self, evt):
        if evt.inaxes != app.frame.mplpanel.ax: return
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        x, y = evt.xdata, evt.ydata
	lo, la = self.lo, self.la         
	button = evt.button
        if button == 1:
            p = ax.plot(x, y, 'ro')
            try:
                self.points.append(p)
                self.area.append( (x, y) )
            except AttributeError:
                self.points = [p]
                self.area = [ (x, y) ]

            mplpanel.canvas.draw()

        elif button == 3:
            grd = self.grd
            path = Path( self.area )
            a, b = lo.shape
            for i in range(a):
                for j in range(b):
                    if path.contains_point( [lo[i, j], 
                                             la[i, j] ] ) == 1:
         
			grd.maskr[i,j] = 0
                        grd.h[i,j] = grd.hmin
			grd.maskcrit[i,j] = 0
            ax.clear()
	    self.plot_grid(grd,self.lon0)	
	    ax.set_xlim([self.lo.min(), self.lo.max()])
	    ax.set_ylim([self.la.min(), self.la.max()])    	
	    ax.plot(self.lo,self.la, 'k', alpha=0.2)
            ax.plot(self.lo.transpose(), self.la.transpose(), 'k', alpha=0.2)
	    mplpanel.canvas.draw()
            del self.points, self.area

    def unmask(self, evt):
        if evt.inaxes != app.frame.mplpanel.ax: return
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        x, y = evt.xdata, evt.ydata
	grd = self.grd
	lo, la = self.lo, self.la
	line, col = find_lower_left_node(lo, la, x, y)
	nf = nc.Dataset(grd.filename, mode='r')
	m = nf.variables['h'][line-4:line+4,col-4:col+4]
	l,c=np.where(m>grd.hmin)
	grd.maskr[line, col] = 1 		# assigning new value
	grd.maskcrit[line,col] = 1
	grd.h[line, col] = np.mean(m[l,c])
	print np.mean(m[l,c])
	self.plot_maskcrit(grd,self.lon0,line,col)
	self.parent.StatusBar.SetStatusText('h(' + str(line)+ ','+str(col)+') was set to '+ str(np.mean(m[l,c]))+'m')	
	mplpanel.canvas.draw()

    def unmask_area(self, evt):
        if evt.inaxes != app.frame.mplpanel.ax: return
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        x, y = evt.xdata, evt.ydata
	button = evt.button
	
        if button == 1:
            p = ax.plot(x, y, 'ro')
	    try:
                self.points.append(p)
                self.area.append( (x, y) )
            except AttributeError:
                self.points = [p]
                self.area = [ (x, y) ]
	    mplpanel.canvas.draw()
        elif button == 3:
            grd = self.grd
            path = Path( self.area )
            lo, la = self.lo, self.la
	    a, b = lo.shape
	    for i in range(a):
                for j in range(b):
                    if path.contains_point( [lo[i, j], 
                                             la[i, j] ] ) == 1:			
			idx = np.where(grd.maskr[i-1:i+1,j-1:j+1]==1)
		        grd.maskr[i,j] = 1
	  	    	#print np.mean(grd.h[idx])
			m = grd.h[i-4:i+4,j-4:j+4]
			l,c = np.where(m>=grd.hmin)
		    	grd.h[i,j] = np.mean(m[l,c])#grd.h[i-1:i+1,j-1:j+1])
			grd.maskcrit[i,j] = 1
			print np.mean(m[l,c])#grd.h[i-1:i+1,j-1:j+1])	     
            ax.clear() 
	    self.plot_grid(grd,self.lon0)
	    ax.set_xlim([self.lo.min(), self.lo.max()])
	    ax.set_ylim([self.la.min(), self.la.max()])        
	    ax.plot(self.lo,self.la, 'k', alpha=0.2)
            ax.plot(self.lo.transpose(), self.la.transpose(), 'k', alpha=0.2)
	    mplpanel.canvas.draw()
            del self.points, self.area
	
def find_lower_left_node(x, y, x0, y0, n=4):
    # need to improve this, not very accurate yet
    dx = np.abs(x - x0); dx = dx / dx.max()
    dy = np.abs(y - y0); dy = dy / dy.max()
    dn = dx + dy    
    line, col, lola = [], [], []

    for k in range(n):
        fn = np.where(dn == dn.min())
        f1, f2 = int(fn[0]), int(fn[1])
        line.append(f1)
        col.append(f2)
        lola.append(x[f1, f2] + y[f1, f2])
        dn[f1, f2] = 1e20

    lola = np.array(lola)
    f = np.where(lola == lola.min())[0][0]

    line = line[f]
    col = col[f]

    return line, col

def near2d(x, y, x0, y0):
    """
    Find the indexes of the grid point that is
    nearest a chosen (x0, y0).
    Usage: line, col = near2d(x, y, x0, y0)
    """
    dx = np.abs(x - x0); dx = dx / dx.max()
    dy = np.abs(y - y0); dy = dy / dy.max()
    dn = dx + dy    
    fn = np.where(dn == dn.min())
    ii = int(fn[0])
    jj  = int(fn[1])
    return ii, jj

def load_bitmap(filename, direc=None):
    """
    Load a bitmap file from the ./icons subdirectory. 
    The filename parameter should not
    contain any path information as this is determined automatically.
    Returns a wx.Bitmap object
    copied from matoplotlib resources
    """

    if not direc:
        basedir = os.path.join(PROJECT_DIR,'icons')
    else:
        basedir = os.path.join(PROJECT_DIR, direc)

    bmpFilename = os.path.normpath(os.path.join(basedir, filename))
    if not os.path.exists(bmpFilename):
        raise IOError('Could not find bitmap file "%s"; dying'%bmpFilename)

    bmp = wx.Bitmap(bmpFilename)
    return bmp


if __name__ == "__main__":
    app = App(False)
    app.MainLoop()

