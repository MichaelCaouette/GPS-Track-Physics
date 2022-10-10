# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 13:48:44 2022

Goal:
    Starting point GUI to select and analysis sets of points
@author: mcaoue2
"""

import gpxpy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from spinmob import egg
import spinmob as sm

from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename

import traceback
_p = traceback.print_last #Very usefull command to use for getting the last-not-printed error


# Debug stuff.
_debug_enabled     = False

def _debug(*a):
    if _debug_enabled: 
        s = []
        for x in a: s.append(str(x))
        print(', '.join(s))
 
        
class TrackSelect(egg.gui.Window):
    """
    GUI
    """
    def __init__(self, name="TrackSelect", size=[1000,500]): 
        """
        Initiate the GUI
        """    
        _debug('TrackSelect:__init__')
        _debug('Make each day your masterpiece. – John Wooden')
        
        # Run the basic stuff for the initialization
        egg.gui.Window.__init__(self, title=name, size=size)
        
        # Button to select the data
        # For load the image
        self.button_load = self.place_object(egg.gui.Button(), row=0, column=0)
        self.button_load.set_text('Load GPX')
        self.button_load.set_style('background-color: rgb(0, 100, 0);')
        self.connect(self.button_load.signal_clicked, self._button_load_clicked)  
        
        # For choping the trajectory
        self.nbBox_chop_min = egg.gui.NumberBox(value=0,  bounds=(0, None),
                                                int=True,
                                                tip='Lower_bound to chop the trajectory.')
        self.nbBox_chop_max = egg.gui.NumberBox(value=10, bounds=(1, None), 
                                                int=True,
                                                tip='Higher_bound to chop the trajectory.')        
        self.button_chop = self.place_object(egg.gui.Button(), row=0, column=1)
        self.button_chop.set_text('Chop trajectory')
        self.button_chop.set_style('background-color: rgb(100, 0, 100);')        
        self.place_object(self.nbBox_chop_min, row=0, column=2, alignment=1)
        self.place_object(self.nbBox_chop_max, row=0, column=3, alignment=1)
        self.connect(self.button_chop.signal_clicked, self._button_chop_clicked) 
        self.connect(self.nbBox_chop_min.signal_changed, self._nbBox_chop_changed) 
        self.connect(self.nbBox_chop_max.signal_changed, self._nbBox_chop_changed)          
                
        # Plots for showing the track in 2D
        # PLEASE DON"T USE DATABOX. WE JUST NEED THE SIMPLE PLOT ! CHECK DIAMOND-LAB-Stuff
        self.databoxplot_xypath = egg.gui.DataboxPlot(autosettings_path='plot_xypath')
        self.databoxplot_tz     = egg.gui.DataboxPlot(autosettings_path='plot_tz')   
        self.place_object(self.databoxplot_xypath, 
                          row=1, column=0,column_span=4, alignment=0)
        self.place_object(self.databoxplot_tz    , 
                          row=2, column=0,column_span=4, alignment=0)   
        
        # For 3D, could use this:
        # egg.pyqtgraph.opengl.GLLinePlotItem()
        
        #Label
        self.label_gen_info = egg.gui.Label(text='Hello', autosettings_path='label_general_info')
        self.place_object(self.label_gen_info, row=3, column=0)  
        
        
        
        
        
        
    def load_GPX(self, gpx_filepath):
        """
        Load the ".GPX" file. 
        It structures the data

        Parameters
        ----------
        gpx_filepath : String
            Path of the GPX file to load.

        Returns
        -------
        None.

        """
        _debug('TrackSelect:load_GPX')
        
        self.gpx_filepath = gpx_filepath
        self.filename = self.gpx_filepath.split(sep='/')[-1]
        
        # Get the t, x, y, z (AKA motion!)
        with open(self.gpx_filepath) as fh:
            self.gpx_file = gpxpy.parse(fh)        
            
        _debug("File has {} track(s).".format(len(self.gpx_file.tracks)))   
        _debug("Track has {} segment(s).".format(len(self.gpx_file.tracks[0].segments))) 
        # So we take the only track and segment. 
        self.segment = self.gpx_file.tracks[0].segments[0]
        _debug('The segment has {} points'.format(len(self.segment.points)))        
        
        # Get some info
        self.nb_total_pts = len(self.segment.points)
        
        # Re-Structure the data points
        self._structure_data(self.segment.points)
        
        # Show the data
        # Initiate some values and plots
        self._nbBox_chop_changed('Blaaaaa')
        self._update_plots()
        
        # Update the info
        self.label_gen_info.set_text('The file is '+ self.filename +
                                     '\n %d points in total (%d after outliers removale and smoothing)'%(self.nb_total_pts, self.nb_used_pts))

    def _button_load_clicked(self):
        """
        Load an image. 
        """
        _debug('TrackSelect:_button_load_clicked')
        # =============================================================================
        # Load the image      
        # =============================================================================
        # we don't want a full GUI, so keep the root window from appearing
        Tk().withdraw() 
        # Get the file directory and name
        filepath = askopenfilename() # show an "Open" dialog box and return the path to the selected file
        self.load_GPX(filepath)        
               
    def _nbBox_chop_changed(self, value):
        """
        Update the value of the minimum and maximum point for the chopping. 
        Also update the plots.
        We don't use directly the input value
        
        """
        _debug('TrackSelect:_nbBox_chop_changed ')
        
        self.chop_min = self.nbBox_chop_min.get_value()
        self.chop_max = self.nbBox_chop_max.get_value()
        _debug('self.chop_min, self.chop_max = ', self.chop_min, self.chop_max)
        
        # Update the plot (should show it)
        # But is works only in certain condition
        if self.chop_min < self.chop_max:
            if self.chop_max < self.nb_used_pts:
                self._update_plots()
        
    
    def _button_chop_clicked(self):
        """
        Chopped a portion of the data
        """
        self.lat_chopped = self.lat[self.chop_min:self.chop_max]
        self.lon_chopped = self.lon[self.chop_min:self.chop_max]
        self.ele_chopped = self.ele[self.chop_min:self.chop_max]
        
        # 3D trajectory        
        self.fig_3D_chop = plt.figure()
        self.ax_3D_chop= self.fig_3D_chop .add_subplot(111, projection='3d')
        self.ax_3D_chop.set_title(self.filename + '  Starting at ' + str(self.t_init_dt)) 
        # The line for the trajectory
        self.ax_3D_chop.plot(self.lat_chopped, 
                            self.lon_chopped , 
                            self.ele_chopped)
        # Initial position
        self.ax_3D_chop.plot(self.lat_chopped[0], self.lon_chopped[0], self.ele_chopped[0], 'or')  
        
        
    def _structure_data(self, list_trackpts, outliers=True):
        """
        Structure the data. 

        Parameters
        ----------
        list_trackpts : List of "gpxpy.gpx.GPXTrackPoint"
            Each "gpxpy.gpx.GPXTrackPoint" is a data structure on which we 
            extract the time and x, y, z position. 
        
        outliers: Bool
            TOBE EXPLAINED LATTER
        
        Returns
        -------
        None.

        """
        _debug('TrackSelect:_structure_data')
        
        self.list_trackpts = list_trackpts
        
        # The structure is that each point has the t, x, y, z information.
        # If we want an array for each degree of freedom, we have to restructure it        
        self.t_init_dt = self.list_trackpts[0].time
        self.t_init_sec = self.t_init_dt.timestamp()        
        
        # Prepare the array of info
        if outliers:
            # Method to remove outliers
            ele_max = 5000 # Maximum elevation
            ele, lat, lon, t = [], [], [], []
            for p in self.list_trackpts:
                if p.elevation <= ele_max:
                    ele.append(p.elevation)
                    lat.append(p.latitude)
                    lon.append(p.longitude)
                    t.append(p.time.timestamp() - self.t_init_sec)
            ele = np.array(ele)
            lat = np.array(lat)
            lon = np.array(lon)
            t   = np.array(t  )              
        else:          
            # Take all the track points, outliers or not. 
            # Pre-allocation
            ele, lat, lon, t = np.zeros((4,len(self.list_trackpts) ))
            for i,p, in enumerate(self.list_trackpts):
                ele[i] = p.elevation
                lat[i] = p.latitude
                lon[i] = p.longitude
                t[i]   = p.time.timestamp() - self.t_init_sec
        
        # Save the info in the object
        self.ele, self.lat, self.lon, self.t = ele, lat, lon, t 
        self.nb_used_pts = len(self.ele)
        
    def _update_plots(self):
        """
        
        Returns
        -------
        None.

        """
        _debug('TrackSelect:_update_plots')
        
        # 2D plot
        self.databoxplot_xypath.append_row([self.lat, self.lon],
                                           ['Latitude','Longitude']).plot()     
        self.databoxplot_tz.append_row([self.t/60, self.ele],
                                       ['Time (min)','Elevation']).plot()               
            
            
        # 3D trajectory        
        self.fig_3D_all = plt.figure()
        self.ax_3D_all = self.fig_3D_all .add_subplot(111, projection='3d')
        self.ax_3D_all.set_title(self.filename + '  Starting at ' + str(self.t_init_dt)) 
        # The line for the trajectory
        self.ax_3D_all.plot(self.lat, self.lon, self.ele)
        # Initial position
        self.ax_3D_all.plot(self.lat[0], self.lon[0], self.ele[0], 'or')
        # Boundary for the chopping
        self.ax_3D_all.plot(self.lat[self.chop_min], 
                            self.lon[self.chop_min], 
                            self.ele[self.chop_min], 'og')
        self.ax_3D_all.plot(self.lat[self.chop_max], 
                            self.lon[self.chop_max], 
                            self.ele[self.chop_max], 'ob')
        
        
        
if __name__ == '__main__':
    # This part of the script runs only when the script file is run alone    
    _debug_enabled     = True

    sm.settings['dark_theme_qt'] = True
    # path = "C://Users/client//Desktop//GPS track physics"
    path = "C://Users//mcaoue2\Desktop\GPS track physics 2022_09_06"
    filename = "TRK4.GPX"    
    
    filepath = path + '/' + filename
    
    self = TrackSelect()
    self.load_GPX(filepath)
    self.show()
    
    
    
    
    
        
        
        
        
        
        
        
        
        