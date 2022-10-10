# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 13:48:44 2022

Goal:
    GUi to select a portion of the traj
@author: mcaoue2
"""

import gpxpy
import numpy as np

from spinmob import egg
import spinmob as sm

from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename

from trajectory_viewer import TrajectoryViewer

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
    def __init__(self, name="TrackSelect", size=[1500,900]): 
        """
        Initiate the GUI
        """    
        _debug('TrackSelect:__init__')
        _debug('Make each day your masterpiece. – John Wooden')
        
        # Run the basic stuff for the initialization
        egg.gui.Window.__init__(self, title=name, size=size)
        
        # Button to select the data
        # For load the image
        self.button_load = self.place_object(egg.gui.Button(), 
                                             row=0, column=0,
                                             alignment=+1)
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
        self.button_chop = self.place_object(egg.gui.Button(), 
                                             row=0, column=1, alignment=+1)
        self.button_chop.set_text('Chop trajectory')
        self.button_chop.set_style('background-color: rgb(100, 0, 100);')        
        self.place_object(self.nbBox_chop_min, row=0, column=2, alignment=0)
        self.place_object(self.nbBox_chop_max, row=0, column=3, alignment=0)
        self.connect(self.button_chop.signal_clicked, self._button_chop_clicked) 
        self.connect(self.nbBox_chop_min.signal_changed, self._nbBox_chop_changed) 
        self.connect(self.nbBox_chop_max.signal_changed, self._nbBox_chop_changed)          
                
        # Trajectory plotter
        
        # Plots for showing the track in 2D
        self.traj_viewer = TrajectoryViewer()  
        self.place_object(self.traj_viewer, 
                          row=1, column=0, column_span=4, alignment=0)          
        
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
        
        # Choose a coordinate system for the viewer
        self.set_coordinate_syst(self, 'cartesian')
        
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
                print('Chop chop ! ',self.chop_min, self.chop_max)
    
                
        
    
    def _button_chop_clicked(self):
        """
        Chopped a portion of the data
        """
        self.lat_chopped = self.lat[self.chop_min:self.chop_max]
        self.lon_chopped = self.lon[self.chop_min:self.chop_max]
        self.ele_chopped = self.ele[self.chop_min:self.chop_max]
        self.tim_chopped = self.t_gps  [self.chop_min:self.chop_max]
        
        # 3D trajectory        
        self.ax_3D_chop.clear()
        self.ax_3D_chop.set_title('Chopped trajectory\n') 
        # The line for the trajectory
        self.ax_3D_chop.plot(self.lat_chopped, 
                             self.lon_chopped , 
                             self.ele_chopped)
        # Initial position
        self.ax_3D_chop.plot([self.lat_chopped[0]],
                             [self.lon_chopped[0]], 
                             [self.ele_chopped[0]], 'or')  
        self.fig_3D_chop.show()
        
        # Some physics (will be in a seperate location latter)
        print()
        print('Info on chopped data:')
        # Elapsed time
        dt = self.tim_chopped[-1]-self.tim_chopped[0]
        print('dt = %.3f s'%dt)
        # Mean drop
        drop = self.ele_chopped[-1]-self.ele_chopped[0]
        print('drop = %.3f m'%drop)
        # Mean vertical speed
        vertical_speed = drop/dt        
        print('vertical_speed = %.3f m/s'%vertical_speed)
        print()
                
        
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
        self.ele, self.lat, self.lon, self.t_gps = ele, lat, lon, t 
        self.nb_used_pts = len(self.ele)
        
    
    def set_coordinate_syst(self, coord_sys):
        """
        Assign a coordinate system to be used and shown. 
        This sets the unit
        Note that the track points must be loaded. 


        Parameters
        ----------
        coord_sys : String
            if 'original':
                Take the original number as given by the GPS.
            elif 'rescaled':
                Scale all the number from 0 to 1
            elif 'cartesian':
                Cartesian projection on the surface of the sphere (ie planet). 
                Flat earthist love this one. 
                However, wanring for deformation near the pole. 

        Returns
        -------
        None.

        """
        _debug('TrackSelect:set_coordinate_syst')
        
        
        self.coord_sys = coord_sys 
        
        # We will define the parameters relevant for the viewers and more
        if self.coord_sys == 'original':
            # Take the number as given by the GPS
            self.ts = self.t_gps/60
            self.xs, self.ys, self.zs =  self.lon, self.lat, self.ele
            self.units = ('min', 'deg', 'deg', 'm')
            self.labels=('time', 'Longitude', 'Latitue', 'Elevation')
            self.rescale=False
            
        elif self.coord_sys == 'rescaled':
            # Scale all the number from 0 to 1
            self.ts = self.t_gps/60
            self.xs, self.ys, self.zs =  self.lon, self.lat, self.ele
            self.units = ('min', 'deg', 'deg', 'm')
            self.labels=('time', 'Longitude', 'Latitue', 'Elevation')
            self.rescale=False            
            
            self.traj_viewer.set_trajectory(self.t_gps/60, self.lat, self.lon, self.ele,
                                            rescale=True,
                                            units=('min', 'rescaled', 'rescaled', 'rescaled'))
        elif self.coord_sys == 'cartesian':
            # Convert the number into cartesian postion
            # The elevation remain, but we approximate the latitude and longitude
            # to span a flat surface, such that we don't bother about the deformation
            # x and y will be the arc lenght
            # This will probably be bad at the poles. 
            # Approximate the planet as a sphere
            R = 6371*1e3 # m, Radius of the earth 
            # Formula for the arc lenght as a function of angle
            self.x = R*(self.lon-np.min(self.lon))*np.pi/180
            self.y = R*(self.lat-np.min(self.lat))*np.pi/180
            self.traj_viewer.set_trajectory(self.t_gps/60, self.x, self.y, self.ele,
                                            units=('min', 'm', 'm', 'm'),
                                            labels=('time', 'Longitude projection', 'Latitue projection', 'Elevation'))  
            
        return
    
    def _update_plots(self):
        """
        
        Returns
        -------
        None.

        """
        _debug('TrackSelect:_update_plots')
        
        
        
        
        # I AM HERE
        # PLEASE ADD THE UNIT
        # AND IMPROVE THE TRAJ VIEWER FOR THE SCALING
        # TODO: have a widget to select the projection to use
        scaling='cartesian'
        if scaling == 'original':
            # Take the number as given by the GPS
            self.traj_viewer.set_trajectory(self.t_gps/60, self.lon, self.lat, self.ele,
                                            units=('min', 'deg', 'deg', 'm'),
                                            labels=('time', 'Longitude', 'Latitue', 'Elevation'))
        elif scaling == 'rescaled':
            # Scale all the number from 0 to 1
            self.traj_viewer.set_trajectory(self.t_gps/60, self.lat, self.lon, self.ele,
                                            rescale=True,
                                            units=('min', 'rescaled', 'rescaled', 'rescaled'))
        elif scaling == 'cartesian':
            # Convert the number into cartesian postion
            # The elevation remain, but we approximate the latitude and longitude
            # to span a flat surface, such that we don't bother about the deformation
            # x and y will be the arc lenght
            # This will probably be bad at the poles. 
            # Approximate the planet as a sphere
            R = 6371*1e3 # m, Radius of the earth 
            # Formula for the arc lenght as a function of angle
            self.x = R*(self.lon-np.min(self.lon))*np.pi/180
            self.y = R*(self.lat-np.min(self.lat))*np.pi/180
            self.traj_viewer.set_trajectory(self.t_gps/60, self.x, self.y, self.ele,
                                            units=('min', 'm', 'm', 'm'),
                                            labels=('time', 'Longitude projection', 'Latitue projection', 'Elevation'))            
            
        
if __name__ == '__main__':
    # This part of the script runs only when the script file is run alone    
    _debug_enabled     = True

    sm.settings['dark_theme_qt'] = True
    sm.settings['ignore_warnings']=True
    
    self = TrackSelect()

    self.show()
    
    
    
    
    
        
        
        
        
        
        
        
        
        