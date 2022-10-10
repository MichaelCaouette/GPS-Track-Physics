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
from traj_physics import TrajPhysics

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
    def __init__(self, name="TrackSelect", size=[1800,900]): 
        """
        Initiate the GUI
        """    
        _debug('TrackSelect:__init__')
        _debug('Make each day your masterpiece. – John Wooden')
        
        # Run the basic stuff for the initialization
        egg.gui.Window.__init__(self, title=name, size=size)
        
        # =============================================================================
        #         # Button to select the data
        # =============================================================================
        # For load the image
        self.button_load = self.place_object(egg.gui.Button(), 
                                             row=0, column=0,
                                             alignment=-1)
        self.button_load.set_text('Load GPX')
        self.button_load.set_style('background-color: rgb(0, 100, 0);')
        self.connect(self.button_load.signal_clicked, self._button_load_clicked)  
        
        # =============================================================================
        #         # Some GUI for choping the trajectory
        # =============================================================================
        self.nbBox_chop_min = egg.gui.NumberBox(value=0,  bounds=(0, None),
                                                int=True,
                                                tip='Lower_bound to chop the trajectory.')
        self.nbBox_chop_len = egg.gui.NumberBox(value=10, bounds=(1, None), 
                                                int=True,
                                                tip='Lenght of the chopping region to chop the trajectory.')        
        self.button_chop = self.place_object(egg.gui.Button(), 
                                             row=0, column=1, alignment=-1)
        self.button_chop.set_text('Chop trajectory')
        self.button_chop.set_style('background-color: rgb(100, 0, 100);')        
        self.place_object(self.nbBox_chop_min, row=0, column=2, alignment=0)
        self.place_object(self.nbBox_chop_len, row=0, column=3, alignment=0)
        self.connect(self.button_chop.signal_clicked, self._button_chop_clicked) 
        self.connect(self.nbBox_chop_min.signal_changed, self._nbBox_chop_changed) 
        self.connect(self.nbBox_chop_len.signal_changed, self._nbBox_chop_changed)          
                
        # =============================================================================
        #         # Trajectory plotter
        # =============================================================================
        self.traj_viewer = TrajectoryViewer()  
        self.place_object(self.traj_viewer, 
                          row=1, column=0, column_span=4, alignment=0)          

        # =============================================================================
        #         # A treeDic to control more parameters
        # =============================================================================
        self.treeDic_settings  = egg.gui.TreeDictionary(autosettings_path='setting_traj')
        self.place_object(self.treeDic_settings, 
                          row=0, column=6, row_span=2, alignment=0)   
        # Choose the coordinate system to display
        self.list_scaling = ['original', 'rescaled', 'cartesian']
        self.treeDic_settings.add_parameter('coord_system', 0, 
                                            type='list', values=self.list_scaling,
                                            tip='Coordinate system to use.\nRescaled is rescaling all the coord from 0 to 1.')   
        self.treeDic_settings.connect_signal_changed('coord_system', self._coord_system_changed)
        # Weither or not to show the sampling points
        self.treeDic_settings.add_parameter('show_sampling', True, 
                                            type='bool',
                                            tip='Weither or not to show the sampling points.')   
        self.treeDic_settings.connect_signal_changed('show_sampling', self._show_sampling_changed)
        
        # =============================================================================
        #         #Label
        # =============================================================================
        self.label_gen_info = egg.gui.Label(text='Hello', autosettings_path='label_general_info')
        self.place_object(self.label_gen_info, row=3, column=0)  
        
        # =============================================================================
        #         # Some stretching/adjustment
        # =============================================================================
        self.set_column_stretch(0)

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
        
        # =============================================================================
        # Update the attributes of the GUI         
        # =============================================================================
        # This update the attributes with the widgets
        self._update_all_attributes(update_plot=True) 
        # Update the printed info
        self.label_gen_info.set_text('The file is '+ self.filename +
                                     '   Flight started on '+self.str_init_t +
                                     '\n %d points in total (%d after outliers removale and smoothing)'%(self.nb_total_pts, self.nb_used_pts))


    def _update_all_attributes(self, update_plot=False):
        """
        Update all the attributre with the values on the widgets. 
        Making the viewer ready to plot. 
        
        It is assumed that the GPS data are already loaded and structured. 
        
        update_plot:
            Bool, optional
            If True, we also update the plot after updating the attributes. 
            The default is False

        Returns
        -------
        None.

        """    
        _debug('TrackSelect:_update_all_attributes ')          
        
        # The coordinate system
        # By setting update_plot to False, we just update the attributes
        self._coord_system_changed(update_plot=False)
        
        # The sampling pts
        # By setting update_plot to False, we just update the attributes
        self._show_sampling_changed(update_plot=False)
        
        # Update the chopper 
        # By setting update_plot to False, we just update the attributes
        self._nbBox_chop_changed(update_plot=False)
       
        if update_plot:
            self._update_plots()

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
               
    def _nbBox_chop_changed(self, value=None, update_plot=True):
        """
        Update the value of the minimum and maximum point for the chopping. 
        Also update the plots.
        We don't use directly the input value of the function, but it is required
        such that the signature of the function matches what the widget wants. 

        update_plot:
            Bool, optional
            If True, we also update the plot after updating the attributes. 
            The default is True.
            
        """
        _debug('TrackSelect:_nbBox_chop_changed ')
        
        self.chop_min = self.nbBox_chop_min.get_value()
        self.chop_max = self.nbBox_chop_len.get_value() + self.chop_min
        _debug('self.chop_min, self.chop_max = ', self.chop_min, self.chop_max)
        
        # Update the plot (should show it)
        # But is works only in certain condition
        if self.chop_min < self.chop_max:
            if self.chop_max < self.nb_used_pts:
                # Create the chopped data
                self.xs_chopped = self.xs[self.chop_min:self.chop_max]
                self.ys_chopped = self.ys[self.chop_min:self.chop_max]
                self.zs_chopped = self.zs[self.chop_min:self.chop_max]
                self.ts_chopped = self.ts[self.chop_min:self.chop_max]                
                # Update the plot if desired
                if update_plot:
                    self._update_plots()
    
    
    def _button_chop_clicked(self):
        """
        Take a portion of the data and show it
        """
        _debug('TrackSelect:_button_chop_clicked ')
        


        # Send the chopped data into the physics calculator
        # MUST BE SI UNIT!
        self.gui_physics = TrajPhysics().show() # Don't forget to show !
        self.gui_physics.set_data(self.ts_chopped*60, # Watch out the unit !
                                  self.xs_chopped, 
                                  self.ys_chopped, 
                                  self.zs_chopped,
                                  units=('sec', 'm', 'm', 'm'),
                                  labels=('time', 'Longitude projection',
                                          'Latitue projection', 'Elevation'))
        
        
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
        self.str_init_t = str(self.t_init_dt) # String telling the starting time
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
    
    def _coord_system_changed(self, update_plot=True):
        """
        Called when the combo boxe for the coordinate system change. 

        update_plot:
            Bool, optional
            If True, we also update the plot after updating the attributes. 
            The default is True.

        """
        self.coord_system = self.treeDic_settings['coord_system']
        _debug('TrackSelect._coord_system_changed --> val= ', self.coord_system)
        
        self.set_coordinate_syst(self.coord_system)
        
        if update_plot:
            self._update_plots()
        
        return
        
    def _show_sampling_changed(self, update_plot=True):
        """
        Called when the checkmark for the show sampling change. 
        
        update_plot:
            Bool, optional
            If True, we also update the plot after updating the attributes. 
            The default is True.        
        
        """
        self.show_sampling = self.treeDic_settings['show_sampling']
        _debug('TrackSelect._show_sampling_changed --> val= ', self.show_sampling)

        # That its. Now the sampling should be updated on the plot when we 
        # update the plot. 
        if update_plot:
            self._update_plots()
        
        return
        
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
        _debug('TrackSelect.set_coordinate_syst')
        
        
        self.coord_sys = coord_sys 
        
        # We will define the parameters relevant for the viewers and more
        if self.coord_sys == 'original':
            # Take the number as given by the GPS
            self.ts = self.t_gps/60
            self.xs, self.ys, self.zs =  self.lon, self.lat, self.ele
            self.units = ('min', 'deg', 'deg', 'm')
            self.labels=('time', 'Longitude', 'Latitue', 'Elevation')
            
        elif self.coord_sys == 'rescaled':
            # Scale all the number from 0 to 1
            # Do the same rescale for each of the 4 variable
            pars_ori = [self.t_gps, self.lon, self.lat, self.ele]
            self.pars_new = ['', '', '', '']
            for i in range(4):
                var = pars_ori[i]
                self.pars_new[i] = (var-np.min(var)) / (np.max(var) - np.min(var))
                
            self.ts, self.xs, self.ys, self.zs = self.pars_new           
            self.units = ('Unitless', 'Unitless',  'Unitless', 'Unitless')
            self.labels=('t', 'x', 'y', 'z')
            
        elif self.coord_sys == 'cartesian':
            # Convert the number into cartesian postion
            # The elevation remain, but we approximate the latitude and longitude
            # to span a flat surface, such that we don't bother about the deformation
            # x and y will be the arc lenght
            # This will probably be bad at the poles. 
            # Approximate the planet as a sphere
            R = 6371*1e3 # m, Radius of the earth 
            # Formula for the arc lenght as a function of angle
            self.xs = R*(self.lon-np.min(self.lon))*np.pi/180
            self.ys = R*(self.lat-np.min(self.lat))*np.pi/180
            self.zs = self.ele
            self.ts = self.t_gps/60
            self.units =('min', 'm', 'm', 'm')
            self.labels=('time', 'Longitude projection', 'Latitue projection', 'Elevation')
        else:
            print('Error TrackSelect.set_coordinate_syst --> coord_sys=',coord_sys, ' is not a proper input.')
            
        return
    
    def _update_plots(self):
        """
        
        Returns
        -------
        None.

        """
        _debug('TrackSelect:_update_plots')
        
        # Update the whole trajectory
        self.traj_viewer.set_trajectory(self.ts, self.xs, self.ys, self.zs,
                                        units=self.units,
                                        labels=self.labels)
        # Add the region for choping 
        self.traj_viewer.set_trajectory(self.ts_chopped,
                                        self.xs_chopped, 
                                        self.ys_chopped,
                                        self.zs_chopped,
                                        overwrite=False,
                                        width=2,
                                        units=self.units,
                                        labels=self.labels,
                                        show_pts=self.show_sampling)
          
            
        
if __name__ == '__main__':
    # This part of the script runs only when the script file is run alone    
    _debug_enabled     = True

    sm.settings['dark_theme_qt'] = True
    sm.settings['ignore_warnings']=True
    
    self = TrackSelect()
    self.show()
    
    
    
    
    
        
        
        
        
        
        
        
        
        