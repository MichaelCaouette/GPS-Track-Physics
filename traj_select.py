# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 13:48:44 2022

Goal:
    GUi to select a portion of the traj
@author: mcaoue2
"""


import numpy as np

import os

from spinmob import egg
import spinmob as sm

from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename

from trajectory_viewer import TrajectoryViewer
from calculate_change_coordinate import Coordinate
from calculate_simple_info import SimplePhysics
from calculate_wind import WindEstimator
from load_gps_data import LoadGPS


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
        self.button_load.set_text('Load Data')
        self.button_load.set_style('background-color: rgb(0, 100, 0);')
        self.connect(self.button_load.signal_clicked, self._button_load_clicked)  
        
        # =============================================================================
        #  Some GUIs for selecting a portion of the trajectory
        # =============================================================================
        self.nbBox_sel_min = egg.gui.NumberBox(value=0,  bounds=(0, None),
                                                int=True,
                                                tip='Lower_bound to select the trajectory.')
        self.nbBox_sel_len = egg.gui.NumberBox(value=10, bounds=(1, None), 
                                                int=True,
                                                tip='Lenght of the selected region of the trajectory.')        
    
        self.place_object(self.nbBox_sel_min, row=0, column=2, alignment=0)
        self.place_object(self.nbBox_sel_len, row=0, column=3, alignment=0)
        self.connect(self.nbBox_sel_min.signal_changed, self._nbBox_sel_changed) 
        self.connect(self.nbBox_sel_len.signal_changed, self._nbBox_sel_changed)    
        
        # =============================================================================
        #   Button to calculate the wind       
        # =============================================================================
        self.button_wind = self.place_object(egg.gui.Button(), 
                                             row=0, column=1, alignment=-1)
        self.button_wind.set_text('Calculate Wind')
        self.button_wind.set_style('background-color: rgb(100, 0, 100);')    
        self.connect(self.button_wind.signal_clicked, self._button_wind_clicked) 
        
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
        self.label_gen_info = egg.gui.Label(text='Waiting to load data')
        self.place_object(self.label_gen_info, row=3, column=0)  
        
        # =============================================================================
        #         # Some stretching/adjustment
        # =============================================================================
        self.set_column_stretch(0)

    def load_data(self, filepath):
        """
        Load the ".GPX" or ".FIT" file. 
        It structures the data.

        Parameters
        ----------
        filepath : String
            Path of the .GPX or .FIT file to load.

        Returns
        -------
        None.

        """
        _debug('TrackSelect:load_data')

        self.filepath = filepath
        self.filename = self.filepath.split(sep='/')[-1]
        
        self.data_load = LoadGPS()
        
        self.file_extension = os.path.splitext(self.filepath)[1]
        
        if   self.file_extension.lower() == '.fit':
            _debug('TrackSelect:load_data: loading .FIT file')
            t, lat, lon, ele = self.data_load.get_FIT(self.filepath )
        elif self.file_extension.lower() == '.gpx':
            _debug('TrackSelect:load_data: loading .GPX file')
            t, lat, lon, ele = self.data_load.get_GPX(self.filepath)
        else:
            raise Exception("Input file must be a .FIT or .GPX file.") 
            
        self.nb_total_pts = len(ele)
        
        
        # =============================================================================
        # Convert the time into seconds
        # =============================================================================

        # Get the initial time       
        self.t_init_dt = t[0]
        self.str_init_t = str(self.t_init_dt) # String telling the starting time
        self.t_init_sec = self.t_init_dt.timestamp()        
        
        # The starting time will be zero sec. 
        self.t_gps = []
        for t_datatime in t:
            self.t_gps.append(t_datatime.timestamp() - self.t_init_sec)
            
        # =============================================================================
        # Make the list numpy array, for math operation      
        # =============================================================================
        ele = np.array(ele)
        lat = np.array(lat)
        lon = np.array(lon)
        self.t_gps   = np.array(self.t_gps )
        
        # =============================================================================
        #         # Remove the outlier
        # =============================================================================
        # Very basic for now. At some point we can have a separeted script with 
        # more fancy methods
        # Remove the common outlier when some elevation are insane
        index = ele<5000 
        self.ele   = ele[index]
        self.lat   = lat[index]
        self.lon   = lon[index]
        self.t_gps = self.t_gps[index] 
        
        self.nb_used_pts = len(self.ele)  
        
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
        self._nbBox_sel_changed(update_plot=False)
       
        if update_plot:
            self._update_plots()

    def _button_load_clicked(self):
        """
        Open a dialog box to select the file to load. 
        And load it. 
        """
        _debug('TrackSelect:_button_load_clicked')

        # we don't want a full GUI, so keep the root window from appearing
        Tk().withdraw() 
        # Get the file directory and name
        # show an "Open" dialog box and return the path to the selected file
        filepath = askopenfilename()         
        self.load_data(filepath)        
               
    def _nbBox_sel_changed(self, value=None, update_plot=True):
        """
        Update the value of the minimum and maximum point for the selected 
        region
        Also update the plots and some info
        We don't use directly the input "v"alue" of the function, but it is required
        such that the signature of the function matches what the widget wants. 

        update_plot:
            Bool, optional
            If True, we also update the plot after updating the attributes. 
            The default is True.
            
        """
        _debug('TrackSelect:_nbBox_sel_changed ')
        
        self.sel_min = self.nbBox_sel_min.get_value()
        self.sel_max = self.nbBox_sel_len.get_value() + self.sel_min
        _debug('self.sel_min, self.sel_max = ', self.sel_min, self.sel_max)
        
        # Update the selected trajectory on certain criteria
        # But is works only in certain condition
        if (self.sel_min < self.sel_max and 
            self.sel_max < self.nb_used_pts):
            # Create the select data
            self.xs_select = self.xs[self.sel_min:self.sel_max]
            self.ys_select = self.ys[self.sel_min:self.sel_max]
            self.zs_select = self.zs[self.sel_min:self.sel_max]
            self.ts_select = self.ts[self.sel_min:self.sel_max]     
            # Compute basic stuff
            self.my_simple = SimplePhysics()
            self.my_simple.set_data(self.ts_select*60, # SI Unit
                                    self.xs_select,
                                    self.ys_select,
                                    self.zs_select,
                                    want_print=True)
            # Update the plot if desired
            if update_plot:
                self._update_plots()
    
    
    def _button_wind_clicked(self):
        """
        Take a portion of the data and show it
        """
        _debug('TrackSelect:_button_wind_clicked ')
        
        self.my_windy = WindEstimator()
        print('TrackSelect:_button_wind_clicked  WATCHOUT THE ASSUMED UNCERTAINY IN X-Y')
        self.my_windy.get_wind_bayes(self.ts_select*60, # SI Unit
                                     self.xs_select,
                                     self.ys_select,
                                     sig_x = 1, sig_y=1)
        
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
            self.my_coor_conv = Coordinate()
            out = self.my_coor_conv.GPS_2_cartesian(self.lon, 
                                                    self.lat, 
                                                    self.ele,
                                                    R_planet = 6371*1e3)
            self.xs, self.ys, self.zs = out
            self.ts = self.t_gps/60
            self.units =('min', 'm', 'm', 'm')
            self.labels=('time', 'x', 'y', 'Elevation')
            
           # Delete when everything works
            # # The elevation remain, but we approximate the latitude and longitude
            # # to span a flat surface, such that we don't bother about the deformation
            # # x and y will be the arc lenght
            # # This will probably be bad at the poles. 
            # # Approximate the planet as a sphere
            # R = 6371*1e3 # m, Radius of the earth 
            # # Formula for the arc lenght as a function of angle
            # self.xs = R*(self.lon-np.min(self.lon))*np.pi/180
            # self.ys = R*(self.lat-np.min(self.lat))*np.pi/180
            # self.zs = self.ele
            # self.ts = self.t_gps/60
            # self.units =('min', 'm', 'm', 'm')
            # self.labels=('time', 'Longitude projection', 'Latitue projection', 'Elevation')
            # print('NEXT IMPROVEMENT: MAKE A REAL CARTESIAN CHANGE AND ROTATE IT WITH THE NORTH POLE')
            
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
        self.traj_viewer.set_trajectory(self.ts_select,
                                        self.xs_select, 
                                        self.ys_select,
                                        self.zs_select,
                                        overwrite=False,
                                        width=2,
                                        units=self.units,
                                        labels=self.labels,
                                        show_pts=self.show_sampling)
          
            
        
if __name__ == '__main__':
    # This part of the script runs only when the script file is run alone    
    _debug_enabled     = True

    sm.settings['dark_theme_qt'] = True
    sm.settings['ignore_warnings']=False
    
    self = TrackSelect()
    self.show()
    

        
        
        
        
        
        
        
        
        