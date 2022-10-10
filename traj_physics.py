# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 11:30:11 2022

@author: mcaoue2
"""


import gpxpy
import numpy as np

from spinmob import egg
import spinmob as sm



from trajectory_viewer import TrajectoryViewer
from wind_estimation import WindEstimator


# Debug stuff.
_debug_enabled     = False

def _debug(*a):
    if _debug_enabled: 
        s = []
        for x in a: s.append(str(x))
        print(', '.join(s))
        
class TrajPhysics(egg.gui.Window):
    """
    GUI
    """
    def __init__(self, name="TrackPhysics", size=[1800,900]): 
        """
        Initiate the GUI
        """    
        _debug('TrajPhysics:__init__')
        _debug('Make each day your masterpiece. – John Wooden')
        
        # Run the basic stuff for the initialization
        egg.gui.Window.__init__(self, title=name, size=size)

        #TODO Add more stuff
        # The viewer of data
        self.traj_viewer = TrajectoryViewer()  
        self.place_object(self.traj_viewer, 
                          row=0, column=0, column_span=4, alignment=0)  
         
        
        
    def set_data(self, 
                 t, x, y, z,
                 labels, units):
        """
        More explanation latter. 
        Watch out for SI unit !
        TODO: remove label and units from the input, because the data are 
        already expected to be in SI unit
        and then use the following input for the trajectory viewer
                                  units=('sec', 'm', 'm', 'm'),
                                  labels=('time', 'Longitude projection',
                                          'Latitue projection', 'Elevation')        
        
        Parameters
        ----------
        t : TYPE
            DESCRIPTION.
        x : TYPE
            DESCRIPTION.
        y : TYPE
            DESCRIPTION.
        z : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.t, self.x, self.y, self.z = t, x, y, z
        
        self.traj_viewer.set_trajectory(self.t, 
                                        self.x, 
                                        self.y, 
                                        self.z,
                                        units=units,
                                        labels=labels)  
        
        # Hack to interact with other GUI
        # Deletet when the coding is better
        # self.t = self.ts_chopped*60        
        # self.x, self.y, self.z = self.xs_chopped, self.xs_chopped, self.zs_chopped
        
        # Some physics (will be in a seperate location latter)
        print()
        print(' --- Physical quantities ---')
        print()
        print('Info on chopped data \nPLEASE VERIFY THE UNITs. COORDINATE SYSTEM MUST BE OKAY:\n')
        # Elapsed time
        dt = (self.t[-1]-self.t[0])
        print('Duration = %.3f s = %.3f min'%(dt, dt/60))
        # Mean drop
        drop = self.z[-1]-self.z[0]
        print('Vertical drop = %.3f m'%drop)
        # Mean vertical speed
        vertical_speed = drop/dt        
        print('vertical_speed = %.3f m/s'%vertical_speed)
        # Total distance travelled
        # Determine a correct size for the points
        # Set by the mean distance between points
        list_dx = np.diff(self.x)
        list_dy = np.diff(self.y)
        list_dz = np.diff(self.z)
        list_grd_dist = (list_dx**2 + list_dy**2)**0.5
        sum_grd_dist = np.sum(list_grd_dist)
        print('Travelled distance projected: %.2f m'%sum_grd_dist)
        list_3D_dist = (list_dx**2 + list_dy**2 + list_dz**2)**0.5
        sum_3D_dist = np.sum(list_3D_dist)   
        print('Travelled distance in 3D:     %.2f m'%sum_3D_dist)
        # Wind speed
        windy = WindEstimator()
        Vaircraft, w, angle =windy.get_wind_rough(self.t, 
                                                  self.x, 
                                                  self.y,
                                                  plot_data=True)
        print('Aircraft speed = %.2f km/h'%(3.6*Vaircraft))
        print('Wind speed     = %.2f km/h'%(3.6*w       ))
        print('Wind direction = %.1f degree'%(angle*180/np.pi))
        print()               
        print('IDEA: Search for the maximum vertical speed (varying the initial and end point.')
        print('IDEA: Estimate the wind. How ? from the turns ?')        
        
        
        
        plt
        
        
        
        
        
        
        
if __name__ == '__main__':
    # This part of the script runs only when the script file is run alone    
            
    # Test some PGX file
    from tkinter import Tk     # from tkinter import Tk for Python 3.x
    from tkinter.filedialog import askopenfilename        
    
    print('Test some stuff please! Import some data')
    
    
    
    
    
    
    
    
    
    
    
    
            
        