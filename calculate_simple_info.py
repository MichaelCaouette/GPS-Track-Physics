# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 11:30:11 2022

@author: mcaoue2
"""


import numpy as np

        
class SimplePhysics():
    """
    Calculate simple information  from the raw GPS data. 
    The goal is to not take too much computation time, so it can be quickly
    updated with new data.
    """
    def __init__(self): 
        """
        Initiate the stuff
        """    
        # Nothing to initiate for now. 
         
        
        
    def set_data(self, 
                 t, x, y, z,
                 want_print=False):
        """
        Calculate (and print if want_print=True) simple physical info from 
        the input data
        Watch out for SI unit !       
        
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
        
        # Some physics 
        # Elapsed time
        dt = (self.t[-1]-self.t[0])
        # Mean drop
        drop = self.z[-1]-self.z[0]
        # Mean vertical speed
        vertical_speed = drop/dt    
        # Total distance travelled
        # Determine a correct size for the points
        # Set by the mean distance between points
        list_dx = np.diff(self.x)
        list_dy = np.diff(self.y)
        list_dz = np.diff(self.z)
        list_grd_dist = (list_dx**2 + list_dy**2)**0.5
        sum_grd_dist = np.sum(list_grd_dist)
        
        list_3D_dist = (list_dx**2 + list_dy**2 + list_dz**2)**0.5
        sum_3D_dist = np.sum(list_3D_dist)     
        
        if want_print:
            print()
            print(' --- Physical quantities ---')
            print()
            print('Info on chopped data \nPLEASE VERIFY THE UNITs. COORDINATE SYSTEM MUST BE OKAY:\n')
            print('Duration = %.3f s = %.3f min'%(dt, dt/60))
            print('Vertical drop = %.3f m'%drop)
            print('vertical_speed = %.3f m/s'%vertical_speed)
            print('Travelled distance projected: %.2f m'%sum_grd_dist)
            print('Travelled distance in 3D:     %.2f m'%sum_3D_dist)
            print('Next update: mean projected speed')
            print('Next update: mean 3D speed')
            print('Next update: bird distance')
            print('Next update: bird speed')
        
        
        
        
if __name__ == '__main__':
    # This part of the script runs only when the script file is run alone    
            
    # Test some PGX file
    from tkinter import Tk     # from tkinter import Tk for Python 3.x
    from tkinter.filedialog import askopenfilename        
    
    print('Test some stuff please! Import some data')
    
    
    
    
    
    
    
    
    
    
    
    
            
        