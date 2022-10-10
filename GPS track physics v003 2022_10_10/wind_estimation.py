# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:44:52 2022

Goal:
    Estimate the wind direction and amplitude from PGX data alone

@author: mcaoue2
"""



import matplotlib.pyplot as plt
import numpy as np


class WindEstimator():
    """
    Estimate the wind (amplitude and direction) from flight GPS data
    """
    def __init__(self): 
        """
        No input for now
        """    
        
    def get_wind_rough(self, t, x, y, 
                       print_data=False,
                       plot_data=False):
        """
        Cleaner description latter
        
        t, x, y:
            Each are list of float with same lenght
            It represents the time and x,y cartesian coordinate of the flight. 
            SI UNIT please !
            
        return Vaircraft, w, angle
            
        """            
        # =============================================================================
        # Rough estimate of the wind from the extrema of the GPS speed
        # It assumes:
        #   - Constant wind (amplitude and direction) over all the data
        #   - The data contain the maximium and minimum speed 
        #     (ie, there is a point when the aircraft points toward the wind and agains it)
        #   - The speed of the wind is lower that the speed of the aircraft
        # =============================================================================
        
        
        # Hack to be deleted
        t, x, y = self.ts_chopped*60, self.xs_chopped, self.ys_chopped
        print_data=True
        plot_data=True
        
        t_ori = t
        
        # Get the speed from naive differentiation (WATCH OUT THE NOISE !)
        # Differential element
        dt = np.diff(t_ori)
        dx = np.diff(x) 
        dy = np.diff(y)
        # The reference time in the differenciated world is shifted with the original
        # and we pop-up the last time
        t_diff = 0.5*dt + t_ori[:-1]
        # Speed
        vx = dx /dt
        vy = dy /dt
        # Magnitude of the speed
        V = np.sqrt(vx**2 + vy**2)
        # Find the extrem speed
        Vmax, Vmin = np.max(V), np.min(V)
        # Estimate the amplitude of the wind and the motor
        w         = 0.5*(Vmax - Vmin)
        Vaircraft = 0.5*(Vmax + Vmin) 
        # Find the direction of the wind. 
        #TODO Improve the accurance by taking the list of speed that ranges within a 
        # certain tolerance near the max. The tolerance should be set by the uncertainty
        # in the speed (coming from uncertainty in the GPS data)
        
        # It is the direction when the speed is maximat
        index = np.argwhere(V==Vmax) # Note that this will generate an array of all occurence
        vx_max = vx[index]
        vy_max = vy[index]
        # My vaforite way to get the angle of a 2D vector
        angle = np.angle(vx_max + 1j*vy_max)
        # Average, in case we have multiple occurence of the max
        angle = np.mean(angle)
        
        if print_data:
            print('Aircraft speed = %.2f km/h'%(3.6*Vaircraft))
            print('Wind speed     = %.2f km/h'%(3.6*w       ))
            print('Wind direction = %.1f degree'%(angle*180/np.pi))

        if plot_data:
            plt.figure(tight_layout=True)
            plt.plot(t_diff/60, 3.6*V)
            plt.plot([t_diff[0]/60, t_diff[-1]/60], 
                 [3.6*Vaircraft, 3.6*Vaircraft], 
                 'k--', 
                 label='Aircraft speed')
            plt.legend()
            plt.xlabel('Time (min)')
            plt.ylabel('Aircraft Speed (km/h)')
            plt.title('Wind Estimation: %.2f km/h, pointing at %.1f degree'%(w, angle*180/np.pi))
           
        return Vaircraft, w, angle
        
    def plot_velocities(self):
        # Show the data
        # Check the data
        plt.figure(tight_layout=True)
        plt.plot(self.t_diff/60, 3.6*self.V)
        plt.plot([self.t_diff[0]/60, self.t_diff[-1]/60], 
                 [3.6*self.Vaircraft, 3.6*self.Vaircraft], 
                 'k--', 
                 label='Aircraft speed')
        plt.legend()
        plt.xlabel('Time (min)')
        plt.ylabel('Aircraft Speed (km/h)')
        plt.title('Wind Estimation: %.2f km/h'+
                  ', pointing at %.1f degree'%(self.w, self.angle*180/np.pi))
        
if __name__ == '__main__':
    # This part of the script runs only when the script file is run alone    

    # Temporary, to used the data from the other GUI
    #MAKE SURE OF THE UNIT
    # my_wind = WindEstimator()
    # my_wind.get_wind_rough(self.ts_chopped*60, 
    #                        self.xs_chopped, 
    #                        self.ys_chopped)

    print('Please test somehting')
    
    
    
    
    # To be included in scripts
    plt.figure(tight_layout=True)
    plt.plot(3.6*vx, 3.6*vy, '.')
    plt.xlabel('Vx (km/h)')
    plt.ylabel('Vy (km/h)')
    plt.title('Flight velocities')
    
    
    
    
    
    
    
    