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

    def get_wind_bayes(self, t, x, y, sig_x, sig_y):
        """
        Sophisticated inference of the wind. Based on Bayes analysis, with 
        a model for the aircraft. 
        See my notebook for more info on the bayesian analysis method
        (hopefully I will update it soon)
        
        TODO: Comment and explain the algorithm. Expecially that I iterate 
        in order to compute faster with a small size of the 3D domain. 
        Which is probably not obvious. 

        TODO: Explain the input parameters
        
        Parameters
        ----------
        t : TYPE
            DESCRIPTION.
        x : TYPE
            DESCRIPTION.
        y : TYPE
            DESCRIPTION.
        print_data : TYPE, optional
            DESCRIPTION. The default is False.
        plot_data : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        # Code to be cleaned later. For now it is the copy paste of a script 
        # that I wrote and tested with fake data
            
        # =============================================================================
        # Prepare the data        
        # =============================================================================
        # to help debugging
        self.t = t
        self.x = x
        self.y = y
        self.sig_x = sig_x
        self.sig_y = sig_y
        
        
        # Get the speed from naive differentiation (WATCH OUT THE NOISE !)
        # Differential element
        dt = np.diff(self.t)
        dx = np.diff(self.x) 
        dy = np.diff(self.y)
        
        # Speed
        meas_vx = dx /dt
        meas_vy = dy /dt 
        # Uncertainty in the speed, assuming not uncertainty in the time and 
        # v = (x2-x1)/dt, x2 and x1 have the same uncertainy and are independent
        # Are they really independent ? I hope so !
        sig_vx = np.sqrt(2)*self.sig_x/dt
        sig_vy = np.sqrt(2)*self.sig_y/dt

        N_pts = len(meas_vx)
        
        # =============================================================================
        # Useful function        
        # =============================================================================
        # Model for the data
        def model_vx(wx, Va, angle):
            return wx + Va*np.cos(angle)
        def model_vy(wy, Va, angle):
            return wy + Va*np.sin(angle)

        # Compute mean and std
        def get_stat(x, P):
            # Input are 1D, same size
            # P is assumed to be normalized
            mean = np.sum(x   *P)
            m2   = np.sum(x**2*P)
            std  = np.sqrt(m2 - mean**2)
            return mean, std
        
        def compute_pdf(list_wx, 
                        list_wy, 
                        list_va):
            
            # The domain of the parameters that we try to infer
            Mwx, Mwy, Mva = np.meshgrid(list_wx, list_wy, list_va)
            
            # Logarithm of the posterior pdf
            L_post = np.zeros(np.shape(Mwx))
            
            N_phi_sum = 20
            phis_sum = np.linspace(0, 2*np.pi, N_phi_sum)
            
            # To avoid numerical problem with logs, we will consider zero to be the smallest number
            smallest_nb = np.nextafter(0, 1)
            
            for i in range(N_pts):
                if i%100 == 0:
                    print('Calculating the Like-lihood of pts %d / %d'%(i, N_pts))
                # Marginalize phi, which means to Integrate on phis
                expo = 0
                for i_phi in range(N_phi_sum):
                    chi_x = (meas_vx[i] - model_vx(Mwx, Mva, phis_sum[i_phi]))**2 / (2*sig_vx[i]**2 )
                    chi_y = (meas_vy[i] - model_vy(Mwy, Mva, phis_sum[i_phi]))**2 / (2*sig_vy[i]**2 )
                    # Add that
                    expo += np.exp(-chi_x - chi_y)
                
                # Don't care about overall factors
                # Replace the zeros by the smallest existing value. This is because of the numberical round-off
                expo = np.where(expo==0, smallest_nb, expo)
                L = np.log(expo) # Logarithm of the likelihood
                L_post += L
                # Reshift the log, to avoid infinities
                L_post -= np.max(L_post)
            
            post = np.exp(L_post)
            # Normalize
            post /= np.sum(post)
            print('Done !')
            
            # Now that we have the post, we can marginalize it on the relevant variable
            # The order of how we sum is choosen such tha it matches.
            post_wx = np.sum(np.sum(post, axis=2), axis=0)
            post_wy = np.sum(np.sum(post, axis=2), axis=1)
            post_va = np.sum(np.sum(post, axis=1), axis=0)
        
            mean_wx, std_wx = get_stat(list_wx, post_wx)
            mean_wy, std_wy = get_stat(list_wy, post_wy)
            mean_va, std_va = get_stat(list_va, post_va)
        
            # Check what we get
            print('mean_wx, std_wx = ', mean_wx*3.6, std_wx*3.6, ' km/h')
            print('mean_wy, std_wy = ', mean_wy*3.6, std_wy*3.6, ' km/h')
            print('mean_va, std_va = ', mean_va*3.6, std_va*3.6, ' km/h')
        
            return (post_wx, post_wy, post_va,
                    mean_wx, std_wx,
                    mean_wy, std_wy,
                    mean_va, std_va)
        
        
        # =============================================================================
        # Iterate with a small prior, to be faster
        # =============================================================================
        N_prior = 10
        min_wx, max_wx = -20/3.6, 20/3.6
        min_wy, max_wy = -20/3.6, 20/3.6
        min_va, max_va =  40/3.6, 60/3.6
        list_wx = np.linspace(min_wx, max_wx, N_prior)
        list_wy = np.linspace(min_wy, max_wy, N_prior)
        list_va = np.linspace(min_va, max_va, N_prior)
        # Compute a first time
        out = compute_pdf(list_wx, 
                          list_wy, 
                          list_va)
        (post_wx, post_wy, post_va,
         mean_wx, std_wx,
         mean_wy, std_wy,
         mean_va, std_va) = out
        
        for iter_pdf in range(4):
            # =============================================================================
            # Set the new range to evaluate the pdf    
            # =============================================================================
            # Update the info
            min_wx, max_wx = np.min(list_wx), np.max(list_wx)
            min_wy, max_wy = np.min(list_wy), np.max(list_wy)
            min_va, max_va = np.min(list_va), np.max(list_va)
            # We will iterate the condition on the list of parameter, in order to 
            # simplify the if statements. Apologize if you find it more confusing !
            list_std  = [std_wx, std_wy, std_va]
            list_mean = [mean_wx, mean_wy, mean_va]
            list_prev_min = [min_wx, min_wy, min_va]
            list_prev_max = [max_wx, max_wy, max_va]
            
            N_prms = len(list_std)
            list_new_min = np.zeros(N_prms)
            list_new_max = np.zeros(N_prms)
            
            # Set the new range for each parameters
            for i_range in range(N_prms):
                std  = list_std [i_range]
                mean = list_mean[i_range]
                # Previous span range of the parameter
                prev_range = list_prev_max[i_range] - list_prev_min[i_range]    
                prev_reso = prev_range/(N_prior-1)
                if std <= prev_reso:
                    # The std is not resolved. We will use the spacing between the points 
                    # as the new range
                    print('new_range = prev_reso')
                    new_range = prev_reso
                elif std >=  prev_range:
                    # The std is too big and we cannot expect the subsequent iteration 
                    # to improve. So we just don't change the range. 
                    print('new_range = prev_range')
                    new_range = prev_range
                else:
                    # The new range will shrink around the pdf
                    print('new_range = std*N_prior/2')
                    new_range = std*N_prior/2# About 2pts in the std
                    
                list_new_min[i_range] = mean - 0.5*new_range
                list_new_max[i_range] = mean + 0.5*new_range  
        
            # Ready to recalculate the pdf on the new range
            list_wx = np.linspace(list_new_min[0], list_new_max[0], N_prior)
            list_wy = np.linspace(list_new_min[1], list_new_max[1], N_prior)
            list_va = np.linspace(list_new_min[2], list_new_max[2], N_prior)
            # Compute a first time
            out = compute_pdf(list_wx, 
                              list_wy, 
                              list_va)
            (post_wx, post_wy, post_va,
             mean_wx, std_wx,
             mean_wy, std_wy,
             mean_va, std_va) = out        
            
            plt.figure(tight_layout=True)
            
            plt.subplot(311)
            plt.plot(list_wx*3.6, post_wx, label='Wind x', color='C0')
            plt.xlabel('Wind x speed (km/h)')
            plt.legend()
            plt.title('Posterior probability\nIteration %d'%iter_pdf)
            
            plt.subplot(312)
            plt.plot(list_wy*3.6, post_wy, label='Wind y', color='C1')
            plt.xlabel('Wind y speed (km/h)')
            plt.legend()
            
            plt.subplot(313)
            plt.plot(list_va*3.6, post_va, label='Aircraft speed', color='C2')
            plt.xlabel('Aircraft speed (km/h)')
            plt.legend()
            
        # =============================================================================
        # Check the data
        # =============================================================================
        
        
        
        # The data and the inference
        plt.figure(tight_layout=True)
        # The data
        plt.errorbar(3.6*meas_vx, 3.6*meas_vy, 
                     xerr= 3.6*sig_vx, yerr=3.6*sig_vy, 
                     fmt='.')
        # The infered model
        phi_fit = np.linspace(0, 2*np.pi, 100)
        vx_fit = model_vx(mean_wx, mean_va, phi_fit)
        vy_fit = model_vy(mean_wy, mean_va, phi_fit)
        plt.plot(vx_fit*3.6, vy_fit*3.6, label='Fit')
        plt.xlabel('Vx (km/h)')
        plt.ylabel('Vy (km/h)')
        plt.title('Flight velocities\n'+
                  '\n mean_wx, std_wx = %.2f +- %.2f km/h'%(mean_wx*3.6, std_wx*3.6) + 
                  '\n mean_wy, std_wy = %.2f +- %.2f km/h'%(mean_wy*3.6, std_wy*3.6) +
                  '\n mean_va, std_va = %.2f +- %.2f km/h'%(mean_va*3.6, std_va*3.6))
        plt.legend()
        plt.axis('equal')
        

        
        
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
        plt.title('Wind Estimation: %.2f km/h'%self.w+
                  ', pointing at %.1f degree'%(self.angle*180/np.pi))
        
if __name__ == '__main__':
    # This part of the script runs only when the script file is run alone   
    
    # Fake some data to test
    (vx, vy) = (5.7, -1.5) # Please keep adding data
    
    self = WindEstimator()
    
    
    

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
    
    
    
    
    
    
    
    