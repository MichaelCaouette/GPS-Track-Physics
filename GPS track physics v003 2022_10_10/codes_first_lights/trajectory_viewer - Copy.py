# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:28:26 2022


Goal:
    Create a convenient viewer for 3D trajectory. 

@author: mcaoue2
"""

import traceback #Very usefull command to use for getting the last-not-printed error
_p = traceback.print_last # Type "_p()" in the console to access the last error

from spinmob import egg # For making the GUI
# Could lso use egg.pyqtgraph.opengl when importing the 3D widgets
import pyqtgraph.opengl as gl # Accessing the 3D GUI
import pyqtgraph as pg

import numpy as np

# Debug stuff.
_debug_enabled     = False
def _debug(*a):
    if _debug_enabled: 
        s = []
        for x in a: s.append(str(x))
        print(', '.join(s))

class TrajectoryViewer(egg.gui.Window):
    """
    GUI to see a 3D trajectory consisting of a set of 3D spatial position VS 
    time.
    """
    def __init__(self, name="Trajecotry viewer", size=[1500,800]): 
        """
        Initiate the GUI
        """    
        _debug('TrajectoryViewer.__init__')
        _debug('Make each day your masterpiece. – John Wooden')
        
        # Run the basic stuff for the initialization
        egg.gui.Window.__init__(self, title=name, size=size)
        
        # We use pyqtgraph for plottinh in 2D and 3D. 
        # https://pyqtgraph.readthedocs.io/en/latest/3dgraphics/index.html
        # The 3D graphics system in pyqtgraph is composed of a view widget and 
        # several graphics items (all subclasses of GLGraphicsItem) which can 
        # be added to a view widget.
        
        # The 3D GUI. It will be filled in the method "_update_3Dviewer". 
        self.view3D = gl.GLViewWidget()
        self.place_object(self.view3D, 
                          row=0, column=0, 
                          column_span=1, row_span=1, alignment=0) # Settting the alignement seems important !
        
        # Plots for showing the track in 2D
        self.win_plot2D = egg.pyqtgraph.GraphicsWindow(title="Trajectory")
        self.plot_groundproj = self.win_plot2D.addPlot(title="Ground Projection")
        self.win_plot2D.nextRow()
        self.plot_elev       = self.win_plot2D.addPlot(title="Elevation VS time")   
        self.place_object(self.win_plot2D, 
                          row=0, column=1,
                          row_span=1, column_span=1, alignment=1)            
        

        
        # Initiate important attributes
        self.n_traj = 0
        
    def _update_2DViewer(self):
        """
        Update the 2D visualization of the trajectory. 

        Returns
        -------
        None.

        """
        _debug('TrajectoryViewer._update_2DViewer')
        # 2D plot
        
        self.plot_groundproj.clear()
        self.plot_elev.clear()
        for i in range(self.n_traj):
            xs, ys, zs = self.list_traj[i]
            ts = self.list_ts[i]
            # The RGB color with pen ranges from 0 to 255
            c = pg.glColor((i, self.n_traj*1.3))
            self.c_tuple = (c[0]*255, c[1]*255, c[2]*255) 
            # The projection of the trajectory on the ground
            self.plot_groundproj.plot(xs, ys, pen=self.c_tuple)#, name="Data") 
            self.plot_groundproj.setLabel('bottom', self.labels[1], units=self.units[1])
            self.plot_groundproj.setLabel('left'  , self.labels[2], units=self.units[2])        
            # The elevation VS time
            self.plot_elev.plot(ts, zs, pen=self.c_tuple)#,  name="Data") 
            self.plot_elev.setLabel('bottom', self.labels[0], units=self.units[0])
            self.plot_elev.setLabel('left'  , self.labels[3], units=self.units[3])
        
        return
    
    def _update_3DViewer(self):
        """
        Update the 3D visualization of the trajectory. 

        Returns
        -------
        None.

        """
        _debug('TrajectoryViewer._update_3DViewer')
        # TODO: Make the color change with altitude :3 (Only if one trajectory!)
        
        # Clear everything that exist
        # The clear function doesn't exist, so I copy pasted the one from the 
        # online code: https://pyqtgraph.readthedocs.io/en/latest/_modules/pyqtgraph/opengl/GLViewWidget.html#GLViewWidget.clear
        # self.view3D.clear()
        for item in self.view3D.items:
            item._setView(None)
        self.view3D.items = []
        self.view3D.update()    
        
        # Show each trajectory
        for i in range(self.n_traj):
            # Extrac the positions
            self.traj = self.list_traj[i]
            self.pos_arg = np.vstack(self.traj).transpose() # Get the proper structure for the Line plotter
            # Plot them
            # How to set the data can be found here: 
            # https://pyqtgraph.readthedocs.io/en/latest/3dgraphics/gllineplotitem.html
            self.my_LinePlotItem = gl.GLLinePlotItem()
            # self.my_LinePlotItem.setData(pos=self.pos_arg)
            self.my_LinePlotItem.setData(pos=self.pos_arg, 
                                         color=pg.glColor((i, self.n_traj*1.3)), #(0.5, 0.9, 0.6)
                                         width=self.list_width[i],
                                         antialias=True) #enables smooth line drawing
            self.view3D.addItem(self.my_LinePlotItem)
            
        # Make the view to span all the trajectories
        # TODO: make this a separated method, which could be called often to
        # reset the view. 
        # Center is the mean center
        self.view3D.opts['center'][0] = self.x0
        self.view3D.opts['center'][1] = self.y0
        self.view3D.opts['center'][2] = self.z0
        # Make the distance from the center such that we see everything
        L = max(self.max_x - self.min_x, self.max_y - self.min_y) # Horizotal distance to cover
        angle = self.view3D.opts['fov'] # horizontal field of view in degrees
        self.view3D.opts['distance'] = 0.5*L / np.tan(angle) # Distance of the camera from the center.

        # Add the ground that span the whole projection of the trajectories on the x-y plane
        # A grid for now. One day it will be the map
        self.ground_plane = gl.GLGridItem()    
        # Set the size
        self.ground_plane.setSize(x=self.dx, 
                                  y=self.dy, 
                                  z=self.dz)     
        # Translate it to the center
        self.ground_plane.translate(dx = self.x0,
                                    dy = self.y0,
                                    dz = 0       ) # The ground will be the 0 atlitude
        # Add a confortable amount of spacing
        N_spacing = 20 
        self.ground_plane.setSpacing(x=self.dx/N_spacing, 
                                     y=self.dy/N_spacing, 
                                     z=self.dz/N_spacing)
        # Add this masterpiece lol
        self.view3D.addItem(self.ground_plane)
        
        
        
        
        return
        
    def set_trajectory(self,
                       ts, xs, ys, zs, 
                       rescale=False,
                       labels=('t', 'x', 'y', 'z'),
                       units=('', '', '', ''),
                       overwrite=True,
                       width=1):
        """
        Set the trajectory in the viewer and update it. 

        Parameters
        ----------
        ts : list of float
            List of the times.           
        xs : list of float
            List of the x positions.
        ys : TYPE
            DESCRIPTION.
        zs : TYPE
            DESCRIPTION.
        rescale : bool, optionanl
            If true, the spatial position will be rescaled to span 0 to 1. 
            The default is False. 
        labels : tuple of 4 strings, optional
            The strings are the labels for the t,x, y and z axis.
            The default is ('t', 'x', 'y', 'z').            
        units : tuple of 4 strings, optional
            Each of the string is the label for the t,x, y, z unit.
            The default is ('', '', '', '').             
        overwrite : Bool, optional
            If true, the previous trajectory will be erased and replaced by the 
            input. Otherwise, it will add on the previous one. The default is True.
        width : float >0, optional
            Width of the line drawn    
            The default is 1. 

        Returns
        -------
        None.

        """
        _debug('TrajectoryViewer.set_trajectory')
        
        self.labels = labels
        self.units = units
        
        if rescale:
            xs = (xs-xs[0]) / (np.max(xs) - np.min(xs))
            ys = (ys-ys[0]) / (np.max(ys) - np.min(ys))
            zs = (zs-zs[0]) / (np.max(zs) - np.min(zs))
        
        # Verify that the size of each axis are the same
        # Instead of verifying each possible combination, we use the logic 
        # that if a=b and b=c, than a=c (and we don't have to verify that a=c)        
        if len(xs) != len(ys) :
            print('ERROR in TrajectoryViewer.set_trajectory --> len(xs) != len(ys)')
            return
        # Now xs and ys has the same size
        if len(xs) != len(zs) :
            print('ERROR in TrajectoryViewer.set_trajectory --> len(xs) != len(zs)')   
            return
        # Now xs and zs has the same size, therefore also same size as ys. 
        if len(xs) != len(ts) :
            print('ERROR in TrajectoryViewer.set_trajectory --> len(xs) != len(ts)')
            return
        # Now we knoen that all the combination have the same size.
        
        # A convienient function that will be used more than once. 
        def init_list_traj(my_class, ts, xs, ys, zs):
            # Initiate the trajectory
            my_class.list_traj = [[xs, ys, zs]]
            my_class.list_ts = [ts]
            my_class.list_width = [width]
            # Initiate the extrema of the plot
            my_class.min_x, my_class.max_x = np.min(xs), np.max(xs)
            my_class.min_y, my_class.max_y = np.min(ys), np.max(ys)
            my_class.min_z, my_class.max_z = np.min(zs), np.max(zs)
            
        # Add or overwrite the trajectory
        if self.n_traj == 0:
            # If it is the first trajectory, we have to initiate the list
            init_list_traj(self, ts, xs, ys, zs)
        else:
            # If it is not the first trajectory, we overwrite or not the 
            # existing trajectory
            if overwrite:
                # Re-initiate the list
                init_list_traj(self, ts, xs, ys, zs)
            else:
                # Add the trajectory to the list
                self.list_traj.append([xs, ys, zs])
                self.list_ts.append(ts)
                self.list_width.append(width)
                # Update the extrema
                if np.min(xs) < self.min_x:
                    self.min_x = np.min(xs)
                if np.min(ys) < self.min_y:
                    self.min_y = np.min(ys)                    
                if np.min(zs) < self.min_z:
                    self.min_z = np.min(zs)   
                if np.max(xs) > self.max_x:
                    self.max_x = np.max(xs)
                if np.max(ys) > self.max_y:
                    self.max_y = np.max(ys)                    
                if np.max(zs) > self.max_z:
                    self.max_z = np.max(zs)   
                
        # Update the number of trajectory
        self.n_traj = len(self.list_traj)
        
        # Update the spatial information
        # The size
        self.dx = self.max_x - self.min_x 
        self.dy = self.max_y - self.min_y 
        self.dz = self.max_z - self.min_z 
        # The center
        self.x0 = 0.5*(self.max_x + self.min_x)
        self.y0 = 0.5*(self.max_y + self.min_y)
        self.z0 = 0.5*(self.max_z + self.min_z)
            
        # Awesome print info
        _debug('TrajectoryViewer.n_traj = ', self.n_traj)
        _debug('TrajectoryViewer.min_x, max_x = ', self.min_x, self.max_x)
        _debug('TrajectoryViewer.max_y, max_y = ', self.min_y, self.max_y)
        _debug('TrajectoryViewer.min_z, max_z = ', self.min_z, self.max_z)
        
        # Update the viewers
        self._update_3DViewer()
        self._update_2DViewer()
        
        return 
        
        
        
if __name__ == '__main__':
    _debug_enabled = True
    
    self = TrajectoryViewer()
    self.show()

    # Set some data
    # A Spiral
    s = np.linspace(0, 6*2*np.pi, 300)
    xs = 20*s*np.cos(s)/np.max(s)
    ys = 10*s*np.sin(s)/np.max(s)
    zs = 40*np.sin(s*0.5*np.pi/np.max(s))    
    self.set_trajectory(s, xs, ys, zs)
    
    # An other trajectory
    s = np.linspace(0, 6*2*np.pi, 300)
    xs = 100 + 10*(s/np.max(s))**2*np.cos(s)
    ys = 30*(1-s/np.max(s))**2*np.sin(s)
    zs = 10+20*np.sin(s*0.5*np.pi/np.max(s)+0.5*np.pi)    
    self.set_trajectory(s, xs, ys, zs, overwrite=False, width=2)












        
        
