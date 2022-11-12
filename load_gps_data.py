# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 18:04:38 2022


Goal:
    Make a clean data extractor
    
    
@author: mcaoue2
"""


import gpxpy
import fitparse


# Debug stuff.
_debug_enabled     = False
def _debug(*a):
    if _debug_enabled: 
        s = []
        for x in a: s.append(str(x))
        print(', '.join(s))
        
        
class LoadGPS():
    """
    Class to load the .FIT or .GPX file and extract the relevant information.
    """
    def __init__(self): 
        """
        Initiate stuff
        """   
    def get_FIT(self, filepath):
        """
        For a .FIT file
        Extract the latitude, longitude, elevation and time.

        Parameters
        ----------        
        filepath : string
            Path of the .FIT file to load
        
        Returns
        -------
        Tuple of array for each dimension of the form
        (list_time, 
         list_lat,
         list_lon,
         list_ele)
        
        The time is given at "datetime.datetime" object.
        The unit for the elevation is meter. 
        The unit for latitude and longitude are degree
        
        """
        # DISCLAIMER and APOLOGIZE
        # The code is not super elegant. I made it to clear up some redondant
        # data and avoid missing data. 
        
        # Initiate the list of what we will record
        self.list_time = []
        self.list_lat  = []
        self.list_lon  = []
        self.list_ele  = []
        
        # The semicircle to degree conversion unit
        # Explained here  https://www.gps-forums.com/threads/explanation-sought-concerning-gps-semicircles.1072/
        sem2deg = 180/2**31
        
        # Load the FIT file
        fitfile = fitparse.FitFile(filepath)
        
        # Iterate over all messages of type "record"
        # (other types include "device_info", "file_creator", "event", etc)
        for record in fitfile.get_messages("record"):
            # Scan each data until we find the type that we care about
            # Note: there might be more efficient way to do this. But this is
            # what I have for now
            
            # Before appending the data, we will first check that none of them
            # are nones. So we first initiate them with Nones. 
            value_t = None
            value_x = None
            value_y = None
            value_z = None
            # The following variable are to avoid recording the double 
            # occurence, which seems to exist with my GARMIN GPS data
            once_t = False
            once_x = False
            once_y = False
            once_z = False
            
            for data in record:
                if data.name == 'timestamp' and not(once_t):
                    value_t = data.value
                    once_t  = True # Avoid double occurence
                if data.name == 'enhanced_altitude' and not(once_z):
                    value_z = data.value
                    once_z  = True # Avoid double occurence
                if data.name == 'position_lat' and not(once_y):
                    value_y = data.value
                    once_y  = True # Avoid double occurence                               
                if data.name == 'position_long' and not(once_x):
                    value_x = data.value 
                    once_x  = True # Avoid double occurence  
                    
            # Append the data point only if there are no None
            if not(None in [value_t, value_x, value_y, value_z]):
                self.list_time.append(value_t)
                self.list_lat .append(value_y*sem2deg)
                self.list_lon .append(value_x*sem2deg)
                self.list_ele .append(value_z)
            
        return  (self.list_time,
                 self.list_lat ,
                 self.list_lon ,
                 self.list_ele  )             

    def get_GPX(self, filepath):
        """
        For a .GPX file
        Extract the latitude, longitude, elevation and time.

        Parameters
        ----------        
        filepath : string
            Path of the .GPX file to load
        
        Returns
        -------
        Tuple of array for each dimension of the form
        (list_time, 
         list_lat,
         list_lon,
         list_ele)
        
        The time is given at "datetime.datetime" object.
        The unit for the elevation is meter. 
        The unit for latitude and longitude are degree
        
        """
        # Initiate the list of what we will record
        self.list_time = []
        self.list_lat  = []
        self.list_lon  = []
        self.list_ele  = []
        
        with open(filepath) as fh:
            self.gpx_file = gpxpy.parse(fh) 
            
        _debug("File has {} track(s).".format(len(self.gpx_file.tracks)))   
        _debug("Track has {} segment(s).".format(len(self.gpx_file.tracks[0].segments))) 
        # So we take the only track and segment. 
        self.segment = self.gpx_file.tracks[0].segments[0]
        _debug('The segment has {} points'.format(len(self.segment.points)))          

        self.list_trackpts = self.segment.points

        # Take all the track points, outliers or not. 
        for i,p, in enumerate(self.list_trackpts):
            self.list_ele .append( p.elevation )
            self.list_lat .append( p.latitude  )
            self.list_lon .append( p.longitude )
            self.list_time.append( p.time      )
            
        return  (self.list_time,
                 self.list_lat ,
                 self.list_lon ,
                 self.list_ele  )         

    def print_all_FIT_data(self, filepath):
        """
        Scan and print all the data in the file. 
        Useful for checking up what is avaiable. 
        
        Parameters
        ----------        
        filepath : string
            Path of the .FIT file to load
            
        Return:
            None
        """
        print('Loading the module')
        # Load the FIT file
        fitfile = fitparse.FitFile(filepath)
        
        # Iterate over all messages of type "record"
        # (other types include "device_info", "file_creator", "event", etc)
        for record in fitfile.get_messages("record"):
            list_cool = []
            list_val  = []
            # Records can contain multiple pieces of data (ex: timestamp, latitude, longitude, etc)
            for data in record:
        
                # Print the name and value of the data (and the units if it has any)
                if data.units:
                    print(" * {}: {} ({})".format(data.name, data.value, data.units))
                else:
                    print(" * {}: {}".format(data.name, data.value))
                # Example of extraction 
                if data.name in ['timestamp',
                                 'enhanced_altitude',
                                 'position_lat',
                                 'position_long']:
                    list_cool.append(data.name)
                    list_val.append(data.value)
                    
            print()
            print('GOT -->')
            print(list_cool)
            print(list_val)
        
            print("---")        

if __name__ == '__main__':
    # This part of the script runs only when the script file is run alone    
    _debug_enabled     = True
    
    self = LoadGPS()
    # Test a .FIT file
    my_file = "C:\\Users\mcaoue2\Desktop\GPS track physics v004 2022_11_02\Raw Data\FIT files\TRK8.FIT"
    t, x,y, z = self.get_FIT(my_file)
    # Test a .GPX file
    # # my_file = "C:\\Users\mcaoue2\Desktop\GPS track physics v004 2022_11_02\Raw Data\TRK8.GPX"
    # my_file = "C:\\Users\mcaoue2\Desktop\GPS track physics v004 2022_11_02\Raw Data\TRK13 testing lay the GPS on top of car.GPX"
    # t, x,y, z = self.get_GPX(my_file)
    
    # Check the data
    print('Verification if all lenght are equal: ',len(t), len(x), len(y), len(z))
    print('Total elapsed time: ', (t[-1].timestamp() - t[0].timestamp())/3600)
    str_init_t = str(t[0]) # String telling the starting time
    print('Starting on ', str_init_t)
    
    import matplotlib.pyplot as plt
    plt.figure()
    plt.subplot(211)
    plt.title('Starting on %s'%str_init_t)
    plt.plot(z,'.-')
    plt.subplot(212)
    plt.plot(x, y, '.-', color='C1')
    
        
    
    
    
    
        







