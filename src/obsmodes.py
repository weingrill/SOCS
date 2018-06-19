#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Joerg Weingrill"
__copyright__ = "Copyright 2018, Leibniz Institute for Astrophysics Potsdam (AIP)"
__credits__ = ["Joerg Weingrill"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Joerg Weingrill"
__email__ = "jweingrill@aip.de"
__status__ = "Development"
__date__="2018-06-18"
obsmodes = {'BVI': {'ExposureTime':     2.0,
                    'ExposureIncrease': '3,30,300,2,20,200,1,10,100',
                    'FilterSequence':   'B,B,B,V,V,V,I,I,I',
                    'pernight':         1
                    },
            'BVR': {'ExposureTime':     3.0,
                    'ExposureIncrease': '2,20,200,1,10,100,1,10,100',
                    'FilterSequence':   'B,B,B,V,V,V,R,R,R',
                    'pernight':         1
                    },
            
            'UBV': {'ExposureTime':    2.0, # 
                    'ExposureIncrease': '3,30,300,2,20,200,1,10,100',
                    'FilterSequence':   'U,U,U,B,B,B,V,V,V',
                    'pernight':         1
                    },
            
            'uvby': {'ExposureTime':    6.0,
                    'ExposureIncrease': '1,1,1,1,10,10,10,10,100,100,100,100',
                    'FilterSequence':   'u,v,b,y,u,v,b,y,u,v,b,y',
                    'MoonDistance.Min': 15,
                    'AirmassTarget.Max':3.0,
                    'pernight':         1
                    },
            
            'Ha': {'ExposureTime':     6.0,
                    'ExposureIncrease': '1,10,100,1,10,100',
                    'FilterSequence':   'haw,haw,haw,han,han,han',
                    'MoonDistance.Min': 15,
                    'AirmassTarget.Max':3.0,
                    'pernight':         1
                    },
            
            'rot': {'ExposureTime':     15.0, 
                    'ExposureIncrease': '1,10,20',
                    'FilterSequence':   'V,V,R',
                    'pernight':         6 
                    },
            
            'frot': {'ExposureTime':     6.0, 
                    'ExposureIncrease': '1,2',
                    'FilterSequence':   'V,V',
                    'pernight':         6 
                    },
            
            'rottest': {'ExposureTime': 120.0, 
                    'ExposureIncrease': '1,1,1',
                    'FilterSequence':   'V,R,rp',
                    'pernight':         6 
                    }
                    
            }