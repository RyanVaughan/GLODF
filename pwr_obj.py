
import numpy as np

class Bus:
    def __init__(self):
        self.no    = np.array([], dtype=int)
        self.btype = np.array([], dtype=int)
        self.vabs  = np.array([], dtype=float)
        self.angle = np.array([], dtype=float)
        self.pgen  = np.array([], dtype=float)
        self.qgen  = np.array([], dtype=float)
        self.pload = np.array([], dtype=float)
        self.qload = np.array([], dtype=float)
        
    def add(self, no, btype, vabs, angle, pgen, qgen, pload, qload):
        
        # A vector can be passed but all must have the same length
        if not all(np.shape(x) == np.shape(no) for x in (btype, vabs, angle, 
                                                         pgen, qgen, pload, 
                                                         qload)):
            raise ValueError("All parameters must have the same length")
            
        # only scalars or vectors 
        if len(np.shape(no)) > 1:
            raise ValueError("All parameters must be scalars or vectors")
                             
        self.no    = np.append(self.no,    no   )
        self.btype = np.append(self.btype, btype)
        self.vabs  = np.append(self.vabs,  vabs )
        self.angle = np.append(self.angle, angle)
        self.pgen  = np.append(self.pgen,  pgen )
        self.qgen  = np.append(self.qgen,  qgen )
        self.pload = np.append(self.pload, pload)
        self.qload = np.append(self.qload, qload)
        
class Line:
    def __init__(self):
        self.no    = np.array([], dtype=int)
        self.bfrom = np.array([], dtype=int)
        self.bto   = np.array([], dtype=int)
        self.R     = np.array([], dtype=float)
        self.X     = np.array([], dtype=float)
        self.B     = np.array([], dtype=float)

    def add(self, no, bfrom, bto, R, X, B):
        
        # A vector can be passed but all must have the same length
        if not all(np.shape(x) == np.shape(no) for x in (bfrom, bto, R, X, B)):
            raise ValueError("All parameters must have the same length")
            
        # only scalars or vectors (can append multiple at a time)
        if len(np.shape(no)) > 1:
            raise ValueError("All parameters must be scalars or vectors")
            
        self.no    = np.append(self.no,    no   )
        self.bfrom = np.append(self.bfrom, bfrom)
        self.bto   = np.append(self.bto,   bto  )
        self.R     = np.append(self.R,     R    )
        self.X     = np.append(self.X,     X    )
        self.B     = np.append(self.B,     B    )
        
class pwr_sys:
    def __init__(self, bus, line):
        self.bus  = bus
        self.line = line
        
def get_sys_problem10():
    # problem 10
    bus = Bus()
    bus.add(1, 0, 1.05, 0,   0, 0,  0.25,  0.1)
    bus.add(2, 1, 1.05, 0, 0.5, 0,  0.15, 0.05)
    bus.add(3, 2, 1.00, 0,   0, 0, 0.275, 0.11)
    bus.add(4, 2, 1.00, 0,   0, 0,     0,    0)
    bus.add(5, 2, 1.00, 0,   0, 0,  0.15, 0.09)
    bus.add(6, 2, 1.00, 0,   0, 0,  0.25, 0.15)

    line = Line()
    line.add(1, 1, 4, 0.020, 0.185, 0.009)
    line.add(2, 1, 6, 0.031, 0.259, 0.010)
    line.add(3, 2, 3, 0.006, 0.025, 0.000)
    line.add(4, 2, 5, 0.071, 0.320, 0.015)
    line.add(5, 4, 6, 0.024, 0.204, 0.010)
    line.add(6, 3, 4, 0.075, 0.067, 0.000)
    line.add(7, 5, 6, 0.025, 0.150, 0.017)
    
    return pwr_sys(bus, line)

def get_sys_problem10_no4():
    # problem 10
    bus = Bus()
    bus.add(1, 0, 1.05, 0,   0, 0,  0.25,  0.1)
    bus.add(2, 1, 1.05, 0, 0.5, 0,  0.15, 0.05)
    bus.add(3, 2, 1.00, 0,   0, 0, 0.275, 0.11)
    bus.add(4, 2, 1.00, 0,   0, 0,     0,    0)
    bus.add(5, 2, 1.00, 0,   0, 0,  0.15, 0.09)
    bus.add(6, 2, 1.00, 0,   0, 0,  0.25, 0.15)

    line = Line()
    line.add(1, 1, 4, 0.020, 0.185, 0.009)
    line.add(2, 1, 6, 0.031, 0.259, 0.010)
    line.add(3, 2, 3, 0.006, 0.025, 0.000)
    line.add(4, 4, 6, 0.024, 0.204, 0.010)
    line.add(5, 3, 4, 0.075, 0.067, 0.000)
    line.add(6, 5, 6, 0.025, 0.150, 0.017)
    
    return pwr_sys(bus, line)


def get_sys_problem10_no45():
    # problem 10
    bus = Bus()
    bus.add(1, 0, 1.05, 0,   0, 0,  0.25,  0.1)
    bus.add(2, 1, 1.05, 0, 0.5, 0,  0.15, 0.05)
    bus.add(3, 2, 1.00, 0,   0, 0, 0.275, 0.11)
    bus.add(4, 2, 1.00, 0,   0, 0,     0,    0)
    bus.add(5, 2, 1.00, 0,   0, 0,  0.15, 0.09)
    bus.add(6, 2, 1.00, 0,   0, 0,  0.25, 0.15)

    line = Line()
    line.add(1, 1, 4, 0.020, 0.185, 0.009)
    line.add(2, 1, 6, 0.031, 0.259, 0.010)
    line.add(3, 2, 3, 0.006, 0.025, 0.000)
    line.add(4, 3, 4, 0.075, 0.067, 0.000)
    line.add(5, 5, 6, 0.025, 0.150, 0.017)
    
    return pwr_sys(bus, line)

def get_sys_example3_11():
    # example 3.11
    bus = Bus()
    bus.add(1, 0, 1.02, 0,   0, 0,   0,   0)
    bus.add(2, 1, 1.00, 0, 0.5, 0,   0,   0)
    bus.add(3, 2, 1.00, 0,   0, 0, 1.2, 0.5)
    
    line = Line()
    line.add(1, 1, 2, 0.02, 0.3, 0.15)
    line.add(2, 1, 3, 0.01, 0.1,  0.1)
    line.add(3, 2, 3, 0.01, 0.1,  0.1)
    
    return pwr_sys(bus, line)

def get_sys_illinois():
    # new example
    # https://aledan.ece.illinois.edu/files/2017/04/TSG_2017_a.pdf
    bus = Bus()
    # no, btype, vabs, angle, pgen, qgen, pload, qload
    bus.add(1, 0, 1.04 , 0, 1.5973, 0,   0,   0)
    bus.add(2, 1, 1.025, 0, 0.7910, 0,   0,   0)
    bus.add(3, 2, 1.00 , 0, 0,      0, 2.35, -0.5)
    
    line = Line()
    #        no, bfrom, bto, R, X, B
    line.add(1, 1, 2, 1.3652, 1/-11.6041, 0.088)
    line.add(2, 2, 3, 0.7598, 1/-6.1168,  0.153)
    line.add(3, 1, 3, 1.1677, 1/-10.7426, 0.079)
    
    return pwr_sys(bus, line)