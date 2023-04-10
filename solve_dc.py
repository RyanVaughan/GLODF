
import numpy as np
from pwr_obj import get_sys_problem10, get_sys_example3_11, get_sys_illinois
from solve_ac import calculate_line_flows, calculate_line_losses

def j_dc(sys):
    num_bus = np.shape(sys.bus.no)[0]
    J = np.zeros((num_bus,num_bus))
    
    for i in range(len(sys.line.no)):
        Y = 1/(sys.line.X[i])
        
        # off diagonals (upper and lower) get -Y
        J[sys.line.bfrom[i]-1, sys.line.bto[i]-1] -= Y
        J[sys.line.bto[i]-1, sys.line.bfrom[i]-1] -= Y
        
        # diagonals of both to and from both get Y
        J[sys.line.bfrom[i]-1, sys.line.bfrom[i]-1] += Y
        J[sys.line.bto[i]-1, sys.line.bto[i]-1] += Y
    
    return J

def solve_dc(sys):
    J = j_dc(sys)
    # it is safe to assume slack bus is #0 so delete #0
    J = J[1:,1:]

    P = sys.bus.pgen[1:] - sys.bus.pload[1:]
    Q = sys.bus.qgen[2:] - sys.bus.qload[2:]

    theta = np.linalg.solve(J, P)
    theta = np.insert(theta, 0, 0)

    v = np.linalg.solve(J[1:,1:], Q)
    v = np.insert(v, 0, [0,0]) + sys.bus.vabs
    v[2:] = v[2:] + 0.05 # because voltages were assumed to be 1 but are 1.05

    J = j_dc(sys)

    P = np.dot(J, theta)
    Q = np.dot(J, v)
    
    return [v, theta, P, Q]

def solve_dc_output(sys, v, theta, P, Q):
    print("angles (rad):\n" + str(theta))
    print("voltages:\n" + str(v))
    print("active powers:\n" + str(P))
    print("reactive powers:\n" + str(Q))
    
    # line loss
    p_loss = np.sum(P)
    print("power loss in system:\n" + str(p_loss))
    
    print("")
    
    f = calculate_line_flows(sys, v, theta)
    print("line real power flows:\n" + str(f[0]))
    #print("line imaginary power flows:\n" + str(f[1]))
    
    l = calculate_line_losses(sys, v, theta)
    print("line power losses:\n" + str(-1*l[0]))
    print("sum losses:\n" + str(np.sum(-1*l[0])))
    

if __name__ == "__main__":
    sys = get_sys_problem10()
    [v, theta, P, Q] = solve_dc(sys)
    solve_dc_output(sys, v, theta, P, Q)