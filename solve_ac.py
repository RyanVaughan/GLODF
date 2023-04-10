
import numpy as np
from pwr_obj import pwr_sys
from pwr_obj import get_sys_problem10, get_sys_example3_11, get_sys_illinois

###############################################################################
# AC Power Flow
def create_Ybus(sys):
    num_bus = np.shape(sys.bus.no)[0]
    Ybus = np.zeros((num_bus,num_bus), dtype=np.complex128)
    
    for i in range(len(sys.line.no)):
        Y = 1/(sys.line.R[i] + 1j*sys.line.X[i])
        
        # off diagonals (upper and lower) get -Y
        Ybus[sys.line.bfrom[i]-1, sys.line.bto[i]-1] -= Y
        Ybus[sys.line.bto[i]-1, sys.line.bfrom[i]-1] -= Y
        
        # diagonals of both to and from both get Y and half B
        Ybus[sys.line.bfrom[i]-1, sys.line.bfrom[i]-1] += Y + 1j*sys.line.B[i]/2
        Ybus[sys.line.bto[i]-1, sys.line.bto[i]-1] += Y + 1j*sys.line.B[i]/2
        
    return Ybus    

# Both inputs should be complex
def residuals(sys, Ybus, V):
    # find appparent power 
    S = (sys.bus.pgen - sys.bus.pload) + (sys.bus.qgen - sys.bus.qload)*1.0j
    S = S - V * np.conj(np.dot(Ybus, V))
    
    # seperate active and reactive power
    residue = np.concatenate((np.real(S), np.imag(S)))
    
    # delete the residues that are not needed
    delRows = del_vars(sys)
    residue = np.delete(residue, delRows)
    return residue
    
# find what vars are not needed
def del_vars(sys):
    num_bus = np.shape(sys.bus.no)[0]
    out = np.arange(2*num_bus)
    for i in range(num_bus):
        if sys.bus.btype[i] == 2:
            d = (np.argwhere(out == i), np.argwhere(out == num_bus+i))
            out = np.delete(out, d)
        elif sys.bus.btype[i] == 1:
            d = np.argwhere(out == i)
            out = np.delete(out, d)
    return out

# both inputs are real
def jacobian(sys, Y, V, angle):
    num_bus = np.shape(sys.bus.no)[0]
    J = np.zeros((2*num_bus,2*num_bus))
    # for parallel computations, repeat V and theta for each row
    V = np.repeat(np.expand_dims(V, axis=1), num_bus, axis=1)
    theta = np.repeat(np.expand_dims(angle, axis=1), num_bus, axis=1)
    
    # fill all off diagonals
    dQ = np.transpose(V) * Y * np.exp((theta - np.transpose(theta))*1.0j)
    dP = V * dQ
    # fix the diagonals
    np.fill_diagonal(dP, np.sum(dP, axis=1) + np.diag(dP))
    np.fill_diagonal(dQ, np.sum(dQ, axis=1) + np.diag(dQ))
    # seperate
    dPdT = np.imag(dP)
    dPdV = np.real(dP)
    dQdT = -1*np.real(dQ) # because j^2 = -1 (Q and theta are imaginary)
    dQdV = np.imag(dQ)
    
    # fill into J
    J[0:num_bus, 0:num_bus] = dPdT
    J[num_bus:2*num_bus, 0:num_bus] = dPdV
    J[0:num_bus, num_bus:2*num_bus] = dQdT
    J[num_bus:2*num_bus, num_bus:2*num_bus] = dQdV
    
    # delete rows and collumns that arent needed
    delRows = del_vars(sys)
    J = np.delete(J, delRows, axis=0)
    J = np.delete(J, delRows, axis=1)
    
    return J

def solve_ac(sys, Y):
    num_bus = np.shape(sys.bus.no)[0]
    x = np.concatenate((sys.bus.angle, sys.bus.vabs))
    # print(str(0) + "th iteration: x = " + str(x))
    
    iters = 0 # iteration number
    for i in range(50): #loop for maximum iterations
        # xprev = x # store previous value in order to determine convergence
        J = jacobian(sys, Y, x[num_bus:2*num_bus], x[0:num_bus])
        residue = residuals(sys, Y, x[num_bus:2*num_bus] * np.exp(x[0:num_bus]*1.0j))
        if (np.allclose(residue, 0, atol=0.000005)): break
        delta = np.linalg.solve(J, residue)
        
        # expand delta by filling zeros for elements that arent being updated
        deltaExpand = np.zeros_like(x)
        keepRows = np.delete(np.arange(num_bus*2), del_vars(sys))
        np.put(deltaExpand, keepRows, delta)
        x = x - deltaExpand
        # print(str(iters+1) + "th iteration: x = " + str(x))
        
        iters += 1
    
    print("took " + str(iters+1) + " iterations")
    V = x[num_bus:2*num_bus] * np.exp(x[0:num_bus]*1.0j)
    S = V * np.conj(np.dot(Y, V))    
    return [V, S]
    
def solve_ac_output(sys, V, S):
    
    v = np.abs(V)
    theta = np.angle(V)
    
    num_bus = np.shape(sys.bus.no)[0]
    print("")
    
    # calculate residuals power at each bus
    r = np.concatenate((np.real(S), np.imag(S)))
    
    print("real power leaving bus:\n" + str(r[0:num_bus]))
    print("reactive power leaving bus:\n" + str(r[num_bus:2*num_bus]))
        
    print("voltages:\n" + str(v))
    print("angles (rad):\n" + str(theta))
        
    # line loss
    p_loss = np.sum(r[0:num_bus])
    print("power loss in system:\n" + str(p_loss))
        
    print("")
    
    f = calculate_line_flows(sys, v, theta)
    print("line real power flows:\n" + str(f[0]))
    print("line imaginary power flows:\n" + str(f[1]))
    
    l = calculate_line_losses(sys, v, theta)
    print("line power losses:\n" + str(-1*l[0]))
    print("sum losses:\n" + str(np.sum(-1*l[0])))
    
def calculate_line_flows(sys, V, theta):
    num_line = np.size(sys.line.no)
    P = np.zeros(num_line)
    Q = np.zeros(num_line)
    
    for i in range(num_line):
        V_from = V[sys.line.bfrom[i]-1]
        V_to = V[sys.line.bto[i]-1]
        t_from = theta[sys.line.bfrom[i]-1]
        t_to = theta[sys.line.bto[i]-1]
        Z = complex(sys.line.R[i], sys.line.X[i])
        Y = 1 / Z
        psi = np.angle(Y)
        Y = np.abs(Y)
        B_half = sys.line.B[i] / 2
        
        P[i] = V_from * V_to * Y * np.cos(t_from - t_to - psi) \
            - V_from**2 * Y * np.cos(psi)
        Q[i] = V_from * V_to * Y * np.sin(t_from - t_to - psi) \
            + V_from**2 * Y * np.sin(psi) - V_from**2 * B_half 
    
    return (P, Q)

def calculate_line_losses(sys, V, theta):
    num_line = np.size(sys.line.no)
    P = np.zeros(num_line)
    Q = np.zeros(num_line)
    
    for i in range(num_line):
        V_from = V[sys.line.bfrom[i]-1]
        V_to = V[sys.line.bto[i]-1]
        t_from = theta[sys.line.bfrom[i]-1]
        t_to = theta[sys.line.bto[i]-1]
        Z = complex(sys.line.R[i], sys.line.X[i])
        Y = 1 / Z
        psi = np.angle(Y)
        Y = np.abs(Y)
        B_half = sys.line.B[i] / 2
        
        P[i] = V_from * V_to * Y * np.cos(t_from - t_to - psi) \
            - V_from**2 * Y * np.cos(psi)
        Q[i] = V_from * V_to * Y * np.sin(t_from - t_to - psi) \
            + V_from**2 * Y * np.sin(psi) - V_from**2 * B_half 
            
    P2 = np.zeros(num_line)
    Q2 = np.zeros(num_line)
    
    for i in range(num_line):
        V_to = V[sys.line.bfrom[i]-1]
        V_from = V[sys.line.bto[i]-1]
        t_to = theta[sys.line.bfrom[i]-1]
        t_from = theta[sys.line.bto[i]-1]
        Z = complex(sys.line.R[i], sys.line.X[i])
        Y = 1 / Z
        psi = np.angle(Y)
        Y = np.abs(Y)
        B_half = sys.line.B[i] / 2
        
        P2[i] = V_from * V_to * Y * np.cos(t_from - t_to - psi) \
            - V_from**2 * Y * np.cos(psi)
        Q2[i] = V_from * V_to * Y * np.sin(t_from - t_to - psi) \
            + V_from**2 * Y * np.sin(psi) - V_from**2 * B_half 
    
    return (P2 + P, Q2 + Q)
    
if __name__ == "__main__":
    sys = get_sys_problem10()
    Y = create_Ybus(sys)
    [V, S] = solve_ac(sys, Y)
    solve_ac_output(sys, V, S)
