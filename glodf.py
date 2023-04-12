
import numpy as np
from pwr_obj import pwr_sys
from pwr_obj import get_sys_problem10, get_sys_example3_11, \
    get_sys_illinois, get_sys_problem10_no4, get_sys_problem10_no45
from solve_ac import solve_ac, create_Ybus, solve_ac_output, calculate_line_flows
from solve_dc import solve_dc, solve_dc_output


# return L x L diagonal matrix
def get_branch_susceptance(sys):
    return np.diag(1/sys.line.X)
    
# return L x (N-1) matrix
def get_reduced_incidence(sys):
    num_line = np.size(sys.line.no)
    num_bus = np.size(sys.bus.no)
    
    A = np.zeros((num_line, num_bus))
    
    for l in range(num_line):
        A[l,sys.line.bfrom[l]-1] = -1
        A[l,sys.line.bto[l] -1]  =  1
    
    # assume slack is first
    Ar = A[:,1:]
    return Ar

# return (N-1) x (N-1) matrix
def get_reduced_nodal_susceptance(B, Ar):
    return -1.0 * np.dot(np.dot(Ar.T, B), Ar)
        
def get_isf(sys):
    Bs = get_branch_susceptance(sys)
    Ar = get_reduced_incidence(sys)
    Bn = get_reduced_nodal_susceptance(Bs, Ar)
    
    psi = np.dot(np.dot(Bs, Ar), np.linalg.inv(Bn))
    
    return psi
    
def ptdf(sys, change_line):
    psi = get_isf(sys)
    bfrom = sys.line.bfrom[change_line - 1] - 2
    bto   = sys.line.bto[change_line   - 1] - 2
    
    phi = psi[:,bfrom] - psi[:,bto]
    
    return phi

def lodf(sys, change_line):
    phi = ptdf(sys, change_line)
    
    sigma = phi / (1 - phi[change_line - 1])
    sigma[change_line - 1] = 0
    
    return sigma #[np.arange(sigma.size) != change_line - 1]
    
def flow_after_lodf(sys, change_line, flow_before):
    sigma = lodf(sys, change_line)
    
    flow_after = flow_before + flow_before[change_line-1] * sigma
    flow_after[change_line-1] = 0
    
    return flow_after

def check_single_outage():
    print("\n\nChecking Single outage code...\n")
    sys = get_sys_problem10()

    # dc or ac
    use_dc = False
    
    if use_dc:
        # Y = create_Ybus(sys)
        [v, theta, P, Q] = solve_dc(sys)
        #solve_dc_output(sys, v, theta, P, Q)
    else:
        Y = create_Ybus(sys)
        [V, S] = solve_ac(sys, Y)
        # solve_ac_output(sys, V, S)
        v = np.abs(V)
        theta = np.angle(V)
    
    
    # isf = get_isf(sys)
    # print("isf:\n" + str(isf))
    # ptdf_slides = ptdf_from_slides(sys)
    # print("ptdf:\n" + str(ptdf_slides))

    f = calculate_line_flows(sys, v, theta)
    print("power flows before:\n" + str(f[0]))
    flow_after = flow_after_lodf(sys, 4, f[0])
    print("take out line for with LODF:\n" + str(flow_after))
    
    sys = get_sys_problem10_no4()
    if use_dc:
        # Y = create_Ybus(sys)
        [v, theta, P, Q] = solve_dc(sys)
        #solve_dc_output(sys, v, theta, P, Q)
    else:
        Y = create_Ybus(sys)
        [V, S] = solve_ac(sys, Y)
        # solve_ac_output(sys, V, S)
        v = np.abs(V)
        theta = np.angle(V)

    f = calculate_line_flows(sys, v, theta)
    print("system solved without line 4:\n" + str(np.insert(f[0], 3, 0)))

def glodf(sys, change_lines):
    change_lines = np.atleast_1d(change_lines)
    
    ########
    # copy from PTDF
    psi = get_isf(sys)
    bfrom = sys.line.bfrom[change_lines - 1] - 2
    bto   = sys.line.bto[change_lines   - 1] - 2
    phi = psi[:,bfrom] - psi[:,bto]
    
    # right side of equation is all lines
    right_side = phi.T
    
    # left side is identity - Phi of change lines
    Phi = right_side[:,change_lines-1]
    left_side = (np.eye(np.shape(Phi)[0]) - Phi)
    
    xi = np.linalg.solve(left_side, right_side)
    
    return xi

def flow_after_glodf(sys, change_lines, flow_before):
    change_lines = np.atleast_1d(change_lines)

    xi = glodf(sys, change_lines)
    
    #GLODFs times flow before
    delta_flow = xi.T @ flow_before[change_lines-1]
    flow_after = flow_before + delta_flow
    
    # ensure lines that are out have no flow
    flow_after[change_lines-1] = 0
    
    return flow_after
    
def check_multi_outage():
    print("\n\nChecking multi outage code...\n")
    sys = get_sys_problem10()

    # dc or ac
    use_dc = False
    
    if use_dc:
        # Y = create_Ybus(sys)
        [v, theta, P, Q] = solve_dc(sys)
        #solve_dc_output(sys, v, theta, P, Q)
    else:
        Y = create_Ybus(sys)
        [V, S] = solve_ac(sys, Y)
        # solve_ac_output(sys, V, S)
        v = np.abs(V)
        theta = np.angle(V)

    f = calculate_line_flows(sys, v, theta)
    print("power flows before:\n" + str(f[0]))
    flow_after = flow_after_glodf(sys, [4,5], f[0])
    print("take out line for with GLODF:\n" + str(flow_after))
    
    sys = get_sys_problem10_no45()
    if use_dc:
        # Y = create_Ybus(sys)
        [v, theta, P, Q] = solve_dc(sys)
        #solve_dc_output(sys, v, theta, P, Q)
    else:
        Y = create_Ybus(sys)
        [V, S] = solve_ac(sys, Y)
        # solve_ac_output(sys, V, S)
        v = np.abs(V)
        theta = np.angle(V)

    f = calculate_line_flows(sys, v, theta)
    print("system solved without line 4 or 5:\n" + str(np.insert(f[0], [3,3], [0,0])))

'''
###############################################################################
# from slides
def ptdf_from_slides(sys):
    num_line = np.size(sys.line.no)
    num_bus = np.size(sys.bus.no)
    
    ptdf = np.zeros((num_line, num_bus), dtype=np.complex128)
    # ptdf = np.zeros((num_line, num_bus))
    
    Z = np.linalg.inv(create_Ybus(sys))
    # X = np.imag(np.linalg.inv(create_Ybus(sys.line, sys.bus)))
    
    for i in range(num_bus):    
        for j in range(num_line):
            bfrom = sys.line.bfrom[j] - 1
            bto   = sys.line.bto[j] - 1
            z = sys.line.R[j] + 1j*sys.line.X[j]
            # b = 1 / sys.line.X[j]
            ptdf[j,i] = (Z[bfrom, i] - Z[bto, i]) / z
            # ptdf[j,i] = (X[bfrom, i] - X[bto, i]) * b

    return np.real(ptdf[:,1:])
    # return ptdf[:,1:]
'''
if __name__ == "__main__":
    np.set_printoptions(precision=5)
    
    check_single_outage()
    
    check_multi_outage()
    