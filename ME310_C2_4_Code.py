import math
import numpy as np

# Data (Fouling)
Rf_i = 0.00035  # m^2.K/W
Rf_o = 0.00088  # m^2.K/W

# Unit Conversion Functions
def m2ft(meters):
    feet = meters * 3.28084
    return feet

def ft2m(feet):
    meters = feet * 0.3048
    return meters

def in2m(inches):
    meters = inches * 0.0254
    return meters

def m2in(meters):
    inches = meters / 0.0254
    return inches

def sqm2sqft(square_meters):
    square_feet = square_meters * 10.7639
    return square_feet

# Assumed Initial Values
Ts_out = 110  # °C
Tt_out = 65 # °C
Uo = 150 # W/m^2.K
L = ft2m(3)  # m 
ODt = in2m(0.5)  # m 
IDt = in2m(0.33)  # m 
PR = 1.25
Pitch_type = "square"
CL = 1
CTP = 0.93
Nt = 10
Nb = 4
Np = 1
B = in2m(9) # m
Ds = in2m(6) # m

# Derived Initial Values
Ki = 1/4 * IDt**2
Ko = 1/4 * ODt**2
percentage = 0.25
a = IDt * percentage
b = Ki / a

# Initialize the variables
iterations = 0
max_iterations = 1000
convergence_limit = 1E-5
red_flag = "off"

## Hydraulic Diameter for Tube
Dht = (4*a*b) / (3*(a+b) - math.sqrt((3*a+b)*(a+3*b)))

## Equivalent Diameter for Tube
Det = 2*math.sqrt(a*b)

## Equivalent for ODt
Det_o = Det + ODt - IDt

## Pitch Length
Pt = PR * Det_o

## Hydraulic Diameter for Shell
if Pitch_type == "square":
    Dhs = (4*Pt**2) / (math.pi * Det_o) - Det_o
elif Pitch_type == "triangular":
    Dhs = (2*math.sqrt(3)*Pt**2) / (math.PI * Det_o) - Det_o

## Equivalent Diameter for Shell
Des = Dhs

## Clearance Length
C = Pt - Det_o

## Bundle Cross-flow Area
As = (Ds * C * B) / Pt

## Heat Transfer Area
Ao = Nt * math.pi * Det_o * L

## Total Tube Inner Surface Area
Ai = Nt * math.pi * Det * L

## Total x-section Area
At = (Nt / Np) * (math.pi / 4) * Det**2


while iterations < max_iterations:

    # Property Tables
    ## Properties of Saturated Water
    T_w = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]  # °C
    rho_w = [999.9, 999.7, 999.1, 998, 997, 996, 994, 992.1, 990.1, 988.1, 985.2, 983.3, 980.4, 977.5]  # Kg/m^3
    Cp_w = [4205, 4194, 4185, 4182, 4180, 4178, 4178, 4179, 4180, 4181, 4183, 4185, 4187, 4190]  # J/KgK
    K_w = [0.571, 0.580, 0.589, 0.598, 0.607, 0.615, 0.623, 0.631, 0.637, 0.644, 0.649, 0.654, 0.659, 0.663]  # W/mK
    u_w = np.array([1.519, 1.307, 1.138, 1.002, 0.891, 0.798, 0.720, 0.653, 0.596, 0.547, 0.504, 0.467, 0.433, 0.404])*1E-3  # Pas
    Pr_w = [11.2, 9.45, 8.09, 7.01, 6.14, 5.42, 4.83, 4.32, 3.91, 3.55, 3.25, 2.99, 2.75, 2.55]

    ## Properties of Engine Oil
    T_o = [77, 87, 97, 107, 117, 127]  # °C
    rho_o = [853.9, 847.8, 841.8, 836, 830.6, 825.1]  # Kg/m^3
    Cp_o = [2118, 2161, 2206, 2250, 2294, 2337]  # J/KgK
    K_o = np.array([138, 138, 137, 136, 135, 134])*1E-3  # W/mK
    u_o = np.array([3.56, 2.52, 1.86, 1.41, 1.1, 0.874])*1E-2  # Pas
    Pr_o = [546, 395, 300, 233, 187, 152]

    ## Properties of Tube Wall (Copper)
    T_wall = [30, 90, 150]
    k_wall = [400.7, 395.9, 391.4]

    # Data (Shell Side - Engine Oil)
    Ts_in = 120  # °C
    Ts_film = np.mean([Ts_in, Ts_out]) # °C
    m_s = 1  # Kg/s
    Ks = np.interp(Ts_film, T_o, K_o)  # W/mK
    rho_s = np.interp(Ts_film, T_o, rho_o)  # Kg/m^3
    Cp_s = np.interp(Ts_film, T_o, Cp_o) # J/KgK
    us = np.interp(Ts_film, T_o, u_o)  # Pas
    Pr_s = np.interp(Ts_film, T_o, Pr_o)
    delta_Ps_max = 70E3  # Pa

    # Data (Tube Side - City Water)
    Tt_in = 60  # °C
    Tt_film = np.mean([Tt_in, Tt_out]) # °C
    Kt = np.interp(Tt_film, T_w, K_w)  # W/mK
    rho_t = np.interp(Tt_film, T_w, rho_w)  # Kg/m^3
    Cp_t = np.interp(Tt_film, T_w, Cp_w) # J/KgK
    ut = np.interp(Tt_film, T_w, u_w)  # Pas
    ub = ut # Pas
    Pr_t = np.interp(Tt_film, T_w, Pr_w)
    delta_Pt_max = 70E3  # Pa

    # Data (Tube Wall)
    Tw = np.mean([Ts_film, Tt_film]) # °C
    K_wall = np.interp(Tw, T_wall, k_wall)  # W/mK
    uw = np.interp(Tw, T_w, u_w)  # Pas

    # Mass Flow Rate in Tube Side
    m_t = (m_s * Cp_s * (Ts_in - Ts_out)) / (Cp_t * (Tt_out - Tt_in))

    ## Fluid Velocity
    Vt = m_t / (rho_t * At)
    
    ## Reynold's Number for Heat Transfer
    Re_t = (rho_t * Vt * Det) / ut
    if Re_t < 2300:
        flow_type_t = "Laminar"
    else:
        flow_type_t = "Turbulent"

    ## Reynold's Number for Pressure Drop
    Re_t_p = (rho_t * Vt * Dht) / ut
    if Re_t_p < 2300:
        flow_type_t = "Laminar"
    else:
        flow_type_t = "Turbulent"

    ## Friction Coefficient
    if 2300 <Re_t_p < 1E6:
        ft = 0.0228 + 16.5 / (Re_t_p**0.84)
    else:
        print(f"The chosen values are not correct. Choose different value.")
        print(f"If you change the temperatures, you will need to change thermal properties accordingly")
        print(f"Or you can change any of the mass flowrate as well, the other mass flowrate will be affected by it")
        print(f"Or you can change the assumed initial values")
        red_flag = "on"
        break

    ## Pressure Drop
    delta_Pt = Np*(ft*(L/Dht) + 4)*(1/2*rho_t*Vt**2)
    if delta_Pt >= delta_Pt_max:
        print(f"Pressure drop for tube side exceeds allowable limit.")
        print(f"The chosen values are not correct. Choose different value.")
        print(f"If you change the temperatures, you will need to change thermal properties accordingly")
        print(f"Or you can change any of the mass flowrate as well, the other mass flowrate will be affected by it")
        print(f"Or you can change the assumed initial values")
        red_flag = "on"
        break

    ## Nusselt Number 
    if 2300 <Re_t < 1E6:
        Nu_t = 0.0391 *(Re_t**0.762) * (Pr_t**0.4)
    else:
        print(f"The chosen values are not correct. Choose different value.")
        print(f"If you change the temperatures, you will need to change thermal properties accordingly")
        print(f"Or you can change any of the mass flowrate as well, the other mass flowrate will be affected by it")
        print(f"Or you can change the assumed initial values")
        red_flag = "on"
        break

    ## Convectional Heat Transfer Coefficient
    hi = (Nu_t * Kt) / Det
    

    # Shell-Side
    ## Velocity of Fluid 
    Vs = m_s / (rho_s*As)

    ## Reynold's Number
    Re_s = (rho_s * Vs * Dhs) / us
    if Re_s < 2300:
        flow_type_s = "Laminar"
    else:
        flow_type_s = "Turbulent"

    ## Friction Coefficient
    fs = math.exp(0.576-0.19*math.log(Re_s))

    ## Pressure Drop
    delta_Ps = fs*(Nb+1)*(Ds/Dhs)*(1/2*rho_s*Vs**2)
    if delta_Ps >= delta_Ps_max:
        print(f"Pressure drop for shell side exceeds allowable limit.")
        print(f"The chosen values are not correct. Choose different value.")
        print(f"If you change the temperatures, you will need to change thermal properties accordingly")
        print(f"Or you can change any of the mass flowrate as well, the other mass flowrate will be affected by it")
        print(f"Or you can change the assumed initial values")
        red_flag = "on"
        break
    
    ## Nusselt Number
    Nu_s = 0.36 * Re_s**0.55 * Pr_s**(1 / 3)

    ## Convectional Heat Transfer Coefficient
    ho = (Nu_s * Ks) / Des
    
    # Overall Heat Transfer Coefficient After Iteration
    Wall_resistance = math.log(ODt / IDt) / (2 * math.pi * K_wall * L)
    Uo_new = 1 / (Ao * (1 / (hi * Ai) + Wall_resistance + 1 / (ho * Ao) + Rf_i / Ai + Rf_o / Ao))

    # Heat Capacity Ratio
    Cs = m_s*Cp_s
    Ct = m_t*Cp_t
    C_max = max(Cs, Ct)
    C_min = min(Cs, Ct)
    c = C_min / C_max

    # Net Transfer Unit
    NTU = Uo*Ao/C_min

    # Effectiveness
    sq = math.sqrt(1 + c**2)
    ex = math.exp(-NTU * sq)
    eff = 2*(1 + c + sq * ((1 + ex)/(1 - ex)))**-1

    # Maximum Temperature Difference
    T_max = max(Ts_in, Ts_out, Tt_in, Tt_out)
    T_min = min(Ts_in, Ts_out, Tt_in, Tt_out)
    delta_T_max = T_max - T_min

    # Maximum Heat Duty
    Q_max = C_min*delta_T_max

    # Actual Heat Duty
    Q = eff * Q_max
    
    # Outlet Temperatures
    Ts_out_new = Ts_in - Q / Cs
    Tt_out_new = Tt_in + Q / Ct
    
    # Check for convergence
    if abs(Ts_out_new - Ts_out) < convergence_limit and abs(Tt_out_new - Tt_out) < convergence_limit and abs(Uo_new - Uo) < convergence_limit:
        break

    # Update Ts_out and Tt_out for the next iteration
    Ts_out = Ts_out_new
    Tt_out = Tt_out_new
    
    # Update Uo for the next iteration
    Uo = Uo_new

    iterations += 1

# Display the results with proper units
if red_flag == "on":
    print(f"Faulty Design. Follow the instructions above.")
else:
    print(f"-----Complete Specification-----")
    print(f"Overall Heat Transfer Coefficient (Uo): {math.ceil(Uo)} W/m^2.K")
    print(f"Heat Duty (Q): {Q/1E3:.2f} KW")
    
    print("\n--Shell Side Specification--")
    print(f"Inlet Temperature (Ts_in): {Ts_in} °C")
    print(f"Outlet Temperature (Ts_out): {Ts_out:.2f} °C")
    print(f"Mass Flow Rate (m_s): {m_s:.2f} Kg/s")
    print(f"Fluid Velocity (Vs): {Vs:.2f} m/s")
    print(f"Pressure Drop (ΔPs): {delta_Ps/1000:.4f} kPa")
    print(f"Inside Diameter (Ds): {m2in(Ds):.2f} inch")
    print(f"Hydraulic Diameter (Dhs): {m2in(Dhs):.2f} inch")
    print(f"Equivalent Diameter (Des): {m2in(Des):.2f} inch")
    print(f"Number of Baffles (Nb): {math.ceil(Nb)}")
    print(f"Baffle Spacing (B): {m2in(B):.2f} inch")
    print(f"Bundle Cross-flow Area (As): {sqm2sqft(As):.4f} ft^2")
    print(f"Pitch Length (Pt): {m2in(Pt):.2f} inch")
    print(f"Clearance Length (C): {m2in(C):.2f} inch")
    print(f"Convectional Heat Transfer Coefficient (ho): {math.ceil(ho)} W/m^2.K")
    print(f"Reynold's Number (Re_s): {math.ceil(Re_s)}")
    print(f"Nusselt Number (Nu_s): {math.ceil(Nu_s)}")
    print(f"Flow Type: {flow_type_s}")
    print(f"Friction Coefficient (fs): {fs:.2f}")

    print("\n--Tube Side Specification--")
    print(f"Inlet Temperature (Tt_in): {Tt_in} °C")
    print(f"Outlet Temperature (Tt_out): {Tt_out:.2f} °C")
    print(f"Mass Flow Rate (m_t): {m_t:.2f} Kg/s")
    print(f"Fluid Velocity (Vt): {Vt:.2f} m/s")
    print(f"Pressure Drop (ΔPt): {delta_Pt/1000:.4f} kPa")
    print(f"Outer Diameter (ODt): {m2in(ODt):.2f} inch")
    print(f"Inner Diameter (IDt): {m2in(IDt):.2f} inch")
    print(f"Hydraulic Diameter (Dht): {m2in(Dht):.2f} inch")
    print(f"Equivalent Diameter (Det): {m2in(Det):.2f} inch")
    print(f"Number of Tubes (Nt): {math.ceil(Nt)}")
    print(f"Tube Passes (Np): {Np}")
    print(f"Tube Length (L): {m2ft(L):.2f} feet")
    print(f"Heat Transfer Area (Ao): {sqm2sqft(Ao):.4f} ft^2")
    print(f"Total Inner Surface Area (Ai): {sqm2sqft(Ai):.4f} ft^2")
    print(f"Total x-section Area (At): {sqm2sqft(At):.4f} ft^2")
    print(f"Pitch Ratio (PR): {PR:.2f}")
    print(f"Tube Layout = 90°")
    print(f"Convectional Heat Transfer Coefficient (hi): {math.ceil(hi)} W/m^2.K")
    print(f"Reynold's Number for Heat Transfer (Re_t): {math.ceil(Re_t)}")
    print(f"Reynold's Number for Pressure Drop (Re_t_p): {math.ceil(Re_t_p)}")
    print(f"Nusselt Number (Nu_t): {math.ceil(Nu_t)}")
    print(f"Flow Type: {flow_type_t}")
    print(f"Friction Coefficient (ft): {ft:.2f}")


