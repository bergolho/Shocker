# ===========================================================================================================
# This script calibrates the diameter of the Purkinje cells in order to adjust the conduction 
# velocity to the reference value used for the cable equation.
# ===========================================================================================================

import sys
import numpy as np

Gi = 7.9
Cf = 3.4
tauf = 0.1

def calc_diameter (v):
    return (v*v*4.0*Cf*tauf) / (Gi) * 100.0

def calc_velocity (d):
    return ((Gi*d)/(4.0*Cf*tauf) )**(0.5) * 0.1

def calc_proportion (d,ref_lat,aprox_lat):
    return d*aprox_lat/ref_lat

def main ():
    if len(sys.argv) != 2:
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <input_reference_cv>" % sys.argv[0])
        print("-------------------------------------------------------------------------")
        print("<input_reference_cv> = Input value for the conduction velocity (CV)")
        print("-------------------------------------------------------------------------")
        return 1

    # Rafa Sebastian paper {m/s}
    #ref_velocity = 1.9
    # T-wave personalisation paper {m/s}
    #ref_velocity = 3.0

    ref_velocity = float(sys.argv[1])

    ref_diameter = calc_diameter(ref_velocity)
    print("Reference diameter = %g" % (ref_diameter))

if __name__ == "__main__":
    main()

