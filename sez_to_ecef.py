# ecef_to_llh.py
#
# Usage: python3 sez_to_ecef.py o_lat_deg o_lon_deg o_hae_km s_km e_km z_km
# Converts ECEF vector components to LLH

# Parameters:
#o_lat_deg SEZ latitude in deg
#o_lon_deg SEZ longitude in deg
#o_hae_km SEZ hae in km
#s_km SEZ s in km
#e_km SEZ e in km
#z_km SEZ z in km
# Output:
# Prints the converged latitude (deg), longitude (deg), and HAE (km)
#
# Written by Samuel Jacobson
# Other contributors:
#
#
# import Python modules
import math
import sys

# Constants
R_E_KM = 6378.137
E_E = 0.081819221456


# Parse script arguments
if len(sys.argv) == 7:
    o_lat_deg = float(sys.argv[1])
    o_lon_deg = float(sys.argv[2])
    o_hae_km = float(sys.argv[3])
    s_km = float(sys.argv[4])
    e_km = float(sys.argv[5])
    z_km = float(sys.argv[6])
else:
    print('Usage: python3 sez_to_ecef.py o_lat_deg o_lon_deg o_hae_km s_km e_km z_km')
    exit()

# Perform SEZ to ECEF conversion
# SEZ to ECEF conversion

# Convert observer's lat/lon to radians
lat_rad = o_lat_deg * math.pi / 180.0
lon_rad = o_lon_deg * math.pi / 180.0

# Trig to keep things clean
sin_lat = math.sin(lat_rad)
cos_lat = math.cos(lat_rad)
sin_lon = math.sin(lon_rad)
cos_lon = math.cos(lon_rad)

#  Transformation matrix

ecef_x_km = (cos_lon * sin_lat * s_km) + (cos_lon * cos_lat * z_km) - (sin_lon * e_km)
ecef_y_km = (sin_lon * sin_lat * s_km) + (sin_lon * cos_lat * z_km) + (cos_lon * e_km)
ecef_z_km =  - (cos_lat * s_km) + (sin_lat * z_km)

# Add the ECEF vector to the SEZ origin

lat_rad = o_lat_deg * math.pi / 180.0
lon_rad = o_lon_deg * math.pi / 180.0

# Calculate N, the prime vertical radius of curvature
def calc_denom(ecc, lat_rad):
    return math.sqrt(1.0 - (ecc ** 2) * (math.sin(lat_rad) ** 2))

denom = calc_denom(E_E, lat_rad)
c_E = R_E_KM / denom

# Calculate ECEF coordinates for the observer
obs_x_ecef = (c_E + o_hae_km) * cos_lat * cos_lon
obs_y_ecef = (c_E + o_hae_km) * cos_lat * sin_lon
obs_z_ecef = (c_E * (1 - E_E ** 2) + o_hae_km) * sin_lat

#Add observer's ECEF coordinates to the transformed SEZ coordinates
ecef_x_km += obs_x_ecef
ecef_y_km += obs_y_ecef
ecef_z_km += obs_z_ecef


# Print results
print(ecef_x_km)
print(ecef_y_km)
print(ecef_z_km)

