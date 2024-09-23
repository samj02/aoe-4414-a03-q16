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
# This work is licensed under CC BY-SA 4.0
# import Python modules
import math
import sys

# Constants
R_E_KM = 6378.137  # Earth's equatorial radius in km
E_E = 0.081819221456  # Earth's eccentricity


# Helper function to convert degrees to radians
def deg_to_rad(deg):
    return deg * math.pi / 180.0


# Helper function to calculate the ECEF coordinates of the observer
def geodetic_to_ecef(lat_deg, lon_deg, hae_km):
    lat_rad = deg_to_rad(lat_deg)
    lon_rad = deg_to_rad(lon_deg)

    # Calculate N, the prime vertical radius of curvature
    N = R_E_KM / math.sqrt(1 - (E_E ** 2) * (math.sin(lat_rad) ** 2))

    # Calculate ECEF coordinates for the observer
    x = (N + hae_km) * math.cos(lat_rad) * math.cos(lon_rad)
    y = (N + hae_km) * math.cos(lat_rad) * math.sin(lon_rad)
    z = (N * (1 - E_E ** 2) + hae_km) * math.sin(lat_rad)

    return x, y, z


# SEZ to ECEF conversion
def sez_to_ecef(o_lat_deg, o_lon_deg, o_hae_km, s_km, e_km, z_km):
    # Convert observer's position to ECEF
    obs_x_ecef, obs_y_ecef, obs_z_ecef = geodetic_to_ecef(o_lat_deg, o_lon_deg, o_hae_km)

    # Convert observer's lat/lon to radians
    lat_rad = deg_to_rad(o_lat_deg)
    lon_rad = deg_to_rad(o_lon_deg)

    # SEZ to ECEF rotation matrix components
    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)

    # Apply SEZ to ECEF transformation
    ecef_x_km = (-sin_lon * s_km) + (cos_lon * e_km)
    ecef_y_km = (-sin_lat * cos_lon * s_km) + (-sin_lat * sin_lon * e_km) + (cos_lat * z_km)
    ecef_z_km = (cos_lat * cos_lon * s_km) + (cos_lat * sin_lon * e_km) + (sin_lat * z_km)

    # Add observer's ECEF coordinates to the transformed SEZ coordinates
    ecef_x_km += obs_x_ecef
    ecef_y_km += obs_y_ecef
    ecef_z_km += obs_z_ecef

    return ecef_x_km, ecef_y_km, ecef_z_km


# Parse script arguments
if len(sys.argv) == 7:
    o_lat_deg = float(sys.argv[1])
    o_lon_deg = float(sys.argv[2])
    o_hae_km = float(sys.argv[3])
    s_km = float(sys.argv[4])
    e_km = float(sys.argv[5])
    z_km = float(sys.argv[6])
else:
    print('Usage: \npython3 sez_to_ecef.py o_lat_deg o_lon_deg o_hae_km s_km e_km z_km')
    exit()

# Perform SEZ to ECEF conversion
ecef_x_km, ecef_y_km, ecef_z_km = sez_to_ecef(o_lat_deg, o_lon_deg, o_hae_km, s_km, e_km, z_km)

# Print results
print(ecef_x_km)
print(ecef_y_km)
print(ecef_z_km)

