import csv
import datetime
import numpy as np


def load_csv_data(csv_file_path):
    """Loads data from the CSV file"""
    try:
        with open(csv_file_path, newline="") as csvfile:
            reader = csv.reader(csvfile, delimiter=",")  # Specify ';' as the delimiter
            data = list(reader)

            if data:
                data_matrix = [row for row in data]  
                return data_matrix
            
    except Exception as e:
                    print(f"Error loading csv")

data = load_csv_data("assets/testset.csv")

    
def getNumberOfPoints(callsign):
    """Function to get the number of data points for a specific flight"""
    count = 0

    for row in data:
        if row[0] == callsign:
              count += 1
    
    return count
       
def getPercentageAboveAltitude(callsign, altitude):
    """Function get the percentage of points for which a specific flight is above a certain altitude"""
    count_rows = 0
    count_alt = 0

    for row in data:
        if row[0] == callsign:
              count_rows += 1
              if float(row[4]) >= altitude:
                   count_alt += 1
    
    return (count_alt/count_rows)*100

      
def findDescentDateTime(callsign):
    """Function to find the date and time when a specific flight starts descending """ 
    altitude = 0
    for row in data:
        if row[0] == callsign:
            if altitude > int(row[4]):
                date = datetime.datetime.fromtimestamp(int(row[1])/1000) # “TIME” units: MILISECONDS since Unix Epoch.
                return date.strftime('%Y-%m-%d %H:%M:%S')
                break
            else:
                 altitude = int(row[4])
    return None

# Constants
A = 6378137.0  # Semi-major axis in meters
E2 = 0.00669437999014  # Eccentricity squared for WGS84
B = 6356752.3142  # Semi-minor axis in meters

def get_translation_vector(lat, lon, alt):
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    nu = A / np.sqrt(1 - E2 * np.sin(lat_rad) ** 2)
    return np.array(
        [
            [(nu + alt) * np.cos(lat_rad) * np.cos(lon_rad)],
            [(nu + alt) * np.cos(lat_rad) * np.sin(lon_rad)],
            [(nu * (1 - E2) + alt) * np.sin(lat_rad)],
        ]
    )

def get_rotation_matrix(lat, lon):
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    return np.array(
        [
            [-np.sin(lon_rad), np.cos(lon_rad), 0],
            [
                -np.sin(lat_rad) * np.cos(lon_rad),
                -np.sin(lat_rad) * np.sin(lon_rad),
                np.cos(lat_rad),
            ],
            [
                np.cos(lat_rad) * np.cos(lon_rad),
                np.cos(lat_rad) * np.sin(lon_rad),
                np.sin(lat_rad),
            ],
        ]
    )

def system_cartesian_to_system_stereographical(c):
    class CoordinatesUVH:
        def __init__(self):
            self.U = 0
            self.V = 0
            self.Height = 0

    res = CoordinatesUVH()
    center = {"Lat": 41.10904, "Lon": 1.226947, "Alt": 3438.954}

    lat_rad = np.radians(center["Lat"])

    R_S = (A * (1.0 - E2)) / (1 - E2 * np.sin(lat_rad) ** 2) ** 1.5

    d_xy2 = c["X"] ** 2 + c["Y"] ** 2
    res.Height = np.sqrt(d_xy2 + (c["Z"] + center["Alt"] + R_S) ** 2) - R_S

    k = (2 * R_S) / (2 * R_S + center["Alt"] + c["Z"] + res.Height)
    res.U = k * c["X"]
    res.V = k * c["Y"]

    return {"U": res.U, "V": res.V, "Height": res.Height}


def geodesic_to_geocentric(lat, lon, alt):
    """
    Converts geodetic coordinates (latitude, longitude, height) to geocentric (x, y, z).
    """

    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    nu = A / np.sqrt(1 - E2 * np.sin(lat_rad) ** 2)

    x = (nu + alt) * np.cos(lat_rad) * np.cos(lon_rad)
    y = (nu + alt) * np.cos(lat_rad) * np.sin(lon_rad)
    z = (nu * (1 - E2) + alt) * np.sin(lat_rad)
    return np.array([x, y, z])


def geocentric_to_system_cartesian(geocentric_coords):
    geo = {
        "X": geocentric_coords[0],
        "Y": geocentric_coords[1],
        "Z": geocentric_coords[2],
    }
    center = {"Lat": 41.10904, "Lon": 1.226947, "Alt": 3438.954}
    R = get_rotation_matrix(center["Lat"], center["Lon"])
    T = get_translation_vector(center["Lat"], center["Lon"], center["Alt"])

    input_vector = np.array([[geo["X"]], [geo["Y"]], [geo["Z"]]])
    result_vector = R @ (input_vector - T)

    return {
        "X": result_vector[0, 0],
        "Y": result_vector[1, 0],
        "Z": result_vector[2, 0],
    }


def get_stereographical_from_lat_lon_alt(lat, lon, alt):
    geocentric_coords = geodesic_to_geocentric(lat, lon, alt)
    cartesian_coords = geocentric_to_system_cartesian(geocentric_coords)

    stereographical_coords = system_cartesian_to_system_stereographical(
        cartesian_coords
    )
    print(stereographical_coords)

    return stereographical_coords

def calculate_distance(U1, V1, U2, V2):
    distance = np.sqrt((U1 - U2) ** 2 + (V1 - V2) ** 2) /1000
    return distance

def getDistanceBetweenPoints(callsign):
    """Function to find the distance (in km) between the first and the last point of a specific flight """
    count = 0
    for row in data:
        if row[0] == callsign:
              if count == 0:
                stereographical_coords1 = get_stereographical_from_lat_lon_alt(float(row[2]), float(row[3]), float(row[4]))
              count +=1
              if count == int(getNumberOfPoints(callsign)):
                stereographical_coords2 = get_stereographical_from_lat_lon_alt(float(row[2]), float(row[3]), float(row[4]))
                break

    distance = calculate_distance(float(stereographical_coords1["U"]), float(stereographical_coords1["V"]), float(stereographical_coords2["U"]), float(stereographical_coords2["V"]))
    return distance

print("The flight RYR7619 has", getNumberOfPoints("RYR7619"), "data points.")
print("A", round(getPercentageAboveAltitude("RYR7619", altitude = 10000),2), "% of the data points of the flight RYR7619 have an altitude above 10000.")
print("The flight RYR7619 started descending on",  findDescentDateTime("RYR7619"))