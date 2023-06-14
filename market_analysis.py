import numpy as np
import json
import os
import math
from flight_performance import FlightPerformance




def mxacos(xyz):
    if xyz > 1:
        xyz = 1
    elif xyz < -1:
        xyz = -1
    return math.acos(xyz)


# Great Circle Distance
def gcdst(gps_from,gps_to):
    r = 6378.39
    conv = (math.pi / 180)
    lat1 = gps_from[0] * conv
    lon1 = gps_from[1] * conv
    lat2 = gps_to[0] * conv
    lon2 = gps_to[1] * conv
    return (r * mxacos((math.cos(lat1)*math.cos(lat2)*math.cos(lon2-lon1)) + (math.sin(lat1)*math.sin(lat2))))

def gcd(g1, g2):
    R = 6373.0
    lat1 = math.radians(g1[0])
    lon1 = math.radians(g1[1])
    lat2 = math.radians(g2[0])
    lon2 = math.radians(g2[1])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c

def read_json(file_name):
    with open(file_name) as json_file:
        dict = json.load(json_file)
    return dict

def get_levels(h_last, h_new):
    start = int((h_last-(h_last%10))+10)
    if h_last%10 == 0:
        start -= 10
    #print(start, h_new)
    return range(start, int(h_new), 10)

def extrapolate2(h1, h2, v1, v2):
    lst = []
    if (v1 < 900 and v2 < 900 and h1 != h2):
        dh = h2-h1
        dv = v2-v1
        grad = dv/dh
        start = int(h1-(h1%10)+10)
        if h1%10 == 0:
            start -= 10
        end = int(h2-(h2%10)+10)
        for i in range(start, end, 10):
            v = v1 + (i-h1)*grad
            lst.append([i, v])
    return lst

def extrapolate(data_h, data_wsp):
    wprof = {}
    max_alt = max(data_h) - (max(data_h)%10)
    alt = np.arange(0, max_alt+10, 10)
    wsp = np.zeros_like(alt, dtype=float)

    for i in range(1, len(data_h)):
        lst = extrapolate2(data_h[i-1], data_h[i], data_wsp[i-1], data_wsp[i])
        for l in lst:
            #print(l)
            idx = np.where(alt==l[0])[0][0]
            wsp[idx] = l[1]
    if data_h[0] >= 10:
        levels = list(np.arange(0, data_h[0], 10))
        for level in levels:
            idx = np.where(alt==level)[0][0]
            wsp[idx] = data_wsp[0]

    for i in range(np.size(alt)):
        wprof[str(alt[i])] = wsp[i]
    
    return wprof, max_alt

def read_RS92(file_txt):
    data = file_txt.split('\n')[194:]
    data = [lst for lst in data if lst]
    alt_lst = []
    wsp_lst = []
    first_gps = []
    last_gps = []
    for line in data:
        x = line.strip().split()
        #print(x)
        if len(x) > 5:
            alt, wsp, lon, lat = float(x[2]), float(x[7]), float(x[9]), float(x[10])
            if alt < 60000 and wsp < 200 and alt >= 0:
                alt_lst.append(alt)
                wsp_lst.append(wsp)
                if lat < 400 and lon < 400:
                    if not first_gps:
                        first_gps = [lat, lon]
                    last_gps = [lat, lon]
    if len(wsp_lst) > 0:
        windprofile, max_alt = extrapolate(alt_lst, wsp_lst)
        dst = gcdst(first_gps, last_gps)
        return windprofile, max_alt, dst, True
    else:
        return 0,0,0,False

def read_iMet(file_txt):
    data = file_txt.split('\n')[117:]
    data = [lst for lst in data if lst]
    alt_lst = []
    wsp_lst = []
    first_gps = []
    last_gps = []
    for line in data:
        x = line.strip().split()
        #print(x)
        if len(x) > 5:
            alt, wsp, lon, lat = float(x[2]), float(x[7]), float(x[9]), float(x[10])
            #print(alt, wsp, lat, lon)
            if alt < 60000 and wsp < 200 and alt >= 0:
                alt_lst.append(alt)
                wsp_lst.append(wsp)
                if lat < 400 and lon < 400:
                    if not first_gps:
                        first_gps = [lat, lon]
                    last_gps = [lat, lon]
    if len(wsp_lst) > 0:
        windprofile, max_alt = extrapolate(alt_lst, wsp_lst)
        dst = gcdst(first_gps, last_gps)
        return windprofile, max_alt, dst, True
    else:
        return 0,0,0,False

stations = ['boulder', 'macquarie', 'debilt', 'paramaribo', 'southpole']
spiral_height = {'boulder':3000, 'macquarie': 1500, 'debilt':1500, 'paramaribo':1500, 'southpole':3500}

AC_params = {"AR": 12, "e": 0.9, "m": 0.75, "S": 0.05, "CD_0": 0.03, "CL_max": 0.666 * 1.55,
             "CL_alpha": 6.1859591823509295 * np.pi / 180}

atmos_dict = read_json('atmospheric_characteristics.json')
res = {}

for station in stations:
    directory = station + '/'
    rng_lst = []
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            with open(f) as file:
                data = file.read()
            if 'RS92' in data:
                wp, h_burst, dst_drift, isnotempty = read_RS92(data)
                if (isnotempty):
                    FP = FlightPerformance(wp, AC_params, h_burst, spiral_height[station], 100)
                    FP.flight_sim()
                    flight_dict = FP.get_flightdict()
                    range_surplus = (flight_dict["distance_travelled"][-1] / 1000) - dst_drift
                    rng_lst.append(range_surplus)
            elif 'iMet-1' in data:
                wp, h_burst, dst_drift, isnotempty = read_iMet(data)
                if isnotempty:
                    FP = FlightPerformance(wp, AC_params, h_burst, 1000, 100)
                    FP.flight_sim()
                    flight_dict = FP.get_flightdict()
                    #print(h_burst, dst_drift)
                    range_surplus = (flight_dict["distance_travelled"][-1] / 1000) - dst_drift
                    rng_lst.append(range_surplus)
                    neg_count = len(list(filter(lambda x: (x < 0), rng_lst)))
                    pos_count = len(list(filter(lambda x: (x >= 0), rng_lst)))
                    print(filename, 'drift=', dst_drift, ' range=', (flight_dict["distance_travelled"][-1] / 1000), 'delta=', range_surplus)
                    print(station, ', pos:', pos_count, ', neg:', neg_count)

    res[station] = rng_lst
