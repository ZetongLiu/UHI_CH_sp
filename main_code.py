# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 17:16:28 2023

@author: Olivier Chavanne


Modified on Sat Jan 20 11:17:47 2024
@author: Zetong Liu

"""

import geopandas as gpd
import pandas as pd
from shapely import box
import os
import matplotlib.pyplot as plt

# Local libraries
from uhiCAD.building import generate_envelope
from uhiCAD.building import generate_buildings
import uhiCAD.xml as xml
# import enerCAD.result as result
# import enerCAD.network as network
# import uhiCAD.production as prod
# import uhiCAD.KPI as KPI

# URL for RegBL API request
GEOADMIN_BASE_URL = "https://api.geo.admin.ch/rest/services/ech/MapServer/ch.bfs.gebaeude_wohnungs_register/"
    
##################################################
# 
#                  Functions
#
##################################################

def create_xml_root(xml_file_to_copy, climate_file, horizon_file):
    '''
    Parameters                                                          
    ----------
    xml_file_to_copy : TYPE
        DESCRIPTION.
    climate_file : TYPE
        DESCRIPTION.
    horizon_file : TYPE
        DESCRIPTION.

    Returns
    -------
    root : TYPE
        DESCRIPTION.
    district : TYPE
        DESCRIPTION.
    '''
    
    # Write XML file for CitySim :
    print("Writing XML file...")    
    # Add Root 
    root = xml.add_root()
    # Add Simulation days
    xml.add_simulation_days(root)
    # Add Climate
    xml.add_climate(root, climate_file)
    # Add District
    district = xml.add_district(root)
    
    # Horizon
    # read in the tab-separated file as a dataframe
    horizon_df = pd.read_csv(horizon_file, sep='\s+', header=None)
    # assign column names to the dataframe
    horizon_df.columns = ['phi', 'theta']
    # Add Far field obstructions
    xml.add_far_field_obstructions(district, horizon_df)
    
    # Add all the composites and profiles, taken from a source XML
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'Composite')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'OccupancyDayProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'OccupancyYearProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DeviceType')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'ActivityType')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DHWDayProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DHWYearProfile')
    
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'Building')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DistrictEnergyCenter')
    
    print("Xml source copied")
    
    return root, district 

def Module_1(gpkg_filepath,XYZfile, GEOADMIN_BASE_URL,
             directory_path, xml_name,
             xml_base_file, climate_file, horizon_file,
             create_geometry_3D=False, calculate_volume_3D=False,
             EGID_column='RegBL_EGID'):
    '''
    Parameters
    ----------
    gpkg_filepath : TYPE
        DESCRIPTION.
    GEOADMIN_BASE_URL : TYPE
        DESCRIPTION.
    directory_path : TYPE
        DESCRIPTION.
    xml_file_to_create : TYPE
        DESCRIPTION.
    xml_base_file : TYPE
        DESCRIPTION.
    climate_file : TYPE
        DESCRIPTION.
    horizon_file : TYPE
        DESCRIPTION.
    create_geometry_3D : TYPE, optional
        DESCRIPTION. The default is False.
    calculate_volume_3D : TYPE, optional
        DESCRIPTION. The default is False.
    EGID_column : TYPE, optional
        DESCRIPTION. The default is 'RegBL_EGID'.

    Returns
    -------
    None.
    '''
    
    ### Exctract geopackage ###
    
    print("Exctracting geopackage layers...")
    
    # MO Cadaster
    MO_all = gpd.read_file(gpkg_filepath, layer = "zone_tout")
    MO_dhn = gpd.read_file(gpkg_filepath, layer = "zone_cad")
    centrale = gpd.read_file(gpkg_filepath, layer = "centrale")
    EGID_column = 'RegBL_EGID'
    
    # Split Multipolygons into Polygons
    zone_all = MO_all.explode(index_parts=False)
    zone_dhn = MO_dhn.explode(index_parts=False)
    
    # List containing EGID of buildings to simulate
    EGID_list = MO_dhn[EGID_column].tolist()
    
    # Save EGID list of buildings connected to CAD
    df_EGID = pd.DataFrame(EGID_list)
    df_EGID.columns = ['EGID']
    EGID_path = os.path.join(directory_path, 'EGID.csv')     
    df_EGID.to_csv(EGID_path, index=False)
    print("EGID.csv created")
    
    # Swissbuildings3D
    print("Swissbuildings3D processing...")
    try:
        floor_data = gpd.read_file(gpkg_filepath, layer = "floor")
        roof_data = gpd.read_file(gpkg_filepath, layer = "roof")
        wall_data = gpd.read_file(gpkg_filepath, layer = "wall")
        green_data = gpd.read_file(gpkg_filepath, layer = 'green')
        ori_street_data = gpd.read_file(gpkg_filepath, layer = 'streets')
        street_data = ori_street_data[~ori_street_data['objektart'].isin( ['Verbindung', 'Platz'])]
        
        # Filter on the zone with 10m buffer around surrounding square box 
        zone_bounds = MO_all.geometry.buffer(10).values.total_bounds
        zone_box = box(zone_bounds[0], zone_bounds[1], zone_bounds[2], zone_bounds[3])
        
        # Cut swissbuildings3D to zone of concern
        floor_data_intersection = floor_data[floor_data.geometry.intersects(zone_box)]
        roof_data_intersection = roof_data[roof_data.geometry.intersects(zone_box)]
        wall_data_intersection = wall_data[wall_data.geometry.intersects(zone_box)]

        # Split Multipolygons into Polygons
        zone_floor = floor_data_intersection.explode(index_parts=True).reset_index()
        zone_roof = roof_data_intersection.explode(index_parts=True).reset_index()
        zone_wall = wall_data_intersection.explode(index_parts=True).reset_index()
        print('Swissbuildings3D cut to zone of interest \n')
    
    except: print('Error : Swissbuildings3D not provided')

    ### Envelope processing ###
    
    try:
        # Get z coordinates of 1st vertex from 1st surface of 1st building's floor polygon as altitude by default for MO footprints
        altitude_default = zone_floor.loc[0].geometry.exterior.coords[0][2]
    except:
        altitude_default = 0
    
    # Create DataFrames containing all necessary information for each building
    print("Creating Buildings GeoDataFrame...")
    footprints, buildings = generate_buildings(zone_all, EGID_list, GEOADMIN_BASE_URL, altitude_default,
                                               create_geometry_3D, calculate_volume_3D, zone_floor, zone_roof, zone_wall)
    print("Buildings GeoDataFrame created \n") 
    
    # Generate the envelope surfaces
    print("Generating Buildings envelope...")
    envelope, buildings_volume_3D, center_coordinates = generate_envelope(footprints, buildings, calculate_volume_3D)
    print("Envelope created \n")
    
    # Merge "volume_3D" and "n_occupants" to main buildings geodataframe according to 'bid'
    merged_buildings = buildings.merge(buildings_volume_3D, left_on='bid', right_on='bid', how='left')    
    if not merged_buildings.empty:
        columns_to_add = ['volume_3D', 'n_occupants']
        for column in columns_to_add:
            buildings[column] = merged_buildings[column]
        print("Buildings 3D volume calculated and merged \n")
    
    ### Buildings XML processing ###
        
    root, district = create_xml_root(xml_base_file, climate_file, horizon_file)
    
    print("Adding buildings...")
    # Add the buildings
    xml.add_all_buildings(district, buildings, envelope, center_coordinates)
    terrain_df = pd.read_table(XYZfile, skiprows=1, delim_whitespace=True, names=['X', 'Y', 'Z'])

    # xml.add_ground_from_XYZ(district, terrain_df, center_coordinates)
    ground_data = xml.add_ground_cut(MO_dhn, district, terrain_df, zone_box, center_coordinates)
    road_index_list, green_index_list, _, _ = xml.modify_type(district, ground_data, green_data, street_data)

    #write xml file of default case
    print('creating xml file \n')

    xml.cut(district, ground_data, MO_dhn, footprints)
    xml_path = os.path.join(directory_path, xml_name+".xml")     
    xml.write_xml_file(root, xml_path) 
    print(f"{xml_name}.xml files created \n")

    # Write XML file of Scenario 1
    sc_id=1
    sc1_data = gpd.read_file(gpkg_filepath, layer = f'scenario{sc_id}')
    road_index_list, green_index_list, _, _  = xml.modify_type(district, ground_data, sc1_data, street_data)

    #approximation comparison
    print('approximation comparison \n') 
    selected_grounds1=ground_data[ground_data['gid'].isin(road_index_list)]
    selected_grounds2=ground_data[ground_data['gid'].isin(green_index_list)]
    # add layers of simulated green areas and streets
    selected_grounds1.to_file(gpkg_filepath, layer='street_grounds_2m', driver='GPKG')
    selected_grounds2.to_file(gpkg_filepath, layer='green_grounds_2m', driver='GPKG')

    scenario_path = os.path.join(directory_path, f"Scenario_{sc_id}")
    os.makedirs(scenario_path, exist_ok=True)
    xml_to_create_path = os.path.join(scenario_path, xml_name+f"_sc_{sc_id}"+".xml")
    xml.write_xml_file(root, xml_to_create_path)

    return envelope, ground_data, buildings, zone_dhn, centrale
   
def simulate_citysim(directory_path, xml_file, citysim_filepath):
    '''
    Parameters
    ----------
    xml_file : TYPE
        DESCRIPTION.

    Returns
    -------
    None.
    '''
    
    import subprocess
    import time
    start = time.time()
    print('Process started')
    print(f'Simulation of {xml_file}.xml...')

    #run CitySim.exe with xml file
    xml_path = os.path.join(directory_path, xml_file+".xml")
    result = subprocess.run([citysim_filepath, '-q', f"{xml_path}"])
    
    end = time.time()
    duration = end - start
    m, s = divmod(duration, 60)
    print('Simulation ended. Time :', "%.0f" %m,'min', "%.0f" %s,'s \n')
    

#------------------Part 2 iterating----------------------------------------------------------

def Module_2(directory_path, xml_name, gpkg_filepath, root, district,
             ground_data, climate_file, horizon_file,
             scenarios_list):   

    for i in range(len(scenarios_list)):
        sc_id = scenarios_list[i]
        sc1_data = gpd.read_file(gpkg_filepath, layer = f'scenario{sc_id}')
        root_copy, district_copy = root, district
        ### Scenarios XML processing ###
        # xml_to_copy_path = os.path.join(directory_path, xml_name+'.xml' )
        # root, district = create_xml_root(xml_to_copy_path, climate_file, horizon_file)
        road_index_list, green_index_list, _, _  = xml.modify_type(district_copy, ground_data, sc1_data)
        # Write XML file
        scenario_path = os.path.join(directory_path, f"Scenario_{sc_id}")
        os.makedirs(scenario_path, exist_ok=True)
        xml_to_create_path = os.path.join(scenario_path, xml_name+f"_sc_{sc_id}"+".xml")
        # xml.cut(district, ground_data, MO_dhn, footprints)

        xml.write_xml_file(root, xml_to_create_path)
        # print(f'{xml_DHN}_sc_{sc_id}.xml file created \n')

#--------------------- KPI calculation

def Module_KPI(ground_data, buffered_streets, itsctd_greens, road_index_list, green_index_list, scenarios, sc_id):

    real_str_area = buffered_streets['geometry'].area.sum()
    sim_str_area = ground_data.loc[ground_data['gid'].isin(road_index_list), 'geometry'].area.sum()
    str_error = abs(real_str_area-sim_str_area)/real_str_area
    real_gr_area = itsctd_greens['geometry'].area.sum()
    sim_gr_area = ground_data.loc[ground_data['gid'].isin(green_index_list), 'geometry'].area.sum()
    gr_error = abs(real_gr_area-sim_gr_area)/real_gr_area

    print('KPI calculated')
    
    return str_error,  gr_error


##################################################
# 
#         Information to provide
#
##################################################

# Geopackage filepath
gpkg_filepath = r"D:\Document\EPFL_Coursework\UHI_CH\dependent_files\lausanne_case.gpkg"                     #TODO

# Create geometry with swissbuildings3D
create_geometry_3D = True                                   #TODO

# Calculate volume from swissbuildings3D
calculate_volume_3D = True                                 #TODO

# CitySim.exe filepath
citysim_filepath = r"D:\Document\SemesterProject\CitySim.exe" #TODO

# XML name to export
directory_path = r"fountaine_Lausanne" #+f"_{Year_of_cli}"                                #TODO
os.makedirs(directory_path, exist_ok=True)                           
                                      
# XML source files
xml_base_file = r"D:\Document\SemesterProject\CAD-O-main\xml_base.xml"    #TODO                   
horizon_file = r"D:\CitySimPro\CitySimPro\Windows\Resources\climateFiles\Lausanne.hor"    #TODO     
XYZfile = r"D:\Document\SemesterProject\New_case\VD\study_area\Lausanne_alt\SWISSALTI3D_2_XYZ_CHLV95_LN02_2538_1152.xyz"  #TODO     
# XYZfile_fine = r"D:\Document\SemesterProject\New_case\VD\study_area\SWISSALTI3D_0.5_XYZ_CHLV95_LN02_2527_1151.xyz"
# Scenarios to simulate
scenarios_list = [1]                                #TODO

do_plot = True

def main(): 
    
    # Generate individual buildings XML
    print('***Module 1*** \n')
    Year_of_cli=['Contemporary', '2030', '2040'] 
    for year in Year_of_cli:
        subdirectory_path = os.path.join(directory_path, f"{year}")
        os.makedirs(subdirectory_path, exist_ok=True)     
        xml_name = directory_path+f'_{year}' 
        climate_file = rf"D:\Document\SemesterProject\New_case\VD\cli\Morges_{year}.cli"       #TODO  
        envelope, ground_data, buildings, zone_dhn, centrale = Module_1(gpkg_filepath, XYZfile, GEOADMIN_BASE_URL, 
                                                subdirectory_path, xml_name,
                                                xml_base_file, climate_file, horizon_file,
                                                create_geometry_3D, calculate_volume_3D,
                                                EGID_column='RegBL_EGID')
    
        # 1st CitySim simulation
        simulate_citysim(subdirectory_path, xml_name, citysim_filepath)
        
        TS_file= os.path.join(subdirectory_path, xml_name+"_TS.out")
        TS_df = pd.read_csv(TS_file, delimiter='\t')
        _, _, all_AST = xml.AST(envelope, TS_df, ground_data)
        all_AST.to_file(gpkg_filepath, layer=f'all_AST_{year}', driver='GPKG')

        # print('***Module 2*** \n')
        # Module_2(directory_path, xml_name, gpkg_filepath,
        #          ground_data, climate_file, horizon_file,
        #          scenarios_list)
        
        KPI_result_list = []
            
        # DHN simulation for each scenario
        for i in range(len(scenarios_list)):
            sc_id = scenarios_list[i]
            print(f'***Scenario {sc_id}*** \n')
            scenario_path = os.path.join(subdirectory_path, f"Scenario_{sc_id}")
            simulate_citysim(scenario_path, f'{xml_name}_sc_{sc_id}', citysim_filepath)
            scenario_TS_path= os.path.join(scenario_path, f'{xml_name}_sc_{sc_id}'+"_TS.out")
            scenario_TS_df = pd.read_csv(scenario_TS_path, delimiter='\t')
            _, _, scenario_all_AST = xml.AST(envelope, scenario_TS_df, ground_data)
            scenario_all_AST.to_file(gpkg_filepath, layer=f'all_AST_{year}_sc_{sc_id}', driver='GPKG')
            merged_df = gpd.GeoDataFrame(pd.merge(all_AST, scenario_all_AST, on='geometry', suffixes=('_dfT', '_s1')))
            merged_df['T_difference'] = merged_df['T_s1'] - merged_df['T_dfT']
            Dif_df = merged_df[['geometry', 'T_difference']]
            Dif_gdf = gpd.GeoDataFrame(Dif_df, geometry='geometry')
            Dif_gdf.to_file(gpkg_filepath, layer=f'all_AST_{year}_sc_{sc_id}_dif', driver='GPKG')       
                
            # KPI calculation
            # df_KPI = Module_KPI(results_production, volume_storage, 
            #                     scenarios, sc_id, scenario_path, do_plot)
            
            # KPI_result_list.append(df_KPI)

            print(f"Scenario {sc_id} at {year}processed \n")
        print(f"year {year} processed \n")
    print("***Overall processing finished***")
    print(f"Find all results and graphs in directory : {directory_path}")

if __name__ == "__main__":
    plt.close("all")
    
    import subprocess
    import time
    start_overall = time.time()
    print('Main code started')
    print('-----------------')
    
    main()    

    print('-----------------')
    print('Main code ended')
    print('-----------------')
    end_overall = time.time()
    duration_overall = end_overall - start_overall
    m, s = divmod(duration_overall, 60)
    print('Overall run time :', "%.0f" %m,'min', "%.0f" %s,'s \n')