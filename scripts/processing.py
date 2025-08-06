import argparse
import pandas as pd
from pathlib import Path
import country_converter as coco

def get_enabled_areas():
    """Reads the areas.csv file and returns a list of enabled area codes."""
    areas_df = pd.read_csv("areas.csv")
    return areas_df[areas_df["enabled"]]["area_code"].tolist()

AREAS = get_enabled_areas()

# dictionary for special cases where coco fails
special_cases = {
    'DE_LU': 'Germany',
    'IE_SEM': 'Ireland',
    'IT_CNOR': 'Italy',
    'IT_CSUD': 'Italy',
    'IT_NORD': 'Italy',
    'IT_SARD': 'Italy',
    'IT_SICI': 'Italy',
    'IT_SUD' : 'Italy',
    'NO_1': 'Norway',
    'NO_2': 'Norway',
    'NO_3': 'Norway',
    'NO_4': 'Norway',
    'NO_5': 'Norway',
    'SE_1': 'Sweden',
    'SE_2': 'Sweden',
    'SE_3': 'Sweden',
    'SE_4': 'Sweden',
    'DK_1': 'Denmark',
    'DK_2': 'Denmark',
}

def area_to_country(area):
    if area in special_cases:
        return special_cases[area]
    return coco.convert(names=area, to="name_short")

area_to_country_dict = {area: area_to_country(area) for area in AREAS}

def process_data(prices_file, load_file, vre_file, output_file):
    """
    Processes raw price, load, and VRE data into a single dataframe.
    """
    print("Processing data...")
    df_price = pd.read_feather(prices_file)
    df_power = pd.read_feather(load_file)
    df_vre = pd.read_feather(vre_file)

    # Resample everything to hourly frequency, as some data might have higher resolution
    df_price = df_price.resample('h').mean()
    df_power = df_power.resample('h').mean()
    df_vre = df_vre.resample('h').mean()

    # Create MultiIndex columns for merging
    df_p = df_price.copy()
    index_p = pd.MultiIndex.from_tuples([(col, 'price') for col in df_p.columns])
    df_p.columns = index_p

    df_d = df_power.copy()
    index_d = pd.MultiIndex.from_tuples([(col, 'demand') for col in df_d.columns])
    df_d.columns = index_d


    df_v = df_vre.copy()

    # Combine the dataframes
    df_EU = pd.concat([df_d, df_p, df_v], axis=1)

    # Process each area
    dfs_EU = {}
    for area in [col for col in df_EU.columns.get_level_values(0).unique() if col in AREAS]:
        # Select data for the current area, resulting in a DataFrame with single-level columns
        df_area_data = df_EU[area].copy()

        vre_new_names = {'Solar': 'solar', 'Wind Onshore': 'wind_onshore', 'Wind Offshore': 'wind_offshore'}
        df_area_data.rename(columns=vre_new_names, inplace=True)
        
        # Perform the assign operations
        df_area_processed = (df_area_data
            .sort_index(axis=1)
            .assign(hour       = lambda x: x.index.hour)
            .assign(doy        = lambda x: x.index.dayofyear)
            .assign(doy_season = lambda x: ((x.doy + 8*30+4)%365)+1)
            .assign(dt         = lambda x: x.index)
            .assign(wind       = lambda x: x.get('wind_onshore', 0) + x.get('wind_offshore', 0))
            .assign(vre        = lambda x: x.get('wind', 0) + x.get('solar', 0))
            .assign(residual   = lambda x: x.get('demand', 0) - x.get('vre', 0))
            .ffill()
            .interpolate()
        )
        
        # Create a new MultiIndex for the columns of the processed area data
        # The first level is the area, the second level are the new and old metrics
        new_columns = pd.MultiIndex.from_product([[area], df_area_processed.columns])
        df_area_processed.columns = new_columns
        
        dfs_EU[area] = df_area_processed

    # Combine all areas into a single dataframe horizontally
    df_EU_all = pd.concat(dfs_EU.values(), axis=1)
    
    # The index is already the timestamp, no need to reset or set names
    # The 'country' column can be derived from the MultiIndex if needed later

    
    # Save the processed data
    df_EU_all.to_feather(output_file)
    print(f"Successfully saved processed data to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process ENTSO-E data.")
    parser.add_argument("prices_file", type=Path, help="Path to the prices feather file.")
    parser.add_argument("load_file", type=Path, help="Path to the load feather file.")
    parser.add_argument("vre_file", type=Path, help="Path to the VRE feather file.")
    parser.add_argument("output_file", type=Path, help="Path to save the processed output feather file.")
    
    args = parser.parse_args()
    
    process_data(args.prices_file, args.load_file, args.vre_file, args.output_file)
