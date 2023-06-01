import argparse

import geopandas
import numpy as np
import pandas
from shapely.geometry.multipolygon import MultiPolygon

def main():

    # Parsing command-line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("path_gpkg_in", help = "File path to input GeoPackage file.")
    parser.add_argument("path_gpkg_out", help = "File path to output GeoPackage file.")
    #
    args = parser.parse_args()
    #
    path_gpkg_in  = args.path_gpkg_in
    path_gpkg_out = args.path_gpkg_out
    
    # Read contours from the GeoPackage.
    contours = geopandas.read_file(path_gpkg_in)

    # Calculate the area enclosed by each contour.
    contours["area"] = contours["geometry"].area

    # Find the unique, sorted list of contour values.
    contour_values = np.sort(contours['ELEV'].unique())
    n_contour_values = len(contour_values)

    # Loop over all the different contour values.
    contours_out_list = []
    #contour_values = [-8200.0]
    for k, contour_value in enumerate(contour_values):

        # Select the contours with a specific elevation value.
        contours_selected = contours.loc[contours['ELEV'] == contour_value]

        # Sort the selected contours in order of increasing area.
        # We are looking for contours which are contained by other
        # contours, which is only possible when the containing contour
        # is larger.
        contours_sorted = contours_selected.sort_values('area')
        
        # Generate a spatial index.
        spatial_index = contours_sorted.sindex
        
        # Count the contours in this group and prepare the output array.
        n_contours = len(contours_sorted)
        #parent_index = np.zeros(n_contours, dtype = int) - 1
        parent_index_list = np.zeros(n_contours, dtype = object)
        has_parent = np.zeros(n_contours, dtype = bool)

        # Loop over all pairs of contours within this group to find
        # contours which are contained by other contours.
        print("Processing contour group {:>4d} of {:>4d} with elevation {:>9.3f}, containing {:>6d} contours".format(
                k + 1, n_contour_values, contour_value, n_contours))
        for i in range(n_contours):

            parent_index_list[i] = []
            
            # The limits of the second loop index are chosen to avoid
            # checking containment by smaller contours, self-comparisons,
            # and duplicate comparisons.
            for j in range(i + 1, n_contours):
                
                # Use the spatial index to quickly identify possible 
                # containments.
                i_is_possible_child_of_j = spatial_index.geometries[i].intersects(spatial_index.geometries[j])
                if i_is_possible_child_of_j:

                    # In the case of possible containment, do a full
                    # calculation to check.
                    i_is_child_of_j = contours_sorted.iloc[j]['geometry'].contains(contours_sorted.iloc[i]['geometry'])
                    if i_is_child_of_j:

                        # Store j as one of the parents of i.
                        parent_index_list[i].append(j)
                        has_parent[i] = True

        # Prepare an output array for the "fixed" contours.
        contours_fixed = contours_sorted.copy()
         
        # "Fix" the parent contours by removing their child contours as holes.
        for i in range(n_contours):
            
            # Contours with no parent do not need to be considered.
            if has_parent[i]:
                
                # Get the index of the parent contour.
                # In the case where [A contains B] and [B contains C], the 
                # contour C will have two parents (A and B). But we only
                # need to consider the first (direct) parent, which is the
                # smallest, i.e. contour B. The indirect parent will be
                # treated in a later iteration of the loop,
                # when we get to contour B.
                j = parent_index_list[i][0]

                # Update the geometry: it should be the larger contour, with
                # the smaller contour removed as hole using the "difference"
                # operator.
                # In the case of multiple containment (discussed above),
                # holes-within-holes become solid again through the repeated
                # action of the the difference() operator.
                new_geometry = contours_fixed.iloc[j]['geometry'].difference(contours_fixed.iloc[i]['geometry'])

                #if not new_geometry.geom_type == 'MultiPolygon':

                #    new_geometry = MultiPolygon([new_geometry])
                
                # Update the geometry of the parent polygon.
                # Here the use of get_loc() avoids "chained indexing" and
                # "SettingWithCopyWarning" assignment of new geometry.
                # Also the new value is cast as a GeoSeries to avoid a bug in
                # GeoPandas.
                # https://github.com/geopandas/geopandas/issues/992#issuecomment-626138400
                contours_fixed.iloc[[j], contours.columns.get_loc('geometry')] = \
                    geopandas.GeoSeries([new_geometry])

        # Remove all of the child contours, which have been turned into holes.
        to_drop = np.where(has_parent)[0]
        contours_fixed.drop(contours_fixed.index[to_drop], inplace = True)

        # Append the contours for this level to the overall list.
        contours_out_list.append(contours_fixed)

    # Merge the contours into a single list.
    contours_out = geopandas.GeoDataFrame(pandas.concat(contours_out_list,
                    ignore_index = True),
                    crs = contours.crs)

    # Save the edited contours.
    print("Writing to {:}".format(path_gpkg_out))
    contours_out.to_file(path_gpkg_out, driver = "GPKG")

    return

if __name__ == '__main__':

    main()
