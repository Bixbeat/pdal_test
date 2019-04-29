import os
import itertools as it
import numpy as np

import pdal
import shapely as shp
from osgeo import ogr, osr, gdal
from shapely.geometry import MultiPoint, mapping
from shapely.ops import cascaded_union

from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit

class PointCloudArrayConverter():
    """Provides functions for converting raw PDAl point cloud clustering outputs"""
    def __init__(self, point_arrays):
      self.arrays = point_arrays
      self.cluster_polys = None

    def clusters_to_polygons(self, method='bbox'):
        array_sorted_by_clusterid = sorted(enumerate(self.arrays[0]), key=self._key_function)

        self.cluster_polygons = []

        for key, group in it.groupby(array_sorted_by_clusterid, key=self._key_function):
            if not key == 0:
                coord_pairs = self._coord_pairs_from_cluster_pts(list(group))
                grouped_multipts = MultiPoint([pair[1] for pair in coord_pairs])

                intensity_stdev = int(np.std([x[0] for x in coord_pairs]))
                intensity_mean = int(np.mean([x[0] for x in coord_pairs]))
                color_mean = int(np.mean([x[0] for x in coord_pairs]))
                num_pts = len(coord_pairs)

                metadata = [key, num_pts, intensity_mean, intensity_stdev]
                if method == 'bbox':
                    geom = grouped_multipts.minimum_rotated_rectangle
                elif method == 'hull':
                    geom = grouped_multipts.convex_hull
                elif method == 'point':
                    geom = grouped_multipts

                self.cluster_polygons.append([metadata, geom])

    def _coord_pairs_from_cluster_pts(self, grouped_cluster_pts):
        return list(map(lambda x: [(x[1][3], (x[1][0], x[1][1])], grouped_cluster_pts))

    def write_output(self, out_fname, crs=None, overwrite=False):
        driver = ogr.GetDriverByName('Esri Shapefile')
        field_names = ['clusterid', 'num_pts', 'avg_int', 'std_int']

        if overwrite or not os.path.isfile(f'{out_fname}.shp'):
            ds = driver.CreateDataSource(f'{out_fname}.shp')
            dest_srs = osr.SpatialReference()
            dest_srs.ImportFromEPSG(crs)
            layer = ds.CreateLayer('', dest_srs)

            for field_name in field_names:
                layer.CreateField(ogr.FieldDefn(field_name, ogr.OFTInteger))

        else:
            ds = driver.Open(f'{out_fname}.shp', 1)
            layer = ds.GetLayer()            

        defn = layer.GetLayerDefn()
        
        for array in self.cluster_polygons:
            feat = ogr.Feature(defn)
            for i, info in enumerate(array[0]):
                feat.SetField(field_names[i], str(info))

            geom = ogr.CreateGeometryFromWkb(array[1].wkb)
            feat.SetGeometry(geom)
            layer.CreateFeature(feat)            

            feat = geom = None 
        ds = layer = feat = geom = None

    def _key_function(self, x):
      return x[1][-1]

if __name__ == '__main__':
    import glob
    filepaths = glob.glob('./data/*.laz')
    # Trim leading dot, Convert double backslashes to forward slash
    filepaths = [f[1:].replace('\\', '/') for f in filepaths]

    for filename in filepaths[-3:]:

        pipeline=f"""{{
        "pipeline": [
          {{
            "type":"readers.las",
            "filename":"C:/Users/d14878/Downloads/scripting/scripts/a4_pcl{filename}"
          }},
          {{
            "type":"filters.cluster",
            "min_points":"100",
            "tolerance":"0.05"
          }}
        ]
        }}"""
        r = pdal.Pipeline(pipeline)
        r.validate()
        r.execute()
        arrays = r.arrays

        pc_converter = PointCloudArrayConverter(arrays)
        pc_converter.clusters_to_polygons(method='point')
        pc_converter.write_output(out_fname='results', crs=28992, overwrite=False)