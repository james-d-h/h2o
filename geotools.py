
## Helper functions

# Define a function to generate random points within 250m of the original coordinates
def maskLand(image): # Update image mask with land mask
    image2 = image.updateMask(LandMask)
    return image2

def getSamplePointsLonLat(feature): # Return Feature Collection of XY Table Data 
    return ee.Feature(ee.Geometry.Point(coords= [feature.getNumber('LONGITUDE'), feature.getNumber('LATITUDE')]), 
                      {'IDn': feature.getString('IDn'), 
                       'SP': feature.getString('SamplePoint'), 
                       'DataSource' : feature.getString('DataSource')});

def WBsAOItoPoints(feature): # Buffer an AOI around the Sample Points
    return feature.setGeometry(feature.get('RandomSamples'))

def getArea(f): # Get area of each polygon 
    return f.set('Area', f.area(1)); 

def getEdge(polygon): 
    edgeCoords = polygon.geometry().coordinates().flatten();
    edge = ee.Geometry.LinearRing(edgeCoords);
    return polygon.setGeometry(edge);

def export_as_csv(collection, folder, asset_id, description):
    """
    Export an ee.FeatureCollection as an Earth Engine asset.
    """
    export_params = {
        'description' : description,
        'driveFolder' : folder,     # Folder in your Google Drive where the CSV file will be saved
        'fileNamePrefix': asset_id,   # Prefix for the CSV file
        'fileFormat'   : 'CSV',      # Format of the exported file
        'selectors'   : []          # Optional list of property names to export. If empty, all properties will be exported.
        }

    task = ee.batch.Export.table.toDrive(collection=collection, **export_params)
    task.start()
    return task

def export_as_asset(collection, asset_id, description):
    """Export an ee.FeatureCollection as an Earth Engine asset."""
    task = ee.batch.Export.table.toAsset(
        collection=collection,
        description=description,
        assetId=asset_id
        )

    task.start()
    return task
