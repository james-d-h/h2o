from .h2o.geotools import maskLand, getSamplePointsLonLat, WBsAOItoPoints, getArea, getEdge, export_as_csv, export_as_asset

#@title #GEE water detector and sampler
def WaterGeoAlgorithm(WaterOccurrence=1, sizeBuffer=300, dJoinWindow=100, bufferReduction=-30, save_directory):

    """
    # Set initial variables
    WaterOccurrence = 50  # [%]
    sizeBuffer       = 500 # [m]
    dJoinWindow     = 250 # [m]
    bufferReduction = -30  # [m]

    1. We flatten the JRC WaterOccurrence Image into polys (WaterOccurrence=1 implying max extent). 
    2. Where breaks less than dJoinWindow exist between jrc polys, we assumme they are connected. 
    3. Using latlon coords we query an area of sizeBuffer, where an intersection with jrc incurrs the spatial join.  
    4. We reduce the polygon of the jrc polys to minimise edge effects. 
    5. Finaally we sample 10 points at random within the union of the query/poly i.e. water.

    *Note the pixel value extractions are left to the user, 
     here we determine the latlon points to sample images downstream.

    """
    RawData  = ee.FeatureCollection((save_directory + 'watergeoPhD'));
    Scotland = ee.Geometry.Rectangle([-7.993835080698091,54.507760429600005, -0.3473507056980907,60.94063272856898]);
    h2o      = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
    LandMask = h2o.clip(Scotland).select('max_extent'); # Set LandMask
    
    # Fixed vars
    scale2use       = 30;  # [m] 
    Nsamples        = 10;  # [-] 
    errorMargin     = 10;  # [m]

    # Define a spatial and or temporal distance filter for join.
    SaveClosest = ee.Join.saveFirst(
            matchKey      = 'WBs',
            ordering      = 'distance',
            measureKey    = 'distance');

    withinRoughAOI = ee.Filter.withinDistance(
            distance      = sizeBuffer,
            leftField     = '.geo',
            rightField    = '.geo',
            maxError      = errorMargin);

    def set_var_optimisers(feature):
        return feature.set("bufferReduction", bufferReduction).set("dJoinWindow", dJoinWindow).set("sizeBuffer", sizeBuffer).set("WaterOccurrence", WaterOccurrence)
    def bufferRoughAOI(feature): # Buffer an AOI around the Sample Points
        return feature.buffer(sizeBuffer, errorMargin)
    def bufferFull(feature): # Buffer an AOI around the Sample Points
        return feature.buffer(ee.Number(sizeBuffer).multiply(2.1), errorMargin)

    def getWaterBodyAsFC(feature): # Gather polygons from spatial/temporal join and recast to Feature Collection 
        Distance2_SP = ee.Feature(feature.get('WBs')).getNumber('distance');
        AreaWB = ee.Feature(feature.get('WBs')).area(10);
        coords = ee.Feature(feature.get('WBs')).geometry();
        randSamples = ee.FeatureCollection.randomPoints(coords, Nsamples, 0, 1).geometry(1);
        return ee.Feature(coords, {'RandomSamples': randSamples, 
                                   'IDn': feature.getString('IDn'), 
                                   'SP': feature.getString('SP'), 
                                   'DataSource' : feature.getString('DataSource'),
                                   'WaterArea': AreaWB, 
                                   'd_SPtoAOI': Distance2_SP})
        
    def reducePolygons(f): # Reduce polygon size by specified input
        return f.buffer(bufferReduction, errorMargin)
    def getCentroid(polygon): 
        centroidCoords = polygon.centroid(errorMargin).geometry().coordinates().flatten();
        centroid = ee.Geometry.Point(centroidCoords);
        return polygon.setGeometry(centroid);
    def addCentroid(polygon): 
        centroidCoords = polygon.centroid(errorMargin).geometry().coordinates().flatten();
        centroid = ee.Geometry.Point(centroidCoords);
        return polygon.set('Centroid', centroid);

    def getMultiPoints(f): # Reduce polygon size by specified input
        pointDict = f.toDictionary();
        listCoords = f.geometry().coordinates(); # Repair the output lists of samples

        def getPoints(coords):
            p = ee.Geometry.Point(coords)
            d_centroid = p.distance(WB_centroid_reduced, errorMargin)
            d_shoreline = p.distance(WB_edges, errorMargin)
            d_AOI = p.distance(WB_edges_reduced, errorMargin)
            d_SP = p.distance(SPsLatLon, errorMargin)
            return ee.Feature(ee.Geometry.Point(coords)).set(pointDict).set({'d_centroid': d_centroid, 'd_shoreline': d_shoreline, 'd_AOI': d_AOI, 'd_SP': d_SP});
        featureCollection = ee.FeatureCollection(listCoords.map(getPoints)) # Repair formatting
        return featureCollection

    def generateRandomPoints(feature):
        buffer        = feature.geometry().buffer(sizeBuffer//4);                   # create a buffer of 1/4 of water other buffer
        randomPoints = ee.FeatureCollection.randomPoints(buffer, Nsamples, 0, 1); # generate 10 random points within the buffer
        feature      = feature.set('DataSource', feature.get('DataSource')).set('WaterArea', 0).set('IDn', feature.get('IDn')).set('SP', feature.geometry()).set('RandomSamples', randomPoints.geometry());
        randomPoints = randomPoints.map(lambda randomPoint: randomPoint.set(feature.toDictionary())); 
        return randomPoints;
        
    def distancesRivers(f):
        fp = f.geometry()
        mp = ee.Geometry(f.get('RandomSamples'))
        cp = ee.Geometry.MultiPoint(mp.coordinates()).centroid()
        bp = ee.Geometry.Point(f.get('SP')).buffer(sizeBuffer//4)
        d_centroid  = fp.distance(cp, errorMargin)
        d_shoreline = sizeBuffer//4
        d_AOI       = d_shoreline + bufferReduction 
        d_SP        = fp.distance(f.get('SP'), errorMargin)
        return f.set({'d_centroid': d_centroid, 'd_shoreline': d_shoreline, 'd_AOI': d_AOI, 'd_SP': d_SP}).set('d_SPtoAOI', 0).set('Centroid', cp);
            
    # Buffer areas of interest
    SPsLatLon = RawData #.map(getSamplePointsLonLat);
    roughAOI  = SPsLatLon.map(bufferRoughAOI).union(dJoinWindow);
    allWB     = SPsLatLon.map(bufferFull).union(dJoinWindow);

    # Get vectors of any recurring water pixels within each Sample Points AOI
    WMask = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').gt(WaterOccurrence);
    RecurringWaterBodies = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').gt(WaterOccurrence).mask(WMask);
    
    RecurringWaterBodies_Shoreline = RecurringWaterBodies.reduceToVectors(
            geometry      = allWB, 
            scale         = scale2use, 
            labelProperty = 'Place');

    RecurringWaterBodies_Vector = RecurringWaterBodies.reduceToVectors(
            geometry      = roughAOI, 
            scale         = scale2use, 
            labelProperty = 'Place').map(reducePolygons);

    # Water bodies joining and cleaning 
    joinedSPsWithWBs = SaveClosest.apply(SPsLatLon, RecurringWaterBodies_Shoreline, withinRoughAOI); # Join SPs and WBs within AOI
    WBs              = joinedSPsWithWBs.map(getWaterBodyAsFC); 
    WB_edges         = WBs.map(getEdge);

    # Reduced Water bodies joining and cleaning 
    joinedSPsWithWBs = SaveClosest.apply(SPsLatLon, RecurringWaterBodies_Vector, withinRoughAOI); # Join SPs and WBs within AOI
    WBs              = joinedSPsWithWBs.map(getWaterBodyAsFC); 
    WB_edges_reduced = WBs.map(getEdge);
    WB_centroid_reduced = WBs.map(getCentroid)

    WBs              = WBs.map(addCentroid);
    WBpoints         = WBs.map(WBsAOItoPoints);
    FC_WBpixels      = WBpoints.map(getMultiPoints).flatten();

    locationFilter   = WBs.geometry(); # Combine WB Geometries to filter

    # Now the inverse list (essentially land)
    pointList        = RawData
    otherCollection  = FC_WBpixels

    # Map the function over the feature list to generate a new feature list of random points
    randomPointList  = pointList.map(generateRandomPoints).flatten();
    joinFilter       = ee.Filter.equals(leftField="IDn", rightField="IDn")
    landFeatures    = ee.Join.inverted().apply(randomPointList, otherCollection, joinFilter)
    return riverFeatures.map(distancesRivers).map(set_var_optimisers), WBpoints.map(set_var_optimisers), FC_WBpixels.map(set_var_optimisers), WBs.map(set_var_optimisers)

