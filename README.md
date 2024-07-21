Water Detection and Sampling Algorithm

This repository contains the code for the Water Detection and Sampling Algorithm, implemented in sampling.py and geotools.py. The algorithm is designed to detect and sample water bodies using Google Earth Engine (GEE).

Features

	•	Water Occurrence Detection: Identifies water bodies based on water occurrence percentage using JRC Global Surface Water dataset.
	•	Spatial Join and Buffering: Connects water bodies within a specified distance and buffers areas of interest.
	•	Polygon Reduction: Minimizes edge effects by reducing polygon sizes.
	•	Random Point Sampling: Generates random sample points within detected water bodies.
	•	Distance Calculations: Calculates distances from sample points to various geographical features.

Requirements

	•	Google Earth Engine (GEE) API
	•	Python 3.x
	•	ee (Earth Engine Python API)

Installation

To use the code, ensure you have the GEE API installed and authenticated. Follow the GEE Python API documentation for installation and authentication.

Usage

	1.	Import the necessary functions from geotools.py
 	2.	Define the Water Detection and Sampling Algorithm in sampling.py
  3.	Set initial parameters
  4.	Run the algorithm
     results = WaterGeoAlgorithm(WaterOccurrence, sizeBuffer, dJoinWindow, bufferReduction, save_directory)
