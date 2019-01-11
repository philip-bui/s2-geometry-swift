//
//  S2AreaCentroid.swift
//  S2Geometry
//
//  Created by Alex Studnicka on 7/1/16.
//  Copyright © 2016 Alex Studnicka. MIT License.
//

/**
	The area of an interior, i.e. the region on the left side of an odd
	number of loops and optionally a centroid.
	The area is between 0 and 4*Pi. If it has a centroid, it is
	the true centroid of the interiord multiplied by the area of the shape.
	Note that the centroid may not be contained by the shape.
*/
public struct S2AreaCentroid {
	
	public let area: Double
	public let centroid: S2Point
	
	public init(area: Double, centroid: S2Point) {
		self.area = area
		self.centroid = centroid
	}
	
}
