//
//  S2Region.swift
//  S2Geometry
//
//  Created by Alex Studnicka on 7/1/16.
//  Copyright © 2016 Alex Studnicka. MIT License.
//

/**
	An S2Region represents a two-dimensional region over the unit sphere. It is
	an abstract interface with various concrete subtypes.

	The main purpose of this interface is to allow complex regions to be
	approximated as simpler regions. So rather than having a wide variety of
	virtual methods that are implemented by all subtypes, the interface is
	restricted to methods that are useful for computing approximations.
*/
public protocol S2Region {
	
	/// Return a bounding spherical cap.
	var capBound: S2Cap { get }
	
	/// Return a bounding latitude-longitude rectangle.
	var rectBound: S2LatLngRect { get }
	
	/**
		If this method returns true, the region completely contains the given cell.
		Otherwise, either the region does not contain the cell or the containment
		relationship could not be determined.
	*/
	func contains(cell: S2Cell) -> Bool
	
	/**
		If this method returns false, the region does not intersect the given cell.
		Otherwise, either region intersects the cell, or the intersection
		relationship could not be determined.
	*/
	func mayIntersect(cell: S2Cell) -> Bool
	
}
