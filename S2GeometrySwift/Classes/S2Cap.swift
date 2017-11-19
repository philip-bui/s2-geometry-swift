//
//  S2Cap.swift
//  S2Geometry
//
//  Created by Alex Studnicka on 7/1/16.
//  Copyright © 2016 Alex Studnicka. MIT License.
//

#if os(Linux)
	import Glibc
#else
	import Darwin.C
#endif

/**
	This class represents a spherical cap, i.e. a portion of a sphere cut off by
	a plane. The cap is defined by its axis and height. This representation has
	good numerical accuracy for very small caps (unlike the (axis,
	min-distance-from-origin) representation), and is also efficient for
	containment tests (unlike the (axis, angle) representation).

	Here are some useful relationships between the cap height (h), the cap
	opening angle (theta), the maximum chord length from the cap's center (d),
	and the radius of cap's base (a). All formulas assume a unit radius.

	h = 1 - cos(theta) = 2 sin^2(theta/2) d^2 = 2 h = a^2 + h^2
*/
public struct S2Cap: S2Region {
	
	/**
		Multiply a positive number by this constant to ensure that the result of a
		floating point operation is at least as large as the true
		infinite-precision result.
	*/
    private static let roundUp: Double = 2 / 4503599627370496
	
	public let axis: S2Point
	public let height: Double
	
	/**
		Create a cap given its axis and the cap height, i.e. the maximum projected
		distance along the cap axis from the cap center. 'axis' should be a
		unit-length vector.
	*/
	public init(axis: S2Point = S2Point(), height: Double = 0) {
		self.axis = axis
		self.height = height
	}
	
	/**
		Create a cap given its axis and the cap opening angle, i.e. maximum angle
		between the axis and a point on the cap. 'axis' should be a unit-length
		vector, and 'angle' should be between 0 and 180 degrees.
	*/
	public init(axis: S2Point, angle: S1Angle) {
		// The height of the cap can be computed as 1-cos(angle), but this isn't
		// very accurate for angles close to zero (where cos(angle) is almost 1).
		// Computing it as 2*(sin(angle/2)**2) gives much better precision.
		
		// assert (S2.isUnitLength(axis));
		let d = sin(0.5 * angle.radians)
		self.init(axis: axis, height: 2 * d * d)
	}
	
	/**
		Create a cap given its axis and its area in steradians. 'axis' should be a
		unit-length vector, and 'area' should be between 0 and 4 * .pi.
	*/
	public init(axis: S2Point, area: Double) {
		// assert (S2.isUnitLength(axis));
		self.init(axis: axis, height: area / (2 * .pi))
	}
	
	/// Return an empty cap, i.e. a cap that contains no points.
	public static let empty = S2Cap(axis: S2Point(x: 1, y: 0, z: 0), height: -1)
	
	/// Return a full cap, i.e. a cap that contains all points.
	public static let full = S2Cap(axis: S2Point(x: 1, y: 0, z: 0), height: 2)
	
	public var area: Double {
		return 2 * .pi * max(0, height)
	}
	
	/// Return the cap opening angle in radians, or a negative number for empty caps.
	public var angle: S1Angle {
		// This could also be computed as acos(1 - height_), but the following
		// formula is much more accurate when the cap height is small. It
		// follows from the relationship h = 1 - cos(theta) = 2 sin^2(theta/2).
		if isEmpty {
			return S1Angle(radians: -1)
		}
		return S1Angle(radians: 2 * asin(sqrt(0.5 * height)))
	}
	
	/// We allow negative heights (to represent empty caps) but not heights greater than 2.
	public var isValid: Bool {
		return S2.isUnitLength(point: axis) && height <= 2
	}
	
	/// Return true if the cap is empty, i.e. it contains no points.
	public var isEmpty: Bool {
		return height < 0
	}
	
	/// Return true if the cap is full, i.e. it contains all points.
	public var isFull: Bool {
		return height >= 2
	}
	
	/**
		Return the complement of the interior of the cap. A cap and its complement
		have the same boundary but do not share any interior points. The complement
		operator is not a bijection, since the complement of a singleton cap
		(containing a single point) is the same as the complement of an empty cap.
	*/
	public var complement: S2Cap {
		// The complement of a full cap is an empty cap, not a singleton.
		// Also make sure that the complement of an empty cap has height 2.
		let cHeight = isFull ? -1 : 2 - max(height, 0.0)
		return S2Cap(axis: -axis, height: cHeight)
	}
	
	/**
	* Return true if and only if this cap contains the given other cap (in a set
	* containment sense, e.g. every cap contains the empty cap).
	*/
	public func contains(other: S2Cap) -> Bool {
		if isFull || other.isEmpty {
			return true
		}
		return angle.radians >= axis.angle(to: other.axis) + other.angle.radians
	}
	
	/**
		Return true if and only if the interior of this cap intersects the given
		other cap. (This relationship is not symmetric, since only the interior of
		this cap is used.)
	*/
	public func interiorIntersects(with other: S2Cap) -> Bool {
		// Interior(X) intersects Y if and only if Complement(Interior(X)) does not contain Y.
		return !complement.contains(other: other)
	}
	
	/**
		Return true if and only if the given point is contained in the interior of
		the region (i.e. the region excluding its boundary). 'p' should be a
		unit-length vector.
	*/
	public func interiorContains(point p: S2Point) -> Bool {
		// assert (S2.isUnitLength(p));
		return isFull || (axis - p).norm2 < 2 * height
	}
	
	/**
		Increase the cap height if necessary to include the given point. If the cap
		is empty the axis is set to the given point, but otherwise it is left
		unchanged. 'p' should be a unit-length vector.
	*/
	public func add(point p: S2Point) -> S2Cap {
		// Compute the squared chord length, then convert it into a height.
		// assert (S2.isUnitLength(p));
		if isEmpty {
			return S2Cap(axis: p, height: 0)
		} else {
			// To make sure that the resulting cap actually includes this point,
			// we need to round up the distance calculation. That is, after
			// calling cap.AddPoint(p), cap.Contains(p) should be true.
			let dist2 = (axis - p).norm2
			let newHeight = max(height, S2Cap.roundUp * 0.5 * dist2)
			return S2Cap(axis: axis, height: newHeight)
		}
	}
	
	// Increase the cap height if necessary to include "other". If the current
	// cap is empty it is set to the given other cap.
	public func add(cap other: S2Cap) -> S2Cap {
		if isEmpty {
			return S2Cap(axis: other.axis, height: other.height)
		} else {
			// See comments for FromAxisAngle() and AddPoint(). This could be
			// optimized by doing the calculation in terms of cap heights rather
			// than cap opening angles.
			let angle = axis.angle(to: other.axis) + other.angle.radians
			if angle >= .pi {
				return S2Cap(axis: axis, height: 2) //Full cap
			} else {
				let d = sin(0.5 * angle)
				let newHeight = max(height, S2Cap.roundUp * 2 * d * d)
				return S2Cap(axis: axis, height: newHeight)
			}
		}
	}
	
	/// Return true if the cap intersects 'cell', given that the cap vertices have alrady been checked.
	public func intersects(with cell: S2Cell, vertices: [S2Point]) -> Bool {
		// Return true if this cap intersects any point of 'cell' excluding its
		// vertices (which are assumed to already have been checked).
		
		// If the cap is a hemisphere or larger, the cell and the complement of the
		// cap are both convex. Therefore since no vertex of the cell is contained,
		// no other interior point of the cell is contained either.
		if height >= 1 {
			return false
		}
		
		// We need to check for empty caps due to the axis check just below.
		if isEmpty {
			return false
		}
		
		// Optimization: return true if the cell contains the cap axis. (This
		// allows half of the edge checks below to be skipped.)
		if cell.contains(point: axis) {
			return true
		}
		
		// At this point we know that the cell does not contain the cap axis,
		// and the cap does not contain any cell vertex. The only way that they
		// can intersect is if the cap intersects the interior of some edge.
		
		let sin2Angle = height * (2 - height) // sin^2(capAngle)
		for k in 0 ..< 4 {
			let edge = cell.getRawEdge(k)
			let dot = axis.dotProd(edge)
			if dot > 0 {
				// The axis is in the interior half-space defined by the edge. We don't
				// need to consider these edges, since if the cap intersects this edge
				// then it also intersects the edge on the opposite side of the cell
				// (because we know the axis is not contained with the cell).
				continue
			}
			// The Norm2() factor is necessary because "edge" is not normalized.
			if dot * dot > sin2Angle * edge.norm2 {
				return false // Entire cap is on the exterior side of this edge.
			}
			// Otherwise, the great circle containing this edge intersects
			// the interior of the cap. We just need to check whether the point
			// of closest approach occurs between the two edge endpoints.
			let dir = edge.crossProd(axis)
			if (dir.dotProd(vertices[k]) < 0 && dir.dotProd(vertices[(k + 1) & 3]) > 0) {
				return true
			}
		}
		return false
	}
	
	public func contains(point p: S2Point) -> Bool {
		// The point 'p' should be a unit-length vector.
		// assert (S2.isUnitLength(p));
		return (axis - p).norm2 <= 2 * height
	}
	
	////////////////////////////////////////////////////////////////////////
	// MARK: S2Region
	////////////////////////////////////////////////////////////////////////
	
	public var capBound: S2Cap {
		return self
	}
	
	public var rectBound: S2LatLngRect {
		if isEmpty {
			return .empty
		}
		
		// Convert the axis to a (lat,lng) pair, and compute the cap angle.
		let axisLatLng = S2LatLng(point: axis)
		let capAngle = angle.radians
		
		var allLongitudes = false
		var lat: (Double, Double) = (0, 0)
		var lng: (Double, Double) = (-.pi, .pi)
		
		// Check whether cap includes the south pole.
		lat.0 = axisLatLng.lat.radians - capAngle
		if lat.0 <= -.pi / 2 {
			lat.0 = -.pi / 2
			allLongitudes = true
		}
		// Check whether cap includes the north pole.
		lat.1 = axisLatLng.lat.radians + capAngle
		if lat.1 >= .pi / 2 {
			lat.1 = .pi / 2
			allLongitudes = true
		}
		if !allLongitudes {
			// Compute the range of longitudes covered by the cap. We use the law
			// of sines for spherical triangles. Consider the triangle ABC where
			// A is the north pole, B is the center of the cap, and C is the point
			// of tangency between the cap boundary and a line of longitude. Then
			// C is a right angle, and letting a,b,c denote the sides opposite A,B,C,
			// we have sin(a)/sin(A) = sin(c)/sin(C), or sin(A) = sin(a)/sin(c).
			// Here "a" is the cap angle, and "c" is the colatitude (90 degrees
			// minus the latitude). This formula also works for negative latitudes.
			//
			// The formula for sin(a) follows from the relationship h = 1 - cos(a).
			
			let sinA = sqrt(height * (2 - height))
			let sinC = cos(axisLatLng.lat.radians)
			if sinA <= sinC {
				let angleA = asin(sinA / sinC)
				lng.0 = remainder(axisLatLng.lng.radians - angleA, 2 * .pi)
				lng.1 = remainder(axisLatLng.lng.radians + angleA, 2 * .pi)
			}
		}
		return S2LatLngRect(lat: R1Interval(lo: lat.0, hi: lat.1), lng: S1Interval(lo: lng.0, hi: lng.1))
	}
	
	public func contains(cell: S2Cell) -> Bool {
		// If the cap does not contain all cell vertices, return false.
		// We check the vertices before taking the Complement() because we can't
		// accurately represent the complement of a very small cap (a height
		// of 2-epsilon is rounded off to 2).
		var vertices: [S2Point] = []
		for k in 0 ..< 4 {
			let vertex = cell.getVertex(k)
			if !contains(point: vertex) {
				return false
			}
			vertices.append(vertex)
		}
		// Otherwise, return true if the complement of the cap does not intersect
		// the cell. (This test is slightly conservative, because technically we
		// want Complement().InteriorIntersects() here.)
		return !complement.intersects(with: cell, vertices: vertices)
	}
	
	public func mayIntersect(cell: S2Cell) -> Bool {
		// If the cap contains any cell vertex, return true.
		var vertices: [S2Point] = []
		for k in 0 ..< 4 {
			let vertex = cell.getVertex(k)
			if contains(point: vertex) {
				return true
			}
			vertices.append(vertex)
		}
		return intersects(with: cell, vertices: vertices)
	}
	
}
