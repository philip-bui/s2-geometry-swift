//
//  S2LatLngRect.swift
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
	An S2LatLngRect represents a latitude-longitude rectangle. It is capable of
	representing the empty and full rectangles as well as single points.
*/
public struct S2LatLngRect: S2Region, Equatable {
	
	public let lat: R1Interval
	public let lng: S1Interval
	
	public var latLo: S1Angle {
		return S1Angle(radians: lat.lo)
	}
	
	public var latHi: S1Angle {
		return S1Angle(radians: lat.hi)
	}
	
	public var lngLo: S1Angle {
		return S1Angle(radians: lng.lo)
	}
	
	public var lngHi: S1Angle {
		return S1Angle(radians: lng.hi)
	}
	
	public var lo: S2LatLng {
		return S2LatLng(lat: S1Angle(radians: lat.lo), lng: S1Angle(radians: lng.lo))
	}
	
	public var hi: S2LatLng {
		return S2LatLng(lat: S1Angle(radians: lat.hi), lng: S1Angle(radians: lng.hi))
	}
	
	/**
		Construct a rectangle from minimum and maximum latitudes and longitudes. If
		lo.lng > hi.lng, the rectangle spans the 180 degree longitude line.
	*/
	public init(lo: S2LatLng, hi: S2LatLng) {
		lat = R1Interval(lo: lo.lat.radians, hi: hi.lat.radians)
		lng = S1Interval(lo: lo.lng.radians, hi: hi.lng.radians)
		// assert (isValid());
	}
	
	/// Construct a rectangle from latitude and longitude intervals.
	public init(lat: R1Interval, lng: S1Interval) {
		self.lat = lat
		self.lng = lng
		// assert (isValid());
	}
	
	/// The canonical empty rectangle
	public static var empty: S2LatLngRect {
		return S2LatLngRect(lat: .empty, lng: .empty)
	}
	
	/// The canonical full rectangle.
	public static var full: S2LatLngRect {
		return S2LatLngRect(lat: fullLat, lng: fullLng)
	}
	
	/// The full allowable range of latitudes.
	public static var fullLat: R1Interval {
		return R1Interval(lo: -(.pi / 2), hi: .pi / 2)
	}
	
	/// The full allowable range of longitudes.
	public static var fullLng: S1Interval {
		return .full
	}
	
	/**
		Return true if the rectangle is valid, which essentially just means that
		the latitude bounds do not exceed Pi/2 in absolute value and the longitude
		bounds do not exceed Pi in absolute value.
	*/
	public var isValid: Bool {
		// The lat/lng ranges must either be both empty or both non-empty.
		return (abs(lat.lo) <= .pi / 2 && abs(lat.hi) <= .pi / 2 && lng.isValid && lat.isEmpty == lng.isEmpty)
	}
	
	/// Return true if the rectangle is empty, i.e. it contains no points at all.
	public var isEmpty: Bool {
		return lat.isEmpty
	}
	
	/// Return true if the rectangle is full, i.e. it contains all points.
	public var isFull: Bool {
		return lat == S2LatLngRect.fullLat && lng.isFull
	}
	
	/// Return true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses the 180 degree latitude line.
	public var isInverted: Bool {
		return lng.isInverted
	}
	
	/// Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order.
	public func getVertex(k: Int) -> S2LatLng {
		// Return the points in CCW order (SW, SE, NE, NW).
		switch (k) {
		case 0:
			return S2LatLng.fromRadians(lat: lat.lo, lng: lng.lo)
		case 1:
			return S2LatLng.fromRadians(lat: lat.lo, lng: lng.hi)
		case 2:
			return S2LatLng.fromRadians(lat: lat.hi, lng: lng.hi)
		default:
			return S2LatLng.fromRadians(lat: lat.hi, lng: lng.lo)
		}
	}
	
	/// Return the center of the rectangle in latitude-longitude space (in general this is not the center of the region on the sphere).
	public var center: S2LatLng {
		return S2LatLng.fromRadians(lat: lat.center, lng: lng.center)
	}
	
	/**
		Return the minimum distance (measured along the surface of the sphere)
		from a given point to the rectangle (both its boundary and its interior).
		The latLng must be valid.
	*/
	public func getDistance(to p: S2LatLng) -> S1Angle {
		// The algorithm here is the same as in getDistance(S2LagLngRect), only
		// with simplified calculations.
		let a = self
		
		precondition(!a.isEmpty)
		precondition(p.isValid)
		
		if lng.contains(point: p.lng.radians) {
			return S1Angle(radians: max(0.0, max(p.lat.radians - a.lat.hi, a.lat.lo - p.lat.radians)))
		}
		
		let interval = S1Interval(lo: a.lng.hi, hi: a.lng.complement.center)
		var aLng = lng.lo
		if interval.contains(point: p.lng.radians) {
			aLng = lng.hi
		}
		
		let lo = S2LatLng.fromRadians(lat: a.lat.lo, lng: aLng).point
		let hi = S2LatLng.fromRadians(lat: a.lat.hi, lng: aLng).point
		let loCrossHi = S2LatLng.fromRadians(lat: 0, lng: aLng - .pi / 2).normalized.point
		return S2EdgeUtil.getDistance(x: p.point, a: lo, b: hi, aCrossB: loCrossHi)
	}
	
	/**
	* Return the minimum distance (measured along the surface of the sphere) to
	* the given S2LatLngRect. Both S2LatLngRects must be non-empty.
	*/
	public func getDistance(other: S2LatLngRect) -> S1Angle {
		let a = self
		let b = other
		
		precondition(!a.isEmpty)
		precondition(!b.isEmpty)
		
		// First, handle the trivial cases where the longitude intervals overlap.
		if a.lng.intersects(with: b.lng) {
			if a.lat.intersects(with: b.lat) {
				return S1Angle(radians: 0)  // Intersection between a and b.
			}
			
			// We found an overlap in the longitude interval, but not in the latitude
			// interval. This means the shortest path travels along some line of
			// longitude connecting the high-latitude of the lower rect with the
			// low-latitude of the higher rect.
			var lo: Double, hi: Double
			if a.lat.lo > b.lat.hi {
				lo = b.lat.hi
				hi = a.lat.lo
			} else {
				lo = a.lat.hi
				hi = b.lat.lo
			}
			return S1Angle(radians: hi - lo)
		}
		
		// The longitude intervals don't overlap. In this case, the closest points
		// occur somewhere on the pair of longitudinal edges which are nearest in
		// longitude-space.
		var aLng: Double, bLng: Double
		let loHi = S1Interval(lo: a.lng.lo, hi: b.lng.hi)
		let hiLo = S1Interval(lo: a.lng.hi, hi: b.lng.lo)
		if loHi.length < hiLo.length {
			aLng = a.lng.lo
			bLng = b.lng.hi
		} else {
			aLng = a.lng.hi
			bLng = b.lng.lo
		}
		
		// The shortest distance between the two longitudinal segments will include
		// at least one segment endpoint. We could probably narrow this down further
		// to a single point-edge distance by comparing the relative latitudes of the
		// endpoints, but for the sake of clarity, we'll do all four point-edge
		// distance tests.
		let aLo = S2LatLng.fromRadians(lat: a.lat.lo, lng: aLng).point
		let aHi = S2LatLng.fromRadians(lat: a.lat.hi, lng: aLng).point
		let aLoCrossHi = S2LatLng.fromRadians(lat: 0, lng: aLng - .pi / 2).normalized.point
		let bLo = S2LatLng.fromRadians(lat: b.lat.lo, lng: bLng).point
		let bHi = S2LatLng.fromRadians(lat: b.lat.hi, lng: bLng).point
		let bLoCrossHi = S2LatLng.fromRadians(lat: 0, lng: bLng - .pi / 2).normalized.point
		
		return min(S2EdgeUtil.getDistance(x: aLo, a: bLo, b: bHi, aCrossB: bLoCrossHi),
			min(S2EdgeUtil.getDistance(x: aHi, a: bLo, b: bHi, aCrossB: bLoCrossHi),
			min(S2EdgeUtil.getDistance(x: bLo, a: aLo, b: aHi, aCrossB: aLoCrossHi),
			S2EdgeUtil.getDistance(x: bHi, a: aLo, b: aHi, aCrossB: aLoCrossHi))))
	}
	
	/**
		Return the width and height of this rectangle in latitude-longitude space.
		Empty rectangles have a negative width and height.
	*/
	public var size: S2LatLng {
		return S2LatLng.fromRadians(lat: lat.length, lng: lng.length)
	}
	
	/// More efficient version of Contains() that accepts a S2LatLng rather than an S2Point.
	public func contains(ll: S2LatLng) -> Bool {
		// assert (ll.isValid());
		return lat.contains(point: ll.lat.radians) && lng.contains(point: ll.lng.radians)
	}
	
	/**
		Return true if and only if the given point is contained in the interior of
		the region (i.e. the region excluding its boundary). The point 'p' does not
		need to be normalized.
	*/
	public func interiorContains(point p: S2Point) -> Bool {
		return interiorContains(ll: S2LatLng(point: p))
	}
	
	/// More efficient version of InteriorContains() that accepts a S2LatLng rather than an S2Point.
	public func interiorContains(ll: S2LatLng)  -> Bool {
		// assert (ll.isValid());
		return lat.interiorContains(point: ll.lat.radians) && lng.interiorContains(point: ll.lng.radians)
	}
	
	/**
	* Return true if and only if the rectangle contains the given other
	* rectangle.
	*/
	public func contains(other: S2LatLngRect) -> Bool {
		return lat.contains(interval: other.lat) && lng.contains(interval: other.lng)
	}
	
	/**
		Return true if and only if the interior of this rectangle contains all
		points of the given other rectangle (including its boundary).
	*/
	public func interiorContains(other: S2LatLngRect) -> Bool {
		return lat.interiorContains(interval: other.lat) && lng.interiorContains(interval: other.lng)
	}
	
	/// Return true if this rectangle and the given other rectangle have any points in common.
	public func intersects(with other: S2LatLngRect) -> Bool {
		
		
		
		print("--- RECT intersects")
		print("--- \(lat) \(other.lat) => \(lat.intersects(with: other.lat))")
		print(latLo.degrees, latHi.degrees, "/", other.latLo.degrees, other.latHi.degrees)
		print("--- \(lng) \(other.lng) (\(lng.isInverted) \(other.lng.isInverted)) => \(lng.intersects(with: other.lng))")
		print(lngLo.degrees, lngHi.degrees, "/", other.lngLo.degrees, other.lngHi.degrees)
		print("---")
		
		return lat.intersects(with: other.lat) && lng.intersects(with: other.lng)
	}
	
	/**
	* Returns true if this rectangle intersects the given cell. (This is an exact
	* test and may be fairly expensive, see also MayIntersect below.)
	*/
	public func intersects(with cell: S2Cell) -> Bool {
		// First we eliminate the cases where one region completely contains the
		// other. Once these are disposed of, then the regions will intersect
		// if and only if their boundaries intersect.
		
		if isEmpty {
			return false
		}
		if contains(point: cell.center) {
			return true
		}
		if cell.contains(point: center.point) {
			return true
		}
		
		// Quick rejection test (not required for correctness).
		if !intersects(with: cell.rectBound) {
			return false
		}
		
		// Now check whether the boundaries intersect. Unfortunately, a
		// latitude-longitude rectangle does not have straight edges -- two edges
		// are curved, and at least one of them is concave.
		
		// Precompute the cell vertices as points and latitude-longitudes.
		var cellV: [S2Point] = []
		var cellLl: [S2LatLng] = []
		
		for i in 0 ..< 4 {
			let vertex = cell.getVertex(i)
			cellV.append(vertex) // Must be normalized.
			let latlng = S2LatLng(point: vertex)
			cellLl.append(latlng)
			if contains(ll: latlng) {
				return true // Quick acceptance test.
			}
		}
		
		for i in 0 ..< 4 {
			let edgeLng = S1Interval(lo: cellLl[i].lng.radians, hi: cellLl[(i + 1) & 3].lng.radians)
			if !lng.intersects(with: edgeLng) {
				continue
			}
			
			let a = cellV[i]
			let b = cellV[(i + 1) & 3]
			if edgeLng.contains(point: lng.lo) {
				if S2LatLngRect.intersectsLngEdge(a: a, b: b, lat: lat, lng: lng.lo) {
					return true
				}
			}
			if edgeLng.contains(point: lng.hi) {
				if S2LatLngRect.intersectsLngEdge(a: a, b: b, lat: lat, lng: lng.hi) {
					return true
				}
			}
			if S2LatLngRect.intersectsLatEdge(a: a, b: b, lat: lat.lo, lng: lng) {
				return true
			}
			if S2LatLngRect.intersectsLatEdge(a: a, b: b, lat: lat.hi, lng: lng) {
				return true
			}
		}
		return false
	}
	
	/**
		Return true if and only if the interior of this rectangle intersects any
		point (including the boundary) of the given other rectangle.
	*/
	public func interiorIntersects(other: S2LatLngRect) -> Bool {
		return lat.interiorIntersects(with: other.lat) && lng.interiorIntersects(with: other.lng)
	}
	
	public func add(point p: S2Point) -> S2LatLngRect {
		return add(point: S2LatLng(point: p))
	}
	
	// Increase the size of the bounding rectangle to include the given point.
	// The rectangle is expanded by the minimum amount possible.
	public func add(point ll: S2LatLng) -> S2LatLngRect {
		// assert (ll.isValid());
		let newLat = lat.add(point: ll.lat.radians)
		let newLng = lng.add(point: ll.lng.radians)
		return S2LatLngRect(lat: newLat, lng: newLng)
	}
	
	/**
		Return a rectangle that contains all points whose latitude distance from
		this rectangle is at most margin.lat(), and whose longitude distance from
		this rectangle is at most margin.lng(). In particular, latitudes are
		clamped while longitudes are wrapped. Note that any expansion of an empty
		interval remains empty, and both components of the given margin must be
		non-negative.
	
		NOTE: If you are trying to grow a rectangle by a certain *distance* on the
		sphere (e.g. 5km), use the ConvolveWithCap() method instead.
	*/
	public func expanded(margin: S2LatLng) -> S2LatLngRect {
		// assert (margin.lat().radians() >= 0 && margin.lng().radians() >= 0);
		if isEmpty {
			return self
		}
		return S2LatLngRect(lat: lat.expanded(radius: margin.lat.radians).intersection(with: S2LatLngRect.fullLat), lng: lng.expanded(radius: margin.lng.radians))
	}
	
	/// Return the smallest rectangle containing the union of this rectangle and the given rectangle.
	public func union(with other: S2LatLngRect) -> S2LatLngRect {
		return S2LatLngRect(lat: lat.union(with: other.lat), lng: lng.union(with: other.lng))
	}
	
	/// The point 'p' does not need to be normalized.
	public func contains(point p: S2Point) -> Bool {
		return contains(ll: S2LatLng(point: p))
	}
	
	/// Return true if the edge AB intersects the given edge of constant longitude.
	private static func intersectsLngEdge(a: S2Point, b: S2Point, lat: R1Interval, lng: Double) -> Bool {
		// Return true if the segment AB intersects the given edge of constant
		// longitude. The nice thing about edges of constant longitude is that
		// they are straight lines on the sphere (geodesics).
		return S2.simpleCrossing(a: a, b: b, c: S2LatLng.fromRadians(lat: lat.lo, lng: lng).point, d: S2LatLng.fromRadians(lat: lat.hi, lng: lng).point)
	}
		
	/// Return true if the edge AB intersects the given edge of constant latitude.
	private static func intersectsLatEdge(a: S2Point, b: S2Point, lat: Double, lng: S1Interval) -> Bool {
		// Return true if the segment AB intersects the given edge of constant
		// latitude. Unfortunately, lines of constant latitude are curves on
		// the sphere. They can intersect a straight edge in 0, 1, or 2 points.
		// assert (S2.isUnitLength(a) && S2.isUnitLength(b));
		
		// First, compute the normal to the plane AB that points vaguely north.
		var z = S2Point.normalize(point: S2.robustCrossProd(a: a, b: b))
		if z.z < 0 {
			z = -z
		}
		
		// Extend this to an orthonormal frame (x,y,z) where x is the direction
		// where the great circle through AB achieves its maximium latitude.
		let y = S2Point.normalize(point: S2.robustCrossProd(a: z, b: S2Point(x: 0, y: 0, z: 1)))
		let x = y.crossProd(z)
		// assert (S2.isUnitLength(x) && x.z >= 0);
		
		// Compute the angle "theta" from the x-axis (in the x-y plane defined
		// above) where the great circle intersects the given line of latitude.
		let sinLat = sin(lat)
		if abs(sinLat) >= x.z {
			return false // The great circle does not reach the given latitude.
		}
		// assert (x.z > 0);
		let cosTheta = sinLat / x.z
		let sinTheta = sqrt(1 - cosTheta * cosTheta)
		let theta = atan2(sinTheta, cosTheta)
		
		// The candidate intersection points are located +/- theta in the x-y
		// plane. For an intersection to be valid, we need to check that the
		// intersection point is contained in the interior of the edge AB and
		// also that it is contained within the given longitude interval "lng".
		
		// Compute the range of theta values spanned by the edge AB.
		let abTheta = S1Interval(lo: atan2(a.dotProd(y), a.dotProd(x)), hi: atan2(b.dotProd(y), b.dotProd(x)))
		
		if abTheta.contains(point: theta) {
			// Check if the intersection point is also in the given "lng" interval.
			let isect = (x * cosTheta) + (y * sinTheta)
			if lng.contains(point: atan2(isect.y, isect.x)) {
				return true;
			}
		}
		if abTheta.contains(point: -theta) {
			// Check if the intersection point is also in the given "lng" interval.
			let intersection = (x * cosTheta) - (y * sinTheta)
			if lng.contains(point: atan2(intersection.y, intersection.x)) {
				return true
			}
		}
		return false
	
	}
	
	////////////////////////////////////////////////////////////////////////
	// MARK: S2Region
	////////////////////////////////////////////////////////////////////////
	
	public var capBound: S2Cap {
		// We consider two possible bounding caps, one whose axis passes
		// through the center of the lat-long rectangle and one whose axis
		// is the north or south pole. We return the smaller of the two caps.
		
		if isEmpty {
			return .empty
		}
		
		var poleZ: Double, poleAngle: Double
		if lat.lo + lat.hi < 0 {
			// South pole axis yields smaller cap.
			poleZ = -1
			poleAngle = .pi / 2 + lat.hi
		} else {
			poleZ = 1
			poleAngle = .pi / 2 - lat.lo
		}
		let poleCap = S2Cap(axis: S2Point(x: 0, y: 0, z: poleZ), angle: S1Angle(radians: poleAngle))
		
		// For bounding rectangles that span 180 degrees or less in longitude, the
		// maximum cap size is achieved at one of the rectangle vertices. For
		// rectangles that are larger than 180 degrees, we punt and always return a
		// bounding cap centered at one of the two poles.
		let lngSpan = lng.hi - lng.lo
		if remainder(lngSpan, 2 * Double.pi) >= 0 {
			if lngSpan < 2 * Double.pi {
				var midCap = S2Cap(axis: center.point, angle: S1Angle(radians: 0))
				for k in 0 ..< 4 {
					midCap = midCap.add(point: getVertex(k: k).point)
				}
				if midCap.height < poleCap.height {
					return midCap
				}
			}
		}
		return poleCap
	}
	
	public var rectBound: S2LatLngRect {
		return self
	}
	
	public func contains(cell: S2Cell) -> Bool {
		// A latitude-longitude rectangle contains a cell if and only if it contains
		// the cell's bounding rectangle. (This is an exact test.)
		return contains(other: cell.rectBound)
	}
	
	/**
		This test is cheap but is NOT exact. Use Intersects() if you want a more
		accurate and more expensive test. Note that when this method is used by an
		S2RegionCoverer, the accuracy isn't all that important since if a cell may
		intersect the region then it is subdivided, and the accuracy of this method
		goes up as the cells get smaller.
	*/
	public func mayIntersect(cell: S2Cell) -> Bool {
		// This test is cheap but is NOT exact (see s2latlngrect.h).
		return intersects(with: cell.rectBound)
	}

	public static func == (lhs: S2LatLngRect, rhs: S2LatLngRect) -> Bool {
		return lhs.lat == rhs.lat && lhs.lng == rhs.lng
	}
	
}
