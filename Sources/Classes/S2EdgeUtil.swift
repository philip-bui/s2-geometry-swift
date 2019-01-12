//
//  S2EdgeUtil.swift
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
	
/// This class allows a vertex chain v0, v1, v2, ... to be efficiently tested for intersection with a given fixed edge AB.
public struct EdgeCrosser {
	
	private let a: S2Point
	private let b: S2Point
	private let aCrossB: S2Point
	
	// The fields below are updated for each vertex in the chain.
	
	// Previous vertex in the vertex chain.
	private var c: S2Point = S2Point()
	// The orientation of the triangle ACB.
	private var acb: Int = 0
	
	/**
	* AB is the given fixed edge, and C is the first vertex of the vertex
	* chain. All parameters must point to fixed storage that persists for the
	* lifetime of the EdgeCrosser object.
	*/
	public init(a: S2Point, b: S2Point, c: S2Point) {
		self.a = a
		self.b = b
		self.aCrossB = a.crossProd(b)
		restartAt(point: c)
	}
	
	/// Call this function when your chain 'jumps' to a new place.
	public mutating func restartAt(point c: S2Point) {
		self.c = c
		acb = -S2.robustCCW(a: a, b: b, c: c, aCrossB: aCrossB)
	}
	
	/**
	* This method is equivalent to calling the S2EdgeUtil.robustCrossing()
	* function (defined below) on the edges AB and CD. It returns +1 if there
	* is a crossing, -1 if there is no crossing, and 0 if two points from
	* different edges are the same. Returns 0 or -1 if either edge is
	* degenerate. As a side effect, it saves vertex D to be used as the next
	* vertex C.
	*/
	public mutating func robustCrossing(point d: S2Point) -> Int {
		// For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC must
		// all be oriented the same way (CW or CCW). We keep the orientation
		// of ACB as part of our state. When each new point D arrives, we
		// compute the orientation of BDA and check whether it matches ACB.
		// This checks whether the points C and D are on opposite sides of the
		// great circle through AB.
		
		// Recall that robustCCW is invariant with respect to rotating its
		// arguments, i.e. ABC has the same orientation as BDA.
		let bda = S2.robustCCW(a: a, b: b, c: d, aCrossB: aCrossB)
		var result: Int
		
		if (bda == -acb && bda != 0) {
			// Most common case -- triangles have opposite orientations.
			result = -1;
		} else if ((bda & acb) == 0) {
			// At least one value is zero -- two vertices are identical.
			result = 0;
		} else {
			// assert (bda == acb && bda != 0);
			result = robustCrossingInternal(point: d) // Slow path.
		}
		// Now save the current vertex D as the next vertex C, and also save the
		// orientation of the new triangle ACB (which is opposite to the current triangle BDA).
		c = d
		acb = -bda
		return result
	}
	
	/**
	* This method is equivalent to the S2EdgeUtil.edgeOrVertexCrossing() method
	* defined below. It is similar to robustCrossing, but handles cases where
	* two vertices are identical in a way that makes it easy to implement
	* point-in-polygon containment tests.
	*/
	public mutating func edgeOrVertexCrossing(point d: S2Point) -> Bool {
		// We need to copy c since it is clobbered by robustCrossing().
		let c2 = S2Point(x: c.x, y: c.y, z: c.z)
		
		let crossing = robustCrossing(point: d)
		if crossing < 0 { return false }
		if crossing > 0 { return true }
		
		return S2EdgeUtil.vertexCrossing(a: a, b: b, c: c2, d: d)
	}
	
	/// This function handles the "slow path" of robustCrossing().
	private func robustCrossingInternal(point d: S2Point) -> Int {
		// ACB and BDA have the appropriate orientations, so now we check the
		// triangles CBD and DAC.
		let cCrossD = c.crossProd(d)
		let cbd = -S2.robustCCW(a: c, b: d, c: b, aCrossB: cCrossD)
		if cbd != acb { return -1 }
		let dac = S2.robustCCW(a: c, b: d, c: a, aCrossB: cCrossD)
		return dac == acb ? 1 : -1
	}
}

/**
	This class computes a bounding rectangle that contains all edges defined by
	a vertex chain v0, v1, v2, ... All vertices must be unit length. Note that
	the bounding rectangle of an edge can be larger than the bounding rectangle
	of its endpoints, e.g. consider an edge that passes through the north pole.
*/
public struct RectBounder {
	// The previous vertex in the chain.
	private var a: S2Point = S2Point()
	
	// The corresponding latitude-longitude.
	private var aLatLng: S2LatLng = S2LatLng()
	
	// The current bounding rectangle.
	/// The bounding rectangle of the edge chain that connects the vertices defined so far.
	public private(set) var bound: S2LatLngRect = .empty
	
	public init() { }
	
	/**
		This method is called to add each vertex to the chain. 'b' must point to
		fixed storage that persists for the lifetime of the RectBounder.
	*/
	public mutating func add(point b: S2Point) {
		// assert (S2.isUnitLength(b));
		
		let bLatLng = S2LatLng(point: b)
		
		if bound.isEmpty {
			bound = bound.add(point: bLatLng)
		} else {
			// We can't just call bound.addPoint(bLatLng) here, since we need to
			// ensure that all the longitudes between "a" and "b" are included.
			bound = bound.union(with: S2LatLngRect(lo: aLatLng, hi: bLatLng))
			
			// Check whether the min/max latitude occurs in the edge interior.
			// We find the normal to the plane containing AB, and then a vector
			// "dir" in this plane that also passes through the equator. We use
			// RobustCrossProd to ensure that the edge normal is accurate even
			// when the two points are very close together.
			let aCrossB = S2.robustCrossProd(a: a, b: b)
			let dir = aCrossB.crossProd(S2Point(x: 0, y: 0, z: 1))
			let da = dir.dotProd(a)
			let db = dir.dotProd(b)
			
			if da * db < 0 {
				// Minimum/maximum latitude occurs in the edge interior. This affects
				// the latitude bounds but not the longitude bounds.
				let absLat = acos(abs(aCrossB.z / aCrossB.norm))
				var lat = bound.lat
				if da < 0 {
					// It's possible that absLat < lat.lo() due to numerical errors.
					lat = R1Interval(lo: lat.lo, hi: max(absLat, bound.lat.hi))
				} else {
					lat = R1Interval(lo: min(-absLat, bound.lat.lo), hi: lat.hi)
				}
				bound = S2LatLngRect(lat: lat, lng: bound.lng)
			}
		}
		a = b
		aLatLng = bLatLng
	}
}

/**
	The purpose of this class is to find edges that intersect a given XYZ
	bounding box. It can be used as an efficient rejection test when attempting to
	find edges that intersect a given region. It accepts a vertex chain v0, v1,
	v2, ... and returns a boolean value indicating whether each edge intersects
	the specified bounding box.

	We use XYZ intervals instead of something like longitude intervals because
	it is cheap to collect from S2Point lists and any slicing strategy should
	give essentially equivalent results.  See S2Loop for an example of use.
*/
public struct XYZPruner {
	
	private var lastVertex: S2Point = S2Point()
	
	// The region to be tested against.
	private var boundSet = false
	private var xmin: Double = 0
	private var ymin: Double = 0
	private var zmin: Double = 0
	private var xmax: Double = 0
	private var ymax: Double = 0
	private var zmax: Double = 0
	private var maxDeformation: Double = 0
	
	public init() { }
	
	/**
		Accumulate a bounding rectangle from provided edges.
	
		- Parameter from: start of edge
		- Parameter to: end of edge.
	*/
	public mutating func addEdgeToBounds(from: S2Point, to: S2Point) {
		if !boundSet {
			boundSet = true
			xmin = from.x
			ymin = from.y
			zmin = from.z
			xmin = from.x
			ymin = from.y
			zmin = from.z
		}
		xmin = min(xmin, min(to.x, from.x))
		ymin = min(ymin, min(to.y, from.y))
		zmin = min(zmin, min(to.z, from.z))
		xmax = max(xmax, max(to.x, from.x))
		ymax = max(ymax, max(to.y, from.y))
		zmax = max(zmax, max(to.z, from.z))
		
		// Because our arcs are really geodesics on the surface of the earth
		// an edge can have intermediate points outside the xyz bounds implicit
		// in the end points.  Based on the length of the arc we compute a
		// generous bound for the maximum amount of deformation.  For small edges
		// it will be very small but for some large arcs (ie. from (1N,90W) to
		// (1N,90E) the path can be wildly deformed.  I did a bunch of
		// experiments with geodesics to get safe bounds for the deformation.
		let absX = abs(from.x - to.x)
		let absY = abs(from.y - to.y)
		let absZ = abs(from.z - to.z)
		let approxArcLen = absX + absY + absZ
		if approxArcLen < 0.025 { // less than 2 degrees
			maxDeformation = max(maxDeformation, approxArcLen * 0.0025)
		} else if approxArcLen < 1.0 { // less than 90 degrees
			maxDeformation = max(maxDeformation, approxArcLen * 0.11)
		} else {
			maxDeformation = approxArcLen * 0.5
		}
	}
	
	public mutating func setFirstIntersectPoint(point v0: S2Point) {
		xmin = xmin - maxDeformation
		ymin = ymin - maxDeformation
		zmin = zmin - maxDeformation
		xmax = xmax + maxDeformation
		ymax = ymax + maxDeformation
		zmax = zmax + maxDeformation
		lastVertex = v0
	}
	
	/**
		Returns true if the edge going from the last point to this point passes
		through the pruner bounding box, otherwise returns false.  So the
		method returns false if we are certain there is no intersection, but it
		may return true when there turns out to be no intersection.
	*/
	public mutating func intersects(with v1: S2Point) -> Bool {
		var result = true
		
		if ((v1.x < xmin && lastVertex.x < xmin) || (v1.x > xmax && lastVertex.x > xmax)) {
			result = false
		} else if ((v1.y < ymin && lastVertex.y < ymin) || (v1.y > ymax && lastVertex.y > ymax)) {
			result = false
		} else if ((v1.z < zmin && lastVertex.z < zmin) || (v1.z > zmax && lastVertex.z > zmax)) {
			result = false
		}
		
		lastVertex = v1
		return result
	}
	
}

/**
	The purpose of this class is to find edges that intersect a given longitude
	interval. It can be used as an efficient rejection test when attempting to
	find edges that intersect a given region. It accepts a vertex chain v0, v1,
	v2, ... and returns a boolean value indicating whether each edge intersects
	the specified longitude interval.

	This class is not currently used as the XYZPruner is preferred for
	S2Loop, but this should be usable in similar circumstances.  Be wary
	of the cost of atan2() in conversions from S2Point to longitude!
*/
public struct LongitudePruner {
	// The interval to be tested against.
	private let interval: S1Interval
	
	// The longitude of the next v0.
	private var lng0: Double
	
	/**
		'interval' is the longitude interval to be tested against,
		and 'v0' is the first vertex of edge chain.
	*/
	public init(interval: S1Interval, point v0: S2Point) {
		self.interval = interval;
		self.lng0 = S2LatLng.longitude(point: v0).radians
	}
	
	/**
		Returns true if the edge (v0, v1) intersects the given longitude
		interval, and then saves 'v1' to be used as the next 'v0'.
	*/
	public mutating func intersects(with v1: S2Point) -> Bool {
		let lng1 = S2LatLng.longitude(point: v1).radians
		let result = interval.intersects(with: S1Interval(p1: lng0, p2: lng1))
		self.lng0 = lng1
		return result
	}
}

/**
	A wedge relation's test method accepts two edge chains A=(a0,a1,a2) and
	B=(b0,b1,b2) where a1==b1, and returns either -1, 0, or 1 to indicate the
	relationship between the region to the left of A and the region to the left
	of B. Wedge relations are used to determine the local relationship between
	two polygons that share a common vertex.

	All wedge relations require that a0 != a2 and b0 != b2. Other degenerate
	cases (such as a0 == b2) are handled as expected. The parameter "ab1"
	denotes the common vertex a1 == b1.
*/
public protocol WedgeRelation {
	static func test(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Int
}

public struct WedgeContains: WedgeRelation {
	/**
		Given two edge chains (see WedgeRelation above), this function returns +1
		if the region to the left of A contains the region to the left of B, and 0 otherwise.
	*/
	public static func test(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Int {
		// For A to contain B (where each loop interior is defined to be its left
		// side), the CCW edge order around ab1 must be a2 b2 b0 a0. We split
		// this test into two parts that test three vertices each.
		return S2.orderedCCW(a: a2, b: b2, c: b0, o: ab1) && S2.orderedCCW(a: b0, b: a0, c: a2, o: ab1) ? 1 : 0
	}
}

public struct WedgeIntersects: WedgeRelation {
	/**
		Given two edge chains (see WedgeRelation above), this function returns -1
		if the region to the left of A intersects the region to the left of B,
		and 0 otherwise. Note that regions are defined such that points along a
		boundary are contained by one side or the other, not both. So for
		example, if A,B,C are distinct points ordered CCW around a vertex O, then
		the wedges BOA, AOC, and COB do not intersect.
	*/
	public static func test(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Int {
		// For A not to intersect B (where each loop interior is defined to be
		// its left side), the CCW edge order around ab1 must be a0 b2 b0 a2.
		// Note that it's important to write these conditions as negatives
		// (!OrderedCCW(a,b,c,o) rather than Ordered(c,b,a,o)) to get correct
		// results when two vertices are the same.
		return S2.orderedCCW(a: a0, b: b2, c: b0, o: ab1) && S2.orderedCCW(a: b0, b: a2, c: a0, o: ab1) ? 0 : -1
	}
}

public struct WedgeContainsOrIntersects: WedgeRelation {
	/**
		Given two edge chains (see WedgeRelation above), this function returns +1
		if A contains B, 0 if A and B are disjoint, and -1 if A intersects but
		does not contain B.
	*/
	public static func test(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Int {
		// This is similar to WedgeContainsOrCrosses, except that we want to
		// distinguish cases (1) [A contains B], (3) [A and B are disjoint],
		// and (2,4,5,6) [A intersects but does not contain B].
		
		if S2.orderedCCW(a: a0, b: a2, c: b2, o: ab1) {
			// We are in case 1, 5, or 6, or case 2 if a2 == b2.
			return S2.orderedCCW(a: b2, b: b0, c: a0, o: ab1) ? 1 : -1 // Case 1 vs. 2,5,6.
		}
		// We are in cases 2, 3, or 4.
		if !S2.orderedCCW(a: a2, b: b0, c: b2, o: ab1) {
			return 0 // Case 3.
		}
		
		// We are in case 2 or 4, or case 3 if a2 == b0.
		return a2 == b0 ? 0 : -1 // Case 3 vs. 2,4.
	}
}

public struct WedgeContainsOrCrosses: WedgeRelation {
	/**
		Given two edge chains (see WedgeRelation above), this function returns +1
		if A contains B, 0 if B contains A or the two wedges do not intersect,
		and -1 if the edge chains A and B cross each other (i.e. if A intersects
		both the interior and exterior of the region to the left of B). In
		degenerate cases where more than one of these conditions is satisfied,
		the maximum possible result is returned. For example, if A == B then the
		result is +1.
	*/
	public static func test(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Int {
		// There are 6 possible edge orderings at a shared vertex (all
		// of these orderings are circular, i.e. abcd == bcda):
		//
		// (1) a2 b2 b0 a0: A contains B
		// (2) a2 a0 b0 b2: B contains A
		// (3) a2 a0 b2 b0: A and B are disjoint
		// (4) a2 b0 a0 b2: A and B intersect in one wedge
		// (5) a2 b2 a0 b0: A and B intersect in one wedge
		// (6) a2 b0 b2 a0: A and B intersect in two wedges
		//
		// In cases (4-6), the boundaries of A and B cross (i.e. the boundary
		// of A intersects the interior and exterior of B and vice versa).
		// Thus we want to distinguish cases (1), (2-3), and (4-6).
		//
		// Note that the vertices may satisfy more than one of the edge
		// orderings above if two or more vertices are the same. The tests
		// below are written so that we take the most favorable
		// interpretation, i.e. preferring (1) over (2-3) over (4-6). In
		// particular note that if orderedCCW(a,b,c,o) returns true, it may be
		// possible that orderedCCW(c,b,a,o) is also true (if a == b or b == c).
		
		if S2.orderedCCW(a: a0, b: a2, c: b2, o: ab1) {
			// The cases with this vertex ordering are 1, 5, and 6,
			// although case 2 is also possible if a2 == b2.
			if S2.orderedCCW(a: b2, b: b0, c: a0, o: ab1) {
				return 1 // Case 1 (A contains B)
			}
			
			// We are in case 5 or 6, or case 2 if a2 == b2.
			return a2 == b2 ? 0 : -1 // Case 2 vs. 5,6.
		}
		// We are in case 2, 3, or 4.
		return S2.orderedCCW(a: a0, b: b0, c: a2, o: ab1) ? 0 : -1 // Case 2,3 vs. 4.
	}
}

/**
	This class contains various utility functions related to edges. It collects
	together common code that is needed to implement polygonal geometry such as
	polylines, loops, and general polygons.
*/
public struct S2EdgeUtil {
	
	/**
		IEEE floating-point operations have a maximum error of 0.5 ULPS (units in
		the last place). For double-precision numbers, this works out to 2**-53
		(about 1.11e-16) times the magnitude of the result. It is possible to
		analyze the calculation done by getIntersection() and work out the
		worst-case rounding error. I have done a rough version of this, and my
		estimate is that the worst case distance from the intersection point X to
		the great circle through (a0, a1) is about 12 ULPS, or about 1.3e-15. This
		needs to be increased by a factor of (1/0.866) to account for the
		edgeSpliceFraction() in S2PolygonBuilder. Note that the maximum error
		measured by the unittest in 1,000,000 trials is less than 3e-16.
	*/
	public static let defaultIntersectionTolerance = S1Angle(radians: 1.5e-15)

	/**
		Return true if edge AB crosses CD at a point that is interior to both
		edges. Properties:
	
		1. simpleCrossing(b,a,c,d) == simpleCrossing(a,b,c,d)
		2. simpleCrossing(c,d,a,b) == simpleCrossing(a,b,c,d)
	*/
	public static func simpleCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Bool {
		// We compute simpleCCW() for triangles ACB, CBD, BDA, and DAC. All
		// of these triangles need to have the same orientation (CW or CCW)
		// for an intersection to exist. Note that this is slightly more
		// restrictive than the corresponding definition for planar edges,
		// since we need to exclude pairs of line segments that would
		// otherwise "intersect" by crossing two antipodal points.
		
		let ab = a.crossProd(b)
		let acb = -(ab.dotProd(c))
		let bda = ab.dotProd(d)
		if acb * bda <= 0 { return false }
		
		let cd = c.crossProd(d)
		let cbd = -(cd.dotProd(b))
		let dac = cd.dotProd(a)
		return (acb * cbd > 0) && (acb * dac > 0)
	}
	
	/**
		Like SimpleCrossing, except that points that lie exactly on a line are
		arbitrarily classified as being on one side or the other (according to the
		rules of S2.robustCCW). It returns +1 if there is a crossing, -1 if there
		is no crossing, and 0 if any two vertices from different edges are the
		same. Returns 0 or -1 if either edge is degenerate. Properties of
		robustCrossing:
	
		1. robustCrossing(b,a,c,d) == robustCrossing(a,b,c,d)
		2. robustCrossing(c,d,a,b) == robustCrossing(a,b,c,d)
		3. robustCrossing(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
		4. robustCrossing(a,b,c,d) <= 0 if a==b or c==d
	
		Note that if you want to check an edge against a *chain* of other edges,
		it is much more efficient to use an EdgeCrosser (above).
	*/
	public static func robustCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Int {
		// For there to be a crossing, the triangles ACB, CBD, BDA, DAC must
		// all have the same orientation (clockwise or counterclockwise).
		//
		// First we compute the orientation of ACB and BDA. We permute the
		// arguments to robustCCW so that we can reuse the cross-product of A and B.
		// Recall that when the arguments to robustCCW are permuted, the sign of the
		// result changes according to the sign of the permutation. Thus ACB and
		// ABC are oppositely oriented, while BDA and ABD are the same.
		let aCrossB = a.crossProd(b)
		let acb = -S2.robustCCW(a: a, b: b, c: c, aCrossB: aCrossB)
		let bda = S2.robustCCW(a: a, b: b, c: d, aCrossB: aCrossB)
		
		// If any two vertices are the same, the result is degenerate.
		if bda & acb == 0 { return 0 }
		
		// If ABC and BDA have opposite orientations (the most common case),
		// there is no crossing.
		if bda != acb { return -1 }
		
		// Otherwise we compute the orientations of CBD and DAC, and check whether
		// their orientations are compatible with the other two triangles.
		let cCrossD = c.crossProd(d)
		let cbd = -S2.robustCCW(a: c, b: d, c: b, aCrossB: cCrossD)
		if cbd != acb { return -1 }
		
		let dac = S2.robustCCW(a: c, b: d, c: a, aCrossB: cCrossD)
		return dac == acb ? 1 : -1
	}
	
	/**
	* Given two edges AB and CD where at least two vertices are identical (i.e.
	* robustCrossing(a,b,c,d) == 0), this function defines whether the two edges
	* "cross" in a such a way that point-in-polygon containment tests can be
	* implemented by counting the number of edge crossings. The basic rule is
	* that a "crossing" occurs if AB is encountered after CD during a CCW sweep
	* around the shared vertex starting from a fixed reference point.
	*
	*  Note that according to this rule, if AB crosses CD then in general CD does
	* not cross AB. However, this leads to the correct result when counting
	* polygon edge crossings. For example, suppose that A,B,C are three
	* consecutive vertices of a CCW polygon. If we now consider the edge
	* crossings of a segment BP as P sweeps around B, the crossing number changes
	* parity exactly when BP crosses BA or BC.
	*
	*  Useful properties of VertexCrossing (VC):
	*
	*  (1) VC(a,a,c,d) == VC(a,b,c,c) == false (2) VC(a,b,a,b) == VC(a,b,b,a) ==
	* true (3) VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c) (3) If
	* exactly one of a,b equals one of c,d, then exactly one of VC(a,b,c,d) and
	* VC(c,d,a,b) is true
	*
	* It is an error to call this method with 4 distinct vertices.
	*/
	public static func vertexCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Bool {
		// If A == B or C == D there is no intersection. We need to check this
		// case first in case 3 or more input points are identical.
		if a == b || c == d { return false }
		
		// If any other pair of vertices is equal, there is a crossing if and only
		// if orderedCCW() indicates that the edge AB is further CCW around the
		// shared vertex than the edge CD.
		if a == d {
			return S2.orderedCCW(a: a.ortho, b: c, c: b, o: a);
		}
		if b == c {
			return S2.orderedCCW(a: b.ortho, b: d, c: a, o: b);
		}
		if a == c {
			return S2.orderedCCW(a: a.ortho, b: d, c: b, o: a);
		}
		if b == d {
			return S2.orderedCCW(a: b.ortho, b: c, c: a, o: b);
		}
		
		// assert (false);
		return false
	}
	
	/**
		A convenience function that calls robustCrossing() to handle cases where
		all four vertices are distinct, and VertexCrossing() to handle cases where
		two or more vertices are the same. This defines a crossing function such
		that point-in-polygon containment tests can be implemented by simply
		counting edge crossings.
	*/
	public static func edgeOrVertexCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Bool {
		let crossing = robustCrossing(a: a, b: b, c: c, d: d)
		if crossing < 0 { return false }
		if crossing > 0 { return true }
		return vertexCrossing(a: a, b: b, c: c, d: d)
	}
	
	public struct CloserResult {
		public private(set) var dmin2: Double
		public private(set) var vmin: S2Point
		
		public init(dmin2: Double, vmin: S2Point) {
			self.dmin2 = dmin2
			self.vmin = vmin
		}
		
		public mutating func replaceIfCloser(x: S2Point, y: S2Point) {
			// If the squared distance from x to y is less than dmin2, then replace
			// vmin by y and update dmin2 accordingly.
			let d2 = (x - y).norm2
			if d2 < dmin2 || (d2 == dmin2 && y < vmin) {
				dmin2 = d2
				vmin = y
			}
		}
	}
	
	/**
		Given two edges AB and CD such that robustCrossing() is true, return their
		intersection point. Useful properties of getIntersection (GI):
	
		1. GI(b,a,c,d) == GI(a,b,d,c) == GI(a,b,c,d)
		2. GI(c,d,a,b) == GI(a,b,c,d)
	
		The returned intersection point X is guaranteed to be close to the edges AB
		and CD, but if the edges intersect at a very small angle then X may not be
		close to the true mathematical intersection point P. See the description of
		"DEFAULT_INTERSECTION_TOLERANCE" below for details.
	*/
	public static func getIntersection(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point) -> S2Point {
		precondition(robustCrossing(a: a0, b: a1, c: b0, d: b1) > 0, "Input edges a0a1 and b0b1 muct have a true robustCrossing.")
		
		// We use robustCrossProd() to get accurate results even when two endpoints
		// are close together, or when the two line segments are nearly parallel.
		let aNorm = S2Point.normalize(point: S2.robustCrossProd(a: a0, b: a1))
		let bNorm = S2Point.normalize(point: S2.robustCrossProd(a: b0, b: b1))
		var x = S2Point.normalize(point: S2.robustCrossProd(a: aNorm, b: bNorm))
		
		// Make sure the intersection point is on the correct side of the sphere.
		// Since all vertices are unit length, and edges are less than 180 degrees,
		// (a0 + a1) and (b0 + b1) both have positive dot product with the
		// intersection point. We use the sum of all vertices to make sure that the
		// result is unchanged when the edges are reversed or exchanged.
		if x.dotProd((a0 + a1) + (b0 + b1)) < 0 {
			x = -x
		}
		
		// The calculation above is sufficient to ensure that "x" is within
		// DEFAULT_INTERSECTION_TOLERANCE of the great circles through (a0,a1) and
		// (b0,b1).
		// However, if these two great circles are very close to parallel, it is
		// possible that "x" does not lie between the endpoints of the given line
		// segments. In other words, "x" might be on the great circle through
		// (a0,a1) but outside the range covered by (a0,a1). In this case we do
		// additional clipping to ensure that it does.
		
		if S2.orderedCCW(a: a0, b: x, c: a1, o: aNorm) && S2.orderedCCW(a: b0, b: x, c: b1, o: bNorm) {
			return x
		}
		
		// Find the acceptable endpoint closest to x and return it. An endpoint is
		// acceptable if it lies between the endpoints of the other line segment.
		var r = CloserResult(dmin2: 10, vmin: x)
		if S2.orderedCCW(a: b0, b: a0, c: b1, o: bNorm) {
			r.replaceIfCloser(x: x, y: a0)
		}
		if S2.orderedCCW(a: b0, b: a1, c: b1, o: bNorm) {
			r.replaceIfCloser(x: x, y: a1)
		}
		if S2.orderedCCW(a: a0, b: b0, c: a1, o: aNorm) {
			r.replaceIfCloser(x: x, y: b0)
		}
		if S2.orderedCCW(a: a0, b: b1, c: a1, o: aNorm) {
			r.replaceIfCloser(x: x, y: b1)
		}
		return r.vmin
	}
	
	/**
		Given a point X and an edge AB, return the distance ratio AX / (AX + BX).
		If X happens to be on the line segment AB, this is the fraction "t" such
		that X == Interpolate(A, B, t). Requires that A and B are distinct.
	*/
	public static func getDistanceFraction(x: S2Point, a0: S2Point, a1: S2Point) -> Double {
		precondition(a0 != a1)
		let d0 = x.angle(to: a0)
		let d1 = x.angle(to: a1)
		return d0 / (d0 + d1)
	}
	
	/**
		Return the minimum distance from X to any point on the edge AB. The result
		is very accurate for small distances but may have some numerical error if
		the distance is large (approximately Pi/2 or greater). The case A == B is
		handled correctly. Note: x, a and b must be of unit length. Throws
		IllegalArgumentException if this is not the case.
	*/
	public static func getDistance(x: S2Point, a: S2Point, b: S2Point) -> S1Angle {
		return getDistance(x: x, a: a, b: b, aCrossB: S2.robustCrossProd(a: a, b: b))
	}
	
	/**
		A slightly more efficient version of getDistance() where the cross product
		of the two endpoints has been precomputed. The cross product does not need
		to be normalized, but should be computed using S2.robustCrossProd() for the
		most accurate results.
	*/
	public static func getDistance(x: S2Point, a: S2Point, b: S2Point, aCrossB: S2Point) -> S1Angle {
		precondition(S2.isUnitLength(point: x))
		precondition(S2.isUnitLength(point: a))
		precondition(S2.isUnitLength(point: b))
		
		// There are three cases. If X is located in the spherical wedge defined by
		// A, B, and the axis A x B, then the closest point is on the segment AB.
		// Otherwise the closest point is either A or B; the dividing line between
		// these two cases is the great circle passing through (A x B) and the
		// midpoint of AB.
		
		if S2.simpleCCW(a: aCrossB, b: a, c: x) && S2.simpleCCW(a: x, b: b, c: aCrossB) {
			// The closest point to X lies on the segment AB. We compute the distance
			// to the corresponding great circle. The result is accurate for small
			// distances but not necessarily for large distances (approaching Pi/2).
			
			let sinDist = abs(x.dotProd(aCrossB)) / aCrossB.norm
			return S1Angle(radians: asin(min(1.0, sinDist)))
		}
		
		// Otherwise, the closest point is either A or B. The cheapest method is
		// just to compute the minimum of the two linear (as opposed to spherical)
		// distances and convert the result to an angle. Again, this method is
		// accurate for small but not large distances (approaching Pi).
		
		let linearDist2 = min((x - a).norm2, (x - b).norm2)
		return S1Angle(radians: 2 * asin(min(1.0, 0.5 * sqrt(linearDist2))))
	}
	
	/**
		Returns the point on edge AB closest to X. x, a and b must be of unit length.
		Throws IllegalArgumentException if this is not the case.
	*/
	public static func getClosestPoint(x: S2Point, a: S2Point, b: S2Point) -> S2Point {
		precondition(S2.isUnitLength(point: x))
		precondition(S2.isUnitLength(point: a))
		precondition(S2.isUnitLength(point: b))
		
		let crossProd = S2.robustCrossProd(a: a, b: b)
		// Find the closest point to X along the great circle through AB.
		let p = (x - (crossProd * (x.dotProd(crossProd) / crossProd.norm2)))
		
		// If p is on the edge AB, then it's the closest point.
		if (S2.simpleCCW(a: crossProd, b: a, c: p) && S2.simpleCCW(a: p, b: b, c: crossProd)) {
			return S2Point.normalize(point: p)
		}
		// Otherwise, the closest point is either A or B.
		return (x - a).norm2 <= (x - b).norm2 ? a : b
	}
	
	// Don't instantiate
	private init() { }
	
}
