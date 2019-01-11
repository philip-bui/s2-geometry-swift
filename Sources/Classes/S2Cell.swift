//
//  S2Cell.swift
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
	An S2Cell is an S2Region object that represents a cell. Unlike S2CellIds, it
	supports efficient containment and intersection tests. However, it is also a
	more expensive representation.
*/
public struct S2Cell: S2Region, Equatable {
	
	private static let maxCellSize = 1 << S2CellId.maxLevel
	
	public let cellId: S2CellId
	public let face: UInt8
	public let level: Int8
	public let orientation: UInt8
	public let uv: [[Double]]
	
	internal init(cellId: S2CellId = S2CellId(), face: UInt8 = 0, level: Int8 = 0, orientation: UInt8 = 0, uv: [[Double]] = [[0, 0], [0, 0]]) {
		self.cellId = cellId
		self.face = face
		self.level = level
		self.orientation = orientation
		self.uv = uv
	}
	
	/// An S2Cell always corresponds to a particular S2CellId. The other constructors are just convenience methods.
	public init(cellId: S2CellId) {
		self.cellId = cellId
		
		var i = 0
		var j = 0
		var mOrientation: Int? = 0
		
		face = UInt8(cellId.toFaceIJOrientation(i: &i, j: &j, orientation: &mOrientation))
		orientation = UInt8(mOrientation!)
		level = Int8(cellId.level)
		
		let cellSize = 1 << (S2CellId.maxLevel - Int(level))
		var _uv: [[Double]] = [[0, 0], [0, 0]]
		for (d, ij) in [i, j].enumerated() {
			// Compute the cell bounds in scaled (i,j) coordinates.
			let sijLo = (ij & -cellSize) * 2 - S2Cell.maxCellSize
			let sijHi = sijLo + cellSize * 2
			_uv[d][0] = S2Projections.stToUV(s: (1.0 / Double(S2Cell.maxCellSize)) * Double(sijLo))
			_uv[d][1] = S2Projections.stToUV(s: (1.0 / Double(S2Cell.maxCellSize)) * Double(sijHi))
		}
		uv = _uv
	}
	
	// This is a static method in order to provide named parameters.
	public init(face: Int, pos: UInt8, level: Int) {
		self.init(cellId: S2CellId(face: face, pos: Int64(pos), level: level))
	}
	
	// Convenience methods.
	public init(point: S2Point) {
		self.init(cellId: S2CellId(point: point))
	}
	
	public init(latlng: S2LatLng) {
		self.init(cellId: S2CellId(latlng: latlng))
	}
	
	public var isLeaf: Bool {
		return Int(level) == S2CellId.maxLevel
	}
	
	public func getVertex(_ k: Int) -> S2Point {
		return S2Point.normalize(point: getRawVertex(k))
	}
	
	/**
		Return the k-th vertex of the cell (k = 0,1,2,3). Vertices are returned in
		CCW order. The points returned by GetVertexRaw are not necessarily unit length.
	*/
	public func getRawVertex(_ k: Int) -> S2Point {
		// Vertices are returned in the order SW, SE, NE, NW.
		return S2Projections.faceUvToXyz(face: Int(face), u: uv[0][(k >> 1) ^ (k & 1)], v: uv[1][k >> 1])
	}
		
	public func getEdge(_ k: Int) -> S2Point {
		return S2Point.normalize(point: getRawEdge(k))
	}
	
	public func getRawEdge(_ k: Int) -> S2Point {
		switch (k) {
		case 0:
			return S2Projections.getVNorm(face: Int(face), v: uv[1][0])		// South
		case 1:
			return S2Projections.getUNorm(face: Int(face), u: uv[0][1])		// East
		case 2:
			return -S2Projections.getVNorm(face: Int(face), v: uv[1][1])	// North
		default:
			return -S2Projections.getUNorm(face: Int(face), u: uv[0][0])	// West
		}
	}
	
	/**
		Return the inward-facing normal of the great circle passing through the
		edge from vertex k to vertex k+1 (mod 4). The normals returned by
		GetEdgeRaw are not necessarily unit length.
	
		If this is not a leaf cell, set children[0..3] to the four children of
		this cell (in traversal order) and return true. Otherwise returns false.
		This method is equivalent to the following:
	
		for (pos=0, id=child_begin(); id != child_end(); id = id.next(), ++pos)
		children[i] = S2Cell(id);
	
		except that it is more than two times faster.
	*/
	public func subdivide() -> [S2Cell] {
		// This function is equivalent to just iterating over the child cell ids
		// and calling the S2Cell constructor, but it is about 2.5 times faster.
		
		guard !cellId.isLeaf else { return [] }
		
		// Compute the cell midpoint in uv-space.
		let uvMid = centerUV
		
		// Create four children with the appropriate bounds.
		var children: [S2Cell] = []
		var id = cellId.childBegin()
		for pos in 0 ..< 4 {
			
			var _uv: [[Double]] = [[0, 0], [0, 0]]
			let ij = S2.posToIJ(orientation: Int(orientation), position: pos)
			
			for d in 0 ..< 2 {
				// The dimension 0 index (i/u) is in bit 1 of ij.
				let m = 1 - ((ij >> (1 - d)) & 1)
				_uv[d][m] = uvMid.get(index: d)
				_uv[d][1 - m] = uv[d][1 - m]
			}
			let child = S2Cell(cellId: id, face: face, level: level + 1, orientation: orientation ^ UInt8(S2.posToOrientation(position: pos)), uv: _uv)
			children.append(child)
			
			id = id.next()
		}
		return children
	}
	
	/**
		Return the direction vector corresponding to the center in (s,t)-space of
		the given cell. This is the point at which the cell is divided into four
		subcells; it is not necessarily the centroid of the cell in (u,v)-space or
		(x,y,z)-space. The point returned by GetCenterRaw is not necessarily unit length.
	*/
	public var center: S2Point {
		return S2Point.normalize(point: rawCenter)
	}
	
	public var rawCenter: S2Point {
		return cellId.rawPoint
	}
	
	/**
		Return the center of the cell in (u,v) coordinates (see `S2Projections`).
		Note that the center of the cell is defined as the point
		at which it is recursively subdivided into four children; in general, it is
		not at the midpoint of the (u,v) rectangle covered by the cell
	*/
	public var centerUV: R2Vector {
		var i = 0
		var j = 0
		var orientation: Int? = nil
		_ = cellId.toFaceIJOrientation(i: &i, j: &j, orientation: &orientation)
		let cellSize = 1 << (S2CellId.maxLevel - Int(level))
		
		// TODO(dbeaumont): Figure out a better naming of the variables here (and elsewhere).
		let si = (i & -cellSize) * 2 + cellSize - S2Cell.maxCellSize
		let x = S2Projections.stToUV(s: (1.0 / Double(S2Cell.maxCellSize)) * Double(si))
		
		let sj = (j & -cellSize) * 2 + cellSize - S2Cell.maxCellSize
		let y = S2Projections.stToUV(s: (1.0 / Double(S2Cell.maxCellSize)) * Double(sj))
		
		return R2Vector(x: x, y: y)
	}
	
	public func contains(point p: S2Point) -> Bool {
		// We can't just call XYZtoFaceUV, because for points that lie on the
		// boundary between two faces (i.e. u or v is +1/-1) we need to return
		// true for both adjacent cells.
		guard let uvPoint = S2Projections.faceXyzToUv(face: Int(face), point: p) else { return false }
		return uvPoint.x >= uv[0][0] && uvPoint.x <= uv[0][1] && uvPoint.y >= uv[1][0] && uvPoint.y <= uv[1][1]
	}
	
	/**
	* Return the average area for cells at the given level.
	*/
	public static func averageArea(level: Int) -> Double {
		return S2Projections.avgArea.getValue(level: level)
	}
	
	/**
		Return the average area of cells at this level. This is accurate to within
		a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely cheap to compute.
	*/
	public var averageArea: Double {
		return S2Cell.averageArea(level: Int(level))
	}
	
	/**
		Return the approximate area of this cell. This method is accurate to within
		3% percent for all cell sizes and accurate to within 0.1% for cells at
		level 5 or higher (i.e. 300km square or smaller). It is moderately cheap to compute.
	*/
	public var approxArea: Double {
		// All cells at the first two levels have the same area.
		if level < 2 { return averageArea }
		
		// First, compute the approximate area of the cell when projected
		// perpendicular to its normal. The cross product of its diagonals gives
		// the normal, and the length of the normal is twice the projected area.
		let flatArea = 0.5 * (getVertex(2) - getVertex(0)).crossProd((getVertex(3) - getVertex(1))).norm
		
		// Now, compensate for the curvature of the cell surface by pretending
		// that the cell is shaped like a spherical cap. The ratio of the
		// area of a spherical cap to the area of its projected disc turns out
		// to be 2 / (1 + sqrt(1 - r*r)) where "r" is the radius of the disc.
		// For example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
		// Here we set Pi*r*r == flat_area to find the equivalent disc.
		return flatArea * 2 / (1 + sqrt(1 - min(M_1_PI * flatArea, 1.0)))
	}
	
	/**
		Return the area of this cell as accurately as possible. This method is more
		expensive but it is accurate to 6 digits of precision even for leaf cells
		(whose area is approximately 1e-18).
	*/
	public var exactArea: Double {
		let v0 = getVertex(0)
		let v1 = getVertex(1)
		let v2 = getVertex(2)
		let v3 = getVertex(3)
		return S2.area(a: v0, b: v1, c: v2) + S2.area(a: v0, b: v2, c: v3)
	}
	
	////////////////////////////////////////////////////////////////////////
	// MARK: S2Region
	////////////////////////////////////////////////////////////////////////
	
	public var capBound: S2Cap {
		// Use the cell center in (u,v)-space as the cap axis. This vector is
		// very close to GetCenter() and faster to compute. Neither one of these
		// vectors yields the bounding cap with minimal surface area, but they
		// are both pretty close.
		//
		// It's possible to show that the two vertices that are furthest from
		// the (u,v)-origin never determine the maximum cap size (this is a
		// possible future optimization).
		
		let u = 0.5 * (uv[0][0] + uv[0][1])
		let v = 0.5 * (uv[1][0] + uv[1][1])
		var cap = S2Cap(axis: S2Point.normalize(point: S2Projections.faceUvToXyz(face: Int(face), u: u, v: v)), height: 0)
		for k in 0 ..< 4 {
			cap = cap.add(point: getVertex(k))
		}
		return cap
	}
	
	// We grow the bounds slightly to make sure that the bounding rectangle
	// also contains the normalized versions of the vertices. Note that the
	// maximum result magnitude is Pi, with a floating-point exponent of 1.
	// Therefore adding or subtracting 2**-51 will always change the result.
	private static let maxError = 1.0 / 2251799813685248
	
	// The 4 cells around the equator extend to +/-45 degrees latitude at the
	// midpoints of their top and bottom edges. The two cells covering the
	// poles extend down to +/-35.26 degrees at their vertices.
	// adding kMaxError (as opposed to the C version) because of asin and atan2
	// roundoff errors
	private static let poleMinLat = asin(sqrt(1.0 / 3.0)) - maxError // 35.26 degrees
	
	public var rectBound: S2LatLngRect {
		if level > 0 {
			// Except for cells at level 0, the latitude and longitude extremes are
			// attained at the vertices. Furthermore, the latitude range is
			// determined by one pair of diagonally opposite vertices and the
			// longitude range is determined by the other pair.
			//
			// We first determine which corner (i,j) of the cell has the largest
			// absolute latitude. To maximize latitude, we want to find the point in
			// the cell that has the largest absolute z-coordinate and the smallest
			// absolute x- and y-coordinates. To do this we look at each coordinate
			// (u and v), and determine whether we want to minimize or maximize that
			// coordinate based on the axis direction and the cell's (u,v) quadrant.
			let u = uv[0][0] + uv[0][1]
			let v = uv[1][0] + uv[1][1]
			let i = S2Projections.getUAxis(face: Int(face)).z == 0 ? (u < 0 ? 1 : 0) : (u > 0 ? 1 : 0)
			let j = S2Projections.getVAxis(face: Int(face)).z == 0 ? (v < 0 ? 1 : 0) : (v > 0 ? 1 : 0)
			
			var lat = R1Interval(p1: getLatitude(i: i, j: j), p2: getLatitude(i: 1 - i, j: 1 - j))
			lat = lat.expanded(radius: S2Cell.maxError).intersection(with: S2LatLngRect.fullLat)
			if (lat.lo == -(Double.pi / 2) || lat.hi == Double.pi / 2) {
				return S2LatLngRect(lat: lat, lng: S1Interval.full)
			}
			let lng = S1Interval(p1: getLongitude(i: i, j: 1 - j), p2: getLongitude(i: 1 - i, j: j))
			return S2LatLngRect(lat: lat, lng: lng.expanded(radius: S2Cell.maxError))
		}
		
		// The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
		// assert (S2Projections.getNorm(face).get(face % 3) == ((face < 3) ? 1 : -1))
		switch face {
		case 0:
			return S2LatLngRect(lat: R1Interval(lo: -(Double.pi / 4), hi: (Double.pi / 4)), lng: S1Interval(lo: -(Double.pi / 4), hi: (Double.pi / 4)))
		case 1:
			return S2LatLngRect(lat: R1Interval(lo: -(Double.pi / 4), hi: (Double.pi / 4)), lng: S1Interval(lo: (Double.pi / 4), hi: 3 * (Double.pi / 4)))
		case 2:
			return S2LatLngRect(lat: R1Interval(lo: S2Cell.poleMinLat, hi: (Double.pi / 2)), lng: S1Interval(lo: -Double.pi, hi: Double.pi))
		case 3:
			return S2LatLngRect(lat: R1Interval(lo: -(Double.pi / 4), hi: (Double.pi / 4)), lng: S1Interval(lo: 3 * (Double.pi / 4), hi: -3 * (Double.pi / 4)))
		case 4:
			return S2LatLngRect(lat: R1Interval(lo: -(Double.pi / 4), hi: (Double.pi / 4)), lng: S1Interval(lo: -3 * (Double.pi / 4), hi: -(Double.pi / 4)))
		default:
			return S2LatLngRect(lat: R1Interval(lo: -(Double.pi / 2), hi: -S2Cell.poleMinLat), lng: S1Interval(lo: -Double.pi, hi: Double.pi))
		}
	}
	
	public func contains(cell: S2Cell) -> Bool {
		return cellId.contains(other: cell.cellId)
	}
	
	public func mayIntersect(cell: S2Cell) -> Bool {
		return cellId.intersects(with: cell.cellId)
	}
	
	// Return the latitude or longitude of the cell vertex given by (i,j),
	// where "i" and "j" are either 0 or 1.
	
	private func getLatitude(i: Int, j: Int) -> Double {
		let p = S2Projections.faceUvToXyz(face: Int(face), u: uv[0][i], v: uv[1][j])
		return atan2(p.z, sqrt(p.x * p.x + p.y * p.y))
	}
	
	private func getLongitude(i: Int, j: Int) -> Double {
		let p = S2Projections.faceUvToXyz(face: Int(face), u: uv[0][i], v: uv[1][j])
		return atan2(p.y, p.x)
	}
	
}

public func ==(lhs: S2Cell, rhs: S2Cell) -> Bool {
	return lhs.face == rhs.face && lhs.level == rhs.level && lhs.orientation == rhs.orientation && lhs.cellId == rhs.cellId
}
