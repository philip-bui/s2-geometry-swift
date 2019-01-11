//
//  S2Projections.swift
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
	This class specifies the details of how the cube faces are projected onto the
	unit sphere. This includes getting the face ordering and orientation correct
	so that sequentially increasing cell ids follow a continuous space-filling
	curve over the entire sphere, and defining the transformation from cell-space
	to cube-space (see s2.h) in order to make the cells more uniform in size.


	We have implemented three different projections from cell-space (s,t) to
	cube-space (u,v): linear, quadratic, and tangent. They have the following
	tradeoffs:

	Linear - This is the fastest transformation, but also produces the least
	uniform cell sizes. Cell areas vary by a factor of about 5.2, with the
	largest cells at the center of each face and the smallest cells in the
	corners.

	Tangent - Transforming the coordinates via atan() makes the cell sizes more
	uniform. The areas vary by a maximum ratio of 1.4 as opposed to a maximum
	ratio of 5.2. However, each call to atan() is about as expensive as all of
	the other calculations combined when converting from points to cell ids, i.e.
	it reduces performance by a factor of 3.

	Quadratic - This is an approximation of the tangent projection that is much
	faster and produces cells that are almost as uniform in size. It is about 3
	times faster than the tangent projection for converting cell ids to points,
	and 2 times faster for converting points to cell ids. Cell areas vary by a
	maximum ratio of about 2.1.

	Here is a table comparing the cell uniformity using each projection. "Area
	ratio" is the maximum ratio over all subdivision levels of the largest cell
	area to the smallest cell area at that level, "edge ratio" is the maximum
	ratio of the longest edge of any cell to the shortest edge of any cell at the
	same level, and "diag ratio" is the ratio of the longest diagonal of any cell
	to the shortest diagonal of any cell at the same level. "ToPoint" and
	"FromPoint" are the times in microseconds required to convert cell ids to and
	from points (unit vectors) respectively.

	Area Edge Diag ToPoint FromPoint Ratio Ratio Ratio (microseconds)
	-------------------------------------------------------
	Linear:		5.200 2.117 2.959 0.103 0.123
	Tangent:	1.414 1.414 1.704 0.290 0.306
	Quadratic:	2.082 1.802 1.932 0.116 0.161

	The worst-case cell aspect ratios are about the same with all three
	projections. The maximum ratio of the longest edge to the shortest edge
	within the same cell is about 1.4 and the maximum ratio of the diagonals
	within the same cell is about 1.7.

	This data was produced using s2cell_unittest and s2cellid_unittest.
*/
public struct S2Projections {
	
	public enum Projection {
		case linear, tan, quadratic
	}
	
	private static let s2Projection: Projection = .quadratic
	
	// All of the values below were obtained by a combination of hand analysis and
	// Mathematica. In general, S2_TAN_PROJECTION produces the most uniform
	// shapes and sizes of cells, S2_LINEAR_PROJECTION is considerably worse, and
	// S2_QUADRATIC_PROJECTION is somewhere in between (but generally closer to
	// the tangent projection than the linear one).
	
	// The minimum area of any cell at level k is at least MIN_AREA.GetValue(k),
	// and the maximum is at most MAX_AREA.GetValue(k). The average area of all
	// cells at level k is exactly AVG_AREA.GetValue(k).
	
	public static let minArea: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 1 / (3 * sqrt(3))				// 0.192
		case .tan:			x = (Double.pi * Double.pi) / (16 * 2.squareRoot())	// 0.436
		case .quadratic:	x = 2 * 2.squareRoot() / 9					// 0.314
		}
		return S2.Metric(dim: 2, deriv: x)
	}()
	
	public static let maxArea: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 1						// 1.000
		case .tan:			x = Double.pi * Double.pi / 16			// 0.617
		case .quadratic:	x = 0.65894981424079037		// 0.659
		}
		return S2.Metric(dim: 2, deriv: x)
	}()
	
	public static let avgArea = S2.Metric(dim: 2, deriv: Double.pi / 6) // 0.524
	
	// Each cell is bounded by four planes passing through its four edges and
	// the center of the sphere. These metrics relate to the angle between each
	// pair of opposite bounding planes, or equivalently, between the planes
	// corresponding to two different s-values or two different t-values. For
	// example, the maximum angle between opposite bounding planes for a cell at
	// level k is MAX_ANGLE_SPAN.GetValue(k), and the average angle span for all
	// cells at level k is approximately AVG_ANGLE_SPAN.GetValue(k).
	
	public static let minAngleSpan: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 0.5						// 0.500
		case .tan:			x = Double.pi / 4					// 0.785
		case .quadratic:	x = 2 / 3					// 0.667
		}
		return S2.Metric(dim: 1, deriv: x)
	}()
	
	public static let maxAngleSpan: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 1						// 1.000
		case .tan:			x = Double.pi / 4					// 0.785
		case .quadratic:	x = 0.85244858959960922		// 0.852
		}
		return S2.Metric(dim: 1, deriv: x)
	}()
	
	public static let avgAngleSpan = S2.Metric(dim: 1, deriv: Double.pi / 4) // 0.785
	
	// The width of geometric figure is defined as the distance between two
	// parallel bounding lines in a given direction. For cells, the minimum
	// width is always attained between two opposite edges, and the maximum
	// width is attained between two opposite vertices. However, for our
	// purposes we redefine the width of a cell as the perpendicular distance
	// between a pair of opposite edges. A cell therefore has two widths, one
	// in each direction. The minimum width according to this definition agrees
	// with the classic geometric one, but the maximum width is different. (The
	// maximum geometric width corresponds to MAX_DIAG defined below.)
	//
	// For a cell at level k, the distance between opposite edges is at least
	// MIN_WIDTH.GetValue(k) and at most MAX_WIDTH.GetValue(k). The average
	// width in both directions for all cells at level k is approximately
	// AVG_WIDTH.GetValue(k).
	//
	// The width is useful for bounding the minimum or maximum distance from a
	// point on one edge of a cell to the closest point on the opposite edge.
	// For example, this is useful when "growing" regions by a fixed distance.
	
	public static let minWidth: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 1 / sqrt(6)				// 0.408
		case .tan:			x = Double.pi / (4 * 2.squareRoot())		// 0.555
		case .quadratic:	x = 2.squareRoot() / 3				// 0.471
		}
		return S2.Metric(dim: 1, deriv: x)
	}()
	
	public static let maxWidth = S2.Metric(dim: 1, deriv: maxAngleSpan.deriv)
	
	public static let avgWidth: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 0.70572967292222848		// 0.706
		case .tan:			x = 0.71865931946258044		// 0.719
		case .quadratic:	x = 0.71726183644304969		// 0.717
		}
		return S2.Metric(dim: 1, deriv: x)
	}()
	
	// The minimum edge length of any cell at level k is at least
	// MIN_EDGE.GetValue(k), and the maximum is at most MAX_EDGE.GetValue(k).
	// The average edge length is approximately AVG_EDGE.GetValue(k).
	//
	// The edge length metrics can also be used to bound the minimum, maximum,
	// or average distance from the center of one cell to the center of one of
	// its edge neighbors. In particular, it can be used to bound the distance
	// between adjacent cell centers along the space-filling Hilbert curve for
	// cells at any given level.
	
	public static let minEdge: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 2.squareRoot() / 3				// 0.471
		case .tan:			x = Double.pi / (4 * 2.squareRoot())		// 0.555
		case .quadratic:	x = 2.squareRoot() / 3				// 0.471
		}
		return S2.Metric(dim: 1, deriv: x)
	}()
	
	public static let maxEdge = S2.Metric(dim: 1, deriv: maxAngleSpan.deriv)
	
	public static let avgEdge: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 0.72001709647780182		// 0.720
		case .tan:			x = 0.73083351627336963		// 0.731
		case .quadratic:	x = 0.72960687319305303		// 0.730
		}
		return S2.Metric(dim: 1, deriv: x)
	}()
	
	// The minimum diagonal length of any cell at level k is at least
	// MIN_DIAG.GetValue(k), and the maximum is at most MAX_DIAG.GetValue(k).
	// The average diagonal length is approximately AVG_DIAG.GetValue(k).
	//
	// The maximum diagonal also happens to be the maximum diameter of any cell,
	// and also the maximum geometric width (see the discussion above). So for
	// example, the distance from an arbitrary point to the closest cell center
	// at a given level is at most half the maximum diagonal length.
	
	public static let minDiag: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 2.squareRoot() / 3				// 0.471
		case .tan:			x = Double.pi / (3 * 2.squareRoot())		// 0.740
		case .quadratic:	x = 4 * 2.squareRoot() / 9			// 0.629
		}
		return S2.Metric(dim: 1, deriv: x)
	}()
	
	public static let maxDiag: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 2.squareRoot()					// 1.414
		case .tan:			x = Double.pi / sqrt(6)			// 1.283
		case .quadratic:	x = 1.2193272972170106		// 1.219
		}
		return S2.Metric(dim: 1, deriv: x)
	}()
	
	public static let avgDiag: S2.Metric = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 1.0159089332094063		// 1.016
		case .tan:			x = 1.0318115985978178		// 1.032
		case .quadratic:	x = 1.03021136949923584		// 1.030
		}
		return S2.Metric(dim: 1, deriv: x)
	}()
	
	// This is the maximum edge aspect ratio over all cells at any level, where
	// the edge aspect ratio of a cell is defined as the ratio of its longest
	// edge length to its shortest edge length.
	public static let maxEdgeAspect: Double = {
		let x: Double
		switch s2Projection {
		case .linear:		x = 2.squareRoot()					// 1.414
		case .tan:			x = 2.squareRoot()					// 1.414
		case .quadratic:	x = 1.44261527445268292		// 1.443
		}
		return x
	}()
	
	// This is the maximum diagonal aspect ratio over all cells at any level,
	// where the diagonal aspect ratio of a cell is defined as the ratio of its
	// longest diagonal length to its shortest diagonal length.
	public static let maxDiagAspect: Double = sqrt(3)	// 1.732
	
	public static func stToUV(s: Double) -> Double {
		switch s2Projection {
		case .linear:
			return s
		case .tan:
			// Unfortunately, tan(Double.pi_4) is slightly less than 1.0. This isn't due
			// to a flaw in the implementation of tan(), it's because the derivative of
			// tan(x) at x=pi/4 is 2, and it happens that the two adjacent floating
			// point numbers on either side of the infinite-precision value of pi/4
			// have tangents that are slightly below and slightly above 1.0 when rounded
			// to the nearest double-precision result.
			let x = tan(Double.pi / 4 * s)
            
			return x + (1.0 / 9007199254740992) * x
		case .quadratic:
			if s >= 0 {
				return (1 / 3) * ((1 + s) * (1 + s) - 1)
			} else {
				return (1 / 3) * (1 - (1 - s) * (1 - s))
			}
		}
	}
	
	public static func uvToST(u: Double) -> Double {
		switch s2Projection {
		case .linear:
			return u
		case .tan:
			return (4 * M_1_PI) * atan(u)
		case .quadratic:
			if u >= 0 {
				return sqrt(1 + 3 * u) - 1
			} else {
				return 1 - sqrt(1 - 3 * u)
			}
		}
	}
	
	/// Convert (face, u, v) coordinates to a direction vector (not necessarily unit length).
	public static func faceUvToXyz(face: Int, u: Double, v: Double) -> S2Point {
		switch face {
		case 0:
			return S2Point(x: 1, y: u, z: v)
		case 1:
			return S2Point(x: -u, y: 1, z: v)
		case 2:
			return S2Point(x: -u, y: -v, z: 1)
		case 3:
			return S2Point(x: -1, y: -v, z: -u)
		case 4:
			return S2Point(x: v, y: -1, z: -u)
		default:
			return S2Point(x: v, y: u, z: -1)
		}
	}
	
	public static func validFaceXyzToUv(face: Int, point p: S2Point) -> R2Vector {
		let pu: Double, pv: Double
		switch face {
		case 0:
			pu = p.y / p.x
			pv = p.z / p.x
		case 1:
			pu = -p.x / p.y
			pv = p.z / p.y
		case 2:
			pu = -p.x / p.z
			pv = -p.y / p.z
		case 3:
			pu = p.z / p.x
			pv = p.y / p.x
		case 4:
			pu = p.z / p.y
			pv = -p.x / p.y
		default:
			pu = -p.y / p.z
			pv = -p.x / p.z
		}
		return R2Vector(x: pu, y: pv)
	}
	
	public static func xyzToFace(point p: S2Point) -> Int {
		var face = p.largestAbsComponent
		if p.get(axis: face) < 0 {
			face += 3
		}
		return face
	}

	public static func faceXyzToUv(face: Int, point p: S2Point) -> R2Vector? {
		if face < 3 {
			if p.get(axis: face) <= 0 {
				return nil
			}
		} else {
			if p.get(axis: face - 3) >= 0 {
				return nil
			}
		}
		return validFaceXyzToUv(face: face, point: p)
	}
	
	public static func getUNorm(face: Int, u: Double) -> S2Point {
		switch face {
		case 0:
			return S2Point(x: u, y: -1, z: 0)
		case 1:
			return S2Point(x: 1, y: u, z: 0)
		case 2:
			return S2Point(x: 1, y: 0, z: u)
		case 3:
			return S2Point(x: -u, y: 0, z: 1)
		case 4:
			return S2Point(x: 0, y: -u, z: 1)
		default:
			return S2Point(x: 0, y: -1, z: -u)
		}
	}
	
	public static func getVNorm(face: Int, v: Double) -> S2Point {
		switch face {
		case 0:
			return S2Point(x: -v, y: 0, z: 1)
		case 1:
			return S2Point(x: 0, y: -v, z: 1)
		case 2:
			return S2Point(x: 0, y: -1, z: -v)
		case 3:
			return S2Point(x: v, y: -1, z: 0)
		case 4:
			return S2Point(x: 1, y: v, z: 0)
		default:
			return S2Point(x: 1, y: 0, z: v)
		}
	}
	
	public static func getNorm(face: Int) -> S2Point {
		return faceUvToXyz(face: face, u: 0, v: 0)
	}
	
	public static func getUAxis(face: Int) -> S2Point {
		switch face {
		case 0:
			return S2Point(x: 0, y: 1, z: 0)
		case 1:
			return S2Point(x: -1, y: 0, z: 0)
		case 2:
			return S2Point(x: -1, y: 0, z: 0)
		case 3:
			return S2Point(x: 0, y: 0, z: -1)
		case 4:
			return S2Point(x: 0, y: 0, z: -1)
		default:
			return S2Point(x: 0, y: 1, z: 0)
		}
	}
	
	public static func getVAxis(face: Int) -> S2Point {
		switch face {
		case 0:
			return S2Point(x: 0, y: 0, z: 1)
		case 1:
			return S2Point(x: 0, y: 0, z: 1)
		case 2:
			return S2Point(x: 0, y: -1, z: 0)
		case 3:
			return S2Point(x: 0, y: -1, z: 0)
		case 4:
			return S2Point(x: 1, y: 0, z: 0)
		default:
			return S2Point(x: 1, y: 0, z: 0)
		}
	}
	
	// Don't instantiate
	private init() { }

}
