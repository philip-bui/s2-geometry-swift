//
//  S1Interval.swift
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
	An S1Interval represents a closed interval on a unit circle (also known as a
	1-dimensional sphere). It is capable of representing the empty interval
	(containing no points), the full interval (containing all points), and
	zero-length intervals (containing a single point).

	Points are represented by the angle they make with the positive x-axis in
	the range [-Pi, Pi]. An interval is represented by its lower and upper bounds
	(both inclusive, since the interval is closed). The lower bound may be
	greater than the upper bound, in which case the interval is "inverted" (i.e.
	it passes through the point (-1, 0)).

	Note that the point (-1, 0) has two valid representations, Pi and -Pi. The
	normalized representation of this point internally is Pi, so that endpoints
	of normal intervals are in the range (-Pi, Pi]. However, we take advantage of
	the point -Pi to construct two special intervals: the Full() interval is
	[-Pi, Pi], and the Empty() interval is [Pi, -Pi].
*/
public struct S1Interval: Equatable {
	
	public let lo: Double
	public let hi: Double
	
	/**
		Both endpoints must be in the range -Pi to Pi inclusive. The value -Pi is
		converted internally to Pi except for the Full() and Empty() intervals.
	*/
	public init(lo: Double, hi: Double) {
		self.init(lo: lo, hi: hi, checked: false)
	}
	
	/**
		Internal constructor that assumes that both arguments are in the correct
		range, i.e. normalization from -Pi to Pi is already done.
	*/
	private init(lo: Double, hi: Double, checked: Bool) {
		var newLo = lo
		var newHi = hi
		if !checked {
			if lo == -.pi && hi != .pi {
				newLo = .pi
			}
			if hi == -.pi && lo != .pi {
				newHi = .pi
			}
		}
		self.lo = newLo
		self.hi = newHi
	}
	
	public static let empty = S1Interval(lo: .pi, hi: -.pi, checked: true)
	
	public static let full = S1Interval(lo: -.pi, hi: .pi, checked: true)
	
	/// Convenience method to construct an interval containing a single point.
	public init(point p: Double) {
		var p = p
		if p == -.pi { p = .pi }
		self.init(lo: p, hi: p)
	}
	
	/**
		Convenience method to construct the minimal interval containing the two
		given points. This is equivalent to starting with an empty interval and
		calling AddPoint() twice, but it is more efficient.
	*/
	public init(p1: Double, p2: Double) {
		var p1 = p1, p2 = p2
		if p1 == -.pi { p1 = .pi }
		if p2 == -.pi { p2 = .pi }
		if S1Interval.positiveDistance(p1, p2) <= .pi {
			self.init(lo: p1, hi: p2, checked: true)
		} else {
			self.init(lo: p2, hi: p1, checked: true)
		}
	}
	
	/**
		An interval is valid if neither bound exceeds Pi in absolute value, and the
		value -Pi appears only in the empty and full intervals.
	*/
	public var isValid: Bool {
		return (abs(lo) <= .pi && abs(hi) <= .pi && !(lo == -.pi && hi != .pi) && !(hi == -.pi && lo != .pi))
	}
	
	/// Return true if the interval contains all points on the unit circle.
	public var isFull: Bool {
		return hi - lo == 2 * .pi
	}
	
	/// Return true if the interval is empty, i.e. it contains no points.
	public var isEmpty: Bool {
		return lo - hi == 2 * .pi
	}
	
	/// Return true if lo() > hi(). (This is true for empty intervals.)
	public var isInverted: Bool {
		return lo > hi
	}
	
	/// Return the midpoint of the interval. For full and empty intervals, the result is arbitrary.
	public var center: Double {
		let center = 0.5 * (lo + hi)
		if !isInverted { return center }
		// Return the center in the range (-Pi, Pi].
		return (center <= 0) ? (center + .pi) : (center - .pi)
	}
	
	/// Return the length of the interval. The length of an empty interval is negative.
	public var length: Double {
		var length = hi - lo
		if length >= 0 { return length }
		length += 2 * .pi
		// Empty intervals have a negative length.
		return (length > 0) ? length : -1
	}
	
	/**
		Return the complement of the interior of the interval. An interval and its
		complement have the same boundary but do not share any interior values. The
		complement operator is not a bijection, since the complement of a singleton
		interval (containing a single value) is the same as the complement of an
		empty interval.
	*/
	public var complement: S1Interval {
		if lo == hi {
			return S1Interval.full // Singleton.
		} else {
			return S1Interval(lo: hi, hi: lo, checked: true) // Handles empty and full.
		}
	}
	
	/// Return true if the interval (which is closed) contains the point 'p'.
	public func contains(point p: Double) -> Bool {
		// Works for empty, full, and singleton intervals.
		// assert (Math.abs(p) <= S2..pi);
		var p = p
		if (p == -.pi) { p = .pi }
		return fastContains(point: p)
	}
	
	/**
		Return true if the interval (which is closed) contains the point 'p'. Skips
		the normalization of 'p' from -Pi to Pi.
	*/
	public func fastContains(point p: Double) -> Bool {
		if isInverted {
			return (p >= lo || p <= hi) && !isEmpty
		} else {
			return p >= lo && p <= hi
		}
	}
	
	/// Return true if the interval contains the given interval 'y'. Works for empty, full, and singleton intervals.
	public func contains(interval y: S1Interval) -> Bool {
		// It might be helpful to compare the structure of these tests to
		// the simpler Contains(double) method above.
		
		if isInverted {
			if y.isInverted {
				return y.lo >= lo && y.hi <= hi
			}
			return (y.lo >= lo || y.hi <= hi) && !isEmpty
		} else {
			if y.isInverted {
				return isFull || y.isEmpty
			}
			return y.lo >= lo && y.hi <= hi
		}
	}
	
	/**
		Returns true if the interior of this interval contains the entire interval
		'y'. Note that x.InteriorContains(x) is true only when x is the empty or
		full interval, and x.InteriorContains(S1Interval(p,p)) is equivalent to
		x.InteriorContains(p).
	*/
	public func interiorContains(interval y: S1Interval) -> Bool {
		if isInverted {
			if !y.isInverted {
				return y.lo > lo || y.hi < hi
			}
			return (y.lo > lo && y.hi < hi) || y.isEmpty
		} else {
			if y.isInverted {
				return isFull || y.isEmpty
			}
			return (y.lo > lo && y.hi < hi) || isFull
		}
	}
	
	/// Return true if the interior of the interval contains the point 'p'.
	public func interiorContains(point p: Double) -> Bool {
		// Works for empty, full, and singleton intervals.
		// assert (Math.abs(p) <= S2..pi);
		var p = p
		if (p == -.pi) { p = .pi }
		
		if isInverted {
			return p > lo || p < hi
		} else {
			return (p > lo && p < hi) || isFull
		}
	}
	
	/**
		Return true if the two intervals contain any points in common. Note that
		the point +/-Pi has two representations, so the intervals [-Pi,-3] and
		[2,Pi] intersect, for example.
	*/
	public func intersects(with y: S1Interval) -> Bool {
		if isEmpty || y.isEmpty {
			return false
		}
		if isInverted {
			// Every non-empty inverted interval contains Pi.
			return y.isInverted || y.lo <= hi || y.hi >= lo
		} else {
			if y.isInverted {
				return y.lo <= hi || y.hi >= lo
			}
			return y.lo <= hi && y.hi >= lo
		}
	}
	
	/**
		Return true if the interior of this interval contains any point of the
		interval 'y' (including its boundary). Works for empty, full, and singleton
		intervals.
	*/
	public func interiorIntersects(with y: S1Interval) -> Bool {
		if isEmpty || y.isEmpty || lo == hi {
			return false
		}
		if isInverted {
			return y.isInverted || y.lo < hi || y.hi > lo
		} else {
			if y.isInverted {
				return y.lo < hi || y.hi > lo
			}
			return (y.lo < hi && y.hi > lo) || isFull
		}
	}
	
	/**
		Expand the interval by the minimum amount necessary so that it contains the
		given point "p" (an angle in the range [-Pi, Pi]).
	*/
	public func add(point p: Double) -> S1Interval {
		// assert (Math.abs(p) <= S2..pi);
		var p = p
		if p == -.pi {
			p = .pi
		}
		
		if fastContains(point: p) {
			return self
		}
		
		if isEmpty {
			return S1Interval(point: p)
		} else {
			// Compute distance from p to each endpoint.
			let dlo = S1Interval.positiveDistance(p, lo)
			let dhi = S1Interval.positiveDistance(hi, p)
			if dlo < dhi {
				return S1Interval(lo: p, hi: hi)
			} else {
				return S1Interval(lo: lo, hi: p)
			}
			// Adding a point can never turn a non-full interval into a full one.
		}
	}
	
	/**
		Return an interval that contains all points within a distance "radius" of
		a point in this interval. Note that the expansion of an empty interval is
		always empty. The radius must be non-negative.
	*/
	public func expanded(radius: Double) -> S1Interval {
		// assert (radius >= 0)
		guard !isEmpty else { return self }
		
		// Check whether this interval will be full after expansion, allowing
		// for a 1-bit rounding error when computing each endpoint.
		if length + 2 * radius >= 2 * .pi - 1e-15 {
			return .full
		}
		
		// NOTE(dbeaumont): Should this remainder be 2 * .pi or just .pi ??
		var lo = remainder(self.lo - radius, 2 * .pi)
		let hi = remainder(self.hi + radius, 2 * .pi)
		if lo == -.pi {
			lo = .pi
		}
		return S1Interval(lo: lo, hi: hi)
	}
	
	/// Return the smallest interval that contains this interval and the given interval "y".
	public func union(with y: S1Interval) -> S1Interval {
		// The y.is_full() case is handled correctly in all cases by the code
		// below, but can follow three separate code paths depending on whether
		// this interval is inverted, is non-inverted but contains Pi, or neither.
		
		if y.isEmpty {
			return self
		}
		if fastContains(point: y.lo) {
			if fastContains(point: y.hi) {
				// Either this interval contains y, or the union of the two
				// intervals is the Full() interval.
				if contains(interval: y) {
					return self // is_full() code path
				}
				return .full
			}
			return S1Interval(lo: lo, hi: y.hi, checked: true)
		}
		if (fastContains(point: y.hi)) {
			return S1Interval(lo: y.lo, hi: hi, checked: true)
		}
		
		// This interval contains neither endpoint of y. This means that either y
		// contains all of this interval, or the two intervals are disjoint.
		if isEmpty || y.fastContains(point: lo) {
			return y
		}
		
		// Check which pair of endpoints are closer together.
		let dlo = S1Interval.positiveDistance(y.hi, lo)
		let dhi = S1Interval.positiveDistance(hi, y.lo)
		if dlo < dhi {
			return S1Interval(lo: y.lo, hi: hi, checked: true)
		} else {
			return S1Interval(lo: lo, hi: y.hi, checked: true)
		}
	}
	
	/**
		Compute the distance from "a" to "b" in the range [0, 2*Pi). This is
		equivalent to (drem(b - a - S2..pi, 2 * S2..pi) + S2..pi), except that
		it is more numerically stable (it does not lose precision for very small
		positive distances).
	*/
	public static func positiveDistance(_ a: Double, _ b: Double) -> Double {
		let d = b - a
		if d >= 0 { return d }
		// We want to ensure that if b == Pi and a == (-Pi + eps),
		// the return result is approximately 2*Pi and not zero.
		return (b + .pi) - (a - .pi)
	}
	
}

public func ==(lhs: S1Interval, rhs: S1Interval ) -> Bool {
	return (lhs.lo == rhs.lo && lhs.hi == rhs.hi) || (lhs.isEmpty && rhs.isEmpty)
}
