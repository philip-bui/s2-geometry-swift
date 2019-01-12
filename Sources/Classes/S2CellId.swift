//
//  S2CellId.swift
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
	An S2CellId is a 64-bit unsigned integer that uniquely identifies a cell in
	the S2 cell decomposition. It has the following format:

	```
	id = [face][face_pos]
	```

	face: a 3-bit number (range 0..5) encoding the cube face.

	face_pos: a 61-bit number encoding the position of the center of this cell
	along the Hilbert curve over this face (see the Wiki pages for details).

	Sequentially increasing cell ids follow a continuous space-filling curve over
	the entire sphere. They have the following properties:
	- The id of a cell at level k consists of a 3-bit face number followed by k
	bit pairs that recursively select one of the four children of each cell. The
	next bit is always 1, and all other bits are 0. Therefore, the level of a
	cell is determined by the position of its lowest-numbered bit that is turned
	on (for a cell at level k, this position is 2 * (MAX_LEVEL - k).)
	- The id of a parent cell is at the midpoint of the range of ids spanned by
	its children (or by its descendants at any level).

	Leaf cells are often used to represent points on the unit sphere, and this
	class provides methods for converting directly between these two
	representations. For cells that represent 2D regions rather than discrete
	point, it is better to use the S2Cell class.
*/
public struct S2CellId: Comparable, Hashable {
	
	// Although only 60 bits are needed to represent the index of a leaf
	// cell, we need an extra bit in order to represent the position of
	// the center of the leaf cell along the Hilbert curve.
	public static let faceBits = 3
	public static let numFaces = 6
	public static let maxLevel = 30 // Valid levels: 0..MAX_LEVEL
	public static let posBits = 2 * maxLevel + 1
	public static let maxSize = 1 << maxLevel
	
	// The following lookup tables are used to convert efficiently between an
	// (i,j) cell index and the corresponding position along the Hilbert curve.
	// "lookup_pos" maps 4 bits of "i", 4 bits of "j", and 2 bits representing the
	// orientation of the current cell into 8 bits representing the order in which
	// that subcell is visited by the Hilbert curve, plus 2 bits indicating the
	// new orientation of the Hilbert curve within that subcell. (Cell
	// orientations are represented as combination of kSwapMask and kInvertMask.)
	//
	// "lookup_ij" is an inverted table used for mapping in the opposite
	// direction.
	//
	// We also experimented with looking up 16 bits at a time (14 bits of position
	// plus 2 of orientation) but found that smaller lookup tables gave better
	// performance. (2KB fits easily in the primary cache.)
	
	// Values for these constants are *declared* in the *.h file. Even though
	// the declaration specifies a value for the constant, that declaration
	// is not a *definition* of storage for the value. Because the values are
	// supplied in the declaration, we don't need the values here. Failing to
	// define storage causes link errors for any code that tries to take the
	// address of one of these values.
	private static let lookupBits = 4
	private static let swapMask = 0x01
	private static let invertMask = 0x02
	
	private static var lookupLoaded = false
	
	private static var lookup: (pos: [Int], ij: [Int]) = {
		return (pos: Array(repeating: 0, count: 1 << (2 * lookupBits + 2)), ij: Array(repeating: 0, count: 1 << (2 * lookupBits + 2)))
	}()
	
	private static var lookupPos: [Int] {
		loadLookup()
		return lookup.pos
	}
	
	private static var lookupIj: [Int]  {
		loadLookup()
		return lookup.ij
	}
	
	private static func loadLookup() {
		guard !lookupLoaded else { return }
		S2CellId.initLookupCell(level: 0, i: 0, j: 0, origOrientation: 0, pos: 0, orientation: 0)
		S2CellId.initLookupCell(level: 0, i: 0, j: 0, origOrientation: swapMask, pos: 0, orientation: swapMask)
		S2CellId.initLookupCell(level: 0, i: 0, j: 0, origOrientation: invertMask, pos: 0, orientation: invertMask)
		S2CellId.initLookupCell(level: 0, i: 0, j: 0, origOrientation: swapMask | invertMask, pos: 0, orientation: swapMask | invertMask)
		lookupLoaded = true
	}
	
	/**
		This is the offset required to wrap around from the beginning of the
		Hilbert curve to the end or vice versa; see next_wrap() and prev_wrap().
	*/
	private static let wrapOffset = Int64((numFaces) << posBits)
	
	public let id: Int64
	
	public var uid: UInt64 {
		return UInt64(bitPattern: id)
	}
	
	public init(id: Int64 = 0) {
		self.id = id
	}
	
	public init(uid: UInt64) {
		self.id = Int64(bitPattern: uid)
	}
	
	public static let none = S2CellId()
	
	public static let sentinel = S2CellId(id: .max)
	
	/**
		Return a cell given its face (range 0..5), 61-bit Hilbert curve position
		within that face, and level (range 0..MAX_LEVEL). The given position will
		be modified to correspond to the Hilbert curve position at the center of
		the returned cell. This is a static function rather than a constructor in
		order to give names to the arguments.
	*/
	public init(face: Int, pos: Int64, level: Int) {
		let id = Int64(face << S2CellId.posBits) + (pos | 1)
		self = S2CellId(id: id).parent(level: level)
	}
	
	/// Return the leaf cell containing the given point (a direction vector, not necessarily unit length).
	public init(point p: S2Point) {
		let face = S2Projections.xyzToFace(point: p)
		let uv = S2Projections.validFaceXyzToUv(face: face, point: p)
		let i = S2CellId.stToIJ(s: S2Projections.uvToST(u: uv.x))
		let j = S2CellId.stToIJ(s: S2Projections.uvToST(u: uv.y))
		self.init(face: face, i: i, j: j)
	}
	
	/// Return the leaf cell containing the given S2LatLng.
	public init(latlng: S2LatLng) {
		self.init(point: latlng.point)
	}
	
	public var point: S2Point {
		return S2Point.normalize(point: rawPoint)
	}
	
	/**
		Return the direction vector corresponding to the center of the given cell.
		The vector returned by ToPointRaw is not necessarily unit length.
	*/
	public var rawPoint: S2Point {
		// First we compute the discrete (i,j) coordinates of a leaf cell contained
		// within the given cell. Given that cells are represented by the Hilbert
		// curve position corresponding at their center, it turns out that the cell
		// returned by ToFaceIJOrientation is always one of two leaf cells closest
		// to the center of the cell (unless the given cell is a leaf cell itself,
		// in which case there is only one possibility).
		//
		// Given a cell of size s >= 2 (i.e. not a leaf cell), and letting (imin,
		// jmin) be the coordinates of its lower left-hand corner, the leaf cell
		// returned by ToFaceIJOrientation() is either (imin + s/2, jmin + s/2)
		// (imin + s/2 - 1, jmin + s/2 - 1). We can distinguish these two cases by
		// looking at the low bit of "i" or "j". In the first case the low bit is
		// zero, unless s == 2 (i.e. the level just above leaf cells) in which case
		// the low bit is one.
		//
		// The following calculation converts (i,j) to the (si,ti) coordinates of
		// the cell center. (We need to multiply the coordinates by a factor of 2
		// so that the center of leaf cells can be represented exactly.)
		
		var i = 0
		var j = 0
		var orientation: Int? = nil
		let face = toFaceIJOrientation(i: &i, j: &j, orientation: &orientation)
		// System.out.println("i= " + i.intValue() + " j = " + j.intValue())
		let delta = isLeaf ? 1 : (((i ^ (Int(id) >> 2)) & 1) != 0) ? 2 : 0
		let si = (i << 1) + delta - S2CellId.maxSize
		let ti = (j << 1) + delta - S2CellId.maxSize
		return S2CellId.faceSiTiToXYZ(face: face, si: si, ti: ti)
	}

	/// Return the S2LatLng corresponding to the center of the given cell.
	public var latLng: S2LatLng {
		return S2LatLng(point: rawPoint)
	}
	
	/// Return true if id() represents a valid cell.
	public var isValid: Bool {
		return face < S2CellId.numFaces && ((lowestOnBit & (0x1555555555555555)) != 0)
	}
	
	/// Which cube face this cell belongs to, in the range 0..5.
	public var face: Int {
		return Int(uid >> UInt64(S2CellId.posBits))
	}
	
	/// The position of the cell center along the Hilbert curve over this face, in the range 0..(2**kPosBits-1).
	public var pos: Int64 {
		return Int64(bitPattern: uid & (UInt64(bitPattern: -1) >> 3))
	}
	
	/// Return the subdivision level of the cell (range 0..MAX_LEVEL).
	public var level: Int {
		// Fast path for leaf cells.
		guard !isLeaf else { return S2CellId.maxLevel }
        var x = Int32(truncatingIfNeeded: id)
		var level = -1
		if x != 0 {
			level += 16
		} else {
			x = Int32(truncatingIfNeeded: id >> Int64(32))
		}
		// We only need to look at even-numbered bits to determine the
		// level of a valid cell id.
		x &= -x // Get lowest bit.
		if ((x & 0x00005555) != 0) {
			level += 8
		}
		if ((x & 0x00550055) != 0) {
			level += 4
		}
		if ((x & 0x05050505) != 0) {
			level += 2
		}
		if ((x & 0x11111111) != 0) {
			level += 1
		}
		// assert (level >= 0 && level <= MAX_LEVEL);
		return level
	}
	
	/// Return true if this is a leaf cell (more efficient than checking whether level() == MAX_LEVEL).
	public var isLeaf: Bool {
		return (id & 1) != 0
	}
	
	/// Return true if this is a top-level face cell (more efficient than checking whether level() == 0).
	public var isFace: Bool {
		return (id & (S2CellId.lowestOnBit(forLevel: 0) - 1)) == 0
	}
	
	/**
		Return the child position (0..3) of this cell's ancestor at the given
		level, relative to its parent. The argument should be in the range
		1..MAX_LEVEL. For example, child_position(1) returns the position of this
		cell's level-1 ancestor within its top-level face cell.
	*/
	public func childPosition(level: Int) -> Int {
		return Int(id >> Int64(2 * (S2CellId.maxLevel - level) + 1)) & 3
	}
	
	// Methods that return the range of cell ids that are contained
	// within this cell (including itself). The range is *inclusive*
	// (i.e. test using >= and <=) and the return values of both
	// methods are valid leaf cell ids.
	//
	// These methods should not be used for iteration. If you want to
	// iterate through all the leaf cells, call child_begin(MAX_LEVEL) and
	// child_end(MAX_LEVEL) instead.
	//
	// It would in fact be error-prone to define a range_end() method,
	// because (range_max().id() + 1) is not always a valid cell id, and the
	// iterator would need to be tested using "<" rather that the usual "!=".
	public var rangeMin: S2CellId {
		return S2CellId(id: id - (lowestOnBit - 1))
	}
	
	public var rangeMax: S2CellId {
		return S2CellId(id: id + (lowestOnBit - 1))
	}
	
	/// Return true if the given cell is contained within this one.
	public func contains(other: S2CellId) -> Bool {
		// assert (isValid() && other.isValid());
		return other >= rangeMin && other <= rangeMax
	}
	
	/// Return true if the given cell intersects this one.
	public func intersects(with other: S2CellId) -> Bool {
		// assert (isValid() && other.isValid());
		return other.rangeMin <= rangeMax && other.rangeMax >= rangeMin
	}
	
	public var parent: S2CellId {
		// assert (isValid() && level() > 0);
		let newLsb = lowestOnBit << 2
		return S2CellId(id: (id & -newLsb) | newLsb);
	}
	
	/**
		Return the cell at the previous level or at the given level (which must be
		less than or equal to the current level).
	*/
	public func parent(level: Int) -> S2CellId {
		// assert (isValid() && level >= 0 && level <= this.level());
		let newLsb = S2CellId.lowestOnBit(forLevel: level)
		return S2CellId(id: (id & -newLsb) | newLsb)
	}
	
	public func childBegin() -> S2CellId {
		// assert (isValid() && level() < MAX_LEVEL);
		let oldLsb = UInt64(bitPattern: lowestOnBit)
		return S2CellId(uid: uid - oldLsb + (oldLsb >> 2))
	}
	
	public func childBegin(level: Int) -> S2CellId {
		// assert (isValid() && level >= this.level() && level <= MAX_LEVEL);
		let lsb = UInt64(bitPattern: lowestOnBit)
		let lsbLevel = UInt64(bitPattern: S2CellId.lowestOnBit(forLevel: level))
		return S2CellId(uid: uid - lsb + lsbLevel)
	}
	
	public func childEnd() -> S2CellId {
		// assert (isValid() && level() < MAX_LEVEL);
		let oldLsb = UInt64(bitPattern: lowestOnBit)
		return S2CellId(uid: uid + oldLsb + (oldLsb >> 2))
	}
	
	public func childEnd(level: Int) -> S2CellId {
		// assert (isValid() && level >= this.level() && level <= MAX_LEVEL);
		let lsb = UInt64(bitPattern: lowestOnBit)
		let lsbLevel = UInt64(bitPattern: S2CellId.lowestOnBit(forLevel: level))
		return S2CellId(uid: uid + lsb + lsbLevel)
	}
	
	// Iterator-style methods for traversing the immediate children of a cell or
	// all of the children at a given level (greater than or equal to the current
	// level). Note that the end value is exclusive, just like standard STL
	// iterators, and may not even be a valid cell id. You should iterate using
	// code like this:
	//
	// for(S2CellId c = id.childBegin(); !c.equals(id.childEnd()); c = c.next())
	// ...
	//
	// The convention for advancing the iterator is "c = c.next()", so be sure
	// to use 'equals()' in the loop guard, or compare 64-bit cell id's,
	// rather than "c != id.childEnd()".
	
	/**
		Return the next cell at the same level along the Hilbert curve. Works
		correctly when advancing from one face to the next, but does *not* wrap
		around from the last face to the first or vice versa.
	*/
	public func next() -> S2CellId {
		
		
		return S2CellId(id: id.addingReportingOverflow(lowestOnBit << 1).partialValue)
	}
	
	/**
		Return the previous cell at the same level along the Hilbert curve. Works
		correctly when advancing from one face to the next, but does *not* wrap
		around from the last face to the first or vice versa.
	*/
	public func prev() -> S2CellId {
        return S2CellId(id: id.subtractingReportingOverflow(lowestOnBit << 1).partialValue)
	}
	
	/**
		Like next(), but wraps around from the last face to the first and vice
		versa. Should *not* be used for iteration in conjunction with
		child_begin(), child_end(), Begin(), or End().
	*/
	public func nextWrap() -> S2CellId {
		let n = next()
		if S2CellId.unsignedLongLessThan(lhs: n.id, rhs: S2CellId.wrapOffset) { return n }
		return S2CellId(id: n.id - S2CellId.wrapOffset)
	}
	
	/**
		Like prev(), but wraps around from the last face to the first and vice
		versa. Should *not* be used for iteration in conjunction with
		child_begin(), child_end(), Begin(), or End().
	*/
	public func prevWrap() -> S2CellId {
		let p = prev()
		if p.id < S2CellId.wrapOffset { return p }
		return S2CellId(id: p.id + S2CellId.wrapOffset)
	}
	
	public static func begin(level: Int) -> S2CellId {
		return S2CellId(face: 0, pos: 0, level: 0).childBegin(level: level)
	}
	
	public static func end(level: Int) -> S2CellId {
		return S2CellId(face: 5, pos: 0, level: 0).childEnd(level: level)
	}

	public struct NumberFormatException: Error {
		public let description: String
		public init(_ description: String = "") {
			self.description = description
		}
	}
	
	/**
		Decodes the cell id from a compact text string suitable for display or
		indexing. Cells at lower levels (i.e. larger cells) are encoded into
		fewer characters. The maximum token length is 16.
	*/
	public init(token: String) throws {
		let chars = Array(token)

		guard chars.count > 0 else {
			throw NumberFormatException("Empty string in S2CellId.fromToken")
		}
	
		guard chars.count <= 16 && token != "X" else {
			self = .none
			return
		}
		
		var value: UInt64 = 0
		for pos in 0 ..< 16 {
			var digit: Int = 0
			if pos < chars.count {
				digit = Int(strtoul(String(chars[pos]), nil, 16))
				if digit == -1 {
					throw NumberFormatException(token)
				}
				if S2CellId.overflowInParse(current: value, digit: digit) {
					throw NumberFormatException("Too large for unsigned long: " + token)
				}
			}
			value = (value * 16) + UInt64(digit)
		}
		self = S2CellId(id: Int64(bitPattern: value))
	}
	
	/**
		Encodes the cell id to compact text strings suitable for display or indexing.
		Cells at lower levels (i.e. larger cells) are encoded into fewer characters.
		The maximum token length is 16.
	
		Simple implementation: convert the id to hex and strip trailing zeros. We
		could use base-32 or base-64, but assuming the cells used for indexing
		regions are at least 100 meters across (level 16 or less), the savings
		would be at most 3 bytes (9 bytes hex vs. 6 bytes base-64).
	*/
	public var token: String {
		guard id != 0 else { return "X" }
		var hex = String(uid, radix: 16).lowercased()
		for _ in 0 ..< max(0, 16 - hex.count) {
			hex.insert("0", at: hex.startIndex)
		}
		
		var len = hex.count
		for char in hex.reversed() {
			guard char == "0" else { break }
			len -= 1
		}
		return String(hex.prefix(len))
	}

	/**
		Returns true if (current * radix) + digit is a number too large to be
		represented by an unsigned long.  This is useful for detecting overflow
		while parsing a string representation of a number.
		Does not verify whether supplied radix is valid, passing an invalid radix
		will give undefined results or an ArrayIndexOutOfBoundsException.
	*/
	private static func overflowInParse(current: UInt64, digit: Int, radix: Int = 10) -> Bool {
		if current >= 0 {
			if current < maxValueDivs[radix] {
				return false
			}
			if current > maxValueDivs[radix] {
				return true
			}
			// current == maxValueDivs[radix]
			return digit > maxValueMods[radix]
		}

		// current < 0: high bit is set
		return true
	}
	
	// calculated as 0xffffffffffffffff / radix
	private static let maxValueDivs: [UInt64] = [0, 0, // 0 and 1 are invalid
		9223372036854775807, 6148914691236517205, 4611686018427387903, // 2-4
		3689348814741910323, 3074457345618258602, 2635249153387078802, // 5-7
		2305843009213693951, 2049638230412172401, 1844674407370955161, // 8-10
		1676976733973595601, 1537228672809129301, 1418980313362273201, // 11-13
		1317624576693539401, 1229782938247303441, 1152921504606846975, // 14-16
		1085102592571150095, 1024819115206086200, 970881267037344821, // 17-19
		922337203685477580, 878416384462359600, 838488366986797800, // 20-22
		802032351030850070, 768614336404564650, 737869762948382064, // 23-25
		709490156681136600, 683212743470724133, 658812288346769700, // 26-28
		636094623231363848, 614891469123651720, 595056260442243600, // 29-31
		576460752303423487, 558992244657865200, 542551296285575047, // 32-34
		527049830677415760, 512409557603043100 ] // 35-36
	
	// calculated as 0xffffffffffffffff % radix
	private static let maxValueMods: [Int] = [0, 0, // 0 and 1 are invalid
		1, 0, 3, 0, 3, 1, 7, 6, 5, 4, 3, 2, 1, 0, 15, 0, 15, 16, 15, 15, // 2-21
		15, 5, 15, 15, 15, 24, 15, 23, 15, 15, 31, 15, 17, 15, 15 ] // 22-36
	
	/**
		Return the four cells that are adjacent across the cell's four edges.
		Neighbors are returned in the order defined by S2Cell::GetEdge. All
		neighbors are guaranteed to be distinct.
	*/
	public func getEdgeNeighbors() -> [S2CellId] {
		let level = self.level
		let size = 1 << (S2CellId.maxLevel - level)
		var i = 0, j = 0, orientation: Int? = nil
		let face = toFaceIJOrientation(i: &i, j: &j, orientation: &orientation)
		
		// Edges 0, 1, 2, 3 are in the S, E, N, W directions.
		return [
			S2CellId(face: face, i: i, j: j - size, sameFace: j - size >= 0).parent(level: level),
			S2CellId(face: face, i: i + size, j: j, sameFace: i + size < S2CellId.maxSize).parent(level: level),
			S2CellId(face: face, i: i, j: j + size, sameFace: j + size < S2CellId.maxSize).parent(level: level),
			S2CellId(face: face, i: i - size, j: j, sameFace: i - size >= 0).parent(level: level)
		]
	}

	/**
		Return the neighbors of closest vertex to this cell at the given level, by
		appending them to "output". Normally there are four neighbors, but the
		closest vertex may only have three neighbors if it is one of the 8 cube
		vertices.
	
		Requires: level < this.evel(), so that we can determine which vertex is
		closest (in particular, level == MAX_LEVEL is not allowed).
	*/
	public func getVertexNeighbors(level: Int) -> [S2CellId] {
		// "level" must be strictly less than this cell's level so that we can
		// determine which vertex this cell is closest to.
		// assert (level < this.level());
		var i = 0, j = 0, orientation: Int? = nil
		let face = toFaceIJOrientation(i: &i, j: &j, orientation: &orientation)
		
		// Determine the i- and j-offsets to the closest neighboring cell in each
		// direction. This involves looking at the next bit of "i" and "j" to
		// determine which quadrant of this->parent(level) this cell lies in.
		let halfsize = 1 << (S2CellId.maxLevel - (level + 1))
		let size = halfsize << 1
		var isame: Bool, jsame: Bool
		var ioffset: Int, joffset: Int
		if ((i & halfsize) != 0) {
			ioffset = size
			isame = (i + size) < S2CellId.maxSize
		} else {
			ioffset = -size
			isame = (i - size) >= 0
		}
		if ((j & halfsize) != 0) {
			joffset = size
			jsame = (j + size) < S2CellId.maxSize
		} else {
			joffset = -size
			jsame = (j - size) >= 0
		}
		
		var output: [S2CellId] = []
		
		output.append(parent(level: level))
		output.append(S2CellId(face: face, i: i + ioffset, j: j, sameFace: isame).parent(level: level))
		output.append(S2CellId(face: face, i: i, j: j + joffset, sameFace: jsame).parent(level: level))
		
		// If i- and j- edge neighbors are *both* on a different face, then this
		// vertex only has three neighbors (it is one of the 8 cube vertices).
		if isame || jsame {
			output.append(S2CellId(face: face, i: i + ioffset, j: j + joffset, sameFace: isame && jsame).parent(level: level))
		}
		
		return output
	}
	
	/**
		Append all neighbors of this cell at the given level to "output". Two cells
		X and Y are neighbors if their boundaries intersect but their interiors do
		not. In particular, two cells that intersect at a single point are neighbors.
	
		Requires: nbr_level >= this->level(). Note that for cells adjacent to a
		face vertex, the same neighbor may be appended more than once.
	*/
	public func getAllNeighbors(level nbrLevel: Int) -> [S2CellId] {
		var i = 0, j = 0, orientation: Int? = nil
		let face = toFaceIJOrientation(i: &i, j: &j, orientation: &orientation)

		// Find the coordinates of the lower left-hand leaf cell. We need to
		// normalize (i,j) to a known position within the cell because nbr_level
		// may be larger than this cell's level.
		let size = 1 << (S2CellId.maxLevel - level)
		i = i & -size
		j = j & -size

		let nbrSize = 1 << (S2CellId.maxLevel - nbrLevel)
		// assert (nbrSize <= size);
		
		var output: [S2CellId] = []

		// We compute the N-S, E-W, and diagonal neighbors in one pass.
		// The loop test is at the end of the loop to avoid 32-bit overflow.
		var k = -nbrSize
		while true {
			let sameFace: Bool
			if (k < 0) {
				sameFace = j + k >= 0
			} else if (k >= size) {
				sameFace = j + k < S2CellId.maxLevel
			} else {
				sameFace = true
				// North and South neighbors.
				output.append(S2CellId(face: face, i: i + k, j: j - nbrSize, sameFace: j - size >= 0).parent(level: nbrLevel))
				output.append(S2CellId(face: face, i: i + k, j: j + size, sameFace: j + size < S2CellId.maxLevel).parent(level: nbrLevel))
			}
			// East, West, and Diagonal neighbors.
			output.append(S2CellId(face: face, i: i - nbrSize, j: j + k, sameFace: sameFace && i - size >= 0).parent(level: nbrLevel))
			output.append(S2CellId(face: face, i: i + size, j: j + k, sameFace: sameFace && i + size < S2CellId.maxLevel).parent(level: nbrLevel))
			if k >= size { break }
			k += nbrSize
		}
		
		return output
	}
	
	// ///////////////////////////////////////////////////////////////////
	// Low-level methods.
	
	/// Return a leaf cell given its cube face (range 0..5) and i- and j-coordinates (see s2.h).
	public init(face: Int, i: Int, j: Int, sameFace: Bool = true) {
		if sameFace {
			// Optimization notes:
			// - Non-overlapping bit fields can be combined with either "+" or "|".
			// Generally "+" seems to produce better code, but not always.
			
			// gcc doesn't have very good code generation for 64-bit operations.
			// We optimize this by computing the result as two 32-bit integers
			// and combining them at the end. Declaring the result as an array
			// rather than local variables helps the compiler to do a better job
			// of register allocation as well. Note that the two 32-bits halves
			// get shifted one bit to the left when they are combined.
			var n: [Int64] = [0, Int64(face << (S2CellId.posBits - 33))]
			
			// Alternating faces have opposite Hilbert curve orientations; this
			// is necessary in order for all faces to have a right-handed
			// coordinate system.
			var bits = (face & S2CellId.swapMask)
			
			// Each iteration maps 4 bits of "i" and "j" into 8 bits of the Hilbert
			// curve position. The lookup table transforms a 10-bit key of the form
			// "iiiijjjjoo" to a 10-bit value of the form "ppppppppoo", where the
			// letters [ijpo] denote bits of "i", "j", Hilbert curve position, and
			// Hilbert curve orientation respectively.
			
			for k in (0 ..< 8).reversed() {
				bits = S2CellId.getBits(n: &n, i: i, j: j, k: k, bits: bits)
			}
			
			self.init(id: (((n[1] << 32) + n[0]) << 1) + 1)
		} else {
			// Given (i, j) coordinates that may be out of bounds, normalize them by
			// returning the corresponding neighbor cell on an adjacent face.
			
			// Convert i and j to the coordinates of a leaf cell just beyond the
			// boundary of this face. This prevents 32-bit overflow in the case
			// of finding the neighbors of a face cell, and also means that we
			// don't need to worry about the distinction between (s,t) and (u,v).
			var i = i, j = j, face = face
			i = max(-1, min(S2CellId.maxSize, i))
			j = max(-1, min(S2CellId.maxSize, j))
			
			// Find the (s,t) coordinates corresponding to (i,j). At least one
			// of these coordinates will be just outside the range [0, 1].
			let kScale = 1.0 / Double(S2CellId.maxSize)
			let s = kScale * Double((i << 1) + 1 - S2CellId.maxSize)
			let t = kScale * Double((j << 1) + 1 - S2CellId.maxSize)
			
			// Find the leaf cell coordinates on the adjacent face, and convert
			// them to a cell id at the appropriate level.
			let p = S2Projections.faceUvToXyz(face: face, u: s, v: t)
			face = S2Projections.xyzToFace(point: p)
			let st = S2Projections.validFaceXyzToUv(face: face, point: p)
			self.init(face: face, i: S2CellId.stToIJ(s: st.x), j: S2CellId.stToIJ(s: st.y))
		}
	}
	
	private static func getBits(n: inout [Int64], i: Int, j: Int, k: Int, bits: Int) -> Int {
		let mask = (1 << lookupBits) - 1
		var bits = bits
		bits += (((i >> (k * lookupBits)) & mask) << (lookupBits + 2))
		bits += (((j >> (k * lookupBits)) & mask) << 2)
		bits = lookupPos[bits]
		n[k >> 2] |= ((Int64(bits) >> 2) << Int64((k & 3) * 2 * lookupBits))
		bits &= (swapMask | invertMask)
		return bits
	}
	
	/**
	* Return the (face, i, j) coordinates for the leaf cell corresponding to this
	* cell id. Since cells are represented by the Hilbert curve position at the
	* center of the cell, the returned (i,j) for non-leaf cells will be a leaf
	* cell adjacent to the cell center. If "orientation" is non-NULL, also return
	* the Hilbert curve orientation for the current cell.
	*/
	public func toFaceIJOrientation(i: inout Int, j: inout Int, orientation: inout Int?) -> Int {
		
		// System.out.println("Entering toFaceIjorientation");
		let face = self.face
		var bits = face & S2CellId.swapMask
		
		// System.out.println("face = " + face + " bits = " + bits);
		
		// Each iteration maps 8 bits of the Hilbert curve position into
		// 4 bits of "i" and "j". The lookup table transforms a key of the
		// form "ppppppppoo" to a value of the form "iiiijjjjoo", where the
		// letters [ijpo] represents bits of "i", "j", the Hilbert curve
		// position, and the Hilbert curve orientation respectively.
		//
		// On the first iteration we need to be careful to clear out the bits
		// representing the cube face.
		for k in (0 ..< 8).reversed() {
			bits = getBits1(i: &i, j: &j, k: k, bits: bits)
			// System.out.println("pi = " + pi + " pj= " + pj + " bits = " + bits);
		}
		
		if orientation != nil {
			// The position of a non-leaf cell at level "n" consists of a prefix of
			// 2*n bits that identifies the cell, followed by a suffix of
			// 2*(MAX_LEVEL-n)+1 bits of the form 10*. If n==MAX_LEVEL, the suffix is
			// just "1" and has no effect. Otherwise, it consists of "10", followed
			// by (MAX_LEVEL-n-1) repetitions of "00", followed by "0". The "10" has
			// no effect, while each occurrence of "00" has the effect of reversing
			// the kSwapMask bit.
			// assert (S2.POS_TO_ORIENTATION[2] == 0);
			// assert (S2.POS_TO_ORIENTATION[0] == S2.SWAP_MASK);
			if (lowestOnBit & 0x1111111111111110) != 0 {
				bits ^= S2.swapMask
			}
			orientation = bits
		}
		
		return face
	}
	
	private func getBits1(i: inout Int, j: inout Int, k: Int, bits: Int) -> Int {
		let nbits = (k == 7) ? (S2CellId.maxLevel - 7 * S2CellId.lookupBits) : S2CellId.lookupBits
		
		var bits = bits
		bits += (Int((id >> Int64(k * 2 * S2CellId.lookupBits + 1))) & ((1 << (2 * nbits)) - 1)) << 2
		
		bits = S2CellId.lookupIj[bits]
		i += ((bits >> (S2CellId.lookupBits + 2)) << (k * S2CellId.lookupBits))
		
		j += ((((bits >> 2) & ((1 << S2CellId.lookupBits) - 1))) << (k * S2CellId.lookupBits))
		bits &= (S2CellId.swapMask | S2CellId.invertMask)
		
		return bits
	}
	
	/// Return the lowest-numbered bit that is on for cells at the given level.
	public var lowestOnBit: Int64 {
		return id & -id
	}
	
	public var hashValue: Int {
		
//		UInt64
		
        return Int(truncatingIfNeeded: (uid >> 32) + uid)
	}
	
	/**
		Return the lowest-numbered bit that is on for this cell id, which is equal
		to (uint64(1) << (2 * (MAX_LEVEL - level))). So for example, a.lsb() <=
		b.lsb() if and only if a.level() >= b.level(), but the first test is more efficient.
	*/
	public static func lowestOnBit(forLevel level: Int) -> Int64 {
		return 1 << Int64(2 * (maxLevel - level))
	}
	
	/// Return the i- or j-index of the leaf cell containing the given s- or t-value.
	private static func stToIJ(s: Double) -> Int {
		// Converting from floating-point to integers via static_cast is very slow
		// on Intel processors because it requires changing the rounding mode.
		// Rounding to the nearest integer using FastIntRound() is much faster.
		let m = Double(maxSize / 2) // scaling multiplier
		let x = round(m * s + (m - 0.5))
		return Int(max(0, min(2 * m - 1, x)))
	}
	
	/// Convert (face, si, ti) coordinates (see s2.h) to a direction vector (not necessarily unit length).
	private static func faceSiTiToXYZ(face: Int, si: Int, ti: Int) -> S2Point {
		let kScale = 1.0 / Double(maxSize)
		let u = S2Projections.stToUV(s: kScale * Double(si))
		let v = S2Projections.stToUV(s: kScale * Double(ti))
		return S2Projections.faceUvToXyz(face: face, u: u, v: v)
	}
	
	/// Returns true if x1 < x2, when both values are treated as unsigned.
	public static func unsignedLongLessThan(lhs: Int64, rhs: Int64) -> Bool {
//		return (lhs + .min) < (rhs + .min)
		return UInt64(bitPattern: lhs) < UInt64(bitPattern: rhs)
	}
	
	/// Returns true if x1 > x2, when both values are treated as unsigned.
	public static func unsignedLongGreaterThan(lhs: Int64, rhs: Int64) -> Bool {
//		return (lhs + .min) > (rhs + .min)
		return UInt64(bitPattern: lhs) > UInt64(bitPattern: rhs)
	}
	
	private static func initLookupCell(level: Int, i: Int, j: Int, origOrientation: Int, pos: Int, orientation: Int) {
		if level == lookupBits {
			let ij = (i << lookupBits) + j
			lookup.pos[(ij << 2) + origOrientation] = (pos << 2) + orientation
			lookup.ij[(pos << 2) + origOrientation] = (ij << 2) + orientation
		} else {
			var level = level, i = i, j = j, pos = pos
			level += 1
			i <<= 1
			j <<= 1
			pos <<= 2
			// Initialize each sub-cell recursively.
			for subPos in 0 ..< 4 {
				let ij = S2.posToIJ(orientation: orientation, position: subPos)
				let orientationMask = S2.posToOrientation(position: subPos)
				initLookupCell(level: level, i: i + (ij >> 1), j: j + (ij & 1), origOrientation: origOrientation, pos: pos + subPos, orientation: orientation ^ orientationMask)
			}
		}
	}

}

public func ==(lhs: S2CellId, rhs: S2CellId) -> Bool {
	return lhs.id == rhs.id
}

public func <(lhs: S2CellId, rhs: S2CellId) -> Bool {
	return S2CellId.unsignedLongLessThan(lhs: lhs.id, rhs: rhs.id)
}
