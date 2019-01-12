//
//  S2EdgeIndex.swift
//  S2Geometry
//
//  Created by Alex Studnicka on 7/1/16.
//  Copyright © 2016 Alex Studnicka. MIT License.
//

public class S2EdgeIndex {
	
	/// Thicken the edge in all directions by roughly 1% of the edge length when thickenEdge is true.
	private static let thickening = 0.01
	
	/// Threshold for small angles, that help lenientCrossing to determine whether two edges are likely to intersect.
	private static let maxDetError = 1e-14
	
	/// The cell containing each edge, as given in the parallel array `edges`.
	private var cells: [Int64] = []
	
	/// The edge contained by each cell, as given in the parallel array `cells`.
	private var edges: [Int] = []
	
	/**
		No cell strictly below this level appears in mapping.
		Initially leaf level, that's the minimum level at which we will ever look for test edges.
	*/
	private var minimumS2LevelUsed: Int = S2CellId.maxLevel
	
	/// Has the index been computed already?
	public private(set) var isIndexComputed: Bool = false
	
	/// Number of queries so far
	private var queryCount: Int = 0
	
	public init() { }
	
	/**
		Overwrite these functions to give access to the underlying data. The
		function getNumEdges() returns the number of edges in the index, while
		edgeFrom(index) and edgeTo(index) return the "from" and "to" endpoints of
		the edge at the given index.
	*/
	public var numEdges: Int {
		fatalError()
	}
	
	public func edgeFrom(index: Int) -> S2Point {
		fatalError()
	}
	
	public func edgeTo(index: Int) -> S2Point {
		fatalError()
	}
	
	/// Empties the index in case it already contained something.
	public func reset() {
		minimumS2LevelUsed = S2CellId.maxLevel
		isIndexComputed = false
		queryCount = 0
		cells = []
		edges = []
	}
	
	/// Computes the index (if it has not been previously done).
	public func computeIndex() {
		if isIndexComputed { return }
		var cellList: [Int64] = []
		var edgeList: [Int] = []
		for i in 0 ..< numEdges {
			let from = edgeFrom(index: i)
			let to = edgeTo(index: i)
			let (level, cover) = getCovering(a: from, b: to, thickenEdge: true)
			minimumS2LevelUsed = min(minimumS2LevelUsed, level);
			for cellId in cover {
				cellList.append(cellId.id)
				edgeList.append(i)
			}
		}
		cells = cellList
		edges = edgeList
		sortIndex()
		isIndexComputed = true
	}
	
	/// Sorts the parallel `cells` and `edges` arrays.
	private func sortIndex() {
		// create an array of indices and sort based on the values in the parallel arrays at each index
		var indices: [Int] = []
		for i in 0 ..< cells.count {
			indices.append(i)
		}
		
		indices.sort(by: { (cells[$0], edges[$0]) < (cells[$1], edges[$1]) })
		
		// copy the cells and edges in the order given by the sorted list of indices
		var newCells: [Int64] = []
		var newEdges: [Int] = []
		for i in 0 ..< indices.count {
			newCells.append(cells[indices[i]])
			newEdges.append(edges[indices[i]])
		}
		
		// replace the cells and edges with the sorted arrays
		cells = newCells
		edges = newEdges
	}
	
	/**
		If the index hasn't been computed yet, looks at how much work has gone into
		iterating using the brute force method, and how much more work is planned
		as defined by 'cost'. If it were to have been cheaper to use a quad tree
		from the beginning, then compute it now. This guarantees that we will never
		use more than twice the time we would have used had we known in advance
		exactly how many edges we would have wanted to test. It is the theoretical
		best.
	
		The value 'n' is the number of iterators we expect to request from this
		edge index.
	
		If we have m data edges and n query edges, then the brute force cost is m
		* n * testCost where testCost is taken to be the cost of
		EdgeCrosser.robustCrossing, measured to be about 30ns at the time of this
		writing.
	
		If we compute the index, the cost becomes: m * costInsert + n *
		costFind(m)
	
		- costInsert can be expected to be reasonably stable, and was measured at
		1200ns with the BM_QuadEdgeInsertionCost benchmark.
	
		- costFind depends on the length of the edge . For m=1000 edges, we got
		timings ranging from 1ms (edge the length of the polygon) to 40ms. The
		latter is for very long query edges, and needs to be optimized. We will
		assume for the rest of the discussion that costFind is roughly 3ms.
	
		When doing one additional query, the differential cost is m * testCost -
		costFind(m) With the numbers above, it is better to use the quad tree (if
		we have it) if m >= 100.
	
		If m = 100, 30 queries will give m*n*testCost = m * costInsert = 100ms,
		while the marginal cost to find is 3ms. Thus, this is a reasonable thing to
		do.
	*/
	public func predictAdditionalCalls(n: Int) {
		if isIndexComputed { return }
		if numEdges > 100 && queryCount + n > 30 {
			computeIndex()
		}
	}
	
	/**
		Appends to "candidateCrossings" all edge references which may cross the
		given edge. This is done by covering the edge and then finding all
		references of edges whose coverings overlap this covering. Parent cells are
		checked level by level. Child cells are checked all at once by taking
		advantage of the natural ordering of S2CellIds.
	*/
	public func findCandidateCrossings(a: S2Point, b: S2Point) -> [Int] {
		precondition(isIndexComputed)
		let cover = getCovering(a: a, b: b, thickenEdge: false).edgeCovering
		
		// Edge references are inserted into the map once for each covering cell, so absorb duplicates here
		var uniqueSet = getEdgesInParentCells(cover: cover)
		
		// TODO(user): An important optimization for long query
		// edges (Contains queries): keep a bounding cap and clip the query
		// edge to the cap before starting the descent.
		getEdgesInChildrenCells(a: a, b: b, cover: cover, candidateCrossings: &uniqueSet)
		
		return Array(uniqueSet)
	}
	
	/**
		Returns the smallest cell containing all four points, or {@link S2CellId#sentinel()} if they are not all on the same face.
		The points don't need to be normalized.
	*/
	private static func containingCell(pa: S2Point, pb: S2Point, pc: S2Point, pd: S2Point) -> S2CellId {
		var a = S2CellId(point: pa)
		var b = S2CellId(point: pb)
		var c = S2CellId(point: pc)
		var d = S2CellId(point: pd)
		
		if a.face != b.face || a.face != c.face || a.face != d.face {
			return S2CellId.sentinel
		}
		
		while a != b || a != c || a != d {
			a = a.parent
			b = b.parent
			c = c.parent
			d = d.parent
		}
		
		return a
	}
	
	/**
		Returns the smallest cell containing both points, or Sentinel if they are
		not all on the same face. The points don't need to be normalized.
	*/
	private static func containingCell(pa: S2Point, pb: S2Point) -> S2CellId {
		var a = S2CellId(point: pa)
		var b = S2CellId(point: pb)
		
		if a.face != b.face {
			return S2CellId.sentinel
		}
		
		while a != b {
			a = a.parent
			b = b.parent
		}
		
		return a
	}
	
	/**
		Computes a cell covering of an edge. Clears edgeCovering and returns the
		level of the s2 cells used in the covering (only one level is ever used for each call).
	
		If thickenEdge is true, the edge is thickened and extended by 1% of its length.
	
		It is guaranteed that no child of a covering cell will fully contain the covered edge.
	*/
	private func getCovering(a: S2Point, b: S2Point, thickenEdge: Bool) -> (level: Int, edgeCovering: [S2CellId]) {
		var edgeCovering: [S2CellId] = []
		
		// Selects the ideal s2 level at which to cover the edge, this will be the
		// level whose S2 cells have a width roughly commensurate to the length of
		// the edge. We multiply the edge length by 2*THICKENING to guarantee the
		// thickening is honored (it's not a big deal if we honor it when we don't
		// request it) when doing the covering-by-cap trick.
		let edgeLength = a.angle(to: b)
		let idealLevel = S2Projections.minWidth.getMaxLevel(value: edgeLength * (1 + 2 * S2EdgeIndex.thickening))
		
		var containingCellId: S2CellId
		if !thickenEdge {
			containingCellId = S2EdgeIndex.containingCell(pa: a, pb: b)
		} else {
			if idealLevel == S2CellId.maxLevel {
				// If the edge is tiny, instabilities are more likely, so we
				// want to limit the number of operations.
				// We pretend we are in a cell much larger so as to trigger the
				// 'needs covering' case, so we won't try to thicken the edge.
				containingCellId = S2CellId(id: 0xFFF0).parent(level: 3)
			} else {
				let pq = (b - a) * S2EdgeIndex.thickening
				let x = edgeLength * S2EdgeIndex.thickening
				let ortho = (S2Point.normalize(point: pq.crossProd(a)) * x)
				let p = a - pq
				let q = b + pq
				// If p and q were antipodal, the edge wouldn't be lengthened,
				// and it could even flip! This is not a problem because
				// idealLevel != 0 here. The farther p and q can be is roughly
				// a quarter Earth away from each other, so we remain
				// Theta(THICKENING).
				containingCellId = S2EdgeIndex.containingCell(pa: (p - ortho), pb: (p + ortho), pc: (q - ortho), pd: (q + ortho))
			}
		}
		
		// Best case: edge is fully contained in a cell that's not too big.
		if containingCellId != S2CellId.sentinel && containingCellId.level >= idealLevel - 2 {
			edgeCovering.append(containingCellId)
			return (containingCellId.level, edgeCovering)
		}
		
		if idealLevel == 0 {
			// Edge is very long, maybe even longer than a face width, so the
			// trick below doesn't work. For now, we will add the whole S2 sphere.
			// TODO(user): Do something a tad smarter (and beware of the
			// antipodal case).
			var cellId = S2CellId.begin(level: 0)
			while cellId != S2CellId.end(level: 0) {
				edgeCovering.append(cellId)
				cellId = cellId.next()
			}
			return (0, edgeCovering)
		}
		
		// TODO(user): Check trick below works even when vertex is at
		// interface
		// between three faces.
		
		// Use trick as in S2PolygonBuilder.PointIndex.findNearbyPoint:
		// Cover the edge by a cap centered at the edge midpoint, then cover
		// the cap by four big-enough cells around the cell vertex closest to the
		// cap center.
		let middle = S2Point.normalize(point: ((a + b) / 2))
		let actualLevel = min(idealLevel, S2CellId.maxLevel - 1);
		edgeCovering = S2CellId(point: middle).getVertexNeighbors(level: actualLevel)
		return (actualLevel, edgeCovering)
	}
	
	/**
		Filters a list of entries down to the inclusive range defined by the given
		cells, in `O(log N)` time.
	
		- Parameter cell1: One side of the inclusive query range.
		- Parameter cell2: The other side of the inclusive query range.
		- Returns: An array of length 2, containing the start/end indices.
	*/
	private func getEdges(cell1: Int64, cell2: Int64) -> (Int, Int) {
		// ensure cell1 <= cell2
		if cell1 > cell2 {
			return getEdges(cell1: cell2, cell2: cell1)
		}
		// The binary search returns -N-1 to indicate an insertion point at index N,
		// if an exact match cannot be found. Since the edge indices queried for are
		// not valid edge indices, we will always get -N-1, so we immediately
		// convert to N.
		return (-1 - binarySearch(cell: cell1, edge: .max), -1 - binarySearch(cell: cell2, edge: .max))
	}
	
	private func binarySearch(cell: Int64, edge: Int) -> Int {
		var low = 0;
		var high = cells.count - 1
		while low <= high {
			let mid = (low + high) >> 1
			let cmp1 = (cells[mid], edges[mid])
			let cmp2 = (cell, edge)
			if cmp1 < cmp2 {
				low = mid + 1
			} else if cmp1 > cmp2 {
				high = mid - 1
			} else {
				return mid
			}
		}
		return -(low + 1)
	}
	
	/**
		Returns candidateCrossings with all the edges present in any ancestor of any
		cell of cover, down to minimumS2LevelUsed. The cell->edge map is in the variable mapping.
	*/
	private func getEdgesInParentCells(cover: [S2CellId]) -> Set<Int> {
		var candidateCrossings: Set<Int> = []
		
		// Find all parent cells of covering cells.
		var parentCells: Set<S2CellId> = []
		for coverCell in cover {
			var parentLevel = coverCell.level - 1
			while parentLevel >= minimumS2LevelUsed {
				if !parentCells.insert(coverCell.parent(level: parentLevel)).inserted { break } // cell is already in => parents are too.
				parentLevel -= 1
			}
		}
		
		// Put parent cell edge references into result.
		for parentCell in parentCells {
			let bounds = getEdges(cell1: parentCell.id, cell2: parentCell.id)
			for i in bounds.0 ..< bounds.1 {
				candidateCrossings.insert(edges[i])
			}
		}
		
		return candidateCrossings
	}
	
	/// Returns true if ab possibly crosses cd, by clipping tiny angles to zero.
	private static func lenientCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Bool {
		// assert (S2.isUnitLength(a));
		// assert (S2.isUnitLength(b));
		// assert (S2.isUnitLength(c));
		
		let acb = a.crossProd(c).dotProd(b)
		let bda = b.crossProd(d).dotProd(a)
		if (abs(acb) < maxDetError || abs(bda) < maxDetError) {
			return true
		}
		if (acb * bda < 0) {
			return false
		}
		let cbd = c.crossProd(b).dotProd(d)
		let dac = c.crossProd(a).dotProd(c)
		if (abs(cbd) < maxDetError || abs(dac) < maxDetError) {
			return true
		}
		return (acb * cbd >= 0) && (acb * dac >= 0)
	}
	
	/// Returns true if the edge and the cell (including boundary) intersect.
	private static func edgeIntersectsCellBoundary(a: S2Point, b: S2Point, cell: S2Cell) -> Bool {
		var vertices: [S2Point] = []
		for i in 0 ..< 4 {
			vertices.append(cell.getVertex(i))
		}
		for i in 0 ..< 4 {
			let fromPoint = vertices[i];
			let toPoint = vertices[(i + 1) % 4]
			if lenientCrossing(a: a, b: b, c: fromPoint, d: toPoint) {
				return true
			}
		}
		return false
	}
	
	/**
		Appends to candidateCrossings the edges that are fully contained in an S2
		covering of edge. The covering of edge used is initially cover, but is
		refined to eliminate quickly subcells that contain many edges but do not
		intersect with edge.
	*/
	private func getEdgesInChildrenCells(a: S2Point, b: S2Point, cover: [S2CellId], candidateCrossings: inout Set<Int>) {
		var cover = cover
		// Put all edge references of (covering cells + descendant cells) into result.
		// This relies on the natural ordering of S2CellIds.
		while !cover.isEmpty {
			let cell = cover.removeLast()
			var bounds = getEdges(cell1: cell.rangeMin.id, cell2: cell.rangeMax.id)
			if bounds.1 - bounds.0 <= 16 {
				for i in bounds.0 ..< bounds.1 {
					candidateCrossings.insert(edges[i])
				}
			} else {
				// Add cells at this level
				bounds = getEdges(cell1: cell.id, cell2: cell.id)
				for i in bounds.0 ..< bounds.1 {
					candidateCrossings.insert(edges[i])
				}
				// Recurse on the children -- hopefully some will be empty.
				let children = S2Cell(cellId: cell).subdivide()
				for child in children {
					// TODO(user): Do the check for the four cells at once,
					// as it is enough to check the four edges between the cells. At
					// this time, we are checking 16 edges, 4 times too many.
					//
					// Note that given the guarantee of AppendCovering, it is enough
					// to check that the edge intersect with the cell boundary as it
					// cannot be fully contained in a cell.
					if S2EdgeIndex.edgeIntersectsCellBoundary(a: a, b: b, cell: child) {
						cover.append(child.cellId)
					}
				}
			}
		}
	}
	
	/*
		An iterator on data edges that may cross a query edge (a,b). Create the
		iterator, call getCandidates(), then hasNext()/next() repeatedly.
	
		The current edge in the iteration has index index(), goes between from() and to().
	*/
	public struct DataEdgeIterator: Sequence, IteratorProtocol {
		
		/// The structure containing the data edges.
		private let edgeIndex: S2EdgeIndex
		
		/// Tells whether getCandidates() obtained the candidates through brute force iteration or using the quad tree structure.
		private var isBruteForce: Bool = false
		
		/// Index of the current edge and of the edge before the last next() call.
		private var currentIndex: Int = 0
		
		/// Cache of edgeIndex.getNumEdges() so that hasNext() doesn't make an extra call
		private var numEdges: Int = 0
		
		/// All the candidates obtained by getCandidates() when we are using a quad-tree (i.e. isBruteForce = false).
		private var candidates: [Int] = []
		
		/// Index within array above. We have: currentIndex = candidates.get(currentIndexInCandidates).
		private var currentIndexInCandidates: Int = 0
		
		public init(edgeIndex: S2EdgeIndex) {
			self.edgeIndex = edgeIndex
		}
		
		/// Initializes the iterator to iterate over a set of candidates that may cross the edge (a,b).
		public mutating func getCandidates(a: S2Point, b: S2Point) {
			edgeIndex.predictAdditionalCalls(n: 1)
			isBruteForce = !edgeIndex.isIndexComputed
			if isBruteForce {
				edgeIndex.queryCount += 1
				currentIndex = 0
				numEdges = edgeIndex.numEdges
			} else {
				candidates = edgeIndex.findCandidateCrossings(a: a, b: b)
				currentIndexInCandidates = 0
				if !candidates.isEmpty {
					currentIndex = candidates[0]
				}
			}
		}
		
		/// False if there are no more candidates; true otherwise.
		private var hasNext: Bool {
			if isBruteForce {
				return currentIndex < numEdges
			} else {
				return currentIndexInCandidates < candidates.count
			}
		}
		
		/// Iterate to the next available candidate.
		public mutating func next() -> Int? {
			guard hasNext else { return nil }
			
			defer {
				if isBruteForce {
					currentIndex += 1
				} else {
					currentIndexInCandidates += 1
					if currentIndexInCandidates < candidates.count {
						currentIndex = candidates[currentIndexInCandidates]
					}
				}
			}
			
			return currentIndex
		}
		
	}
	
}
