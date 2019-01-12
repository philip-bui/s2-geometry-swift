//
//  S2RegionCoverer.swift
//  S2Geometry
//
//  Created by Alex Studnicka on 7/1/16.
//  Copyright © 2016 Alex Studnicka. MIT License.
//

/**
	An S2RegionCoverer is a class that allows arbitrary regions to be
	approximated as unions of cells (S2CellUnion). This is useful for
	implementing various sorts of search and precomputation operations.

	Typical usage: {@code S2RegionCoverer coverer; coverer.setMaxCells(5); S2Cap
	cap = S2Cap.fromAxisAngle(...); S2CellUnion covering;
	coverer.getCovering(cap, covering); * }

	This yields a cell union of at most 5 cells that is guaranteed to cover the
	given cap (a disc-shaped region on the sphere).

	The approximation algorithm is not optimal but does a pretty good job in
	practice. The output does not always use the maximum number of cells allowed,
	both because this would not always yield a better approximation, and because
	max_cells() is a limit on how much work is done exploring the possible
	covering as well as a limit on the final output size.

	One can also generate interior coverings, which are sets of cells which are
	entirely contained within a region. Interior coverings can be empty, even for
	non-empty regions, if there are no cells that satisfy the provided
	constraints and are contained by the region. Note that for performance
	reasons, it is wise to specify a max_level when computing interior coverings -
	otherwise for regions with small or zero area, the algorithm may spend a
	lot of time subdividing cells all the way to leaf level to try to find
	contained cells.

	This class is thread-unsafe. Simultaneous calls to any of the getCovering
	methods will conflict and produce unpredictable results.
*/
public class S2RegionCoverer {
	
	/**
		By default, the covering uses at most 8 cells at any level. This gives a
		reasonable tradeoff between the number of cells used and the accuracy of
		the approximation (see table below).
	*/
	public static let defaultMaxCells = 8
	
	private static let faceCells: [S2Cell] = [
		S2Cell(face: 0, pos: 0, level: 0),
		S2Cell(face: 1, pos: 0, level: 0),
		S2Cell(face: 2, pos: 0, level: 0),
		S2Cell(face: 3, pos: 0, level: 0),
		S2Cell(face: 4, pos: 0, level: 0),
		S2Cell(face: 5, pos: 0, level: 0)
	]
	
	private var _minLevel: Int = 0
	public var minLevel: Int {
		get {
			return _minLevel
		}
		set {
			_minLevel = max(0, min(S2CellId.maxLevel, newValue))
		}
	}
	
	private var _maxLevel: Int = S2CellId.maxLevel
	public var maxLevel: Int {
		get {
			return _maxLevel
		}
		set {
			_maxLevel = max(0, min(S2CellId.maxLevel, newValue))
		}
	}
	
	private var _levelMod: Int = 1
	public var levelMod: Int {
		get {
			return _levelMod
		}
		set {
			_levelMod = max(1, min(3, newValue))
		}
	}
	
	private var _maxCells: Int = S2RegionCoverer.defaultMaxCells
	/**
		Sets the maximum desired number of cells in the approximation (defaults to
		kDefaultMaxCells). Note the following:
	
		- For any setting of max_cells(), up to 6 cells may be returned if that
		is the minimum number of cells required (e.g. if the region intersects all
		six face cells). Up to 3 cells may be returned even for very tiny convex
		regions if they happen to be located at the intersection of three cube
		faces.
	
		- For any setting of max_cells(), an arbitrary number of cells may be
		returned if min_level() is too high for the region being approximated.
	
		- If max_cells() is less than 4, the area of the covering may be
		arbitrarily large compared to the area of the original region even if the
		region is convex (e.g. an S2Cap or S2LatLngRect).
	
		Accuracy is measured by dividing the area of the covering by the area of
		the original region. The following table shows the median and worst case
		values for this area ratio on a test case consisting of 100,000 spherical
		caps of random size (generated using s2regioncoverer_unittest):
	
		```
		max_cells: 3 4 5 6 8 12 20 100 1000
		median ratio: 5.33 3.32 2.73 2.34 1.98 1.66 1.42 1.11 1.01
		worst case: 215518 14.41 9.72 5.26 3.91 2.75 1.92 1.20 1.02
		```
	*/
	public var maxCells: Int {
		get {
			return _maxCells
		}
		set {
			_maxCells = newValue
		}
	}
	
	// True if we're computing an interior covering.
	private var interiorCovering: Bool = false
	
	/**
		We save a temporary copy of the pointer passed to GetCovering() in order to
		avoid passing this parameter around internally. It is only used (and only
		valid) for the duration of a single GetCovering() call.
	*/
	private var region: S2Region? = nil
	
	/**
		A temporary variable used by GetCovering() that holds the cell ids that
		have been added to the covering so far.
	*/
	private var result: [S2CellId] = []
	
	struct Candidate {
		let cell: S2Cell
		var isTerminal: Bool // Cell should not be expanded further.
		var children: [Candidate] // Actual size may be 0, 4, 16, or 64 elements.
		
		init(cell: S2Cell, isTerminal: Bool, children: [Candidate] = []) {
			self.cell = cell
			self.isTerminal = isTerminal
			self.children = children
		}
	}
	
	struct QueueEntry: Comparable {
		let id: Int
		let candidate: Candidate
		
		init(id: Int, candidate: Candidate) {
			self.id = id
			self.candidate = candidate
		}
	}
	
	/**
		We keep the candidates in a priority queue. We specify a vector to hold the
		queue entries since for some reason priority_queue<> uses a deque by default.
	*/
	private var candidateQueue: [QueueEntry] = []
	
	/// Default constructor, sets all fields to default values.
	public init() { }
	
	/**
		Computes a list of cell ids that covers the given region and satisfies the
		various restrictions specified above.
	
		- Parameter region: The region to cover
	
		- Returns: The list filled in by this method
	*/
	public func getCovering(region: S2Region) -> [S2CellId] {
		// Rather than just returning the raw list of cell ids generated by
		// GetCoveringInternal(), we construct a cell union and then denormalize it.
		// This has the effect of replacing four child cells with their parent
		// whenever this does not violate the covering parameters specified
		// (min_level, level_mod, etc). This strategy significantly reduces the
		// number of cells returned in many cases, and it is cheap compared to
		// computing the covering in the first place.
		
		let tmp = getCoveringUnion(region: region)
		return tmp.denormalize(minLevel: minLevel, levelMod: levelMod)
	}
	
	/**
		Computes a list of cell ids that is contained within the given region and
		satisfies the various restrictions specified above.
	
		- Parameter region: The region to fill
		
		- Returns: The list filled in by this method
	*/
	public func getInteriorCovering(region: S2Region) -> [S2CellId] {
		let tmp =  getInteriorCoveringUnion(region: region)
		return tmp.denormalize(minLevel: minLevel, levelMod: levelMod)
	}
	
	/**
		Return a normalized cell union that covers the given region and satisfies
		the restrictions *EXCEPT* for min_level() and level_mod(). These criteria
		cannot be satisfied using a cell union because cell unions are
		automatically normalized by replacing four child cells with their parent
		whenever possible. (Note that the list of cell ids passed to the cell union
		constructor does in fact satisfy all the given restrictions.)
	*/
	public func getCoveringUnion(region: S2Region) -> S2CellUnion {
		interiorCovering = false
		getCoveringInternal(region: region)
		return S2CellUnion(cellIds: result)
	}

	/**
		Return a normalized cell union that is contained within the given region
		and satisfies the restrictions *EXCEPT* for min_level and level_mod.
	*/
	public func getInteriorCoveringUnion(region: S2Region) -> S2CellUnion {
		interiorCovering = true
		getCoveringInternal(region: region)
		return S2CellUnion(cellIds: result)
	}
	
	/**
		If the cell intersects the given region, return a new candidate with no
		children, otherwise return null. Also marks the candidate as "terminal" if
		it should not be expanded further.
	*/
	private func newCandidate(cell: S2Cell) -> Candidate? {
		guard let region = region, region.mayIntersect(cell: cell) else { return nil }
		
		let cellLevel = Int(cell.level)
		var isTerminal = false
		if cellLevel >= minLevel {
			if interiorCovering {
				if region.contains(cell: cell) {
					isTerminal = true
				} else if cellLevel + levelMod > maxLevel {
					return nil
				}
			} else {
				if cellLevel + levelMod > maxLevel || region.contains(cell: cell) {
					isTerminal = true
				}
			}
		}
		return Candidate(cell: cell, isTerminal: isTerminal)
	}
	
	/// Return the log base 2 of the maximum number of children of a candidate.
	private var maxChildrenShift: Int {
		return 2 * levelMod
	}
	
	/**
		Process a candidate by either adding it to the result list or expanding its
		children and inserting it into the priority queue. Passing an argument of
		NULL does nothing.
	*/
	private func addCandidate(candidate: Candidate?) {
		guard var candidate = candidate else { return }
		
		print("addCandidate:", candidate)
		
		if candidate.isTerminal {
//			print("isTerminal")
			result.append(candidate.cell.cellId)
			return
		}
		// assert (candidate.numChildren == 0);
		
		// Expand one level at a time until we hit min_level_ to ensure that
		// we don't skip over it.
		var numLevels = (Int(candidate.cell.level) < minLevel) ? 1 : levelMod
		let numTerminals = expandChildren(candidate: &candidate, cell: candidate.cell, numLevels: &numLevels)
		
		if candidate.children.isEmpty {
			// Do nothing
		} else if !interiorCovering && numTerminals == 1 << maxChildrenShift && Int(candidate.cell.level) >= minLevel {
			// Optimization: add the parent cell rather than all of its children.
			// We can't do this for interior coverings, since the children just
			// intersect the region, but may not be contained by it - we need to
			// subdivide them further.
			candidate.isTerminal = true
			addCandidate(candidate: candidate)
		} else {
			// We negate the priority so that smaller absolute priorities are returned
			// first. The heuristic is designed to refine the largest cells first,
			// since those are where we have the largest potential gain. Among cells
			// at the same level, we prefer the cells with the smallest number of
			// intersecting children. Finally, we prefer cells that have the smallest
			// number of children that cannot be refined any further.
			let priority = -((((Int(candidate.cell.level) << maxChildrenShift) + candidate.children.count) << maxChildrenShift) + numTerminals)
			candidateQueue.append(QueueEntry(id: priority, candidate: candidate))
			// logger.info("Push: " + candidate.cell.id() + " (" + priority + ") ");
		}
	}
	
	/**
		Populate the children of "candidate" by expanding the given number of levels from the given cell.
		Returns the number of children that were marked "terminal".
	*/
	private func expandChildren(candidate: inout Candidate, cell: S2Cell, numLevels: inout Int) -> Int {
		guard let region = region else { return 0 }
		
		print("expandChildren:", candidate)
		
		numLevels -= 1
		let childCells = cell.subdivide()
		var numTerminals = 0
		for childCell in childCells {
			if numLevels > 0 {
				if region.mayIntersect(cell: childCell) {
					numTerminals += expandChildren(candidate: &candidate, cell: childCell, numLevels: &numLevels)
				}
				continue
			}
			if let child = newCandidate(cell: childCell) {
				candidate.children.append(child)
				if child.isTerminal {
					numTerminals += 1
				}
			}
		}
		return numTerminals
	}
	
	/// Computes a set of initial candidates that cover the given region.
	private func getInitialCandidates() {
		guard let region = region else { return }
		
		// Optimization: if at least 4 cells are desired (the normal case),
		// start with a 4-cell covering of the region's bounding cap. This
		// lets us skip quite a few levels of refinement when the region to
		// be covered is relatively small.
		if false && maxCells >= 4 {
			// Find the maximum level such that the bounding cap contains at most one cell vertex at that level.
			let cap = region.capBound
			var level = min(S2Projections.minWidth.getMaxLevel(value: 2 * cap.angle.radians), min(maxLevel, S2CellId.maxLevel - 1))
			if levelMod > 1 && level > minLevel {
				level -= (level - minLevel) % levelMod
			}
			// We don't bother trying to optimize the level == 0 case, since more than four face cells may be required.
			if level > 0 {
				// Find the leaf cell containing the cap axis, and determine which subcell of the parent cell contains it.
				let base = S2CellId(point: cap.axis).getVertexNeighbors(level: level)
				for id in base {
					addCandidate(candidate: newCandidate(cell: S2Cell(cellId: id)))
				}
				return
			}
		}
		
		// Default: start with all six cube faces.
		for face in 0 ..< 6 {
			addCandidate(candidate: newCandidate(cell: S2RegionCoverer.faceCells[face]))
		}
	}
	
	/// Generates a covering and stores it in result.
	private func getCoveringInternal(region: S2Region) {
		// Strategy: Start with the 6 faces of the cube. Discard any
		// that do not intersect the shape. Then repeatedly choose the
		// largest cell that intersects the shape and subdivide it.
		//
		// result contains the cells that will be part of the output, while the
		// priority queue contains cells that we may still subdivide further. Cells
		// that are entirely contained within the region are immediately added to
		// the output, while cells that do not intersect the region are immediately
		// discarded.
		// Therefore pq_ only contains cells that partially intersect the region.
		// Candidates are prioritized first according to cell size (larger cells
		// first), then by the number of intersecting children they have (fewest
		// children first), and then by the number of fully contained children
		// (fewest children first).
		
		precondition(candidateQueue.isEmpty && result.isEmpty)
		
		self.region = region
		
		getInitialCandidates()
		candidateQueue.sort()
		
		print("---------------------------------------------")
		print("candidateQueue:", candidateQueue)
		print("---------------------------------------------")
		
		while !candidateQueue.isEmpty && (!interiorCovering || result.count < maxCells) {
			var candidate = candidateQueue.removeFirst().candidate
			if (Int(candidate.cell.level) < minLevel || candidate.children.count == 1
				|| result.count + (interiorCovering ? 0 : candidateQueue.count) + candidate.children.count <= maxCells) {
				// Expand this candidate into its children.
				for child in candidate.children {
					addCandidate(candidate: child)
				}
			} else if interiorCovering {
				// Do nothing
			} else {
				candidate.isTerminal = true
				addCandidate(candidate: candidate)
			}
		}
		
		print("result:", result)
		
		candidateQueue.removeAll()
		self.region = nil
	}
	
}

func ==(lhs: S2RegionCoverer.QueueEntry, rhs: S2RegionCoverer.QueueEntry) -> Bool {
	return lhs.id == rhs.id
}

func <(lhs: S2RegionCoverer.QueueEntry, rhs: S2RegionCoverer.QueueEntry) -> Bool {
	return lhs.id > rhs.id
}
