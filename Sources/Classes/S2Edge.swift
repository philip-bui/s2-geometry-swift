//
//  S2Edge.swift
//  S2Geometry
//
//  Created by Alex Studnicka on 7/1/16.
//  Copyright © 2016 Alex Studnicka. MIT License.
//

/**
	An abstract directed edge from one S2Point to another S2Point.
*/
public struct S2Edge: Equatable {
	
	public let start: S2Point
	public let end: S2Point
	
	public init(start: S2Point, end: S2Point) {
		self.start = start
		self.end = end
	}
	
}

public func ==(lhs: S2Edge, rhs: S2Edge) -> Bool {
	return lhs.start == rhs.start && lhs.end == rhs.end
}
