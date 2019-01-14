//
//  GeometryTestCase.swift
//  S2GeometrySwift
//
//  Created by Philip on 13/1/19.
//

import Darwin
@testable import S2GeometrySwift
import XCTest

class GeometryTestCase: XCTestCase {

    /**
     * Return a random cell id at the given level or at a randomly chosen level.
     * The distribution is uniform over the space of cell ids, but only
     * approximately uniform over the surface of the sphere.
     */
    public func getRandomCellId(level: Int = Int.random(in: 0..<S2CellId.maxLevel)) -> S2CellId {
        let face = Int.random(in: 0..<S2CellId.numFaces)
        let pos = Int64.random(in: 0...Int64.max) & ((Int64(1) << (2 * S2CellId.maxLevel)) - Int64(1))
        return S2CellId(face: face, pos: pos, level: level)
    }
}
