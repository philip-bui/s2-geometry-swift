//
//  S2PointTests.swift
//  S2GeometrySwift iOS Tests
//
//  Created by Philip on 11/1/19.
//

@testable import S2GeometrySwift
import XCTest

class S2PointTests: XCTestCase {
    
    func testEquality() {
        XCTAssertEqual(S2Point(x: 1, y: 1, z: 1), S2Point(x: 1, y: 1, z: 1))
        XCTAssertNotEqual(S2Point(x: 1, y: 1, z: 1), S2Point(x: 1, y: 1, z: 2))
    }
}
