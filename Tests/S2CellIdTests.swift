//
//  S2CellIdTests.swift
//  S2GeometrySwift
//
//  Created by Philip on 13/1/19.
//

import XCTest
@testable import S2GeometrySwift

class S2CellIdTests: GeometryTestCase {

    private func getCellId(_ latDegrees: Double, _ lngDegrees: Double) -> S2CellId{
        let id = S2CellId(latlng: S2LatLng.fromDegrees(lat: latDegrees, lng: lngDegrees))
        return id
    }

    func testBasic() {
        NSLog("TestBasic")
        
        var id = S2CellId()
        XCTAssertEqual(id.id, 0)
        XCTAssertTrue(!id.isValid)
        
        id = S2CellId.init(face: 3, pos: 0x12345678, level: S2CellId.maxLevel - 4)
        XCTAssertTrue(id.isValid)
        XCTAssertEqual(id.face, 3)
        // XCTAssertEqual(id.pos, 0x12345700)
        XCTAssertEqual(id.level, S2CellId.maxLevel - 4)
        XCTAssertTrue(!id.isLeaf)
        
        // Check face definitions
        XCTAssertEqual(getCellId(0, 0).face, 0)
        XCTAssertEqual(getCellId(0, 90).face, 1)
        XCTAssertEqual(getCellId(90, 0).face, 2)
        XCTAssertEqual(getCellId(0, 180).face, 3)
        XCTAssertEqual(getCellId(0, -90).face, 4)
        XCTAssertEqual(getCellId(-90, 0).face, 5)
        
        // Check parent/child relationships.
        XCTAssertEqual(id.childBegin(level: id.level + 2).pos, 0x12345610)
        XCTAssertEqual(id.childBegin().pos, 0x12345640)
        XCTAssertEqual(id.parent.pos, 0x12345400)
        XCTAssertEqual(id.parent(level: id.level - 2).pos, 0x12345000)
        
        // Check ordering of children relative to parents.
        XCTAssertTrue(id.childBegin() < id)
        XCTAssertTrue(id.childEnd() > id)
        XCTAssertEqual(id.childBegin().next().next().next().next(), id.childEnd())
        XCTAssertEqual(id.childBegin(level: S2CellId.maxLevel), id.rangeMin)
        XCTAssertEqual(id.childEnd(level: S2CellId.maxLevel), id.rangeMax.next())
        
        // Check wrapping from beginning of Hilbert curve to end and vice versa.
        XCTAssertEqual(S2CellId.begin(level: 0).prevWrap(), S2CellId.end(level: 0).prev())
        // The original Java test used >>> which bit shifts all including the sign bit.
        // We use >> on UInt64 which has no sign bit then convert to Int64.
        XCTAssertEqual(S2CellId.begin(level: S2CellId.maxLevel).prevWrap(),
            S2CellId(face: 5, pos: Int64(bitPattern: UInt64(bitPattern: ~0) >> S2CellId.faceBits), level: S2CellId.maxLevel))
        
        XCTAssertEqual(S2CellId.end(level: 4).prev().nextWrap(), S2CellId.begin(level: 4))
        XCTAssertEqual(S2CellId.end(level: S2CellId.maxLevel).prev().nextWrap(),
                     S2CellId(face: 0, pos: 0, level: S2CellId.maxLevel))
        
        // Check that cells are represented by the position of their center
        // along the Hilbert curve.
        XCTAssertEqual(id.rangeMin.uid + id.rangeMax.uid, 2 * id.uid)
    }
}
