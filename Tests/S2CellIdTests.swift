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
        XCTAssertFalse(id.isValid)
        
        id = S2CellId(face: 3, pos: 0x12345678, level: S2CellId.maxLevel - 4)
        XCTAssertTrue(id.isValid)
        XCTAssertEqual(id.face, 3)
        XCTAssertEqual(id.pos, 0x12345700)
        XCTAssertEqual(id.level, S2CellId.maxLevel - 4)
        XCTAssertFalse(id.isLeaf)
        
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
        XCTAssertLessThan(id.childBegin(), id)
        XCTAssertGreaterThan(id.childEnd(), id)
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
    
    func testInverses() {
        NSLog("TestInverses")
        
        // Check the conversion of random leaf cells to S2LatLngs and back.
        for _ in 0...20 { // 20000 to 20 for quicker testing.
            let id = getRandomCellId(level: S2CellId.maxLevel)
            XCTAssertTrue(id.isLeaf && id.level == S2CellId.maxLevel)
            let center = id.latLng
            XCTAssertEqual(S2CellId(latlng: center).id, id.id)
        }
    }
    
    func testToToken() {
        XCTAssertEqual("000000000000010a", S2CellId(id: 266).token)
        XCTAssertEqual("80855c", S2CellId(id: -9185834709882503168).token)
    }
    
    func testTokens() throws {
        NSLog("TestTokens")
        
        // Test random cell ids at all levels.
        for _ in 0...20 { // 20000 to 20 for quicker testing.
            let id = getRandomCellId(level: S2CellId.maxLevel)
            if !id.isValid {
                continue
            }
            let token = id.token
            XCTAssertLessThanOrEqual(token.count, 16)
            XCTAssertEqual(try S2CellId(token: token), id)
        }
        // Check that invalid cell ids can be encoded.
        let token = S2CellId.none.token
        XCTAssertEqual(try S2CellId(token: token), S2CellId.none)
    }
    
    static let kMaxExpandedLevel = 3
    
    private func expandCell(parent: S2CellId, cells: [S2CellId], parentMap: [S2CellId: S2CellId]) {
        var parentMap2 = parentMap
        guard parent.level != S2CellIdTests.kMaxExpandedLevel else {
            return
        }
        var i = 0
        var j = 0
        var orientation: Int? = 0
        let face = parent.toFaceIJOrientation(i: &i, j: &j, orientation: &orientation)
        XCTAssertEqual(face, parent.face)
        
        var pos = 0
        var child = parent.childBegin()
        while child != parent.childEnd() {
            XCTAssertEqual(child.level, parent.level + 1)
            XCTAssertFalse(child.isLeaf)
            var childOrientation: Int? = 0
            XCTAssertEqual(child.toFaceIJOrientation(i: &i, j: &j, orientation: &childOrientation), face)
            XCTAssertEqual(childOrientation!, orientation! ^ S2.posToOrientation(position: pos))
            
            parentMap2[child] = parent
            expandCell(parent: child, cells: cells, parentMap: parentMap2)
            pos = pos + 1
            child = child.next()
        }
    }
    
    func testContainment() {
        NSLog("TestContainment")
        
        var parentMap = [S2CellId: S2CellId]()
        var cells = [S2CellId]()
        for face in 0..<6 {
            expandCell(parent: S2CellId(face: face, pos: 0, level: 0), cells: cells, parentMap: parentMap)
        }
        for i in 0..<cells.count {
            for j in 0..<cells.count {
                var contained = true
                var id = cells[j]
                while id != cells[i] {
                    guard parentMap[id] != nil else {
                        contained = false
                        break
                    }
                    id = parentMap[id]!
                }
                XCTAssertEqual(cells[i].contains(other: cells[j]), contained)
                XCTAssertEqual(cells[j] >= cells[i].rangeMin && cells[j] <= cells[i].rangeMax, contained)
                XCTAssertEqual(cells[i].intersects(with: cells[j]) || cells[i].contains(other: cells[j]), cells[j].contains(other: cells[i]))
            }
        }
    }
    
    static let maxWalkLevel = 4 // 8 to 4 for quicker testing.
    
    func testContinuity() {
        NSLog("TestContinuity")
        
        // Make sure that sequentially increasing cell ids form a continuous
        // path over the surface of the sphere, i.e. there are no
        // discontinuous jumps from one region to another.
        let maxDist = S2Projections.maxEdge.getValue(level: S2CellIdTests.maxWalkLevel)
        let end = S2CellId.end(level: S2CellIdTests.maxWalkLevel)
        var id = S2CellId.begin(level: S2CellIdTests.maxWalkLevel)
        while id != end {
            let p = id.rawPoint
            XCTAssertLessThanOrEqual(p.angle(to: id.nextWrap().rawPoint), maxDist)
            
            // Check that the ToPointRaw() returns the center of each cell
            // in (s,t) coordinates.
            let face = S2Projections.xyzToFace(point: p)
            let uv = S2Projections.validFaceXyzToUv(face: face, point: p)
            let d = (1.0 / Double(1 << S2CellIdTests.maxWalkLevel))
            XCTAssertEqual(S2Projections.uvToST(u: uv.x).remainder(dividingBy: d), 0, accuracy: 1.0)
            XCTAssertEqual(S2Projections.uvToST(u: uv.y).remainder(dividingBy: d), 0, accuracy: 1.0)
            id = id.next()
        }
    }
    
    func testCoverage() {
        NSLog("TestCoverage")
        
        // Make sure that random points on the sphere can be represented to the
        // expected level of accuracy, which in the worst case is sqrt(2/3) times
        // the maximum arc length between the points on the sphere associated with
        // adjacent values of "i" or "j". (It is sqrt(2/3) rather than 1/2 because
        // the cells at the corners of each face are stretched -- they have 60 and
        // 120 degree angles.)
        let maxDist = 0.5 * S2Projections.maxDiag.getValue(level: S2CellId.maxLevel)
        for _ in 0..<2 {
            let p = S2Point(x: 0.37861576725894824, y: 0.2772406863275093, z: 0.8830558887338725)
            let q = S2CellId(point: p).rawPoint
            XCTAssertLessThanOrEqual(p.angle(to: q), maxDist)
        }
    }
    
    func testAllNeighbours(id: S2CellId, level: Int) {
        XCTAssertTrue(level >= id.level && level < S2CellId.maxLevel)
        
        // We compute GetAllNeighbors, and then add in all the children of "id"
        // at the given level. We then compare this against the result of finding
        // all the vertex neighbors of all the vertices of children of "id" at the
        // given level. These should give the same result.
        var all: Set<S2CellId> = Set()
        for id in id.getAllNeighbors(level: level) {
            all.insert(id)
        }
        var expected: Set<S2CellId> = Set()
        
        let end = id.childEnd(level: level + 1)
        var c = id.childBegin(level: level + 1)
        while c != end {
            all.insert(c.parent)
            for id in c.getVertexNeighbors(level: level) {
                expected.insert(id)
            }
            c = c.next()
        }
        XCTAssertEqual(all, expected)
    }
    
    func testNeighbours() {
        NSLog("TestNeighbours")
        
        // Check the edge neighbors of face 1.
        let outFaces = [5, 3, 2, 0]
        let faceNmbrs = S2CellId(face: 1, pos: 0, level: 0).getEdgeNeighbors()
        for i in 0..<4 {
            XCTAssertTrue(faceNmbrs[i].isFace)
            XCTAssertEqual(faceNmbrs[i].face, outFaces[i])
        }
        
        // Check the vertex neighbors of the center of face 2 at level 5.
        var nbrs = S2CellId(point: S2Point(x: 0, y: 0, z: 1)).getVertexNeighbors(level: 5);
        nbrs.sort()
        for i in 0..<4 {
            XCTAssertEqual(nbrs[i], S2CellId(face: 2, i: (1 << 29) - (i < 2 ? 1 : 0), j: (1 << 29) - ((i == 0 || i == 3) ? 1 : 0)).parent(level: 5))
        }
        
        // Check the vertex neighbors of the corner of faces 0, 4, and 5.
        let id = S2CellId(face: 0, pos: 0, level: S2CellId.maxLevel)
        nbrs = id.getVertexNeighbors(level: 0);
        nbrs.sort()
        XCTAssertEqual(nbrs.count, 3);
        XCTAssertEqual(nbrs[0], S2CellId(face: 0, pos: 0, level: 0))
        XCTAssertEqual(nbrs[1], S2CellId(face: 4, pos: 0, level: 0))
        XCTAssertEqual(nbrs[2], S2CellId(face: 5, pos: 0, level: 0))
        
        // Check that GetAllNeighbors produces results that are consistent
        // with GetVertexNeighbors for a bunch of random cells.
        for _ in 0..<10 { // 10000 to 10 for quicker testing.
            var id1 = getRandomCellId()
            if id1.isLeaf {
                id1 = id1.parent
            }
            
            // TestAllNeighbors computes approximately 2**(2*(diff+1)) cell id1s,
            // so it's not reasonable to use large values of "diff".
            let maxDiff = min(6, S2CellId.maxLevel - id1.level - 1)
            let level = id1.level + Int.random(in: 0...maxDiff)
            testAllNeighbours(id: id1, level: level)
        }
    }
}
