// swift-tools-version:4.0
import PackageDescription

let package = Package(
    name: "S2GeometrySwift",
    products: [
        .library(name: "S2GeometrySwift", targets: ["S2GeometrySwift"]),
    ],
    targets: [
        .target(name: "S2GeometrySwift", path: "Sources"),
        .testTarget(name: "S2GeometrySwiftTests", dependencies: ["S2GeometrySwift"], path: "Tests")
    ]
)
