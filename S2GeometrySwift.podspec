Pod::Spec.new do |s|
  s.name = 'S2GeometrySwift'
  s.version = '1.0.2'
  s.license= { :type => 'MIT', :file => 'LICENSE' }
  s.summary = 'S2 Geometry library in Swift.'
  s.description = 'Swift port of S2 Geometry by Google.'
  s.homepage = 'https://github.com/philip-bui/s2-geometry-swift'
  s.author = { 'Philip Bui' => 'philip.bui.developer@gmail.com' }
  s.source = { :git => 'https://github.com/philip-bui/s2-geometry-swift.git', :tag => s.version }
  s.documentation_url = 'https://github.com/philip-bui/s2-geometry-swift'

  s.ios.deployment_target = '8.0'
  s.osx.deployment_target = '10.9'
  s.tvos.deployment_target = '9.0'
  s.watchos.deployment_target = '2.0'
 
  s.source_files = 'Sources/Classes/**/*'
end
