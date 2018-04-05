#
# Be sure to run `pod lib lint S2GeometrySwift.podspec' to ensure this is a
# valid spec before submitting.
#
# Any lines starting with a # are optional, but their use is encouraged
# To learn more about a Podspec see http://guides.cocoapods.org/syntax/podspec.html
#

Pod::Spec.new do |s|
  s.name             = 'S2GeometrySwift'
  s.version          = '1.0.0'
  s.summary          = 'Geometry library in Swift.'

# This description is used to generate tags and improve search results.
#   * Think: What does it do? Why did you write it? What is the focus?
#   * Try to keep it short, snappy and to the point.
#   * Write the description between the DESC delimiters below.
#   * Finally, don't worry about the indent, CocoaPods strips it!

  s.description      = 'Contains a Swift import of Googles S2 Geometry Library'

  s.homepage         = 'https://github.com/philip-bui/s2-geometry-swift'
  # s.screenshots     = 'www.example.com/screenshots_1', 'www.example.com/screenshots_2'
  s.license          = { :type => 'MIT', :file => 'LICENSE' }
  s.author           = { 'Philip' => 'philip.bui.developer@gmail.com' }
  s.source           = { :git => 'https://github.com/philip-bui/s2-geometry-swift.git', :tag => s.version.to_s }
  # s.social_media_url = 'https://twitter.com/<TWITTER_USERNAME>'

  s.ios.deployment_target = '8.0'

  s.source_files = 'S2GeometrySwift/Classes/**/*'
  
  # s.resource_bundles = {
  #   'S2GeometrySwift' => ['S2GeometrySwift/Assets/*.png']
  # }

  # s.public_header_files = 'Pod/Classes/**/*.h'
  # s.frameworks = 'UIKit', 'MapKit'
  # s.dependency 'AFNetworking', '~> 2.3'
end
